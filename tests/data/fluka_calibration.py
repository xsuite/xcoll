import json
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


def __main__():
    result_file = Path('data/fluka_calibration_data.json')
    data_file = Path('data/fluka_calibration_losses.npz')
    if result_file.exists() and data_file.exists():
        with open(result_file, 'r') as fp:
            final_result = json.load(fp)
        npzfile = np.load(data_file)
        E_diff_ionisation = npzfile['ionisation']
        E_diff_full = npzfile['full']
        print("Loaded existing calibration data from files.")

    else:
        final_result, E_diff_ionisation, E_diff_full = create_calibration_data()
        with result_file.open('w') as fp:
            json.dump(final_result, fp, cls=xo.JEncoder, indent=4)
        np.savez_compressed(data_file, ionisation=E_diff_ionisation, full=E_diff_full)

    plot_calibration_data(E_diff_ionisation, E_diff_full)


def create_calibration_data():
    num_part = 20_000
    batches = 500
    MASTER_SEED = 123456789
    rng = np.random.default_rng(MASTER_SEED)
    coll = xc.FlukaCollimator(length=0.6, angle=0, jaw=0.001, material='MG6403Fc')

    particle_ref = xt.Particles('proton', p0c=6.8e12)
    xc.fluka.engine.particle_ref = particle_ref
    xc.fluka.engine.capacity = 2*num_part
    xc.fluka.engine.include_elastic = False
    xc.fluka.engine.include_inelastic = False
    xc.fluka.engine.include_showers = False
    xc.fluka.engine.include_single_coulomb = False
    xc.fluka.engine.include_multiple_coulomb = False

    final_result = {
        "num_particles": num_part*batches,
        "num_batches": batches,
        "batch_size": num_part,
        "seed": MASTER_SEED
    }

    # Only collisional ionisation loss  => mean stopping power
    xc.fluka.engine.seed = 44444444
    xc.fluka.engine.include_ionisation_fluctuations = False
    xc.fluka.engine.include_pair_production = False
    xc.fluka.engine.include_bremsstrahlung = False
    xc.fluka.engine.start(elements=coll, clean=True, verbose=False, fortran_debug_level=0)
    part = xp.build_particles(
        x=np.random.uniform(0.002-1e-6, 0.002+1e-6, num_part),
        px=np.random.uniform(-1e-6, 1e-6, num_part),
        y=np.random.uniform(-1e-6, 1e-6, num_part),
        py=np.random.uniform(-1e-6, 1e-6, num_part),
        particle_ref=xc.fluka.engine.particle_ref,
        _capacity=xc.fluka.engine.capacity)
    coll.track(part)
    xc.fluka.engine.stop(clean=True)

    mean_stopping_power = np.unique(part.energy0[part.state > 0] - part.energy[part.state > 0])
    assert len(mean_stopping_power) == 1, "Mean stopping power is not unique"
    # We expect 410.12 MeV, which corresponds to 1/rho dE/dx of 2.68 MeV cm2/g
    final_result['mean_stopping_power'] = mean_stopping_power[0]
    print()
    print(f"RESULT only mean")
    print(f"    Mean stopping power: {mean_stopping_power[0]/1e6} MeV")
    print()

    # Regular ionisation losses
    xc.fluka.engine.include_ionisation_fluctuations = True
    xc.fluka.engine.include_pair_production = False
    xc.fluka.engine.include_bremsstrahlung = False
    E_diff_ionisation = _loop_tracking(coll, num_part, batches, rng)
    _analyse(E_diff_ionisation, num_part, batches, final_result, 'ionisation')

    # Full energy losses
    xc.fluka.engine.include_ionisation_fluctuations = True
    xc.fluka.engine.include_pair_production = True
    xc.fluka.engine.include_bremsstrahlung = True
    E_diff_full = _loop_tracking(coll, num_part, batches, rng)
    _analyse(E_diff_full, num_part, batches, final_result, 'full')

    return final_result, E_diff_ionisation, E_diff_full


def plot_calibration_data(E_diff_ionisation, E_diff_full):
    nbins = 1000
    E_min = min(E_diff_ionisation.min(), E_diff_full.min())/1e9
    E_max = max(E_diff_ionisation.max(), E_diff_full.max())/1e9
    bins = np.logspace(np.log10(E_min), np.log10(E_max), nbins + 1)
    bin_centres = np.sqrt(bins[:-1] * bins[1:])
    dlog = np.diff(np.log10(bins))

    _, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].hist(E_diff_ionisation/1e6, bins=nbins, range=(0, 1000), density=True, label='Ionisation only')
    ax[0].hist(E_diff_full/1e6, bins=nbins, range=(0, 1000), density=True, label='Full energy loss')
    ax[0].set_xlabel("Energy loss [MeV]")
    ax[0].set_ylabel("Probability density")
    ax[0].legend()

    counts, _ = np.histogram(E_diff_ionisation/1e9, bins=bins)
    dNdlogE = counts / (len(E_diff_ionisation) * dlog)
    ax[1].step(bin_centres, dNdlogE, where='mid', label='Ionisation only')
    counts, _ = np.histogram(E_diff_full/1e9, bins=bins)
    dNdlogE = counts / (len(E_diff_full) * dlog)
    ax[1].step(bin_centres, dNdlogE, where='mid', label='Full energy loss')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('Energy [GeV]')
    ax[1].set_ylabel(r'Normalised frequency $\frac{dN}{d\log E}$')
    ax[1].legend()

    plt.tight_layout()
    plt.show()


def _loop_tracking(coll, num_part, batches, rng):
    # Higher precision ionisation loss
    extra_card = "IONFLUCT         1.0       0.0       4.0  BLCKHOLE  @LASTMAT"
    data = np.zeros(batches*num_part)
    for bb in range(batches):
        if bb % 25 == 0:
            print(f"Tracking batch {bb} of {batches}")
        xc.fluka.engine.seed = int(rng.integers(1, 100_000_000))
        xc.fluka.engine.start(elements=coll, clean=True, verbose=False, fortran_debug_level=0, extra_cards=[extra_card])
        part = xp.build_particles(
            x=rng.uniform(0.002-1e-6, 0.002+1e-6, num_part),
            px=rng.uniform(-1e-6, 1e-6, num_part),
            y=rng.uniform(-1e-6, 1e-6, num_part),
            py=rng.uniform(-1e-6, 1e-6, num_part),
            particle_ref=xc.fluka.engine.particle_ref,
            _capacity=xc.fluka.engine.capacity)
        coll.track(part)
        xc.fluka.engine.stop(clean=True)
        E_diff = part.energy0[part.state > 0] - part.energy[part.state > 0]
        data[bb*num_part:(bb+1)*num_part] = E_diff
    return data

def _analyse(losses, num_part, batches, final_result, name):
    median = np.median(losses)
    mean = np.mean(losses)
    std = np.std(losses)
    q_lo, q_hi = np.quantile(losses, [0.01, 0.90])
    counts, edges = np.histogram(losses, bins=1000, range=(q_lo, q_hi))
    centres = 0.5 * (edges[:-1] + edges[1:])
    mpv = centres[np.argmax(counts)]
    tail_probs = np.array([1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5])
    quantiles = {str(pp): qq for pp, qq in zip(tail_probs, np.quantile(losses, 1 - tail_probs))}

    batch_mean = [np.mean(losses[bb*num_part:(bb+1)*num_part]) for bb in range(batches)]
    batch_std  = [np.std(losses[bb*num_part:(bb+1)*num_part]) for bb in range(batches)]
    batch_median = [np.median(losses[bb*num_part:(bb+1)*num_part]) for bb in range(batches)]
    batch_mpv = []
    for bb in range(batches):
        counts, edges = np.histogram(losses[bb*num_part:(bb+1)*num_part], bins=100, range=(q_lo, q_hi))
        centres = 0.5 * (edges[:-1] + edges[1:])
        batch_mpv.append(centres[np.argmax(counts)])
    batch_q_1 = [np.count_nonzero(losses[bb*num_part:(bb+1)*num_part] > quantiles['5e-05'])
                 for bb in range(batches)]
    batch_q_10 = [np.count_nonzero(losses[bb*num_part:(bb+1)*num_part] > quantiles['0.0005'])
                  for bb in range(batches)]
    batch_q_100 = [np.count_nonzero(losses[bb*num_part:(bb+1)*num_part] > quantiles['0.005'])
                   for bb in range(batches)]
    final_result[name] = {
        "mean": mean,
        "std": std,
        "median": median,
        "mpv": mpv,
        "1st percentile": q_lo,
        "90th percentile": q_hi,
        "thresholds": {
            1: quantiles['5e-05'],
            10: quantiles['0.0005'],
            100: quantiles['0.005']
        },
        "upper_tail_quantiles": quantiles,
        "batches": {
            "mean": batch_mean,
            "std": batch_std,
            "median": batch_median,
            "mpv": batch_mpv,
            "q_1": batch_q_1,
            "q_10": batch_q_10,
            "q_100": batch_q_100,
        }
    }
    print()
    print(f"RESULT {name}")
    print(f"    Mean: {mean/1e6} MeV")
    print(f"    Std: {std/1e6} MeV")
    print(f"    Median: {median/1e6} MeV")
    print(f"    Most probable value: {mpv/1e6} MeV")
    print(f"    1% quantile: {q_lo/1e6} MeV")
    print(f"    90% quantile: {q_hi/1e6} MeV")
    print(f"    Particle thresholds:")
    print(f"        1: {quantiles['5e-05']/1e6} MeV")
    print(f"        10: {quantiles['0.0005']/1e6} MeV")
    print(f"        100: {quantiles['0.005']/1e6} MeV")
    print(f"    Per-batch results:")
    print(f"        Mean: {np.mean(batch_mean)/1e6} MeV  [std: {np.std(batch_mean)/1e6} MeV, min: {np.min(batch_mean)/1e6} MeV, max: {np.max(batch_mean)/1e6} MeV]")
    print(f"        Std: {np.mean(batch_std)/1e6} MeV  [std: {np.std(batch_std)/1e6} MeV, min: {np.min(batch_std)/1e6} MeV, max: {np.max(batch_std)/1e6} MeV]")
    print(f"        Median: {np.mean(batch_median)/1e6} MeV  [std: {np.std(batch_median)/1e6} MeV, min: {np.min(batch_median)/1e6} MeV, max: {np.max(batch_median)/1e6} MeV]")
    print(f"        Most probable value: {np.mean(batch_mpv)/1e6} MeV  [std: {np.std(batch_mpv)/1e6} MeV, min: {np.min(batch_mpv)/1e6} MeV, max: {np.max(batch_mpv)/1e6} MeV]")
    print(f"        1-particle outlier: {np.mean(batch_q_1)}  [std: {np.std(batch_q_1)}, min: {np.min(batch_q_1)}, max: {np.max(batch_q_1)}]")
    print(f"        10-particle outlier: {np.mean(batch_q_10)}  [std: {np.std(batch_q_10)}, min: {np.min(batch_q_10)}, max: {np.max(batch_q_10)}]")
    print(f"        100-particle outlier: {np.mean(batch_q_100)}  [std: {np.std(batch_q_100)}, min: {np.min(batch_q_100)}, max: {np.max(batch_q_100)}]")
    print()



if __name__ == "__main__":
    __main__()
