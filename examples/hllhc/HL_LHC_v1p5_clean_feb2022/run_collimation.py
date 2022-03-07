import os
import sys
import time
import glob
import gzip
import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import namedtuple
from itertools import cycle, islice, dropwhile
from multiprocessing import Pool
from pathlib import Path

import yaml
import pandas as pd
from IPython import embed

import xobjects as xo
import xtrack as xt
import xpart as xp
import collimasim as cs

from contextlib import redirect_stdout, redirect_stderr

PART_DATA_VARS = ['s', 'x', 'y', 'px', 'py', 'zeta', 'delta',
                  'charge_ratio', 'weight', 'particle_id',
                  'at_element', 'at_turn', 'state',
                  'parent_particle_id']


HLLHC_WARM_REGIONS = np.array([
    [0.00000000e+00, 2.25000000e+01],
    [8.31530000e+01, 1.36689000e+02],
    [1.82965500e+02, 2.01900000e+02],
    [2.10584700e+02, 2.24300000e+02],
    [3.09545428e+03, 3.15562858e+03],
    [3.16774008e+03, 3.18843308e+03],
    [3.21144458e+03, 3.26386758e+03],
    [3.30990008e+03, 3.35497408e+03],
    [3.40100558e+03, 3.45342858e+03],
    [3.47644008e+03, 3.49406558e+03],
    [3.50588528e+03, 3.56831858e+03],
    [6.40540880e+03, 6.45791380e+03],
    [6.46877850e+03, 6.85951380e+03],
    [6.87037850e+03, 6.92353380e+03],
    [9.73590702e+03, 9.82473052e+03],
    [9.83083202e+03, 9.86173052e+03],
    [9.87873202e+03, 9.93998552e+03],
    [9.95054802e+03, 1.00434620e+04],
    [1.00540245e+04, 1.01152780e+04],
    [1.01322795e+04, 1.01639705e+04],
    [1.01700720e+04, 1.02576030e+04],
    [1.31036000e+04, 1.31200300e+04],
    [1.31238892e+04, 1.31471237e+04],
    [1.31918002e+04, 1.32476472e+04],
    [1.33067940e+04, 1.33520892e+04],
    [1.34110312e+04, 1.34670082e+04],
    [1.35114547e+04, 1.35357845e+04],
    [1.35388592e+04, 1.35552845e+04],
    [1.63946378e+04, 1.64508713e+04],
    [1.64569728e+04, 1.64872713e+04],
    [1.64933728e+04, 1.68308713e+04],
    [1.68369728e+04, 1.68672713e+04],
    [1.68733728e+04, 1.69282948e+04],
    [1.97348504e+04, 1.97606997e+04],
    [1.97715644e+04, 2.02179087e+04],
    [2.02287734e+04, 2.02529744e+04],
    [2.30899797e+04, 2.31385770e+04],
    [2.31503967e+04, 2.31713755e+04],
    [2.31943870e+04, 2.32468100e+04],
    [2.32928425e+04, 2.33379155e+04],
    [2.33839480e+04, 2.34363710e+04],
    [2.34593825e+04, 2.34800825e+04],
    [2.34921940e+04, 2.35531160e+04],
    [2.64334879e+04, 2.64483032e+04],
    [2.64569832e+04, 2.64759232e+04],
    [2.65221932e+04, 2.65757332e+04],
    [2.66363832e+04, 2.66588832e+04],
])


C_LIGHT = 299792458  # m / s

# A thin container to hold some information about a particle
Particle = namedtuple('Particle', ['mass', 'momentum', 'kinetic_energy',
                                   'total_energy', 'charge', 'pdgid'])

# Some utility functions
_kinetic_energy = lambda mass, energy: energy - mass
_momentum = lambda mass, energy: np.sqrt(energy ** 2 - mass ** 2)
_energy = lambda mass, momentum: np.sqrt(momentum ** 2 + mass ** 2)


def cycle_line(line, name):
    names = line.element_names
    elements = line.elements

    if name not in names:
        raise ValueError(f'No element name {name} found in the line.')

    # Store the original s of the element
    s0 = line.get_s_elements(mode='upstream')[names.index(name)]

    zipper = zip(names, elements)
    cycler = cycle(zipper) # get an infinite cyclic iterator
    skipper = dropwhile(lambda n: n[0] != name, cycler) # move to the desired name
    cycled = list(islice(skipper, None, len(names))) # the cycled list
    names_cyc, elements_cyc = zip(*cycled)

    line.element_names = list(names_cyc)
    line.elements = list(elements_cyc)

    return line, s0


def _save_particles_hdf(particles=None, lossmap_data=None, filename='part'):
    if not filename.endswith('.hdf'):
        filename += '.hdf'

    fpath = Path(filename)
    # Remove a potential old file as the file is open in append mode
    if fpath.exists():
        fpath.unlink()

    if particles is not None:
        df = particles.to_pandas(compact=True)
        df.to_hdf(fpath, key='particles', format='table', mode='a',
                  complevel=9, complib='blosc')


    if lossmap_data is not None:
        for key, lm_df in lossmap_data.items():
            lm_df.to_hdf(fpath, key=key, mode='a', format='table',
                         complevel=9, complib='blosc')


def _load_particles_hdf(filename):
    return xp.Particles.from_pandas(pd.read_hdf(filename, key='particles'))


def _read_particles_hdf(filename):
    return pd.read_hdf(filename, key='particles')


def _load_lossmap_hdf(filename):
    keys = ('lossmap_scalar', 'lossmap_aper', 'lossmap_coll')

    lm_dict = {}
    for key in keys:
        # Pandas HDF file table format doesn't save empty dataframes
        try:
            lm_dict[key] = pd.read_hdf(filename, key=key)
        except KeyError:
            lm_dict[key] = None
    return lm_dict


def check_warm_loss(s, warm_regions):
    return np.any((warm_regions.T[0] < s) & (warm_regions.T[1] > s))


def load_config(config_file):
    with open(config_file, 'r') as stream:
        config_dict = yaml.safe_load(stream)
    return config_dict


def load_and_process_line(config_dict):

    beam = config_dict['beam']
    inp = config_dict['input']
    run = config_dict['run']

    # Convert the parameters that are not strings
    try:
        emittance = (float(beam['emittance']['x']), float(beam['emittance']['y']))
    except TypeError:
        emittance = float(beam['emittance'])

    pdg_id = int(beam['pdg_id'])
    mass = float(beam['mass'])
    p0 = float(beam['momentum'])
    start_element = beam.get('start_element')

    energy_cut = float(run['energy_cut'])
    seed = int(float(run['seed']))
    batch_mode = bool(run['batch_mode'])

    mat_rename = inp.get('material_rename_map', {})

    ke0 = _kinetic_energy(mass, p0)
    e0 = _energy(mass, p0)

    g4man = cs.Geant4CollimationManager(collimator_file=inp['collimator_file'],
                                        bdsim_config_file=inp['bdsim_config'],
                                        tfs_file=inp['tfs_file'],
                                        reference_pdg_id=pdg_id,
                                        reference_kinetic_energy=ke0,
                                        emittance_norm=emittance,
                                        relative_energy_cut=energy_cut,
                                        seed=seed,
                                        material_rename_map=mat_rename,
                                        batchMode=batch_mode,
                                        )

    # load the line
    with open(inp['xtrack_line'], 'r') as fid:
        line = xt.Line.from_dict(json.load(fid))

    g4man.place_all_collimators(line)

    # Get the S of the starting element for transform purposes later
    #idx0 = line.element_names.index(start_element)
    #s0 = np.array(line.get_s_elements(mode='downstream'))[idx0]

    if start_element:
        line, s0 = cycle_line(line, start_element)

    ref_part = Particle(mass, p0, ke0, e0, 1, beam['pdg_id'])

    return line, ref_part, start_element, s0


def load_gpdist_distr(dist_file, p0c, capacity=5e4):
    # Load a file with initial coordinates from gpdist and convert it to MADX inrays format

    # Be careful of the header, it may change
    names = ['pid', 'genid', 'weight', 'x', 'y', 'z', 'xp', 'yp', 'zp', 'A', 'Z', 'm', 'p', 't']
    #names = ['pid', 'genid', 'weight', 'x', 'y', 'xp', 'yp', 'm', 'p', 't', 'q', 'g4pid']

    coords = pd.read_csv(dist_file, delim_whitespace=True, index_col=False, names=names)

    loadpart = None
    if loadpart is not None:
        print(f'Running only {loadpart} particles.')

    coords= coords.iloc[0:loadpart]

    # transform the coorindates
    coords['delta'] = (coords['p'] - p0c)/ p0c
    coords['zeta'] = coords['t'] * C_LIGHT # TODO: check that zeta is ct, minus sign needed?

    # TODO: Are the px, py coordinates the same as the ones in xtrack?
    coords.rename({'xp': 'px', 'yp': 'py'}, axis=1, inplace=True)

    xtrack_columns = ['x', 'px', 'y', 'py', 'delta', 'zeta']
    coords = coords[xtrack_columns]

    particles = xp.Particles(
        _capacity=capacity,
		p0c=p0c * 1e3, # In MeV here
        x=coords['x'],
        px=coords['px'],
        y=coords['y'],
        py=coords['py'],
        zeta=coords['zeta'],
        #delta=coords['delta'], # TODO: change back after fix
    )

    particles.delta[:len(coords['delta'])] = coords['delta']

    return particles


def prepare_particles(config_dict):
    dist = config_dict['dist']
    capacity = int(float(config_dict['run']['max_particles']))
    p0c = float(config_dict['beam']['momentum'])

    _supported_dist = ['gpdist']

    if dist['kind'] == 'gpdist':
        particles = load_gpdist_distr(dist['file'], p0c=p0c, capacity=capacity)
    else:
        raise ValueError('Unsupported distribution kind: {}. Supported ones are: {}'
                         .format(dist['kind'], ','.join(_supported_dist)))

    return particles


def run(config_dict, line, particles, s0):
    ## Chose a context
    context = xo.ContextCpu() # Only support CPU for Geant4 coupling

    ## Transfer lattice on context and compile tracking code
    tracker = xt.Tracker(_context=context, line=line)

    nturns = config_dict['run'].get('turns', 1)

    t0 = time.time()

    ## Track (saving turn-by-turn data)
    print(f'{particles._num_active_particles=}')
    print(f'{particles._num_lost_particles=}')
    print('turn=0')

    for turn in range(nturns):
        tracker.track(particles, num_turns=1)
        print(f'{particles._num_active_particles=}')
        print(f'{particles._num_lost_particles=}')
        print(f'{turn=}')

        if particles._num_active_particles==0:
            break

    print(f'Tracking {nturns} turns done in: {time.time() -t0} s')

    mask_allocated = particles.state > -100

    s_ele = np.array(line.get_s_elements(mode='downstream'))

    s_range = (0, max(s_ele)) # Get the end point as assume the start point is zero
    particles.s[mask_allocated] = np.take(s_ele, particles.at_element[mask_allocated])

    particles_before = particles.copy()

    #import pdb; pdb.set_trace()
    loss_refiner = xt.LossLocationRefinement(tracker, n_theta=360, r_max=1, dr=50e-6, ds=0.1)
    loss_refiner.refine_loss_location(particles)

    if ('lossmap' in config_dict
        and config_dict['lossmap'].get('make_lossmap', False)):
        binwidth = config_dict['lossmap']['aperture_binwidth']
        lossmap_data = prepare_lossmap(particles, line, s0, binwidth=binwidth)
    else:
        lossmap_data = None

    output_file = config_dict['run'].get('outputfile', 'part.hdf')
    _save_particles_hdf(particles, lossmap_data=lossmap_data, filename=output_file)

    # Save another file with uninterpolated losses
    # DEBUG
    binwidth = config_dict['lossmap']['aperture_binwidth']
    lossmap_data = prepare_lossmap(particles_before, line, s0, binwidth=binwidth)
    _save_particles_hdf(particles, lossmap_data=lossmap_data, filename="Outputdata/part_nointerp.hdf")


def load_output(directory, output_file, match_pattern='*part.hdf*',
                imax=None, load_lossmap=True, load_particles=False):

    t0 = time.time()

    job_dirs = glob.glob(os.path.join(directory, 'Job.*')) # find directories to loop over

    job_dirs_sorted = []
    for i in range(len(job_dirs)):
        # Very inefficient, but it sorts the directories by their numerical index
        job_dir_idx = job_dirs.index(os.path.join(directory, 'Job.{}'.format(i)))
        job_dirs_sorted.append(job_dirs[job_dir_idx])

    part_hdf_files = []
    part_dataframes = []
    lossmap_dicts = []

    print(f'Parsing directories...')
    dirs_visited = 0
    files_loaded = 0
    for i, d in enumerate(job_dirs_sorted):
        if imax is not None and i > imax:
            break

        #print(f'Processing {d}')
        dirs_visited += 1
        output_dir = os.path.join(d, 'Outputdata')
        output_files = glob.glob(os.path.join(output_dir, match_pattern))
        if output_files:
            of = output_files[0]
            part_hdf_files.append(of)
            files_loaded += 1
        else:
            print(f'No output found in {d}')

    part_merged = None
    if load_particles:
        print(f'Loading particles...')
        p = Pool()
        part_dataframes = p.map(_read_particles_hdf, part_hdf_files)
        part_objects = [xp.Particles.from_pandas(pdf) for pdf in part_dataframes]
        print('Particles load finished, merging...')
        part_merged = xp.Particles.merge(part_objects)

    # Load the loss maps
    lmd_merged = None
    if load_lossmap:
        print(f'Loading loss map data...')
        p = Pool()
        lossmap_dicts = p.map(_load_lossmap_hdf, part_hdf_files)

        print('Loss map load finished, merging..')

        num_tol = 1e-9
        lmd_merged = lossmap_dicts[0]
        for lmd in lossmap_dicts[1:]:
            # Scalar parameters
            # Ensure consistency
            identical_params = ('s_min', 's_max', 'binwidth', 'nbins')
            for vv in identical_params:
                assert np.isclose(lmd_merged['lossmap_scalar'][vv],
                                  lmd['lossmap_scalar'][vv],
                                  num_tol)

            lmd_merged['lossmap_scalar']['n_primaries'] += lmd['lossmap_scalar']['n_primaries']

            # Collimator losses
            # These cannot be empty dataframes even if there is no losses
            assert np.isclose(lmd_merged['lossmap_coll']['coll_start'],
                              lmd['lossmap_coll']['coll_start'],
                              num_tol).all()

            assert np.isclose(lmd_merged['lossmap_coll']['coll_end'],
                              lmd['lossmap_coll']['coll_end'],
                              num_tol).all()

            lmd_merged['lossmap_coll']['coll_loss'] += lmd['lossmap_coll']['coll_loss']

            # Aperture losses
            alm = lmd_merged['lossmap_aper']
            al = lmd['lossmap_aper']

            # If the aperture loss dataframe is empty, it is not stored on HDF
            if al is not None:
                if alm is None:
                    lmd_merged['lossmap_aper'] = al
                else:
                    lm = alm.aper_loss.add(al.aper_loss, fill_value=0)
                    lmd_merged['lossmap_aper'] = pd.DataFrame({'aper_loss': lm})

    _save_particles_hdf(particles=part_merged, lossmap_data=lmd_merged, filename=output_file)

    print('Directories visited: {}, files loaded: {}'.format(dirs_visited, files_loaded))
    print(f'Processing done in {time.time() -t0} s')



def prepare_lossmap(particles, line, s0, binwidth):
    s_ele = np.array(line.get_s_elements(mode='downstream'))
    max_s = max(s_ele)

    # Transform particles back to the orginal S frame
    particles = particles.copy() # Copy the particles to avoid changing the coordinates

    #print("$#########!@Q#!@#!@#!@#@#!@#!@#HACKYFIX")    # DEBUG
    particles.s = (particles.s + s0) % max_s 

    # Transform elements back to the orginal S frame
    s_ele = (s_ele + s0) % max_s

    s_range = (0, max_s) # Get the end point as assume the start point is zero

    coll_idx = []
    for idx, elem in enumerate(line.elements):
        if isinstance(elem, xt.BeamInteraction) and isinstance(elem.interaction_process, cs.Geant4Collimator):
            coll_idx.append(idx)
    coll_idx = np.array(coll_idx)

    # The masks return array views, so no extra copying is performed
    mask_allocated = particles.state > -100
    particles = particles.filter(mask_allocated)

    mask_prim = particles.parent_particle_id == particles.particle_id
    n_prim = len(particles.filter(mask_prim).x)

    mask_part_type = abs(particles.chi - 1) < 1e-7 # Select particles same as the beam particle
    mask_lost = particles.state == 0
    particles = particles.filter(mask_part_type & mask_lost)

    mask_losses_coll = np.in1d(particles.at_element, coll_idx) # Get a mask for the collimator losses
    #mask_warm = np.array([check_warm_loss(s, HLLHC_WARM_REGIONS) for s in particles.s])

    # Collimator losses binned per element
    h_coll, edges_coll= np.histogram(particles.at_element[mask_losses_coll], bins=range(max(coll_idx)+2))

    # Process the collimator per element histogram for plotting
    coll_lengths = np.array([line.elements[ci].length for ci in coll_idx])
    coll_values = np.take(h_coll, coll_idx) # reduce the empty bars in the histogram

    coll_end = np.take(s_ele, coll_idx)
    coll_start = coll_end - coll_lengths

    # Aperture losses binned in S
    nbins_ap = int(np.ceil((s_range[1] - s_range[0])/binwidth))
    bins_ap = np.linspace(s_range[0], s_range[1], nbins_ap)

    aper_loss, _ = np.histogram(particles.s[~mask_losses_coll], bins=bins_ap)

    # Prepare structures for optimal storage
    aper_loss_series = pd.Series(aper_loss)

    # Scalar variables go in their own DF to avoid replication
    # The bin edges can be re-generated with linspace, no need to store
    scalar_dict = {
        'binwidth': binwidth,
        'nbins': nbins_ap,
        's_min': s_range[0],
        's_max': s_range[1],
        'n_primaries': n_prim,
    }

    # Drop the zeros while preserving the index
    aperloss_dict = {
        'aper_loss': aper_loss_series[aper_loss_series > 0],
    }

    coll_dict = {
        'coll_start': coll_start,
        'coll_end': coll_end,
        'coll_loss': coll_values
        }

    scalar_df = pd.DataFrame(scalar_dict, index=[0])
    coll_df = pd.DataFrame(coll_dict)
    aper_df = pd.DataFrame(aperloss_dict)

    lm_dict = {'lossmap_scalar': scalar_df,
               'lossmap_aper': aper_df,
               'lossmap_coll': coll_df}

    return lm_dict

def plot_lossmap(lossmap_file, extra_ranges=None):
    # A dict of three dataframes
    lossmap_data = _load_lossmap_hdf(lossmap_file)

    lms = lossmap_data['lossmap_scalar']
    lma = lossmap_data['lossmap_aper']
    lmc = lossmap_data['lossmap_coll']

    s_min = lms['s_min'][0]
    s_max = lms['s_max'][0]
    nbins = lms['nbins'][0]
    binwidth = lms['binwidth'][0]
    s_range = (s_min, s_max)

    # Collimator losses
    coll_start = lmc['coll_start']
    coll_end = lmc['coll_end']
    coll_values = lmc['coll_loss']

    coll_lengths = coll_end - coll_start
    norm = sum(coll_values)

    coll_values /= (norm * coll_lengths)

    # There can be an alternative way of plotting using a bar plot
    # Make the correct edges to get the correct width of step plot
    # The dstack and flatten merges the arrays one set of values at a time
    zeros = np.full_like(lmc.index, 0) # Zeros to pad the bars
    coll_edges = np.dstack([coll_start, coll_start, coll_end, coll_end]).flatten()
    coll_loss = np.dstack([zeros, coll_values, coll_values, zeros]).flatten()

    # Aperture losses
    aper_edges =  np.linspace(s_min, s_max, nbins)

    aper_loss = lma['aper_loss'].reindex(range(0, nbins-1), fill_value=0)
    #warm_loss = lma['warm_loss'].reindex(range(0, nbins-1), fill_value=0)

    aper_loss /= (norm * binwidth)

    # Check if the start of the bin is in a warm region
    mask_warm = np.array([check_warm_loss(s, HLLHC_WARM_REGIONS) for s in aper_edges[:-1]])

    warm_loss = aper_loss * mask_warm
    cold_loss = aper_loss * ~mask_warm

    # Make the plots
    fig, ax = plt.subplots(figsize=(12, 4))

    # The zorder determines the plotting order = warm -> collimator-> cold (on top)
    # The edge lines on the plots provide the dynamic scaling feature of the plots
    # e.g a 10 cm aperture loss spike is still resolvable for a full ring view
    lw=1
    ax.stairs(warm_loss, aper_edges, color='r', lw=lw, ec='r', fill=True, zorder=20)
    ax.stairs(cold_loss, aper_edges, color='b', lw=lw, ec='b', fill=True, zorder=30)

    ax.fill_between(coll_edges, coll_loss, step='pre', color='k', zorder=9)
    ax.step(coll_edges, coll_loss, color='k', lw=lw, zorder=10)

    ax.set_xlim(s_range[0], s_range[1])
    ax.set_ylim(1e-7, 2)

    ax.yaxis.grid(b=True, which='major', zorder=0)
    ax.yaxis.grid(b=True, which='minor', zorder=0)

    ax.set_yscale('log', nonpositive='clip')
    ax.set_xlabel(r'S [$m$]')
    ax.set_ylabel(r'Cleaning inefficiency [$m^{-1}$]')

    #plt.show()
    #return

    plt.savefig('loss_map_full.pdf', bbox_inches='tight')

    if extra_ranges is not None:
        for srange in extra_ranges:
            ax.set_xlim(srange)
            plt.savefig("lossmap_{}_{}.pdf".format(srange[0], srange[1]))


def main():
    if len(sys.argv) != 2:
        raise ValueError('The script takes only one input - the configuration file')

    config_file = sys.argv[1]
    config_dict = load_config(config_file)

    line, ref_part, start_elem, s0 = load_and_process_line(config_dict)
    particles = prepare_particles(config_dict)

    output_file = config_dict['run'].get('outputfile', 'part.hdf')

    run(config_dict, line, particles, s0) # modifies the Particles object in place

    #plot_lossmap(output_file)


if __name__ == '__main__':
    with open('output.log', 'w') as of:
        with redirect_stdout(of):
            main()
