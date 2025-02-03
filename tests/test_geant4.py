import numpy as np
import pytest
import time
# from scipy.stats import ks_2samp

import xobjects as xo
import xpart as xp
import xcoll as xc

from xobjects.test_helpers import for_all_test_contexts

# try the import here and skip tests if missing
# also need the import here in case of pytest --forked
try:
    import collimasim as cs
except ImportError:
    cs = None

path = xc._pkg_root.parent / 'tests' / 'data'


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_black_absorbers(test_context):
    n_part = 10000
    angles=[0,45,90]
    jaws = [0.03, -0.02]
    co = [-0.01, 0.01]
    L = 0.873

    g4_collimators = []
    ba_collimators = []
    for ii, angle in enumerate(angles):
        shift = co[0]*np.cos(angle) + co[1]*np.sin(angle)
        g4coll = xc.Geant4Collimator(length=L, angle=angle, jaw=jaws+shift,
                                     _context=test_context, material='cu',
                                     geant4_id=f'g4coll_{ii}')
        g4_collimators.append(g4coll)
        bacoll = xc.BlackAbsorber(length=L, angle=angle, jaw=jaws+shift,
                                  _context=test_context)
        ba_collimators.append(bacoll)

    xc.Geant4Engine.start(elements=g4_collimators, seed=1993, particle_ref='proton', p0c=7e12,
                          bdsim_config_file=str(path / f'geant4_ba_protons.gmad'))

    x = np.random.uniform(-0.1, 0.1, n_part)
    y = np.random.uniform(-0.1, 0.1, n_part)
    px = np.random.uniform(-0.1, 0.1, n_part)
    py = np.random.uniform(-0.1, 0.1, n_part)
    part = xp.build_particles(x=x, y=y, px=px, py=py, _context=test_context,
                              particle_ref=xc.Geant4Engine().particle_ref)
    part_ba = part.copy()
    # xc.Geant4Engine().particle_ref

    for coll in g4_collimators:
        coll.track(part)

    for coll in ba_collimators:
        coll.track(part_ba)

    part.sort(interleave_lost_particles=True)
    part_ba.sort(interleave_lost_particles=True)

    assert np.all(part.filter(part.state==1).particle_id 
                  == part_ba.filter(part_ba.state==1).particle_id)

    # Stop the Geant4 connection
    xc.Geant4Engine.stop(clean=True)  

@pytest.mark.parametrize('num_part', [1000, 5000])
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_simple_track(num_part):
    print(f"Running test_simple_track with {num_part} particles")
    _capacity = num_part*2

    # If a previous test failed, stop the server manually
    if xc.Geant4Engine.is_running():
        xc.Geant4Engine.stop(clean=True)

    # Define collimator and start the FLUKA server
    coll = xc.Geant4Collimator(length=0.6, jaw=0.001, material='cu')
    coll_name = 'tcp.c6l7.b1'
    xc.Geant4Engine.start(elements=coll, names=coll_name, seed=1993, particle_ref='proton', p0c=7e12,
                          bdsim_config_file=path / 'geant4_protons.gmad')

    # Particle distribution
    x_init   = np.random.normal(loc=0.002, scale=1e-3, size=num_part)
    px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                   particle_ref=particle_ref, _capacity=_capacity)
    part_geant4 = part_init.copy()
    part_drift = part_init.copy()

    # Geant4 tracking
    start = time.time()
    coll.track(part_geant4)
    geant4_time = round(time.time()-start, 3)
    expected_time = round(0.025 * num_part)    # To be adapted to reasonable Geant4 time
    if geant4_time > expected_time:
        raise RuntimeError(f"Geant4 tracking took {geant4_time}s, exceeding threshold of {expected_time}s")
    else:
        print(f"Geant4 tracking took {geant4_time}s")

    # Drift tracking
    coll._equivalent_drift.length = coll.length
    coll._equivalent_drift.track(part_drift)

    mask_geant4 = part_geant4.state > 0
    mask_drift = part_drift.state > 0
    # ks_stat, p_value = ks_2samp(part_geant4.x[mask_geant4], part_drift.x[mask_drift])
    # assert p_value <= 0.05, f"Distributions are not significantly different (p = {p_value})"
    # print(f"KS test passed with p = {p_value}")
    # mask_geant4 = part_geant4.state > 0
    # mask_drift = part_drift.state > 0
    # ks_stat, p_value = ks_2samp(part_geant4.x[mask_geant4], part_drift.x[mask_drift])
    # assert p_value <= 0.05, f"Distributions are not significantly different (p = {p_value})"
    # print(f"KS test passed with p = {p_value}")

    # Stop the Geant4 connection
    xc.Geant4Engine.stop(clean=True)


# @pytest.mark.parametrize('tilt', [0, [2.2e-6, 1.3e-6], [1.9e-6, -2.7e-6]],
#                          ids=['no_tilt', 'positive_tilt', 'pos_neg_tilt'])
@pytest.mark.parametrize('angle', [0, 90, 130.5])
@pytest.mark.parametrize('jaw', [0.001, [0.0013, -0.002789], [-1.2e-6, -3.2e-3], [3.789e-3, 4.678e-7]],
                         ids=['symmetric', 'asymmetric', 'negative', 'positive'])
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_jaw(jaw, angle): #, tilt):
    tilt = 0    # For now, need to implement test for tilted jaws
    _ACCURACY = 1.5e-9  # Anything in this region around the jaw might or might not hit; we can't be sure
    num_part = 5000
    _capacity = num_part*2
    jaw_band = 1.e-6

    # If a previous test failed, stop the server manually
    if xc.Geant4Engine.is_running():
        xc.Geant4Engine.stop(clean=True)

    # Define collimator and start the FLUKA server
    coll = xc.Geant4Collimator(length=0.6, jaw=jaw, angle=angle, tilt=tilt, material='Ti', geant4_id=f'g4coll_0')
    coll_name = 'tcp.c6l7.b1'
    xc.Geant4Engine.start(elements=coll, names=coll_name, seed=1993, particle_ref='proton', p0c=7.e12,
                          bdsim_config_file=path / 'geant4_protons.gmad')
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=7.e12)

    # Particle distribution (x and y are in the frame of the collimator)
    num_part_step = num_part//5
    x = np.random.uniform(-0.02, 0.02, num_part_step)
    x = np.concatenate([x, np.random.uniform(coll.jaw_L - jaw_band, coll.jaw_L -_ACCURACY, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_L +_ACCURACY, coll.jaw_L + jaw_band, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_R - jaw_band, coll.jaw_R -_ACCURACY, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_R +_ACCURACY, coll.jaw_R + jaw_band, num_part_step)])
    y = np.random.uniform(-0.02, 0.02, 5*num_part_step)
    x_new = np.cos(np.deg2rad(angle))*x - np.sin(np.deg2rad(angle))*y
    y_new = np.sin(np.deg2rad(angle))*x + np.cos(np.deg2rad(angle))*y
    part_init = xp.build_particles(x=x_new, y=y_new, particle_ref=xc.Geant4Engine().particle_ref,
                                   _capacity=_capacity)
    mask = np.concatenate([(x >= min(coll.jaw_LU, coll.jaw_LD)) | (x <= max(coll.jaw_RU, coll.jaw_RD)),
                          np.full(5*num_part_step, False)])
    hit_ids = part_init.particle_id[mask & (part_init.state > 0)]
    not_hit_ids = part_init.particle_id[~mask & (part_init.state > 0)]

    # TODO: jaw tilts, and particle angles

    # Track
    part = part_init.copy()
    coll.track(part)
    mask_hit = np.isin(part.parent_particle_id, hit_ids)
    mask_not_hit = np.isin(part.parent_particle_id, not_hit_ids)

    # mask = mask_not_hit & (part.state < 1)
    # wrong_ids = part.parent_particle_id[mask]
    # mask_wrong_init = np.isin(part_init.particle_id, wrong_ids)

    # mask_hit_init = np.isin(part_init.particle_id, hit_ids)
    # mask_not_hit_init = np.isin(part_init.particle_id, not_hit_ids)
    # plt.scatter(part_init.x[mask_not_hit_init], part_init.y[mask_not_hit_init], c='b', s=2)
    # plt.scatter(part_init.x[mask_hit_init], part_init.y[mask_hit_init], c='g', s=2)
    # plt.scatter(part_init.x[mask_wrong_init], part_init.y[mask_wrong_init], c='r', s=2)
    # plt.axvline(coll.jaw_L, c='k', linestyle='--')
    # plt.axvline(coll.jaw_R, c='k', linestyle='--')
    # plt.show()


    # Particles that are supposed to not have hit the collimator, but have a kick or are dead, are considered faulty
    assert not np.any(abs(part.px[mask_not_hit]) > _ACCURACY)
    assert not np.any(abs(part.py[mask_not_hit]) > _ACCURACY)
    assert not np.any(part.state[mask_not_hit] < 1)

    # Particles that are supposed to have hit the collimator, but are alive and have no kick, are considered faulty
    faulty =  mask_hit & (abs(part.px) < _ACCURACY) & (abs(part.py) < _ACCURACY)
    faulty &= (part.state > 0)
    assert len(part.x[faulty]) <= 1  # We allow for a small margin of error

    # Stop the Geant4 connection
    xc.Geant4Engine.stop(clean=True)






#         'CONDA_BACKUP_HOST': 'x86_64-conda-linux-gnu',
#         'HOST': 'x86_64-conda-linux-gnu',
#         'CONDA_BACKUP_CXX': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-c++',
#         'CONDA_BACKUP_GXX': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-g++',
#         'CONDA_BACKUP_CXXFLAGS': '-fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_DEBUG_CXXFLAGS': '-fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_CXX_FOR_BUILD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-c++',
#         'CONDA_BACKUP_GFORTRAN': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'CONDA_BACKUP_F95': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-f95',
#         'CONDA_BACKUP_FC_FOR_BUILD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'CONDA_BACKUP_FFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_FORTRANFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_DEBUG_FFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fcheck=all -fbacktrace -fimplicit-none -fvar-tracking-assignments -ffunction-sections -pipe -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_DEBUG_FORTRANFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fcheck=all -fbacktrace -fimplicit-none -fvar-tracking-assignments -ffunction-sections -pipe -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_FC': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'CONDA_BACKUP_F77': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'CONDA_BACKUP_F90': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'CONDA_BACKUP_BUILD': 'x86_64-conda-linux-gnu',
#         'CONDA_BACKUP_CONDA_TOOLCHAIN_HOST': 'x86_64-conda-linux-gnu',
#         'CONDA_BACKUP_CONDA_TOOLCHAIN_BUILD': 'x86_64-conda-linux-gnu',
#         'CONDA_BACKUP_CC': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-cc',
#         'CONDA_BACKUP_CPP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-cpp',
#         'CONDA_BACKUP_GCC': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc',
#         'CONDA_BACKUP_GCC_AR': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ar',
#         'CONDA_BACKUP_GCC_NM': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-nm',
#         'CONDA_BACKUP_GCC_RANLIB': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ranlib',
#         'CONDA_BACKUP_CPPFLAGS': '-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_CFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_LDFLAGS': '-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib',
#         'CONDA_BACKUP_DEBUG_CPPFLAGS': '-D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include -D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include -D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include -D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_DEBUG_CFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CONDA_BACKUP_CMAKE_PREFIX_PATH': '/apps/miniforge3/envs/Geant4:/apps/miniforge3/envs/Geant4/x86_64-conda-linux-gnu/sysroot/usr',
#         'CONDA_BACKUP__CONDA_PYTHON_SYSCONFIGDATA_NAME': '_sysconfigdata_x86_64_conda_cos6_linux_gnu',
#         'CONDA_BACKUP_CONDA_BUILD_SYSROOT': '/apps/miniforge3/envs/Geant4/x86_64-conda-linux-gnu/sysroot',
#         'CONDA_BACKUP_CC_FOR_BUILD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-cc',
#         'CONDA_BACKUP_build_alias': 'x86_64-conda-linux-gnu',
#         'CONDA_BACKUP_host_alias': 'x86_64-conda-linux-gnu',
#         'CONDA_BACKUP_MESON_ARGS': '-Dbuildtype=release',
#         'CONDA_BACKUP_ADDR2LINE': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-addr2line',
#         'CONDA_BACKUP_AR': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ar',
#         'CONDA_BACKUP_AS': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-as',
#         'CONDA_BACKUP_CXXFILT': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-c++filt',
#         'CONDA_BACKUP_ELFEDIT': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-elfedit',
#         'CONDA_BACKUP_GPROF': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gprof',
#         'CONDA_BACKUP_LD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ld',
#         'CONDA_BACKUP_LD_GOLD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ld.gold',
#         'CONDA_BACKUP_NM': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-nm',
#         'CONDA_BACKUP_OBJCOPY': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-objcopy',
#         'CONDA_BACKUP_OBJDUMP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-objdump',
#         'CONDA_BACKUP_RANLIB': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ranlib',
#         'CONDA_BACKUP_READELF': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-readelf',
#         'CONDA_BACKUP_SIZE': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-size',
#         'CONDA_BACKUP_STRINGS': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-strings',
#         'CONDA_BACKUP_STRIP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-strip',
#         'CXX': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-c++',
#         'GXX': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-g++',
#         'CXXFLAGS': '-fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'DEBUG_CXXFLAGS': '-fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CXX_FOR_BUILD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-c++',
#         'GFORTRAN': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'F95': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-f95',
#         'FC_FOR_BUILD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'FFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'FORTRANFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'DEBUG_FFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fcheck=all -fbacktrace -fimplicit-none -fvar-tracking-assignments -ffunction-sections -pipe -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'DEBUG_FORTRANFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fcheck=all -fbacktrace -fimplicit-none -fvar-tracking-assignments -ffunction-sections -pipe -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -I/apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'FC': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'F77': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'F90': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gfortran',
#         'BUILD': 'x86_64-conda-linux-gnu',
#         'CONDA_TOOLCHAIN_HOST': 'x86_64-conda-linux-gnu',
#         'CONDA_TOOLCHAIN_BUILD': 'x86_64-conda-linux-gnu',
#         'CC': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-cc',
#         'CPP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-cpp',
#         'GCC': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc',
#         'GCC_AR': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ar',
#         'GCC_NM': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-nm',
#         'GCC_RANLIB': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ranlib',
#         'CPPFLAGS': '-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /apps/miniforge3/envs/Geant4/include',
#         'CFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'LDFLAGS': '-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/apps/miniforge3/envs/Geant4/lib -Wl,-rpath-link,/apps/miniforge3/envs/Geant4/lib -L/apps/miniforge3/envs/Geant4/lib',
#         'DEBUG_CPPFLAGS': '-D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include -D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include -D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include -D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include -D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem /apps/miniforge3/envs/Geant4/include',
#         'DEBUG_CFLAGS': '-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem /apps/miniforge3/envs/Geant4/include',
#         'CMAKE_PREFIX_PATH': '/apps/miniforge3/envs/Geant4:/apps/miniforge3/envs/Geant4/x86_64-conda-linux-gnu/sysroot/usr',
#         '_CONDA_PYTHON_SYSCONFIGDATA_NAME': '_sysconfigdata_x86_64_conda_cos6_linux_gnu',
#         'CONDA_BUILD_SYSROOT': '/apps/miniforge3/envs/Geant4/x86_64-conda-linux-gnu/sysroot',
#         'CC_FOR_BUILD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-cc',
#         'build_alias': 'x86_64-conda-linux-gnu',
#         'host_alias': 'x86_64-conda-linux-gnu',
#         'MESON_ARGS': '-Dbuildtype=release',
#         'ADDR2LINE': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-addr2line',
#         'AR': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ar',
#         'AS': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-as',
#         'CXXFILT': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-c++filt',
#         'ELFEDIT': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-elfedit',
#         'GPROF': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gprof',
#         'LD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ld',
#         'NM': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-nm',
#         'OBJCOPY': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-objcopy',
#         'OBJDUMP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-objdump',
#         'RANLIB': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ranlib',
#         'READELF': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-readelf',
#         'SIZE': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-size',
#         'STRINGS': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-strings',
#         'STRIP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-strip',
#         'LD_GOLD': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ld.gold',
#         'G4NEUTRONHPDATA': '/newhome/fvanderv/pythondev/G4NDL4.5',
#         'G4LEDATA': '/newhome/fvanderv/pythondev/G4EMLOW7.3',
#         'G4LEVELGAMMADATA': '/newhome/fvanderv/pythondev/PhotonEvaporation5.2',
#         'G4RADIOACTIVEDATA': '/newhome/fvanderv/pythondev/RadioactiveDecay5.2',
#         'G4NEUTRONXSDATA': '/newhome/fvanderv/pythondev/G4NEUTRONXS1.4',
#         'G4PIIDATA': '/newhome/fvanderv/pythondev/G4PII1.3',
#         'G4REALSURFACEDATA': '/newhome/fvanderv/pythondev/RealSurface2.1.1',
#         'G4SAIDXSDATA': '/newhome/fvanderv/pythondev/G4SAIDDATA1.1',
#         'G4ABLADATA': '/newhome/fvanderv/pythondev/G4ABLA3.1',
#         'G4ENSDFSTATEDATA': '/newhome/fvanderv/pythondev/G4ENSDFSTATE2.2',
#         'LD_LIBRARY_PATH': '/newhome/fvanderv/pythondev/bdsim/bin/../lib',
#         'ROOT_INCLUDE_PATH': '/newhome/fvanderv/pythondev/bdsim/bin/../include/bdsim/:/newhome/fvanderv/pythondev/bdsim/bin/../include/bdsim/analysis/:/newhome/fvanderv/pythondev/bdsim/bin/../include/bdsim/parser/',
#         'CONDA_BACKUP_DWP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-dwp',
#         'CONDA_PREFIX_2': '/apps/miniforge3/envs/xcoll',
#         'DWP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-dwp',
#         'CMAKE_ARGS': '-DCMAKE_AR=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ar -DCMAKE_CXX_COMPILER_AR=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ar -DCMAKE_C_COMPILER_AR=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ar -DCMAKE_RANLIB=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ranlib -DCMAKE_CXX_COMPILER_RANLIB=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ranlib -DCMAKE_C_COMPILER_RANLIB=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ranlib -DCMAKE_LINKER=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ld -DCMAKE_STRIP=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-strip -DCMAKE_BUILD_TYPE=Release',
#         'PYTHIA8': '/apps/miniforge3/envs/Geant4',
#         'PYTHIA8_DIR': '/apps/miniforge3/envs/Geant4',
#         'PYTHIA8DATA': '/apps/miniforge3/envs/Geant4/share/Pythia8/xmldoc',
#         'ROOTSYS': '/apps/miniforge3/envs/Geant4',
#         'GSETTINGS_SCHEMA_DIR_CONDA_BACKUP': '',
#         'GSETTINGS_SCHEMA_DIR': '/apps/miniforge3/envs/Geant4/share/glib-2.0/schemas',
#         'XML_CATALOG_FILES': 'file:///apps/miniforge3/envs/Geant4/etc/xml/catalog file:///etc/xml/catalog',


# import os
# _old_os_environ = dict(os.environ)
# geant4_dir = '/newhome/fvanderv/pythondev/'
# geant4_env = {
#         'G4NEUTRONHPDATA': geant4_dir + 'G4NDL4.5',
#         'G4LEDATA': geant4_dir + 'G4EMLOW7.3',
#         'G4LEVELGAMMADATA': geant4_dir + 'PhotonEvaporation5.2',
#         'G4RADIOACTIVEDATA': geant4_dir + 'RadioactiveDecay5.2',
#         'G4NEUTRONXSDATA': geant4_dir + 'G4NEUTRONXS1.4',
#         'G4PIIDATA': geant4_dir + 'G4PII1.3',
#         'G4REALSURFACEDATA': geant4_dir + 'RealSurface2.1.1',
#         'G4SAIDXSDATA': geant4_dir + 'G4SAIDDATA1.1',
#         'G4ABLADATA': geant4_dir + 'G4ABLA3.1',
#         'G4ENSDFSTATEDATA': geant4_dir + 'G4ENSDFSTATE2.2',
#         'LD_LIBRARY_PATH': geant4_dir + 'bdsim/lib',
#         'ROOT_INCLUDE_PATH': geant4_dir + 'bdsim//include/bdsim/:' + geant4_dir + 'bdsim/include/bdsim/analysis/:' + geant4_dir + 'bdsim/include/bdsim/parser/',
#         # 'CONDA_BACKUP_DWP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-dwp',
#         # 'CONDA_PREFIX_2': '/apps/miniforge3/envs/xcoll',
#         # 'DWP': '/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-dwp',
#         # 'CMAKE_ARGS': '-DCMAKE_AR=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ar -DCMAKE_CXX_COMPILER_AR=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ar -DCMAKE_C_COMPILER_AR=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ar -DCMAKE_RANLIB=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ranlib -DCMAKE_CXX_COMPILER_RANLIB=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ranlib -DCMAKE_C_COMPILER_RANLIB=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-gcc-ranlib -DCMAKE_LINKER=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-ld -DCMAKE_STRIP=/apps/miniforge3/envs/Geant4/bin/x86_64-conda-linux-gnu-strip -DCMAKE_BUILD_TYPE=Release',
#         'PYTHIA8': '/apps/miniforge3/envs/Geant4',
#         'PYTHIA8_DIR': '/apps/miniforge3/envs/Geant4',
#         'PYTHIA8DATA': '/apps/miniforge3/envs/Geant4/share/Pythia8/xmldoc',
#         'ROOTSYS': '/apps/miniforge3/envs/Geant4',
# }
# for kk, val in geant4_env.items():
#     os.environ[kk] = val
# # This is not enough: need to pip install collimasim to get pybind etc
# # sys.path.append((xc._pkg_root / 'scattering_routines' / 'geant4').as_posix())


# os.environ.clear()
# os.environ.update(_old_os_environ)






# import numpy as np
# import pytest
# import time

# import xobjects as xo
# import xpart as xp
# import xcoll as xc

# #import pydevd
# #pydevd.settrace(suspend=True)
# #pydevd.settrace('127.0.0.1', port=5678, stdoutToServer=True, stderrToServer=True, suspend=True)

# import collimasim as cs

# path = xc._pkg_root.parent / 'tests' / 'data'

# jaw = 0.001
# #jaw = [0.0013, -0.002789]
# #jaw = [-1.2e-6, -3.2e-3]
# angle = 0

# tilt = 0    # For now, need to implement test for tilted jaws
# _ACCURACY = 1.5e-12  # Anything in this region around the jaw might or might not hit; we can't be sure
# num_part = 50
# _capacity = num_part*2
# jaw_band = 2.e-9

# # If a previous test failed, stop the server manually
# if xc.Geant4Engine.is_running():
#     xc.Geant4Engine.stop(clean=True)

# # Define collimator and start the Geant4 server
# coll = xc.Geant4Collimator(length=0.6, jaw=0.001, angle=0, tilt=0, material='Ti', geant4_id=f'g4coll_0')
# xc.Geant4Engine.start(elements=coll, seed=1993, particle_ref='proton', p0c=7.e12,
#                       bdsim_config_file=str(path / 'geant4_protons.gmad'))

# # Particle distribution (x and y are in the frame of the collimator)
# num_part_step = num_part//5
# context = xo.ContextCpu()
# particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=7e12, _context=context)
# x = np.random.uniform(-0.02, 0.02, num_part_step)
# x = np.concatenate([x, np.random.uniform(coll.jaw_L - jaw_band, coll.jaw_L -_ACCURACY, num_part_step)])
# x = np.concatenate([x, np.random.uniform(coll.jaw_L +_ACCURACY, coll.jaw_L + jaw_band, num_part_step)])
# x = np.concatenate([x, np.random.uniform(coll.jaw_R - jaw_band, coll.jaw_R -_ACCURACY, num_part_step)])
# x = np.concatenate([x, np.random.uniform(coll.jaw_R +_ACCURACY, coll.jaw_R + jaw_band, num_part_step)])
# y = np.random.uniform(-0.02, 0.02, 5*num_part_step)
# x_new = np.cos(np.deg2rad(angle))*x - np.sin(np.deg2rad(angle))*y
# y_new = np.sin(np.deg2rad(angle))*x + np.cos(np.deg2rad(angle))*y
# part_init = xp.build_particles(x=x_new, y=y_new, particle_ref=xc.Geant4Engine().particle_ref,
#                                _capacity=_capacity)

# mask = np.concatenate([(x >= min(coll.jaw_LU, coll.jaw_LD)) | (x <= max(coll.jaw_RU, coll.jaw_RD)),
#                       np.full(5*num_part_step, False)])
# hit_ids = part_init.particle_id[mask & (part_init.state > 0)]
# not_hit_ids = part_init.particle_id[~mask & (part_init.state > 0)]

# # TODO: jaw tilts, and particle angles

# # Track
# part = part_init.copy()
# coll.track(part)
# mask_hit = np.isin(part.parent_particle_id, hit_ids)
# mask_not_hit = np.isin(part.parent_particle_id, not_hit_ids)
# # Particles that are supposed to not have hit the collimator, but have a kick or are dead, are considered faulty
# #assert not np.any(abs(part.px[mask_not_hit]) > _ACCURACY)
# #assert not np.any(abs(part.py[mask_not_hit]) > _ACCURACY)
# #assert not np.any(part.state[mask_not_hit] < 1)

# # Particles that are supposed to have hit the collimator, but are alive and have no kick, are considered faulty
# ##faulty =  mask_hit & (abs(part.px) < _ACCURACY) & (abs(part.py) < _ACCURACY)
# #f#aulty &= (part.state > 0)
# #assert len(part.x[faulty]) <= 1  # We allow for a small margin of error

# xc.Geant4Engine.stop()
# print('hej')

# xc.Geant4Engine.start(elements=coll, seed=1993, particle_ref='proton', p0c=7.e12,
#                       bdsim_config_file=str(path / 'geant4_protons.gmad'))
# part2 = part_init.copy()
# coll.track(part2)