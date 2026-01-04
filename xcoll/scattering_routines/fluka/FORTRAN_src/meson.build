project('pyflukaf', ['c', 'fortran'])

py = import('python').find_installation()

np_inc = run_command(py, '-c', 'import numpy; print(numpy.get_include())', check: true).stdout().strip()
f2py_inc = run_command(py, '-c', 'import numpy.f2py as f2py; print(f2py.get_include())', check: true).stdout().strip()

inc = include_directories('.', np_inc, f2py_inc)

# ---- your Fortran sources (in the correct order you already have) ----
deps_fortran = static_library(
  'pyfluka_deps',
  [
    'core_tools.f90',
    'constants.f90',
    'strings.f90',
    'mod_alloc.f90',
    'common_modules.f90',
    'string_tools.f90',
    'mod_units.f90',
    'pdgid.f90',
    'mod_fluka.f90',
    # add any other .f90 you need here
  ],
  # where to FIND modules produced by those sources (adjust if needed)
  fortran_args: ['-I.'],
)

# ---- the FLUKA IO static archive you want embedded ----
flukaio_a = files('libFlukaIO.a')

whole_archive_args = []
if host_machine.system() == 'darwin'
  whole_archive_args = ['-Wl,-force_load,' + flukaio_a[0].full_path()]
elif host_machine.system() == 'linux'
  whole_archive_args = ['-Wl,--whole-archive', flukaio_a[0].full_path(), '-Wl,--no-whole-archive']
endif

# ---- build the Python extension ----
py.extension_module(
  'pyflukaf',
  sources: [
    'pyfluka.f90',
    'pyflukaf-f2pywrappers.f',
    'pyflukafmodule.c',
    join_paths(f2py_inc, 'fortranobject.c'),
  ],
  include_directories: inc,
  link_with: deps_fortran,
  link_args: whole_archive_args,
  install: true,
)
