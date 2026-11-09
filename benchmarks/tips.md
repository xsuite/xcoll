# Setup
```bash
mkdir -p .benchmarks
```

# Example benchmarks
All example benchmarks must be ran from within the Xcoll root folder

Run all:
```bash
python -m pytest -c benchmarks/pytest.ini benchmarks --benchmark-json=.benchmarks/result.json
```

Run quick benchmarks:
```bash
python -m pytest -c benchmarks/pytest.ini benchmarks -m quick --benchmark-json=.benchmarks/quick.json
```

More expensive cases; includes a full LHC Run 3 crystal loss map:
```bash
python -m pytest -c benchmarks/pytest.ini benchmarks -m full --benchmark-json=.benchmarks/full.json
```

An individual case (or a group):
```bash
python -m pytest -c benchmarks/pytest.ini benchmarks -k 'collimator_cfc_grazing' --benchmark-json=.benchmarks/grazing.json
```

Save measured round data as well, for later statistical inspection:
```bash
python -m pytest -c benchmarks/pytest.ini benchmarks -m quick --benchmark-save-data --benchmark-json=.benchmarks/quick_full_data.json
```

# Compare different branches
```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BENCHMARK_TESTBED="frederik-local-macbook-m3-cpu-serial"
```

In branch 1:
```bash
python -m pytest -c benchmarks/pytest.ini benchmarks --benchmark-save-data --benchmark-json=.benchmarks/branch1.json
```

In branch 2 (with the IDENTICAL benchmark files). Need to prebuild the kernels potentially
again (note that this cannot be done in the root folder, as it will wrongly read the local xcoll folder):
```bash
cd tests; xsuite-prebuild r --kind='serial'; cd ..
python -m pytest -c benchmarks/pytest.ini benchmarks --benchmark-save-data --benchmark-json=.benchmarks/branch2.json
```

Compare the results with `pytest-benchmark`:
```bash
pytest-benchmark compare .benchmarks/branch1.json .benchmarks/branch2.json
```

Simple offline visual dashboard:
```bash
python benchmarks/report.py .benchmarks/branch1.json .benchmarks/branch2.json -o .benchmarks/dashboard.html
open .benchmarks/dashboard.html
```

