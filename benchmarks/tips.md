# Example benchmarks
All example benchmarks must be ran from within the xcoll/benchmarks folder

Run all:
```bash
pytest
```

Run quick benchmarks:
```bash
pytest -m quick
```

More expensive cases; includes a full LHC Run 3 crystal loss map:
```bash
pytest -m full
```

An individual case (could be a benchmark group as well) and specify an output file:
```bash
python -m pytest -k 'collimator_cfc_grazing' --benchmark-json=grazing.json
```

# Compare different branches
```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BENCHMARK_TESTBED="frederik-local-macbook-m3-cpu-serial"
```

In branch 1:
```bash
python -m pytest
```

In branch 2 (with the IDENTICAL benchmark files). Need to prebuild the kernels potentially
again (note that this cannot be done in the root folder, as it will wrongly read the local xcoll folder):
```bash
cd tests; xsuite-prebuild r --kind='serial'; cd ..
python -m pytest
```

Compare the results with `pytest-benchmark`:
```bash
pytest-benchmark compare results/branch1.json results/branch2.json
```

Simple offline visual dashboard:
```bash
python report.py results/branch1.json results/branch2.json -o dashboard.html
open dashboard.html
```

