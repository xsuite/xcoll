"""Material initialisation and conversion, separate from particle tracking."""

import pytest
import xcoll as xc


def _material():
    return xc.Material(Z=74, A=183.84, density=19.3, name="W-benchmark")


@pytest.mark.quick
@pytest.mark.benchmark(group="Materials | construction")
def bench_material_elemental_construction(benchmark):
    material = benchmark.pedantic(_material, rounds=10, iterations=100)
    assert material.Z == 74
    benchmark.extra_info.update(protocol=2, physics_case="material-elemental", material="W")


@pytest.mark.quick
@pytest.mark.benchmark(group="Materials | serialisation")
def bench_material_to_dict(benchmark):
    material = xc.materials.CarbonFibreComposite
    data = benchmark.pedantic(material.to_dict, rounds=10, iterations=100)
    assert isinstance(data, dict) and data
    benchmark.extra_info.update(protocol=2, physics_case="material-to-dict", material="CFC")


@pytest.mark.full
@pytest.mark.benchmark(group="Materials | serialisation")
def bench_material_from_dict(benchmark):
    data = _material().to_dict()
    material = benchmark.pedantic(xc.Material.from_dict, args=(data,), rounds=10, iterations=100)
    assert material.Z == 74
    benchmark.extra_info.update(protocol=2, physics_case="material-from-dict", material="W")
