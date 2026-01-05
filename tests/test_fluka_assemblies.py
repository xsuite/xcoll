# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import re
import pytest
from pathlib import Path

import xcoll as xc

# Do NOT mark these as FLUKA tests, as they have to be ran in serial

@pytest.mark.fluka
@pytest.mark.serial
def test_registry_initialisation(request):
    worker_id = getattr(request.config, "workerinput", None)
    if worker_id is not None:
        pytest.skip("Running test in parallel - FLUKA assembly initialisation test cannot run.")
    assert xc.FlukaPrototype._registry == []
    assert xc.FlukaAssembly._registry == []
    xc.fluka.environment # Force initialization
    assert xc.FlukaPrototype._registry is xc.FlukaAssembly._registry
    assert len(xc.FlukaPrototype._registry) >= 7 # Has at least Maria crystals and Donadon test assembly
    assert len([pro for pro in xc.FlukaPrototype._registry if isinstance(pro, xc.FlukaAssembly)]) >= 1
    assert len([pro for pro in xc.FlukaPrototype._registry if not isinstance(pro, xc.FlukaAssembly)]) >= 6
    assert xc.FlukaPrototype._assigned_registry == {}
    assert xc.FlukaAssembly._assigned_registry == {}
    assert xc.FlukaPrototype._assigned_registry is not xc.FlukaAssembly._assigned_registry


@pytest.mark.fluka
@pytest.mark.serial
def test_new_null():
    xc.fluka.environment # Force initialization
    new_pro = xc.FlukaPrototype()
    assert new_pro._is_null
    assert new_pro not in xc.FlukaPrototype._registry
    assert new_pro not in xc.FlukaPrototype._assigned_registry
    assert new_pro not in xc.FlukaAssembly._assigned_registry

    new_pro = xc.FlukaAssembly()
    assert new_pro._is_null
    assert new_pro not in xc.FlukaPrototype._registry
    assert new_pro not in xc.FlukaPrototype._assigned_registry
    assert new_pro not in xc.FlukaAssembly._assigned_registry


@pytest.mark.fluka
@pytest.mark.serial
def test_new_prototype():
    xc.fluka.environment # Force initialization
    prototypes_before = xc.FlukaPrototype._registry.copy()
    with pytest.raises(ValueError, match="Both 'fedb_series' and 'fedb_tag' must be provided."):
        new_pro = xc.FlukaPrototype(fedb_series='test')
    with pytest.raises(ValueError, match="Both 'fedb_series' and 'fedb_tag' must be provided."):
        new_pro = xc.FlukaPrototype(fedb_tag='proto')
    new_pro = xc.FlukaPrototype(fedb_series='test', fedb_tag='proto')
    assert not new_pro._is_null
    assert new_pro in xc.FlukaPrototype._registry
    assert new_pro not in xc.FlukaPrototype._assigned_registry
    assert new_pro not in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaPrototype._registry == prototypes_before + [new_pro]

    prototypes_before = xc.FlukaPrototype._registry.copy()
    new_pro_copy = xc.FlukaPrototype(fedb_series='test', fedb_tag='proto')
    assert new_pro_copy is new_pro
    assert xc.FlukaPrototype._registry == prototypes_before
    assert new_pro.name == 'proto'
    assert str(new_pro) == "FlukaPrototype 'proto' (unassigned): tag proto in test series"
    assert new_pro.body_file is not None
    assert new_pro.material_file is not None
    assert new_pro.region_file is not None
    assert not hasattr(new_pro, 'assembly_file')
    assert new_pro.files == [new_pro.body_file, new_pro.material_file, new_pro.region_file]
    assert new_pro.body_file.as_posix() == (xc._pkg_root / 'lib' / 'fedb' / 'bodies' / 'test_proto.bodies').as_posix()
    assert new_pro.material_file.as_posix() == (xc._pkg_root / 'lib' / 'fedb' / 'materials' / 'test_proto.assignmat').as_posix()
    assert new_pro.region_file.as_posix() == (xc._pkg_root / 'lib' / 'fedb' / 'regions' / 'test_proto.regions').as_posix()
    assert not new_pro.exists()
    assert not new_pro.assigned
    assert new_pro.elements == []
    assert new_pro._id is None

    assert new_pro.to_dict() == {'__class__': 'FlukaPrototype',
                                 'name': 'proto',
                                 'fedb_series': 'test',
                                 'fedb_tag': 'proto',
                                 'side': None,
                                 'angle': 0,
                                 'length': None,
                                 'width': None,
                                 'height': None,
                                 'material': None,
                                 'is_crystal': False,
                                 'bending_radius': None,
                                 'info': None,
                                 'extra_commands': None,
                                 'is_broken': False}
    with pytest.raises(ValueError, match="Prototype 'proto' does not exist in the FEDB!"):
        coll = xc.FlukaCollimator(assembly=new_pro)
    new_pro.delete()


@pytest.mark.fluka
@pytest.mark.serial
def test_prototype_with_files():
    new_pro1 = xc.FlukaPrototype(fedb_series='test', fedb_tag='proto1')
    new_pro2 = xc.FlukaPrototype(fedb_series='test', fedb_tag='proto2')
    with open("pro_body.txt", "w") as f:
        f.write("RPP PRO_BODY   0.0 9.0 -4.42 4.42 -60. 60.\n")
    with open("pro_region.txt", "w") as f:
        f.write("PRO_Jawl    25 | +PRO_BODY\n")
    with open("pro_material.txt", "w") as f:
        f.write("ASSIGNMA    AC150GPH  PRO_Jawl\n")
    new_pro1.body_file = "pro_body.txt"
    new_pro1.region_file = "pro_region.txt"
    new_pro1.material_file = "pro_material.txt"
    new_pro2.body_file = "pro_body.txt"
    new_pro2.region_file = "pro_region.txt"
    new_pro2.material_file = "pro_material.txt"
    Path("pro_body.txt").unlink()
    Path("pro_region.txt").unlink()
    Path("pro_material.txt").unlink()
    assert new_pro1.body_file.exists()
    assert new_pro1.region_file.exists()
    assert new_pro1.material_file.exists()
    assert new_pro1.exists()
    assert new_pro2.body_file.exists()
    assert new_pro2.region_file.exists()
    assert new_pro2.material_file.exists()
    assert new_pro2.exists()
    assert not new_pro1.assigned
    assert not new_pro2.assigned
    assert new_pro1.elements == []
    assert new_pro2.elements == []

    coll1 = xc.FlukaCollimator(assembly=new_pro1)
    assert coll1.assembly is new_pro1
    assert new_pro1.assigned
    assert not new_pro2.assigned
    assert not new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == [coll1]
    assert new_pro2.elements == []
    assert new_pro1.active_elements == []
    assert new_pro2.active_elements == []
    assert new_pro1._id == 0
    assert new_pro2._id is None
    assert new_pro1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro2.fluka_position is None
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto1'] is new_pro1
    assert 'proto2' not in xc.FlukaPrototype._assigned_registry
    assert str(new_pro1) == "FlukaPrototype 'proto1' (1 inactive element): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (unassigned): tag proto2 in test series"

    coll2 = xc.FlukaCollimator(assembly=new_pro1)
    assert coll1.assembly is new_pro1
    assert coll2.assembly is new_pro1
    assert new_pro1.assigned
    assert not new_pro2.assigned
    assert not new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == [coll1, coll2]
    assert new_pro2.elements == []
    assert new_pro1.active_elements == []
    assert new_pro2.active_elements == []
    assert new_pro1._id == 0
    assert new_pro2._id is None
    assert new_pro1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro2.fluka_position is None
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto1'] is new_pro1
    assert 'proto2' not in xc.FlukaPrototype._assigned_registry
    assert str(new_pro1) == "FlukaPrototype 'proto1' (2 inactive elements): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (unassigned): tag proto2 in test series"

    coll2.jaw = 2e-3
    assert coll1.assembly is new_pro1
    assert coll2.assembly is new_pro1
    assert new_pro1.assigned
    assert not new_pro2.assigned
    assert new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == [coll1, coll2]
    assert new_pro2.elements == []
    assert new_pro1.active_elements == [coll2]
    assert new_pro2.active_elements == []
    assert new_pro1._id == 0
    assert new_pro2._id is None
    assert new_pro1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro2.fluka_position is None
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto1'] is new_pro1
    assert 'proto2' not in xc.FlukaPrototype._assigned_registry
    assert str(new_pro1) == "FlukaPrototype 'proto1' (1 active and 1 inactive elements): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (unassigned): tag proto2 in test series"

    coll2.assembly = new_pro1
    assert coll1.assembly is new_pro1
    assert coll2.assembly is new_pro1
    assert new_pro1.assigned
    assert not new_pro2.assigned
    assert new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == [coll1, coll2]
    assert new_pro2.elements == []
    assert new_pro1.active_elements == [coll2]
    assert new_pro2.active_elements == []
    assert new_pro1._id == 0
    assert new_pro2._id is None
    assert new_pro1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro2.fluka_position is None
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto1'] is new_pro1
    assert 'proto2' not in xc.FlukaPrototype._assigned_registry
    assert str(new_pro1) == "FlukaPrototype 'proto1' (1 active and 1 inactive elements): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (unassigned): tag proto2 in test series"

    coll3 = xc.FlukaCollimator(assembly=new_pro2)
    coll3.jaw = 1.2e-3
    assert coll1.assembly is new_pro1
    assert coll2.assembly is new_pro1
    assert coll3.assembly is new_pro2
    assert new_pro1.assigned
    assert new_pro2.assigned
    assert new_pro1.active
    assert new_pro2.active
    assert new_pro1.elements == [coll1, coll2]
    assert new_pro2.elements == [coll3]
    assert new_pro1.active_elements == [coll2]
    assert new_pro2.active_elements == [coll3]
    assert new_pro1._id == 0
    assert new_pro2._id == 1
    assert new_pro1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro2.fluka_position == (0.0, 0.0, 0.0, -500.0, -3000.0, 1000.0)
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto1'] is new_pro1
    assert 'proto2' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto2'] is new_pro2
    assert str(new_pro1) == "FlukaPrototype 'proto1' (1 active and 1 inactive elements): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (1 active element): tag proto2 in test series"

    coll1.assembly = new_pro2
    assert coll1.assembly is new_pro2
    assert coll2.assembly is new_pro1
    assert coll3.assembly is new_pro2
    assert new_pro1.assigned
    assert new_pro2.assigned
    assert new_pro1.active
    assert new_pro2.active
    assert new_pro1.elements == [coll2]
    assert new_pro2.elements == [coll3, coll1]
    assert new_pro1.active_elements == [coll2]
    assert new_pro2.active_elements == [coll3]
    assert new_pro1._id == 0
    assert new_pro2._id == 1
    assert new_pro1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro2.fluka_position == (0.0, 0.0, 0.0, -500.0, -3000.0, 1000.0)
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto1'] is new_pro1
    assert 'proto2' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto2'] is new_pro2
    assert str(new_pro1) == "FlukaPrototype 'proto1' (1 active element): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (1 active and 1 inactive elements): tag proto2 in test series"

    with pytest.raises(ValueError, match=re.escape("Element 'FlukaCollimator(length=0, active=np.int8(1), fluka_id=np.int16(0), length_front=0, length_back=0)' already assigned proto2 prototype!")):
        new_pro1.add_element(coll3)
    new_pro2.remove_element(coll3)
    assert coll1.assembly is new_pro2
    assert coll2.assembly is new_pro1
    assert coll3.assembly is None
    assert new_pro1.assigned
    assert new_pro2.assigned
    assert new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == [coll2]
    assert new_pro2.elements == [coll1]
    assert new_pro1.active_elements == [coll2]
    assert new_pro2.active_elements == []
    assert new_pro1._id == 0
    assert new_pro2._id == 1
    assert new_pro1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro2.fluka_position == (0.0, 0.0, 0.0, -500.0, -3000.0, 1000.0)
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto1'] is new_pro1
    assert 'proto2' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto2'] is new_pro2
    assert str(new_pro1) == "FlukaPrototype 'proto1' (1 active element): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (1 inactive element): tag proto2 in test series"

    new_pro1.add_element(coll3)
    assert coll1.assembly is new_pro2
    assert coll2.assembly is new_pro1
    assert coll3.assembly is new_pro1
    assert new_pro1.assigned
    assert new_pro2.assigned
    assert new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == [coll2, coll3]
    assert new_pro2.elements == [coll1]
    assert new_pro1.active_elements == [coll2, coll3]
    assert new_pro2.active_elements == []
    assert new_pro1._id == 0
    assert new_pro2._id == 1
    assert new_pro1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro2.fluka_position == (0.0, 0.0, 0.0, -500.0, -3000.0, 1000.0)
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto1'] is new_pro1
    assert 'proto2' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto2'] is new_pro2
    assert str(new_pro1) == "FlukaPrototype 'proto1' (2 active elements): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (1 inactive element): tag proto2 in test series"

    with pytest.raises(ValueError, match="Cannot delete prototype 'proto1' while it has elements assigned!"):
        new_pro1.delete()
    new_pro1.clear()
    assert coll1.assembly is new_pro2
    assert coll2.assembly is None
    assert coll3.assembly is None
    assert not new_pro1.assigned
    assert new_pro2.assigned
    assert not new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == []
    assert new_pro2.elements == [coll1]
    assert new_pro1.active_elements == []
    assert new_pro2.active_elements == []
    assert new_pro1._id is None
    assert new_pro2._id == 0
    assert new_pro1.fluka_position is None
    assert new_pro2.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_pro1.body_file.exists()
    assert new_pro1.region_file.exists()
    assert new_pro1.material_file.exists()
    assert new_pro1.exists()
    assert new_pro2.body_file.exists()
    assert new_pro2.region_file.exists()
    assert new_pro2.material_file.exists()
    assert new_pro2.exists()
    assert new_pro1 in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' not in xc.FlukaPrototype._assigned_registry
    assert 'proto2' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto2'] is new_pro2
    assert str(new_pro1) == "FlukaPrototype 'proto1' (unassigned): tag proto1 in test series"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (1 inactive element): tag proto2 in test series"

    new_pro1.delete()
    assert coll1.assembly is new_pro2
    assert coll2.assembly is None
    assert coll3.assembly is None
    assert not new_pro1.assigned
    assert new_pro2.assigned
    assert not new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == []
    assert new_pro2.elements == [coll1]
    assert new_pro1.active_elements == []
    assert new_pro2.active_elements == []
    assert new_pro1._id is None
    assert new_pro2._id == 0
    assert new_pro1.fluka_position is None
    assert new_pro2.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert not new_pro1.body_file.exists()
    assert not new_pro1.region_file.exists()
    assert not new_pro1.material_file.exists()
    assert not new_pro1.exists()
    assert new_pro2.body_file.exists()
    assert new_pro2.region_file.exists()
    assert new_pro2.material_file.exists()
    assert new_pro2.exists()
    assert new_pro1 not in xc.FlukaPrototype._registry
    assert new_pro2 in xc.FlukaPrototype._registry
    assert 'proto1' not in xc.FlukaPrototype._assigned_registry
    assert 'proto2' in xc.FlukaPrototype._assigned_registry
    assert xc.FlukaPrototype._assigned_registry['proto2'] is new_pro2
    assert str(new_pro1) == "FlukaPrototype 'proto1' (unassigned): tag proto1 in test series <defunct>"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (1 inactive element): tag proto2 in test series"

    new_pro2.clear()
    new_pro2.delete()
    assert coll1.assembly is None
    assert coll2.assembly is None
    assert coll3.assembly is None
    assert not new_pro1.assigned
    assert not new_pro2.assigned
    assert not new_pro1.active
    assert not new_pro2.active
    assert new_pro1.elements == []
    assert new_pro2.elements == []
    assert new_pro1.active_elements == []
    assert new_pro2.active_elements == []
    assert new_pro1._id is None
    assert new_pro2._id is None
    assert new_pro1.fluka_position is None
    assert new_pro2.fluka_position is None
    assert not new_pro1.body_file.exists()
    assert not new_pro1.region_file.exists()
    assert not new_pro1.material_file.exists()
    assert not new_pro1.exists()
    assert not new_pro2.body_file.exists()
    assert not new_pro2.region_file.exists()
    assert not new_pro2.material_file.exists()
    assert not new_pro2.exists()
    assert new_pro1 not in xc.FlukaPrototype._registry
    assert new_pro2 not in xc.FlukaPrototype._registry
    assert 'proto1' not in xc.FlukaPrototype._assigned_registry
    assert 'proto2' not in xc.FlukaPrototype._assigned_registry
    assert str(new_pro1) == "FlukaPrototype 'proto1' (unassigned): tag proto1 in test series <defunct>"
    assert str(new_pro2) == "FlukaPrototype 'proto2' (unassigned): tag proto2 in test series <defunct>"


@pytest.mark.fluka
@pytest.mark.serial
def test_new_assembly():
    xc.fluka.environment # Force initialization
    assemblies_before = xc.FlukaAssembly._registry.copy()
    with pytest.raises(ValueError, match="Both 'fedb_series' and 'fedb_tag' must be provided."):
        new_assm = xc.FlukaAssembly(fedb_series='test')
    with pytest.raises(ValueError, match="Both 'fedb_series' and 'fedb_tag' must be provided."):
        new_assm = xc.FlukaAssembly(fedb_tag='assm')
    new_assm = xc.FlukaAssembly(fedb_series='test', fedb_tag='assm')
    assert not new_assm._is_null
    assert new_assm in xc.FlukaPrototype._registry
    assert new_assm not in xc.FlukaAssembly._assigned_registry
    assert new_assm not in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaPrototype._registry == assemblies_before + [new_assm]
    assemblies_before = xc.FlukaAssembly._registry.copy()
    new_assm_copy = xc.FlukaAssembly(fedb_series='test', fedb_tag='assm')
    assert new_assm_copy is new_assm
    assert xc.FlukaPrototype._registry == assemblies_before
    assert new_assm.name == 'assm'
    assert str(new_assm) == "FlukaAssembly 'assm' (unassigned): tag assm in test series"
    assert new_assm.body_file is None
    assert new_assm.material_file is None
    assert new_assm.region_file is None
    assert hasattr(new_assm, 'assembly_file')
    assert new_assm.assembly_file is not None
    assert new_assm.files is None
    assert new_assm.assembly_file.as_posix() == (xc._pkg_root / 'lib' / 'fedb' / 'assemblies' / 'test_assm.lbp').as_posix()
    assert not new_assm.exists()
    assert not new_assm.assigned
    assert new_assm.elements == []
    assert new_assm._id is None
    assert new_assm.to_dict() == {'__class__': 'FlukaAssembly',
                                 'name': 'assm',
                                 'fedb_series': 'test',
                                 'fedb_tag': 'assm',
                                 'side': None,
                                 'angle': 0,
                                 'length': None,
                                 'width': None,
                                 'height': None,
                                 'material': None,
                                 'is_crystal': False,
                                 'bending_radius': None,
                                 'info': None,
                                 'extra_commands': None,
                                 'is_broken': False}
    with pytest.raises(ValueError, match="Assembly 'assm' does not exist in the FEDB!"):
        coll = xc.FlukaCollimator(assembly=new_assm)
    new_assm.delete(_ignore_files=True)


@pytest.mark.fluka
@pytest.mark.serial
def test_assembly_with_files():
    xc.fluka.environment # Force initialization
    # Create prototypes needed for the assembly
    new_pro_jaw = xc.FlukaPrototype(fedb_series='test', fedb_tag='PROTO_B')
    with open("pro_body.txt", "w") as f:
        f.write("RPP PROTO_B   0.0 9.0 -4.42 4.42 -60. 60.\n")
    with open("pro_region.txt", "w") as f:
        f.write("PROTO_B    5 +PROTO_B\n")
    with open("pro_material.txt", "w") as f:
        f.write("ASSIGNMA    AC150GPH  PROTO_B\n")
    new_pro_jaw.body_file = "pro_body.txt"
    new_pro_jaw.region_file = "pro_region.txt"
    new_pro_jaw.material_file = "pro_material.txt"
    Path("pro_body.txt").unlink()
    Path("pro_region.txt").unlink()
    Path("pro_material.txt").unlink()
    new_pro_tank = xc.FlukaPrototype(fedb_series='test', fedb_tag='PROTO_T')
    with open("pro_body.txt", "w") as f:
        f.write("RPP PROTO_T  -28 28 -28 28 -65 65\n")
        f.write("RPP PROTO_I  -28 28 -28 28 -65 65\n")
    with open("pro_region.txt", "w") as f:
        f.write("PROTO_T     5 +PROTO_T -PROTO_I\n")
        f.write("PROTO_I     5 +PROTO_I\n")
    with open("pro_material.txt", "w") as f:
        f.write("ASSIGNMA    VACUUM  PROTO_T\n")
        f.write("ASSIGNMA    VACUUM  PROTO_I\n")
    new_pro_tank.body_file = "pro_body.txt"
    new_pro_tank.region_file = "pro_region.txt"
    new_pro_tank.material_file = "pro_material.txt"
    Path("pro_body.txt").unlink()
    Path("pro_region.txt").unlink()
    Path("pro_material.txt").unlink()

    # Create assemblies
    new_assm1 = xc.FlukaAssembly(fedb_series='test', fedb_tag='ASSM1')
    with open("assembly.txt", "w") as f:
        f.write("PROTOTYPE       PROTO_T\n")
        f.write("FEDB_TAG        PROTO_T\n")
        f.write("FEDB_SERIES     test\n")
        f.write("ROT-DEFI         0.0       0.0       0.0     100.0   -3000.0    1000.0 proto\n")
        f.write("PROTOTYPE       PROTO_B\n")
        f.write("FEDB_TAG        PROTO_B\n")
        f.write("FEDB_SERIES     test\n")
        f.write("ROT-DEFI         0.0       0.0       0.0       0.0   -3000.0    1000.0 proto\n")
        f.write("ASSEMBLY        ASSM1\n")
        f.write("BODY        CONTAINO    CONTAINO    PROTO_T     test     PROTO_T    1\n")
        f.write("BODY        CONTAINI    CONTAINO    PROTO_I     test     PROTO_T    1\n")
        f.write("BODY        JAW_POS     JAW_POS     PROTO_B     test     PROTO_B    1\n")
        f.write("BODY        JAW_NEG     JAW_NEG     PROTO_B     test     PROTO_B    1\n")
        f.write("REGION      *           EXTERNAL    *           -CONTAINO\n")
        f.write("REGION      TANK        LATTICE     CONTAINO    +CONTAINO -CONTAINI\n")
        f.write("REGION      INNERVAC    VACUUM      *           +CONTAINI -JAW_POS -JAW_NEG\n")
        f.write("REGION      JAW_POS     LATTICE     JAW_POS     +JAW_POS\n")
        f.write("REGION      JAW_NEG     LATTICE     JAW_NEG     +JAW_NEG\n")
        f.write("ROT-DEFI             0.0         0.0         0.0         0.0         0.0         0.0 CONTAINO\n")
        f.write("ROT-DEFI             0.0         0.0         0.0         0.0         0.0         0.0 JAW_POS\n")
        f.write("ROT-DEFI           300.0         0.0       180.0         0.0         0.0         0.0 JAW_NEG\n")
    new_assm1.assembly_file = "assembly.txt"
    Path("assembly.txt").unlink()
    new_assm2 = xc.FlukaAssembly(fedb_series='test', fedb_tag='ASSM2')
    with open("assembly.txt", "w") as f:
        f.write("PROTOTYPE       PROTO_T\n")
        f.write("FEDB_TAG        PROTO_T\n")
        f.write("FEDB_SERIES     test\n")
        f.write("ROT-DEFI         0.0       0.0       0.0     100.0   -3000.0    1000.0 proto\n")
        f.write("PROTOTYPE       PROTO_B\n")
        f.write("FEDB_TAG        PROTO_B\n")
        f.write("FEDB_SERIES     test\n")
        f.write("ROT-DEFI         0.0       0.0       0.0       0.0   -3000.0    1000.0 proto\n")
        f.write("ASSEMBLY        ASSM2\n")
        f.write("BODY        CONTAINO    CONTAINO    PROTO_T     test     PROTO_T    1\n")
        f.write("BODY        CONTAINI    CONTAINO    PROTO_I     test     PROTO_T    1\n")
        f.write("BODY        JAW_POS     JAW_POS     PROTO_B     test     PROTO_B    1\n")
        f.write("BODY        JAW_NEG     JAW_NEG     PROTO_B     test     PROTO_B    1\n")
        f.write("REGION      *           EXTERNAL    *           -CONTAINO\n")
        f.write("REGION      TANK        LATTICE     CONTAINO    +CONTAINO -CONTAINI\n")
        f.write("REGION      INNERVAC    VACUUM      *           +CONTAINI -JAW_POS -JAW_NEG\n")
        f.write("REGION      JAW_POS     LATTICE     JAW_POS     +JAW_POS\n")
        f.write("REGION      JAW_NEG     LATTICE     JAW_NEG     +JAW_NEG\n")
        f.write("ROT-DEFI             0.0         0.0         0.0         0.0         0.0         0.0 CONTAINO\n")
        f.write("ROT-DEFI             0.0         0.0         0.0         0.0         0.0         0.0 JAW_POS\n")
        f.write("ROT-DEFI           300.0         0.0       180.0         0.0         0.0         0.0 JAW_NEG\n")
    new_assm2.assembly_file = "assembly.txt"
    Path("assembly.txt").unlink()

    assert new_assm1.assembly_file.exists()
    assert new_assm2.assembly_file.exists()
    assert new_assm1.exists()
    assert new_assm2.exists()
    assert not new_assm1.assigned
    assert not new_assm2.assigned
    assert new_assm1.elements == []
    assert new_assm2.elements == []
    assert new_assm1.prototypes == [new_pro_tank, new_pro_jaw]
    assert new_assm2.prototypes == [new_pro_tank, new_pro_jaw]
    files  = set([ff.as_posix() for ff in new_pro_jaw.files])
    files |= set([ff.as_posix() for ff in new_pro_tank.files])
    files |= set([new_assm1.assembly_file.as_posix()])
    assert set([ff.as_posix() for ff in new_assm1.files]) == files
    files  = set([ff.as_posix() for ff in new_pro_jaw.files])
    files |= set([ff.as_posix() for ff in new_pro_tank.files])
    files |= set([new_assm2.assembly_file.as_posix()])
    assert set([ff.as_posix() for ff in new_assm2.files]) == files
    assert new_assm1.check_file_valid()
    assert new_assm2.check_file_valid()

    coll1 = xc.FlukaCollimator(assembly=new_assm1)
    assert coll1.assembly is new_assm1
    assert new_assm1.assigned
    assert not new_assm2.assigned
    assert not new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == [coll1]
    assert new_assm2.elements == []
    assert new_assm1.active_elements == []
    assert new_assm2.active_elements == []
    assert new_assm1._id == 0
    assert new_assm2._id is None
    assert new_assm1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm2.fluka_position is None
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM1'] is new_assm1
    assert 'ASSM2' not in xc.FlukaAssembly._assigned_registry
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (1 inactive element): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (unassigned): tag ASSM2 in test series"

    coll2 = xc.FlukaCollimator(assembly=new_assm1)
    assert coll1.assembly is new_assm1
    assert coll2.assembly is new_assm1
    assert new_assm1.assigned
    assert not new_assm2.assigned
    assert not new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == [coll1, coll2]
    assert new_assm2.elements == []
    assert new_assm1.active_elements == []
    assert new_assm2.active_elements == []
    assert new_assm1._id == 0
    assert new_assm2._id is None
    assert new_assm1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm2.fluka_position is None
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM1'] is new_assm1
    assert 'ASSM2' not in xc.FlukaAssembly._assigned_registry
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (2 inactive elements): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (unassigned): tag ASSM2 in test series"

    coll2.jaw = 2e-3
    assert coll1.assembly is new_assm1
    assert coll2.assembly is new_assm1
    assert new_assm1.assigned
    assert not new_assm2.assigned
    assert new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == [coll1, coll2]
    assert new_assm2.elements == []
    assert new_assm1.active_elements == [coll2]
    assert new_assm2.active_elements == []
    assert new_assm1._id == 0
    assert new_assm2._id is None
    assert new_assm1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm2.fluka_position is None
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM1'] is new_assm1
    assert 'ASSM2' not in xc.FlukaAssembly._assigned_registry
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (1 active and 1 inactive elements): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (unassigned): tag ASSM2 in test series"

    coll2.assembly = new_assm1
    assert coll1.assembly is new_assm1
    assert coll2.assembly is new_assm1
    assert new_assm1.assigned
    assert not new_assm2.assigned
    assert new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == [coll1, coll2]
    assert new_assm2.elements == []
    assert new_assm1.active_elements == [coll2]
    assert new_assm2.active_elements == []
    assert new_assm1._id == 0
    assert new_assm2._id is None
    assert new_assm1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm2.fluka_position is None
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM1'] is new_assm1
    assert 'ASSM2' not in xc.FlukaAssembly._assigned_registry
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (1 active and 1 inactive elements): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (unassigned): tag ASSM2 in test series"

    coll3 = xc.FlukaCollimator(assembly=new_assm2)
    coll3.jaw = 1.2e-3
    assert coll1.assembly is new_assm1
    assert coll2.assembly is new_assm1
    assert coll3.assembly is new_assm2
    assert new_assm1.assigned
    assert new_assm2.assigned
    assert new_assm1.active
    assert new_assm2.active
    assert new_assm1.elements == [coll1, coll2]
    assert new_assm2.elements == [coll3]
    assert new_assm1.active_elements == [coll2]
    assert new_assm2.active_elements == [coll3]
    assert new_assm1._id == 0
    assert new_assm2._id == 1
    assert new_assm1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm2.fluka_position == (0.0, 0.0, 0.0, -500.0, -3000.0, 1000.0)
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM1'] is new_assm1
    assert 'ASSM2' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM2'] is new_assm2
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (1 active and 1 inactive elements): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (1 active element): tag ASSM2 in test series"

    coll1.assembly = new_assm2
    assert coll1.assembly is new_assm2
    assert coll2.assembly is new_assm1
    assert coll3.assembly is new_assm2
    assert new_assm1.assigned
    assert new_assm2.assigned
    assert new_assm1.active
    assert new_assm2.active
    assert new_assm1.elements == [coll2]
    assert new_assm2.elements == [coll3, coll1]
    assert new_assm1.active_elements == [coll2]
    assert new_assm2.active_elements == [coll3]
    assert new_assm1._id == 0
    assert new_assm2._id == 1
    assert new_assm1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm2.fluka_position == (0.0, 0.0, 0.0, -500.0, -3000.0, 1000.0)
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM1'] is new_assm1
    assert 'ASSM2' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM2'] is new_assm2
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (1 active element): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (1 active and 1 inactive elements): tag ASSM2 in test series"

    with pytest.raises(ValueError, match=re.escape("Element 'FlukaCollimator(length=0, active=np.int8(1), fluka_id=np.int16(0), length_front=0, length_back=0)' already assigned ASSM2 assembly!")):
        new_assm1.add_element(coll3)
    new_assm2.remove_element(coll3)
    assert coll1.assembly is new_assm2
    assert coll2.assembly is new_assm1
    assert coll3.assembly is None
    assert new_assm1.assigned
    assert new_assm2.assigned
    assert new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == [coll2]
    assert new_assm2.elements == [coll1]
    assert new_assm1.active_elements == [coll2]
    assert new_assm2.active_elements == []
    assert new_assm1._id == 0
    assert new_assm2._id == 1
    assert new_assm1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm2.fluka_position == (0.0, 0.0, 0.0, -500.0, -3000.0, 1000.0)
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM1'] is new_assm1
    assert 'ASSM2' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM2'] is new_assm2
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (1 active element): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (1 inactive element): tag ASSM2 in test series"

    new_assm1.add_element(coll3)
    assert coll1.assembly is new_assm2
    assert coll2.assembly is new_assm1
    assert coll3.assembly is new_assm1
    assert new_assm1.assigned
    assert new_assm2.assigned
    assert new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == [coll2, coll3]
    assert new_assm2.elements == [coll1]
    assert new_assm1.active_elements == [coll2, coll3]
    assert new_assm2.active_elements == []
    assert new_assm1._id == 0
    assert new_assm2._id == 1
    assert new_assm1.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm2.fluka_position == (0.0, 0.0, 0.0, -500.0, -3000.0, 1000.0)
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM1'] is new_assm1
    assert 'ASSM2' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM2'] is new_assm2
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (2 active elements): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (1 inactive element): tag ASSM2 in test series"

    with pytest.raises(ValueError, match="Cannot delete assembly 'ASSM1' while it has elements assigned!"):
        new_assm1.delete()
    new_assm1.clear()
    assert coll1.assembly is new_assm2
    assert coll2.assembly is None
    assert coll3.assembly is None
    assert not new_assm1.assigned
    assert new_assm2.assigned
    assert not new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == []
    assert new_assm2.elements == [coll1]
    assert new_assm1.active_elements == []
    assert new_assm2.active_elements == []
    assert new_assm1._id is None
    assert new_assm2._id == 0
    assert new_assm1.fluka_position is None
    assert new_assm2.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert new_assm1.assembly_file.exists()
    assert all(ff.exists() for ff in new_assm1.files)
    assert new_assm1.exists()
    assert new_assm2.assembly_file.exists()
    assert all(ff.exists() for ff in new_assm2.files)
    assert new_assm2.exists()
    assert new_pro_jaw.body_file.exists()
    assert new_pro_jaw.region_file.exists()
    assert new_pro_jaw.material_file.exists()
    assert new_pro_jaw.exists()
    assert new_pro_tank.body_file.exists()
    assert new_pro_tank.region_file.exists()
    assert new_pro_tank.material_file.exists()
    assert new_pro_tank.exists()
    assert new_pro_jaw in xc.FlukaPrototype._registry
    assert new_pro_tank in xc.FlukaPrototype._registry
    assert new_assm1 in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' not in xc.FlukaAssembly._assigned_registry
    assert 'ASSM2' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM2'] is new_assm2
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (unassigned): tag ASSM1 in test series"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (1 inactive element): tag ASSM2 in test series"

    new_assm1.delete()
    assert coll1.assembly is new_assm2
    assert coll2.assembly is None
    assert coll3.assembly is None
    assert not new_assm1.assigned
    assert new_assm2.assigned
    assert not new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == []
    assert new_assm2.elements == [coll1]
    assert new_assm1.active_elements == []
    assert new_assm2.active_elements == []
    assert new_assm1._id is None
    assert new_assm2._id == 0
    assert new_assm1.fluka_position is None
    assert new_assm2.fluka_position == (0.0, 0.0, 0.0, -1000.0, -3000.0, 1000.0)
    assert not new_assm1.assembly_file.exists()
    assert new_assm1.files is None
    assert not new_assm1.exists()
    assert new_assm2.assembly_file.exists()
    assert all(ff.exists() for ff in new_assm2.files)
    assert new_assm2.exists()
    assert new_pro_jaw.body_file.exists()
    assert new_pro_jaw.region_file.exists()
    assert new_pro_jaw.material_file.exists()
    assert new_pro_jaw.exists()
    assert new_pro_tank.body_file.exists()
    assert new_pro_tank.region_file.exists()
    assert new_pro_tank.material_file.exists()
    assert new_pro_tank.exists()
    assert new_pro_jaw in xc.FlukaPrototype._registry
    assert new_pro_tank in xc.FlukaPrototype._registry
    assert new_assm1 not in xc.FlukaPrototype._registry
    assert new_assm2 in xc.FlukaPrototype._registry
    assert 'ASSM1' not in xc.FlukaAssembly._assigned_registry
    assert 'ASSM2' in xc.FlukaAssembly._assigned_registry
    assert xc.FlukaAssembly._assigned_registry['ASSM2'] is new_assm2
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (unassigned): tag ASSM1 in test series <defunct>"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (1 inactive element): tag ASSM2 in test series"

    new_assm2.clear()
    new_assm2.delete()
    assert coll1.assembly is None
    assert coll2.assembly is None
    assert coll3.assembly is None
    assert not new_assm1.assigned
    assert not new_assm2.assigned
    assert not new_assm1.active
    assert not new_assm2.active
    assert new_assm1.elements == []
    assert new_assm2.elements == []
    assert new_assm1.active_elements == []
    assert new_assm2.active_elements == []
    assert new_assm1._id is None
    assert new_assm2._id is None
    assert new_assm1.fluka_position is None
    assert new_assm2.fluka_position is None
    assert not new_assm1.assembly_file.exists()
    assert new_assm1.files is None
    assert not new_assm1.exists()
    assert not new_assm2.assembly_file.exists()
    assert new_assm2.files is None
    assert not new_assm2.exists()
    assert not new_pro_jaw.body_file.exists()
    assert not new_pro_jaw.region_file.exists()
    assert not new_pro_jaw.material_file.exists()
    assert not new_pro_jaw.exists()
    assert not new_pro_tank.body_file.exists()
    assert not new_pro_tank.region_file.exists()
    assert not new_pro_tank.material_file.exists()
    assert not new_pro_tank.exists()
    assert new_pro_jaw not in xc.FlukaPrototype._registry
    assert new_pro_tank not in xc.FlukaPrototype._registry
    assert new_assm1 not in xc.FlukaPrototype._registry
    assert new_assm2 not in xc.FlukaPrototype._registry
    assert 'ASSM1' not in xc.FlukaAssembly._assigned_registry
    assert 'ASSM2' not in xc.FlukaAssembly._assigned_registry
    assert str(new_assm1) == "FlukaAssembly 'ASSM1' (unassigned): tag ASSM1 in test series <defunct>"
    assert str(new_assm2) == "FlukaAssembly 'ASSM2' (unassigned): tag ASSM2 in test series <defunct>"
