# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import pytest
import numpy as np

import xobjects as xo
import xtrack as xt

from xcoll.scattering_routines.engine import BaseEngine


class DummyInterface:
    compiled = True

    def __init__(self):
        self.assert_ready_calls = 0

    def assert_environment_ready(self):
        self.assert_ready_calls += 1


class DummyElement:
    def __init__(
        self,
        name=None,
        jaw=1e-3,
        active=True,
        tracking=True,
    ):
        self.name = name
        self.jaw = jaw
        self.active = active
        self._tracking = tracking


class DummyEngine(BaseEngine):
    _element_classes = (DummyElement,)
    _uses_input_file = False
    _uses_run_folder = False

    def __init__(self, **kwargs):
        self.started_kwargs = None
        self.stop_calls = 0
        self.removed = []
        self.restored = []
        self._running = False

        super().__init__(**kwargs)

        self._interface = DummyInterface()

    def _start_engine(self, **kwargs):
        self.started_kwargs = kwargs
        self._running = True

    def _stop_engine(self, **kwargs):
        self.stop_calls += 1
        self._running = False
        return kwargs

    def _is_running(self):
        return self._running

    def _remove_element(self, el):
        self.removed.append(el)

    def _restore_element(self, el):
        self.restored.append(el)


def make_engine():
    engine = DummyEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    return engine


def test_engine_exposes_physics_settings():
    engine = make_engine()
    assert engine.return_protons is True
    assert engine.return_pions is False
    engine.return_none = True
    engine.return_pions = True
    assert engine._physics_settings.return_pions is True
    assert engine.return_pions is True
    assert engine.return_protons is False
    engine.reset_physics_settings()
    assert engine.return_protons is True
    assert engine.return_pions is False


def test_engine_temporary_physics_settings():
    engine = make_engine()
    element = DummyElement()
    # Defaults before starting.
    assert engine.return_protons is True
    assert engine.return_ions is False
    engine.start(elements=element, return_none=True, return_ions=True)
    assert engine.return_protons is False
    assert engine.return_ions is True
    engine.stop(clean=True)
    # Temporary start() settings must have disappeared.
    assert engine.return_protons is True
    assert engine.return_ions is False


def test_engine_persistent_physics_settings():
    engine = make_engine()
    element = DummyElement()
    engine.return_none = True
    engine.return_pions = True
    engine.start(elements=element)
    assert engine.return_pions is True
    assert engine.return_protons is False
    engine.stop(clean=True)
    # Attribute settings persist.
    assert engine.return_pions is True
    assert engine.return_protons is False


def test_engine_dynamic_defaults_follow_reference():
    engine = make_engine()
    assert engine.return_protons is True
    assert engine.return_pions is False
    engine.particle_ref = xt.Particles("pi+", p0c=450e9)
    assert engine.return_protons is False
    assert engine.return_pions is True


def test_engine_explicit_setting_survives_reference_change():
    engine = make_engine()
    engine.return_pions = False
    engine.particle_ref = xt.Particles("pi+", p0c=450e9)
    assert engine.return_pions is False


def test_engine_pdg_id_api():
    engine = make_engine()
    engine.return_none = True
    engine.return_pdg_id(411)
    assert engine.pdg_id_is_returned(411)
    assert not engine.pdg_id_is_returned(-411)
    engine.dont_return_pdg_id(411)
    assert not engine.pdg_id_is_returned(411)
    engine.reset_physics_settings()
    assert not engine.pdg_id_is_returned(411)


def test_temporary_settings_preserve_dynamic_defaults():
    engine = make_engine()
    element = DummyElement()
    assert engine._physics_settings._return_protons is None
    engine.start(elements=element, return_none=True, return_ions=True)
    engine.stop(clean=True)
    assert engine._physics_settings._return_protons is None
    assert engine.return_protons is True


def test_engine_defaults():
    engine = DummyEngine()
    assert engine.name == "dummy"
    assert engine.particle_ref is None
    assert engine.seed is None
    assert engine.line is None
    assert engine.cwd is None
    assert engine.element_dict == {}
    assert engine.verbose is False
    assert engine.is_running() is False


def test_particle_ref():
    engine = DummyEngine()
    ref = xt.Particles("proton", p0c=7e12)
    engine.particle_ref = ref
    assert engine.particle_ref is not None
    assert engine.particle_ref.pdg_id[0] == 2212
    assert engine.particle_ref.p0c[0] == 7e12
    del engine.particle_ref
    assert engine.particle_ref is None


def test_temporary_particle_ref_is_updated_and_restored():
    engine = make_engine()
    original_engine_ref = engine.particle_ref.copy()
    temporary = xt.Particles("proton", p0c=6e12, mass0=900e6)
    temporary_before = temporary.copy()

    # Pretend that the external code uses a different proton mass.
    engine._masses = {2212: 938.2720813e6}
    engine.start(elements=DummyElement("coll"), particle_ref=temporary)
    # Engine was normalised.
    assert np.isclose(engine.particle_ref.mass0, 938.2720813e6)
    # Caller-owned object was normalised as well.
    assert np.isclose(temporary.mass0, 938.2720813e6)
    engine.stop(clean=True)

    # Persistent engine reference restored.
    assert np.isclose(engine.particle_ref.p0c[0], original_engine_ref.p0c[0])
    assert np.isclose(engine.particle_ref.mass0, original_engine_ref.mass0)
    # Caller-owned temporary reference restored in place.
    assert np.isclose(temporary.p0c[0], temporary_before.p0c[0])
    assert np.isclose(temporary.mass0, temporary_before.mass0)


def test_particle_ref_validation():
    engine = DummyEngine()
    with pytest.raises(ValueError):
        engine.particle_ref = "proton"
    with pytest.raises(ValueError):
        engine.particle_ref = xt.Particles(
            pdg_id=[2212, 2212],
            p0c=[7e12, 7e12],
        )
    with pytest.raises(ValueError, match="valid pdg_id"):
        engine.particle_ref = xt.Particles(p0c=7e12)


def test_seed():
    engine = DummyEngine()
    assert engine.seed is None
    engine.seed = 1234
    assert engine.seed == 1234
    engine.seed = None
    assert engine.seed is None
    with pytest.raises(ValueError):
        engine.seed = -1


class DummyInt32Engine(DummyEngine):
    _int32 = True


def test_seed_int32_casting():
    engine = DummyInt32Engine()
    engine.seed = 123
    assert engine.seed == 123


def test_start_with_single_element():
    engine = make_engine()
    element = DummyElement(name="coll")
    engine.start(elements=element)
    assert engine.is_running()
    assert engine.element_dict == {"coll": element}
    engine.stop(clean=True)


def test_element_gets_generated_name():
    engine = make_engine()
    element = DummyElement(name=None)
    engine.start(elements=element)
    assert element.name == "dummy_el_0"
    assert engine.element_dict == {"dummy_el_0": element}
    engine.stop(clean=True)


def test_explicit_name():
    engine = make_engine()
    element = DummyElement(name="old")
    engine.start(elements=element, names="new")
    assert element.name == "new"
    assert engine.element_dict == {"new": element}
    engine.stop(clean=True)


def test_names_length_mismatch():
    engine = make_engine()
    with pytest.raises(ValueError, match="Length of `elements` and `names`"):
        engine.start(
            elements=[DummyElement(), DummyElement()],
            names=["one"],
        )


def test_duplicate_names():
    engine = make_engine()
    with pytest.raises(ValueError, match="Duplicate names"):
        engine.start(
            elements=[
                DummyElement(name="same"),
                DummyElement(name="same"),
            ]
        )


def test_wrong_element_type():
    engine = make_engine()
    with pytest.raises(ValueError, match="is not a DummyElement"):
        engine.start(elements=object())


def test_inactive_elements_are_ignored_and_restored():
    engine = make_engine()
    active = DummyElement(name="active")
    inactive = DummyElement(name="inactive", active=False)
    engine.start(elements=[active, inactive])
    assert engine.element_dict == {"active": active}
    assert inactive in engine.removed
    engine.stop(clean=True)
    assert inactive in engine.restored
    assert inactive.active is False

def test_element_without_jaw_is_ignored():
    engine = make_engine()
    good = DummyElement(name="good")
    bad = DummyElement(name="bad", jaw=None)
    engine.start(elements=[good, bad])
    assert engine.element_dict == {"good": good}
    engine.stop(clean=True)


class DummyBeamElement(xt.BeamElement):
    _xofields = {
        "jaw": xo.Float64,
        "active": xo.Int8,
        "_tracking": xo.Int8,
    }

class DummyLineEngine(DummyEngine):
    _element_classes = (DummyBeamElement,)


def test_start_from_line():
    engine = DummyLineEngine()
    e1 = DummyBeamElement(jaw=1e-3, active=True, _tracking=True)
    e2 = DummyBeamElement(jaw=2e-3, active=True, _tracking=True)
    line = xt.Line(elements=[e1, e2], element_names=["a", "b"])
    line.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.start(line=line)
    assert set(engine.element_dict) == {"a", "b"}
    engine.stop(clean=True)


def test_start_from_line_selected_names():
    engine = DummyLineEngine()
    e1 = DummyBeamElement(jaw=1e-3, active=True, _tracking=True)
    e2 = DummyBeamElement(jaw=2e-3, active=True, _tracking=True)
    line = xt.Line(elements=[e1, e2], element_names=["a", "b"])
    line.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.start(line=line, names=["b"])
    assert list(engine.element_dict) == ["b"]
    engine.stop(clean=True)


def test_line_and_elements_are_mutually_exclusive():
    engine = DummyLineEngine()
    e1 = DummyBeamElement(jaw=1e-3, active=True, _tracking=True)
    e2 = DummyBeamElement(jaw=2e-3, active=True, _tracking=True)
    line = xt.Line(elements=[e1, e2], element_names=["a", "b"])
    line.particle_ref = xt.Particles("proton", p0c=7e12)
    with pytest.raises(ValueError, match="Cannot provide both"):
        engine.start(line=line, elements=e1)


def test_particle_ref_taken_from_line():
    engine = DummyLineEngine()
    e1 = DummyBeamElement(jaw=1e-3, active=True, _tracking=True)
    e2 = DummyBeamElement(jaw=2e-3, active=True, _tracking=True)
    line = xt.Line(elements=[e1, e2], element_names=["a", "b"])
    line.particle_ref = xt.Particles("proton", p0c=7e12)
    assert engine.particle_ref is None
    engine.start(line=line)
    assert engine.particle_ref.p0c[0] == 7e12
    engine.stop(clean=True)
    assert engine.particle_ref is None


def test_existing_engine_particle_ref_overrides_line_temporarily():
    engine = DummyLineEngine()
    engine.particle_ref = xt.Particles("proton", p0c=6e12)
    e1 = DummyBeamElement(jaw=1e-3, active=True, _tracking=True)
    e2 = DummyBeamElement(jaw=2e-3, active=True, _tracking=True)
    line = xt.Line(elements=[e1, e2], element_names=["a", "b"])
    line.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.start(line=line)
    assert engine.particle_ref.p0c[0] == 6e12
    assert line.particle_ref.p0c[0] == 6e12
    engine.stop(clean=True)
    assert engine.particle_ref.p0c[0] == 6e12
    assert line.particle_ref.p0c[0] == 7e12


def test_temporary_particle_ref():
    engine = make_engine()
    original = engine.particle_ref.copy()
    temporary = xt.Particles("proton", p0c=6e12)
    engine.start(elements=DummyElement("coll"), particle_ref=temporary)
    assert engine.particle_ref.p0c[0] == 6e12
    engine.stop(clean=True)
    assert engine.particle_ref.p0c[0] == original.p0c[0]


def test_temporary_verbose():
    engine = make_engine()
    element = DummyElement("coll")
    assert engine.verbose is False
    engine.start(elements=element, verbose=True)
    assert engine.verbose is True
    engine.stop(clean=True)
    assert engine.verbose is False


def test_temporary_seed_is_restored():
    engine = make_engine()
    engine.seed = 123
    engine.start(elements=DummyElement("coll"), seed=456)
    assert engine.seed == 456
    engine.stop(clean=True)
    assert engine.seed == 123


def test_seed_generated_if_missing():
    engine = make_engine()
    assert engine.seed is None
    engine.start(elements=DummyElement("coll"))
    assert engine.seed is not None
    engine.stop(clean=True)
    assert engine.seed is None


def test_start_and_stop():
    engine = make_engine()
    element = DummyElement("coll")
    assert not engine.is_running()
    engine.start(elements=element)
    assert engine.is_running()
    assert engine._tracking_initialised is False
    assert engine.interface.assert_ready_calls == 1
    engine.stop(clean=True)
    assert not engine.is_running()
    assert engine.element_dict == {}


def test_start_when_already_running_is_noop():
    engine = make_engine()
    element = DummyElement("coll")
    engine.start(elements=element)
    first_kwargs = engine.started_kwargs
    engine.start(elements=element)
    assert engine.started_kwargs is first_kwargs
    engine.stop(clean=True)


def test_backend_specific_kwargs_are_forwarded():
    engine = make_engine()
    engine.start(elements=DummyElement("coll"), custom_backend_option=123)
    assert engine.started_kwargs == {"custom_backend_option": 123}
    engine.stop(clean=True)


class DummyFolderEngine(DummyEngine):
    _uses_run_folder = True


def test_cwd_created_and_reset(tmp_path):
    engine = DummyFolderEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    requested = tmp_path / "run"
    engine.start(elements=DummyElement("coll"), cwd=requested)
    assert engine.cwd == requested
    assert requested.exists()
    engine.stop(clean=True)
    assert engine.cwd is None
    assert not requested.exists()


def test_existing_cwd_gets_unique_suffix(tmp_path):
    requested = tmp_path / "run"
    requested.mkdir()
    engine = DummyFolderEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.start(elements=DummyElement("coll"), cwd=requested)
    assert engine.cwd.name == "run_0000"
    engine.stop(clean=True)


class DummyInputEngine(DummyEngine):
    _uses_input_file = True
    _uses_run_folder = True
    _multiple_input_files = False

    def _generate_input_file(self, **kwargs):
        path = self.cwd / "dummy.in"
        path.write_text("generated")
        return path, kwargs

    def _match_input_file(self):
        self.match_calls = getattr(self, "match_calls", 0) + 1

    def _get_input_files_to_clean(self, input_file=None, cwd=None, **kwargs):
        if cwd is None:
            return []
        return [input_file, cwd / "auxiliary.in", cwd / "input_aux"]

    def _get_output_files_to_clean(self, input_file=None, cwd=None, **kwargs):
        if cwd is None:
            return []
        return [cwd / "dummy.out", cwd / "output_aux"]

    def _all_input_files(self, input_file=None):
        if input_file is None:
            input_file = self.input_file
        return [input_file]


def test_start_generates_input_file(tmp_path):
    engine = DummyInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.start(
        elements=DummyElement("coll"),
        cwd=tmp_path / "run",
        clean=False,
    )
    assert engine.input_file.exists()
    assert engine.input_file.read_text() == "generated"
    assert engine.match_calls == 1
    engine.stop(clean=True)


def test_start_with_existing_input_file(tmp_path):
    source = tmp_path / "input.in"
    source.write_text("custom")
    engine = DummyInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.start(
        elements=DummyElement("coll"),
        input_file=source,
        cwd=tmp_path / "run",
        clean=False,
    )
    assert engine.input_file.read_text() == "custom"
    engine.stop(clean=True)


def test_missing_input_file():
    engine = DummyInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    with pytest.raises(ValueError, match="does not exist"):
        engine.start(
            elements=DummyElement("coll"),
            input_file="does_not_exist.in",
        )


def test_generate_input_file_does_not_start_engine(tmp_path):
    engine = DummyInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    path = engine.generate_input_file(
        elements=DummyElement("coll"),
        filename=tmp_path / "saved.in",
    )
    assert path == tmp_path / "saved.in"
    assert path.exists()
    assert not engine.is_running()
    assert engine.element_dict == {}


def test_ready_to_track():
    engine = make_engine()
    coll = DummyElement("coll")
    engine.start(elements=coll)
    particles = xt.Particles(
        "proton",
        p0c=7e12,
        x=[0, 0],
        _capacity=10,
    )
    assert engine.assert_ready_to_track_or_skip(coll, particles)
    engine.stop(clean=True)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"active": False},
        {"tracking": False},
        {"jaw": None},
    ],
)
def test_ready_to_track_skip_element(kwargs):
    engine = make_engine()
    coll = DummyElement("coll", **kwargs)
    particles = xt.Particles(
        "proton",
        p0c=7e12,
        x=[0, 0],
        _capacity=10,
    )
    assert not engine.assert_ready_to_track_or_skip(coll, particles)


def test_ready_to_track_skips_empty_particles():
    engine = make_engine()
    coll = DummyElement("coll")
    engine.start(elements=coll)
    particles = xt.Particles(
        "proton",
        p0c=7e12,
        x=[],
        _capacity=10,
    )
    assert not engine.assert_ready_to_track_or_skip(coll, particles)
    engine.stop(clean=True)


def test_ready_to_track_requires_running_engine():
    engine = make_engine()
    coll = DummyElement("coll")
    particles = xt.Particles(
        "proton",
        p0c=7e12,
        x=[0, 0],
        _capacity=10,
    )
    with pytest.raises(RuntimeError, match="not yet running"):
        engine.assert_ready_to_track_or_skip(coll, particles)


def test_ready_to_track_requires_secondary_capacity():
    engine = make_engine()
    coll = DummyElement("coll")
    engine.start(elements=coll)
    particles = xt.Particles(
        "proton",
        p0c=7e12,
        x=[0, 0],
        _capacity=2,
    )
    with pytest.raises(ValueError, match="capacity equal to size"):
        engine.assert_ready_to_track_or_skip(coll, particles)
    engine.stop(clean=True)


def test_ready_to_track_requires_pdg_ids():
    engine = make_engine()
    coll = DummyElement("coll")
    engine.start(elements=coll)
    particles = xt.Particles(
        "proton",
        p0c=7e12,
        x=[0, 0],
        _capacity=10,
    )
    particles.pdg_id = 0
    with pytest.raises(ValueError, match="pdg_id"):
        engine.assert_ready_to_track_or_skip(coll, particles)
    engine.stop(clean=True)


def test_ready_to_track_rejects_mass_mismatch():
    engine = make_engine()
    coll = DummyElement("coll")
    engine.start(elements=coll)
    particles = xt.Particles(
        p0c=7e12,
        mass0=0.9 * engine.particle_ref.mass0,
        q0=engine.particle_ref.q0,
        pdg_id=2212,
        x=[0, 0],
        _capacity=10,
    )
    with pytest.raises(ValueError, match="reference mass"):
        engine.assert_ready_to_track_or_skip(coll, particles)
    assert not engine.is_running()


def test_ready_to_track_rejects_charge_mismatch():
    engine = make_engine()
    coll = DummyElement("coll")
    engine.start(elements=coll)
    particles = xt.Particles(
        p0c=7e12,
        mass0=engine.particle_ref.mass0,
        q0=2,
        pdg_id=2212,
        x=[0, 0],
        _capacity=10,
    )
    with pytest.raises(ValueError, match="reference charge"):
        engine.assert_ready_to_track_or_skip(coll, particles)
    assert not engine.is_running()


def test_ready_to_track_resynchronises_known_mass():
    engine = make_engine()
    coll = DummyElement("coll")
    engine.start(elements=coll)
    old_mass = engine.particle_ref.mass0
    new_mass = old_mass * 1.001
    engine._masses = {2212: new_mass}
    # Mimic the situation after the backend has corrected its
    # reference-particle mass.
    engine.particle_ref.mass0 = new_mass
    engine.particle_ref._update_refs(p0c=engine.particle_ref.p0c[0])
    particles = xt.Particles(
        p0c=7e12,
        mass0=old_mass,
        q0=1,
        pdg_id=2212,
        x=[0, 0],
        _capacity=10,
    )
    assert engine.assert_ready_to_track_or_skip(coll, particles)
    assert np.isclose(particles.mass0, new_mass)
    engine.stop(clean=True)


def test_engine_cleaning(tmp_path):
    engine = DummyInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.start(
        elements=DummyElement("coll"),
        cwd=tmp_path / "run",
        clean=False,
    )
    cwd = engine.cwd
    input_file = engine.input_file
    auxiliary_input = cwd / "auxiliary.in"
    output_file = cwd / "dummy.out"
    input_dir = cwd / "input_aux"
    output_dir = cwd / "output_aux"
    auxiliary_input.write_text("aux")
    output_file.write_text("output")
    input_dir.mkdir()
    (input_dir / "foo").write_text("foo")
    output_dir.mkdir()
    (output_dir / "foo").write_text("foo")

    # Default cleaning preserves the actual input file.
    engine.clean(clean_all=False)
    assert input_file.exists()
    assert not auxiliary_input.exists()
    assert not output_file.exists()
    assert not input_dir.exists()
    assert not output_dir.exists()

    # Recreate some auxiliary/output files.
    auxiliary_input.write_text("aux")
    output_file.write_text("output")

    # clean_all also removes the actual input.
    engine.clean(clean_all=True)
    assert not input_file.exists()
    assert not auxiliary_input.exists()
    assert not output_file.exists()
    engine.stop(clean=True)


class DummyMultiInputEngine(DummyInputEngine):
    _multiple_input_files = True

    def _generate_input_file(self, **kwargs):
        first = self.cwd / "dummy.in"
        second = self.cwd / "dummy_aux.in"
        first.write_text("generated primary")
        second.write_text("generated auxiliary")
        return [first, second], kwargs

    def _all_input_files(self, input_file=None):
        if input_file is None:
            input_file = self.input_file
        if not isinstance(input_file, (list, tuple)):
            input_file = [input_file]
        return list(input_file)


def test_start_with_multiple_input_files(tmp_path):
    first = tmp_path / "first.in"
    second = tmp_path / "second.in"
    first.write_text("first")
    second.write_text("second")

    engine = DummyMultiInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)

    engine.start(
        elements=DummyElement("coll"),
        input_file=[first, second],
        cwd=tmp_path / "run",
        clean=False,
    )

    assert isinstance(engine.input_file, list)
    assert len(engine.input_file) == 2
    assert [path.name for path in engine.input_file] == [
        "first.in",
        "second.in",
    ]
    assert engine.input_file[0].read_text() == "first"
    assert engine.input_file[1].read_text() == "second"

    engine.stop(clean=True)


def test_start_multi_input_accepts_single_path(tmp_path):
    source = tmp_path / "input.in"
    source.write_text("single")

    engine = DummyMultiInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)

    engine.start(
        elements=DummyElement("coll"),
        input_file=source,
        cwd=tmp_path / "run",
        clean=False,
    )

    assert isinstance(engine.input_file, list)
    assert len(engine.input_file) == 1
    assert engine.input_file[0].read_text() == "single"

    engine.stop(clean=True)


def test_generate_multiple_input_files(tmp_path):
    engine = DummyMultiInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)

    paths = engine.generate_input_file(
        elements=DummyElement("coll"),
        cwd=tmp_path / "run",
        filename=tmp_path / "saved.in",
    )

    assert isinstance(paths, list)
    assert len(paths) == 2
    assert paths[0] == tmp_path / "saved.in"
    assert paths[0].read_text() == "generated primary"
    assert paths[1].parent == tmp_path
    assert paths[1].read_text() == "generated auxiliary"

    assert not engine.is_running()
    assert engine.element_dict == {}


class DummyResetEngine(DummyEngine):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.backend_setting = 123

    def _reset_engine_settings(self):
        self.backend_setting = 0


def test_engine_reset():
    engine = DummyResetEngine()

    line = xt.Line()
    engine.line = line
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.seed = 12345
    engine.verbose = True
    engine.return_none = True
    engine.return_pions = True
    engine.return_pdg_id(411)
    engine._element_index = 17

    engine.reset(clean=True)

    assert engine.line is None
    assert engine.particle_ref is None
    assert engine.seed is None
    assert engine.verbose is False
    assert engine.backend_setting == 0

    assert not engine._physics_settings._extra_pdg_ids_to_return
    assert not engine._physics_settings._extra_pdg_ids_to_kill

    # The name generator is intentionally process-lifetime state.
    assert engine._element_index == 17


def test_ready_to_track_requires_compiled_interface():
    engine = make_engine()
    coll = DummyElement("coll")
    engine.start(elements=coll)

    engine.interface.compiled = False

    particles = xt.Particles(
        "proton",
        p0c=7e12,
        x=[0, 0],
        _capacity=10,
    )

    with pytest.raises(RuntimeError, match="interface not compiled"):
        engine.assert_ready_to_track_or_skip(coll, particles)

    assert not engine.is_running()


class FaultInjectionEngine(DummyEngine):
    def __init__(self, fail_at, **kwargs):
        self.fail_at = fail_at
        super().__init__(**kwargs)

    def _pre_start(self, **kwargs):
        if self.fail_at == "pre_start":
            raise RuntimeError("Injected failure in _pre_start")
        return kwargs

    def _pre_input(self, **kwargs):
        if self.fail_at == "pre_input":
            raise RuntimeError("Injected failure in _pre_input")
        return kwargs

    def _start_engine(self, **kwargs):
        if self.fail_at == "start":
            # Deliberately mimic a partially started backend.
            self._running = True
            raise RuntimeError("Injected failure in _start_engine")
        return super()._start_engine(**kwargs)


class FaultInjectionInputEngine(DummyInputEngine):
    def _generate_input_file(self, **kwargs):
        raise RuntimeError("Injected failure in _generate_input_file")


@pytest.mark.parametrize("fail_at", ["pre_input", "start"])
def test_failed_start_restores_temporary_state(fail_at):
    engine = FaultInjectionEngine(fail_at)

    persistent_ref = xt.Particles("proton", p0c=7e12)
    engine.particle_ref = persistent_ref
    engine.seed = 123

    temporary_ref = xt.Particles(
        "proton",
        p0c=6e12,
        mass0=900e6,
    )
    temporary_before = temporary_ref.copy()

    engine._masses = {
        2212: 938.2720813e6,
    }

    engine.return_none = True
    engine.return_pions = True

    with pytest.raises(RuntimeError, match="Injected failure"):
        engine.start(
            elements=DummyElement("coll"),
            particle_ref=temporary_ref,
            seed=456,
            verbose=True,
            return_none=True,
            return_ions=True,
        )

    assert not engine.is_running()
    assert engine.seed == 123
    assert engine.verbose is False
    assert engine.element_dict == {}

    assert np.isclose(engine.particle_ref.p0c[0], persistent_ref.p0c[0])

    assert np.isclose(
        temporary_ref.p0c[0],
        temporary_before.p0c[0],
    )
    assert np.isclose(
        temporary_ref.mass0,
        temporary_before.mass0,
    )

    # Persistent physics settings restored.
    assert engine.return_pions is True
    assert engine.return_ions is False


def test_failed_generate_input_file_restores_state(tmp_path):
    engine = FaultInjectionInputEngine()
    engine.particle_ref = xt.Particles("proton", p0c=7e12)
    engine.seed = 123

    with pytest.raises(
        RuntimeError,
        match="Injected failure in _generate_input_file",
    ):
        engine.generate_input_file(
            elements=DummyElement("coll"),
            cwd=tmp_path / "run",
            seed=456,
            verbose=True,
            return_none=True,
            return_ions=True,
        )

    assert engine.seed == 123
    assert engine.verbose is False
    assert engine.cwd is None
    assert engine.element_dict == {}
    assert not engine.is_running()


def test_line_particle_ref_is_updated_and_restored():
    engine = DummyLineEngine()
    element = DummyBeamElement(jaw=1e-3, active=True, _tracking=True)
    line = xt.Line(elements=[element], element_names=["coll"])
    line.particle_ref = xt.Particles("proton", p0c=7e12, mass0=900e6)
    original_ref = line.particle_ref._resolved
    original_snapshot = original_ref.copy()
    engine._masses = {2212: 938.2720813e6}
    engine.start(line=line)
    assert np.isclose(engine.particle_ref.mass0, 938.2720813e6)
    assert np.isclose(line.particle_ref.mass0, 938.2720813e6)
    engine.stop(clean=True)
    assert engine.particle_ref is None
    assert np.isclose(original_ref.mass0, original_snapshot.mass0)
    assert np.isclose(original_ref.p0c[0], original_snapshot.p0c[0])


class FailingLineEngine(DummyLineEngine):
    def _start_engine(self, **kwargs):
        self._running = True
        raise RuntimeError(
            "Injected start failure"
        )


def test_line_particle_ref_restored_after_failed_start():
    engine = FailingLineEngine()
    element = DummyBeamElement(jaw=1e-3, active=True, _tracking=True)
    line = xt.Line(elements=[element], element_names=["coll"])
    line.particle_ref = xt.Particles("proton", p0c=7e12, mass0=900e6)
    original_ref = line.particle_ref._resolved
    snapshot = original_ref.copy()
    engine._masses = {2212: 938.2720813e6}
    with pytest.raises(RuntimeError, match="Injected start failure"):
        engine.start(line=line)
    assert not engine.is_running()
    assert engine.particle_ref is None
    assert np.isclose(original_ref.mass0, snapshot.mass0)
    assert np.isclose(original_ref.p0c[0], snapshot.p0c[0])


def test_duplicate_names_checked_before_deactivation():
    engine = make_engine()
    active = DummyElement(name="same", active=True)
    inactive = DummyElement(name="same", active=False)
    with pytest.raises(ValueError, match="Duplicate names"):
        engine.start(elements=[active, inactive])
    # Duplicate validation happened before element mutation.
    assert inactive not in engine.removed
    assert inactive.active is False
