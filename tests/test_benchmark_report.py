"""Zero-Xsuite tests for the optional reporting script."""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]



def _load_report():
    spec = importlib.util.spec_from_file_location("report", Path(__file__).resolve().parents[1] / "benchmarks" / "report.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_static_dashboard_round_trip(tmp_path):
    mod = _load_report()
    paths = []
    for i, val in enumerate((0.01, 0.011)):
        doc = {
            "datetime": f"2026-11-09T0{i}:00:00Z", "machine_info": {"node": "dev", "system": "Linux", "processor": "CPU"},
            "commit_info": {"id": str(i) * 40, "branch": "main" if i == 0 else "feature"},
            "testbed": "dev", "xsuite_versions": {"xcoll": "0.13.0rc0"},
            "benchmarks": [{"fullname": "benchmarks/bench_everest.py::bench_x", "group": "Everest",
                            "stats": {"median": val, "stddev": 0.0001, "rounds": 5},
                            "extra_info": {"protocol": 2, "particles": 1_000_000}}],
        }
        path = tmp_path / f"{i}.json"
        path.write_text(json.dumps(doc))
        paths.append(path)
    out = mod.create_dashboard(paths, tmp_path / "dashboard.html")
    html = out.read_text()
    assert "Xsuite performance dashboard" in html
    assert "1000000" in html
    assert "benchmarks/bench_everest.py::bench_x" in html
    assert "__RUN_DATA__" not in html


def test_data_is_safely_embedded(tmp_path):
    mod = _load_report()
    payload = {"benchmarks": [{"fullname": "</script><script>alert(1)</script>",
                               "stats": {"median": 1.0}}]}
    inp = tmp_path / "input.json"
    inp.write_text(json.dumps(payload))
    output = mod.create_dashboard([inp], tmp_path / "out.html")
    assert "</script><script>alert(1)</script>" not in output.read_text()


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "benchmarks" / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    # Dataclasses inspect sys.modules for their defining module.
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def record(name="bench_test", lost=.2, particles=1_000_000, **kw):
    ex = {"protocol": 2, "particles": particles, "lost_fraction": lost,
          "survivors_with_kick_gt_0.1urad": .25}
    ex.update(kw)
    return {"name": name, "fullname": f"bench_everest.py::{name}",
            "stats": {"median": .10, "stddev": .001, "rounds": 5},
            "extra_info": ex}


def data(records=None, node="runner-a", testbed="github-ubuntu-24.04-cpu-informational", sha="a"*40):
    return {"datetime": "2026-11-09T12:00:00Z",
            "machine_info": {"node": node, "system": "Linux", "processor": "Intel Example"},
            "commit_info": {"id": sha, "branch": "main"}, "testbed": testbed,
            "benchmarks": records if records is not None else [record()]}


def test_physical_ok():
    m = module("compare_physics")
    r = m.compare_documents(data(), data([record(lost=.203)]))
    assert r[0].status == "OK"


def test_physical_fail_change():
    m = module("compare_physics")
    r = m.compare_documents(data(), data([record(lost=.26)]))
    assert r[0].status == "CHANGED"
    assert "lost_fraction" in r[0].differences[0]


def test_missing_observables_are_not_silently_accepted():
    m = module("compare_physics")
    r = m.compare_documents(data(), data([{"name": "bench_test", "extra_info": {"protocol": 2,"particles":1_000_000}}]))
    assert r[0].status == "UNVERIFIED"


def test_changed_workload_not_comparable():
    m = module("compare_physics")
    r = m.compare_documents(data(), data([record(particles=10_000)]))
    assert r[0].status == "INCOMPATIBLE"


def test_missing_benchmark_not_ignored():
    m = module("compare_physics")
    r = m.compare_documents(data(), data([]))
    assert r[0].status == "MISSING"


def test_lhc_counts_are_normalised():
    m = module("compare_physics")
    a = record(lost=.2, particles=3000, reference_collimator_losses=180, reference_aperture_losses=30)
    b = record(lost=.2, particles=3000, reference_collimator_losses=185, reference_aperture_losses=29)
    assert m.compare_entry(a,b).status == "OK"
    b["extra_info"]["reference_collimator_losses"] = 300
    assert m.compare_entry(a,b).status == "CHANGED"


def test_history_package_single_record(tmp_path):
    m = module("history")
    inp = tmp_path / "result"
    inp.mkdir()
    (inp / "data.json").write_text(json.dumps(data()))
    out = m.pack(inp, tmp_path / "history", "1234", "1", "b"*40, "release/v0.13.0", "quick")
    payload = json.loads(out.read_text())
    assert payload["ci"]["source_branch"] == "release/v0.13.0"
    assert payload["ci"]["source_sha"] == "b"*40
    assert payload["benchmarks"][0]["extra_info"]["lost_fraction"] == .2
    with pytest.raises(FileExistsError):
        m.pack(inp, tmp_path / "history", "1234", "1", "b"*40, "release/v0.13.0", "quick")


def test_history_refuses_multiple_or_empty_runs(tmp_path):
    m = module("history")
    inp = tmp_path / "result"
    inp.mkdir()
    with pytest.raises(ValueError):
        m.pack(inp, tmp_path / "history", "123", "1", "b"*40, "main", "quick")
    for i in (1,2):
        (inp / f"{i}.json").write_text(json.dumps(data()))
    with pytest.raises(ValueError):
        m.pack(inp, tmp_path / "history", "123", "1", "b"*40, "main", "quick")


def test_dashboard_ci_hosts_group_same_machine(tmp_path):
    m = module("report")
    paths=[]
    for i,host in enumerate(("runner-111", "runner-222")):
        result = data(node=host, sha=str(i+1)*40)
        result["datetime"] = f"2026-11-0{i+8}T12:00:00Z"
        f = tmp_path / f"run{i}.json"
        f.write_text(json.dumps(result))
        paths.append(f)
    runs = m.read_runs(paths)
    assert runs[0]["machine"] == runs[1]["machine"]
    assert runs[0]["entries"][0]["name"] == "bench_everest.py::bench_test"
    output = m.create_dashboard(paths, tmp_path / "site" / "index.html")
    html = output.read_text()
    assert "Benchmark history" in html
    assert "physicalStatus" in html
    assert "runner-111" not in html
    assert "run0.json" in html and "run1.json" in html


def test_cli_fails_when_physics_changes(tmp_path):
    first=tmp_path/"a.json";second=tmp_path/"b.json"
    first.write_text(json.dumps(data()))
    second.write_text(json.dumps(data([record(lost=.3)])))
    proc = subprocess.run([sys.executable,str(ROOT/"benchmarks"/"compare.py"),
                           str(first),str(second),"--no-timings"],capture_output=True,text=True)
    assert proc.returncode == 1
    assert "CHANGED" in proc.stdout
