"""Zero-Xsuite tests for the optional reporting script."""
import importlib.util
import json
from pathlib import Path


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
