#pytest -n 14  --html="pytest_fluka_results.html" --self-contained-html test_fluka.py test_lossmap.py::test_lossmap_fluka test_jaw_position.py::test_fluka test_fluka_many_particles.py
pytest  --html="pytest_fluka_results.html" --self-contained-html test_fluka.py test_lossmap.py::test_lossmap_fluka test_jaw_position.py::test_fluka test_fluka_many_particles.py
