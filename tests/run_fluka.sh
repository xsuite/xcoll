python ../examples/fluka_init.py
if [ $? -ne 0 ]; then
    echo "FLUKA initialization failed"
    exit 1
fi
pytest test_fluka_assemblies.py  # Cannot be ran in parallel as it creates/deletes files and adapts the prototype registry
if [ $? -ne 0 ]; then
    echo "FLUKA assembly failed"
    exit 1
fi
pytest --html="pytest_fluka_results.html" -vv -n 12 -m fluka
