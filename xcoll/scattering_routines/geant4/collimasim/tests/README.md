# Test suite for collimasim

To run all the tests successfully, both `Xsuite` and `pyAT` must be installed

To run with pytest, the extension `pytest-forked` must be installed

```
  pip install pytest-forked
```

and then 

```
  pytest --forked --continue-on-collection-errors
```
in the `collimasim/tests` directory to run the tests.

This is because collimasim is currently not re-entry safe. 
Once an instance of BDSIM has been launched in a process, 
no further instances can be launched in that process. 
This appears to be due to the global state not being reset 
properly after the BDISM instance terminates. 
To be investigated in the future. 
The flag --continue-on-collection-errors ensures that 
the Xsuite tests run if pyAT isn't isntalled and vice versa.