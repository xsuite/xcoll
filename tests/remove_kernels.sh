#!/bin/bash
path=$(python -c 'import xsuite; print(xsuite.__file__)')
path=${path%__init__.py}

rm ${path}lib/*.{so,c,json}
