#!/bin/bash

cat pytest_results.xml | python -c 'import sys, xml.dom.minidom; print(xml.dom.minidom.parseString(sys.stdin.read()).toprettyxml())'

