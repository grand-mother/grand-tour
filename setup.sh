#!/bin/bash

# Script base directory.
basedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Python module.
lib_dir=$basedir/lib/python
[[ "$PYTHONPATH" =~ "${lib_dir}" ]] || export PYTHONPATH=${lib_dir}:$PYTHONPATH
