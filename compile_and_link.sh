#!/bin/sh
# Small script to compile and link the fortran module
# and to do other small tasks.
# Tested with gfortran (Sorry!)
ln -s multipole/multipole.py multipole.py
cd f2py_modules
f2py --fcompiler="gfortran" --f90flags="-ffree-form"  extmod_spf_mpole.pyf -c spf_mpole.f90
cd ..
ln -s f2py_modules/extmod_spf_mpole.so extmod_spf_mpole.so
