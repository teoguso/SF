#!/bin/sh
# Small script that removes logic links
# and does other small tasks.
# Basically it is the 'undo' 
# of 'compile_and_link.sh'.
rm multipole.py
rm extmod_spf_mpole.so
cd f2py_modules
rm extmod_spf_mpole.so
cd ..
