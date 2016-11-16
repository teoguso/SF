#!/bin/sh
### Written by Matteo Guzzo ###
### A.D. MMXIV (2014)       ###
# Small script that removes logic links
# and does other small tasks.
# Basically it is the 'undo' 
# of 'compile_and_link.sh'.
echo " Removing link 'multipole.py'..."
rm multipole.py
echo " Removing link 'outread.py'..."
rm outread.py
echo " Removing link 'broad.py'..."
rm broad.py
echo " Removing links  to compiled f2py objects..."
rm ./*.so
echo " Changing directory to 'f2py_modules'..."
cd f2py_modules
echo " Removing compiled fortran object..."
rm ./*.so
cd ..
echo " Removing test spectral functions..."
for i in `ls -1 test` 
do rm -r test/$i/Spfunctions
done
