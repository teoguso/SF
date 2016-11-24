### Written by Matteo Guzzo ###
### A.D. MMXIV (2014)       ###
# Small script to compile and link the fortran module
# and to do other small tasks.
# Tested with gfortran (Sorry!)
cd f2py_modules
f2py --fcompiler="gfortran" --f90flags="-fopenmp -ffree-form -ffast-math -funroll-loops -O3 -march=corei7-avx" -lgomp extmod_spf_mpole.pyf -c spf_mpole.f90
f2py -c --fcompiler="gfortran" --f90flags="-fopenmp -ffree-form -ffast-math -funroll-loops -O3 -march=corei7-avx" -lgomp calc_ct_f2.pyf calc_ct_fort.f90
#f2py --fcompiler="ifort" extmod_spf_mpole.pyf -c spf_mpole.f90
#f2py --f90exec=/usr/global/intel/bin/ifort extmod_spf_mpole.pyf -c spf_mpole.f90
cd ..
