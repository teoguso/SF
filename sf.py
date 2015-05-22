#!/usr/bin/env python
"""
### Written by Matteo Guzzo ###
### A.D. MMXV (2015)       ###
New version, trying to give it a more c++/object-orented feel.
List of files needed:
- invar.in with input variables.
- _SIG for the self-energy.
- s.dat, p_even.dat, p_odd.dat, d_even.dat, etc. 
for the orbital character and symmetries.
- cs*.dat for the photon cross sections.

NOT ANYMORE - hartree.dat or elda.dat and vxc.dat for the hartree energies.
NOT ANYMORE - wtk.dat for the k-points weights.
TODO - a_wp.dat for the extrinsic/interference effects and additional lifetime.
"""
from sf_modules import *
from outread import *
import numpy as np;
import matplotlib.pylab as plt;
#from scipy.interpolate import interp1d
#from scipy import optimize
import sys
from os.path import isfile, join, isdir
from os import getcwd, pardir, mkdir, chdir

### ============================= ###
###  ==  PROGRAM BEGINS HERE  ==  ###
### ============================= ###

# ======== READING INPUT VARIABLES ======= #
print 52*"="
print " SF :: START"
print 52*"="
invar_dict = read_invar()
#print ('%12s, %9.4f' % invar_dict.keys(), invar_dict.values())
print " "+4*"="+" Input variables "+4*"="
for key in  invar_dict.keys():
    #print key, invar_dict[key]
    print ('%12s :: %12s' % (key, invar_dict[key]))
print 52*"="
print " SF :: END"
print 52*"="

gwout = CodeOutReader(invar_dict['gwcode'])

print gwout.fname
print gwout.hartree
hartree = gwout.hartree
#TODO implement a __str__ method for CodeOutReader 
# that prints out a well-formatted set of info.

# ====== READING HARTREE ===== #
#hartree = read_hartree()
#hartree = hartree # - efermi
# ======== READING WTK ======= #
wtk = read_wtk()
# ======== READING _SIG FILE ======= #
#en, res, ims = read_sigfile(nkpt,nband,sigfilename)
#print " enmin, enmax"
#print  enmin, enmax
enmit = float(invar_dict['enmin']) + float(invar_dict['efermi'])
enmat = float(invar_dict['enmax']) + float(invar_dict['efermi'])
minkpt = int(invar_dict['minkpt']) 
maxkpt = int(invar_dict['maxkpt']) 
minband = int(invar_dict['minband']) 
maxband = int(invar_dict['maxband']) 
sigfilename = invar_dict['sigmafile']
en, res, ims = read_sigfile2(sigfilename,enmit,enmat,minkpt,maxkpt,minband,maxband)
# Rescale energy if in hartree
if enhartree is not None:
    print " ### Converting energies from Hartree to eV ###"
    print " ### 1 Hartree = 27.2116 eV ###"
    en = 2.0*13.6058*en
# Reset wrt efermi
en = en - efermi
res[:,:] = res[:,:] - efermi
print "en[0], en[-1], enmin, enmax"
print en[0], en[-1], enmin, enmax
print " ### nkpt, nband:", nkpt, nband
print " # ------------------------------------------------ # ";
if penergy != 0:
    # ======== CROSS SECTIONS ======= #
    cs = read_cross_sections(penergy)
    # ====== BAND TYPE AND SYMMETRY ==== #
    sp = read_band_type_sym(sfac,pfac,nband)
    # ===== EFFECTIVE STATE-DEPENDENT PREFACTOR ==== #
    pdos = 10000.*np.dot(cs,sp)
else:
    pdos=np.ones((nband))
print " pdos:", pdos
print " Size(pdos):",np.size(pdos)

