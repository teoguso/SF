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
from __future__ import print_function
from threading import Thread
from sf_modules import *
from outread import *
import numpy as np;
import matplotlib.pylab as plt;
#from scipy.interpolate import interp1d
#from scipy import optimize
import sys
from os.path import isfile, join, isdir
from os import getcwd, pardir, mkdir, chdir
import time

### ============================= ###
###  ==  PROGRAM BEGINS HERE  ==  ###
### ============================= ###
start_time = time.time()
# ======== READING INPUT VARIABLES ======= #
for i in range(52): print('=',end='')
print()
print( " SF :: START")
for i in range(52): print('=',end='')
print()
invar_dict = read_invar()
#print ('%12s, %9.4f' % invar_dict.keys(), invar_dict.values())
print(" "+"===="+" Input variables "+"====")
print()
for key in  invar_dict.keys():
    #print key, invar_dict[key]
    print ('%12s :: %12s' % (key, invar_dict[key]))
for i in range(52): print('=',end='')
print()
#print 52*"="
print (" SF :: END")
for i in range(52): print('=',end='')
print()
#print 52*"="
if len(sys.argv) > 1:
    out_file = sys.argv[1]
else:
    out_file = None
gwout = CodeOutReader(invar_dict['gwcode'],out_file)
print(gwout)
#print(gwout.fname)
"""for x in gwout.hartree: 
    s = ' '.join(['{:10}']*len(x))
    #print(s)
    print(s.format(*x), end='\n')
"""
# ====== READING HARTREE ===== #
hartree = gwout.hartree
# ======== READING WTK ======= #
if invar_dict['gwcode']=='abinit' and gwout.nversion <= 5: # FOR OLDER ABINIT VERSIONS
    wtk = read_wtk()
    invar_dict['wtk'] = wtk
# ======== READING _SIG FILE ======= #
efermi =  float(invar_dict['efermi'])
enmin = float(invar_dict['enmin'])
enmax = float(invar_dict['enmax'])
en, res, ims = read_sigfile(invar_dict)
# Rescale energy if in hartree
print(invar_dict['enhartree'])
enhartree = invar_dict['enhartree']
if enhartree and enhartree != 0:
    print(" ### Converting energies from Hartree to eV ###")
    print(" ### 1 Hartree = 27.2116 eV ###")
    en = 2.0*13.6058*en
#TODO: enmin and emax are unchanged. Check if this is consistent!
# Reset wrt efermi
en = en - efermi
res[:,:] = res[:,:] - efermi
print(" en[0], en[-1], enmin, enmax \n", en[0], en[-1], enmin, enmax)
nkpt =  int(invar_dict['nkpt']) 
minband = int(invar_dict['minband']) 
maxband = int(invar_dict['maxband']) 
nband = maxband - minband +1
invar_dict['nband'] = nband
print(" ### nkpt, nband:", nkpt, nband)
print(" # ------------------------------------------------ # ")
pdos = calc_pdos(invar_dict)
print(" pdos:", pdos)
print(" Size(pdos):",np.size(pdos))
#TODO: Check if consistent use of numpy arrays. 
### ===================================================== ###
print(" # ------------------------------------------------ # ")
# Here we move to a subdirectory to avoid flooding-up the current directory
newdirname = "Spfunctions"
origdir = getcwd() # remember where we are
newdir = join(origdir, newdirname) # Complete path of the new directory
print(" Moving into output directory:\n ", newdir)
if not isdir(newdir) :
    mkdir(newdir)
chdir(newdir)
### ================================= ###
### ===== GW SPECTRAL FUNCTION ====== ###
t_part1 = time.time() - start_time
print(" --- Time spent so far: {} seconds. ---".format(t_part1))
# GW spectral function part (if requested)
if int(invar_dict['calc_gw']) == 1:
    #newen, spftot = calc_spf_gw(minkpt,maxkpt,minband,maxband,wtk,pdos,en,enmin,enmax,res,ims,hartree)
    newen, spftot, allkb = calc_sf_gw(invar_dict,hartree,pdos,en,res,ims)
        ### ==== WRITING OUT GW SPECTRAL FUNCTION === ###
    #newen = newen-efermi
    print(" ### Writing out A(\omega)_GW...  ")
    # Start a new thread
    thread = Thread(target = write_spfkb, args = (invar_dict, newen, allkb))
    thread.start()
    #write_spfkb(invar_dict,newen,allkb)
    sfac = invar_dict['sfactor']
    pfac = invar_dict['pfactor']
    penergy = invar_dict['penergy']
    outname = "spftot_gw"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev"+".dat"
    outfile = open(outname,'w')
    for i in xrange(np.size(newen)):
        outfile.write("%7.4f %15.10e\n"% (newen[i],spftot[i])) # Dump string representations of arrays
    outfile.close()
    print(" A(\omega)_GW written in", outname)
    plt.plot(newen,spftot,label="ftot_gw");

print(" --- Time spent for GW: {} seconds. ---".format(time.time() - t_part1))
print(" --- Time spent so far: {} seconds. ---".format(time.time() - start_time))
# ============================= ###


