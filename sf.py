#!/usr/bin/env python
"""
### Written by Matteo Guzzo ###
### A.D. MMXV (2015)       ###
New version, trying to give it a more c++/object-oriented feel.
List of files needed:
- invar.in with input variables.
- Output file from GW calculation for wtk, hartree, eqp and whatnot.
- _SIG/SelfXC.dat for the self-energy.
- s.dat, p_even.dat, p_odd.dat, d_even.dat, etc. 
for the orbital character and symmetries.
- cs*.dat for the photon cross sections.

TODO: Fix the inconsistency with ik,ib/ikeff,ibeff: Easiest way to do it:
    Read all kpt and bd, and then read just between minkpt,maxkpt/minbd,maxbd
ONLY FOR LEGACY - hartree.dat or elda.dat and vxc.dat for the hartree energies.
ONLY FOR LEGACY - wtk.dat for the k-points weights.
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
if len(sys.argv) > 1: 
    fname = sys.argv[1]
    invar_dict = read_invar(fname)
else:
    invar_dict = read_invar()
dict_c = invar_dict # So as not to change invar_dict
#print(('%12s, %9.4f' % invar_dict.keys(), invar_dict.values()))
print(" "+"===="+" Input variables "+"====")
print()
for key in sorted(invar_dict.keys()):
    #print(key, invar_dict[key])
    print(('%12s :: %12s' % (key, invar_dict[key])))
for i in range(52): print('=',end='')
print()
#print(52*"=")
print((" SF :: END"))
for i in range(52): print('=',end='')
print()
#print(52*"=")
if len(sys.argv) > 2:
    out_file = sys.argv[2]
else:
    out_file = None
# TODO: remove the debugging prints here below
if invar_dict['gwcode'] == 'abinit':
   #print("1")
    gwout = AbinitOutReader(filename = out_file,is_sc = invar_dict['is_sc']) 
elif invar_dict['gwcode'] == 'exciting':
   #print("2")
    gwout = ExcitingOutReader() 
else: 
   #print("3")
    gwout = CodeOutReader(invar_dict['gwcode'],filename=out_file,is_sc=invar_dict['is_sc'])
print(gwout)
#print(gwout.fname)
"""for x in gwout.hartree: 
    s = ' '.join(['{:10}']*len(x))
    #print(s)
    print(s.format(*x), end='\n')
"""
# ---------------------------------- #

# ======== READING _SIG FILE ======= #
en, res, ims, sig_bdgw = read_sigfile(dict_c)
dict_c['sig_bdgw'] = sig_bdgw
# ---------------------------------- #

# ====== READING HARTREE ===== #
hartree = gwout.hartree
# ---------------------------------- #

# ======== READING WTK ======= #
if 'add_wtk' in dict_c and int(dict_c['add_wtk']) == 0:
    print("K-point weights are neglected, i.e. all equal to 1.")
    dict_c['wtk'] = [1 for i in range(len(hartree[:][0]))]
elif dict_c['gwcode']=='abinit' and gwout.nversion <= 5 \
        and not 'wtk' in dict_c: # FOR OLDER ABINIT VERSIONS
    wtk = read_wtk()
    dict_c['wtk'] = wtk
elif not 'wtk' in dict_c:
    dict_c['wtk'] = gwout.var_dict['wtk']
# ---------------------------------- #

# Reset wrt efermi
efermi =  float(dict_c['efermi'])
enmin = float(dict_c['enmin'])
enmax = float(dict_c['enmax'])
en = en - efermi
res[:,:] = res[:,:] - efermi
print(" en[0], en[-1], enmin, enmax \n", en[0], en[-1], enmin, enmax)
# ---------------------------------- #

# TODO: This part is broken (it would work if only en was in hartree, but not all the other quantities)
# Rescale energy if in hartree
#print(invar_dict['enhartree'])
#enhartree = int(invar_dict['enhartree'])
#if enhartree and enhartree != 0:
#    print(" ### Converting energies from Hartree to eV ###")
#    print(" ### 1 Hartree = 27.2116 eV ###")
#    en = 2.0*13.6058*en
##TODO: enmin and emax are unchanged. Check if this is consistent!
nkpt =  int(invar_dict['nkpt']) 
minband = int(invar_dict['minband']) 
maxband = int(invar_dict['maxband']) 
nband = maxband - minband +1
invar_dict['nband'] = nband
print(" ### nkpt, nband:", nkpt, nband)
print(" # ------------------------------------------------ # ")
pdos = calc_pdos(dict_c,res)
print(" pdos:", pdos)
print(" Size(pdos):",np.size(pdos))
#TODO: Check if consistent use of numpy arrays. 
### ===================================================== ###
print(" # ------------------------------------------------ # ")
# ---------------------------------- #

# Here we move to a subdirectory to avoid flooding-up the current directory
newdirname = "Spfunctions"
origdir = getcwd() # remember where we are
newdir = join(origdir, newdirname) # Complete path of the new directory
print(" Moving into output directory:\n ", newdir)
if not isdir(newdir) :
    mkdir(newdir)
chdir(newdir)
# ---------------------------------- #

### WRITING OUT HARTREE ###
with open('hartree.dat','w') as of:
    for ik in range(len(hartree)):
        for ib in range(len(hartree[0])):
            of.write("%14.5f" % (hartree[ik][ib]))
        of.write("\n")
# ---------------------------------- #

### ================================= ###
### ===== GW SPECTRAL FUNCTION ====== ###
t_pregw = time.time() 
#print(" --- Time spent so far: {} seconds. ---".format(t_part1))
sfac = dict_c['sfactor']
pfac = dict_c['pfactor']
penergy = dict_c['penergy']
minkpt = int(dict_c['minkpt'])
maxkpt = int(dict_c['maxkpt'])
npoles = int(dict_c['npoles'])
extinf = int(dict_c['extinf'])
print("penergy = ", penergy)
print("minkpt = ", minkpt)
print("maxkpt = ", maxkpt)
print("minband = ", minband)
print("maxband = ", maxband)
print("npoles = ", npoles)
print("extinf = ", extinf)
# GW spectral function part 
    #newen, spftot = calc_spf_gw(minkpt,maxkpt,minband,maxband,wtk,pdos,en,enmin,enmax,res,ims,hartree)
# allkb contains A_GW, Re(Sigma), w-e_H-Re(Sigma), Im(Sigma)
# for each ik,ib on the 'newen' array of energies. 
# only A_GW is multiplied by wtk*pdos
if minband < sig_bdgw[0]: 
    print("WARNING: requested first band", minband, " is below the lowest available GW band.")
    dict_c['minband'] = sig_bdgw[0]
    minband = int(dict_c['minband']) 
    print("The first available band (", minband, ") will be taken.")
    if maxband < minband:
        maxband = minband
        print("WARNING: maxband was below minband and was changed accordingly.")
if maxband > sig_bdgw[1]:
    print("WARNING: requested last band", maxband, " is above the highest available GW band.")
    dict_c['maxband'] = sig_bdgw[1]
    maxband = int(dict_c['maxband']) 
    print("The last available band", maxband, " will be taken.")
    if minband > maxband:
        minband = maxband
        print("WARNING: minband was above maxband and was changed accordingly.")
# Maybe a good global variable?
dict_c['bdrange'] = range(minband - 1,maxband)
dict_c['kptrange'] = range(minkpt - 1, maxkpt)

### RESTART??? ###
plt.figure(1) # This is the figure for the plot at the end of the calculation
if int(dict_c['restart']) == 1:  
    print("### RESTART: READING SPECTRAL FUNCTIONS ###")
#if int(dict_c['restart']) == 1:
    if int(dict_c['calc_gw']) == 1:
        print("### GW SPECTRAL FUNCTIONS ###")
        en, ftot_a = read_sf(dict_c, pdos, 'gw') 
        dict_c['calc_gw'] = 0 
        plt.plot(en, ftot_a, label="ftot_gw "+str(minkpt)+" "+str(maxkpt)+" "+str(minband)+" "+str(maxband))
        if extinf == 1: 
            str_exi = "_extinf" 
        else: 
            str_exi = ""
        outname = "sftot_gw"+"_kpt_"+str(minkpt)+"_"+str(maxkpt)+"_bd_"+str(minband)+"_"+str(maxband)+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev"+".dat"
        print(" Output file: ", outname)
        with open(outname,'w') as fout:
            fout.write("# kpt "+str(minkpt)+" "+str(maxkpt)+"\n") 
            fout.write("# bd  "+str(minband)+" "+str(maxband)+"\n")
           #fout.write("# minkpt maxkpt %8.4f %12.8f\n" % (minkpt, maxkpt))
           #fout.write("# minband maxband %8.4f %12.8f\n" % (minband, maxband))
            for i in range(en.size): 
                fout.write("%8.4f %12.8f\n" % (en[i], ftot_a[i]))
    if int(dict_c['calc_exp']) == 1:
        print("### EXP SPECTRAL FUNCTIONS ###")
        enexp, ftot_b = read_sf(dict_c, pdos, 'exp') 
        dict_c['calc_exp'] = 0
        plt.plot(enexp, ftot_b, label="ftot_exp "+str(minkpt)+" "+str(maxkpt)+" "+str(minband)+" "+str(maxband))
        if extinf == 1: 
            str_exi = "_extinf" 
        else: 
            str_exi = ""
        outname = "sftot_exp"+"_kpt_"+str(minkpt)+"_"+str(maxkpt)+"_bd_"+str(minband)+"_"+str(maxband)+"_"+str(penergy)+"ev_np"+str(npoles)+str_exi+".dat"
        print(" Output file: ", outname)
        with open(outname,'w') as fout:
            fout.write("# kpt "+str(minkpt)+" "+str(maxkpt)+"\n") 
            fout.write("# bd  "+str(minband)+" "+str(maxband)+"\n")
            for i in range(enexp.size): 
                fout.write("%8.4f %12.8f\n" % (enexp[i], ftot_b[i]))
    ### ==== WRITING OUT GW SPECTRAL FUNCTION === ###
    #newen = newen-efermi
# ---------------------------------- #

else:
#if int(dict_c['restart']) == 0:  
    newen, spftot, allkb = calc_sf_gw(dict_c,hartree,pdos,en,res,ims)
    print("="*40)
    print("MEMORY USAGE TEST")
    print("="*40)
    dummy = 0
    for item in allkb:
        dummy += item.nbytes
    print(" GB size: ", float((newen.nbytes + spftot.nbytes + dummy)) / 1000000000)
    print("="*40)
    if int(dict_c['calc_gw']) == 1:
        print()
        print(" ### Writing out A(\omega)_GW...  ")
        # Start a new thread
        #write_spfkb(invar_dict,newen,allkb)
        thread = Thread(target = write_spfkb, args = (dict_c, newen, allkb))
        thread.start()
        outname = "spftot_gw"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev"+".dat"
        #outfile = open(outname,'w')
        with open(outname,'w') as outfile: 
            for i in xrange(np.size(newen)): 
                outfile.write("%7.4f %15.10e\n"% (newen[i], spftot[i])) # Dump string representations of arrays
        #outfile.close()
        print(" A(\omega)_GW written in", outname)
        plt.plot(newen,spftot,label="ftot_gw");
    
   #print(" --- Time spent for GW: {} seconds. ---".format(time.time() - t_part1))
   #print(" --- Time spent so far: {} seconds. ---".format(time.time() - start_time))
    # ============================= ###
    
    ### ===================================== ###
    ### ===== EXPONENTIAL SPECTRAL FUNCTION ====== ###
    if int(dict_c['calc_exp']) == 1:
        # Time section
        #import time
       #e0=time.time()
       #c0=time.clock()
       #elaps1=time.time() - e0
       #cpu1=time.clock() - c0
       #print(str(" Starting time (elaps, cpu): %10.6e %10.6e"% (elaps1, cpu1)))
        print()
        print(" ### Calculation of exponential A ### ")
        ### ==== Finding zero in res --> Eqp ===== ###
        print(" Finding zeros in real parts...")
        eqp, imeqp = calc_eqp_imeqp(en, res, ims, hartree, 0)
        print(" Test imeqp:\n", imeqp)
        # Writing out eqp
        # Writing out imeqp
        thread = Thread(target = write_eqp_imeqp, args = (eqp, imeqp))
        thread.start()
        dict_c['origdir'] = origdir
        enexp, ftot, sfkb = calc_sf_c(dict_c, hartree, pdos, eqp, imeqp, newen, allkb)
        if int(dict_c['calc_crc']) == 1:
            print()
            print(" ### Calculation of constrained retarded cumulant ### ")
            B_crc_kb = calc_B_crc(dict_c, eqp, newen, allkb)
            sftot_crc, sfkb_crc = calc_sf_crc(dict_c, B_crc_kb, hartree, newen, allkb)
       ## TEST!!!!
       #dict_c['npoles'] = 1  
       #enexp1, ftot1, sfkb = calc_sf_c_para(dict_c, hartree, pdos, eqp, imeqp, newen, allkb)
       #enexp2, ftot2, sfkb = calc_sf_c_serial(dict_c, hartree, pdos, eqp, imeqp, newen, allkb)
       #plt.plot(enexp1, ftot1, label="ftot_para")
       #plt.plot(enexp2, ftot2, label="ftot_serial")
       #plt.legend(loc=2)
       #plt.show()
        print("Integral test, spftot: ", np.trapz(ftot,enexp))
        print("="*40)
        print("MEMORY USAGE TEST")
        print("="*40)
        print(" GB size: ", float((enexp.nbytes+ftot.nbytes+sfkb.nbytes))/1000000000)
        # Writing out sfkb
        print("TEST: sfkb.shape:", sfkb.shape)
        thread = Thread(target = write_sfkb_c, args = (invar_dict, enexp, sfkb))
        thread.start()
    
    
        ### TODO: fix all below ###
        mylabel = "ftot np "+str(npoles)+" mnb "+str(minband)+" mxb "+str(maxband)
        plt.plot(enexp, ftot, label = mylabel)
    
# Now go back to original directory
print(" Moving back to parent directory:\n", origdir)
chdir(newdir)
end_time = time.time()
print()
t2 = t_pregw - start_time
t3 = end_time - start_time
#print(t2, t2/60, t2/3600)
#print(" --- Time spent for GW: {} seconds. ---".format(time.time() - t_part1))
#print(" --- Time spent so far: {} seconds. ---".format(time.time() - start_time))
#title = 'Spectral function '+ 'A (' + r'$\omega $ ' + ') - '+r'$ h\nu = $'+str(penergy)+' eV'
#plt.title(title)
plt.legend(loc=2)
plt.show()
t4 = t3
t5 = t3
if int(t3/3600) >= 1: 
    t4 = t3 - int(t3/3600)*3600
if int(t4/60) >= 1: 
    t5 = t4 - int(t4/60)*60
print(" Calculation lasted "+str(int(t3/3600))+" Hours, "+str(int(t4/60))+" Minutes and "+str(int(t5))+" Seconds")
print(" "+30*"-")
print(" End of program reached.")
print(" "+30*"-")
