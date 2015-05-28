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
#print(('%12s, %9.4f' % invar_dict.keys(), invar_dict.values()))
print(" "+"===="+" Input variables "+"====")
print()
for key in  invar_dict.keys():
    #print(key, invar_dict[key])
    print(('%12s :: %12s' % (key, invar_dict[key])))
for i in range(52): print('=',end='')
print()
#print(52*"=")
print((" SF :: END"))
for i in range(52): print('=',end='')
print()
#print(52*"=")
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

### ===================================== ###
### ===== EXPONENTIAL SPECTRAL FUNCTION ====== ###
if int(invar_dict['calc_exp']) == 1:
    # Time section
    #import time
    e0=time.time()
    c0=time.clock()
    elaps1=time.time() - e0
    cpu1=time.clock() - c0
    print(str(" Starting time (elaps, cpu): %10.6e %10.6e"% (elaps1, cpu1)))
    print(" ### Calculation of exponential A...  ")
    ### ==== Finding zero in res --> Eqp ===== ###
    print(" Finding zeros in real parts...")
    eqp, imeqp = calc_eqp_imeqp(nkpt,nband,en,res,ims,hartree,0)
    print(" Test imeqp:", imeqp)
    # Writing out eqp
    # Writing out imeqp
    thread = Thread(target = write_eqp_imeqp, args = (eqp,imeqp))
    thread.start()
    npoles = int(invar_dict['npoles'])
    if npoles==999:
        omegampole = np.ones((nkpt,nband))*omega_p
        ampole =  np.zeros((nkpt,nband))
        for ik in xrange(nkpt):
            for ib in xrange(nband):
                print(" ik, ib", ik, ib)
                #interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
                #if eqp[ik,ib]<=efermi:
                if eqp[ik,ib]<=0:
                    tmpen = en[ims[ik,ib]>=0]
                    tmpim = ims[ik,ib,ims[ik,ib]>=0]
                else:
                    tmpen = en[ims[ik,ib]<0]
                    tmpim = ims[ik,ib,ims[ik,ib]<0]
                ampole[ik,ib] = abs(np.trapz(tmpim,tmpen))/np.pi
                print(" 1/pi*\int\Sigma   =", ampole[ik,ib])
                # Workaround correction for small energy plasmons
                ampole[ik,ib] = ampole[ik,ib]/(abs(tmpen[-1]-tmpen[0]))*omega_p
#                # Workaround for small energy plasmons
#                if eqp[ik,ib]<=efermi:
#                    tmpim = tmpim[tmpen>=eqp[ik,ib]-2.5]
#                    tmpen = tmpen[tmpen>=eqp[ik,ib]-2.5]
#                else:
#                    tmpim = tmpim[tmpen <eqp[ik,ib]+2.5]
#                    tmpen = tmpen[tmpen <eqp[ik,ib]+2.5]
#                ampole[ik,ib] = np.trapz(tmpim,tmpen)/np.pi
        #ampole = ampole/omega_p**2
                #ampole[ik,ib] = np.trapz(en[ims[ik,ib]>=0],ims[ik,ib,ims[ik,ib]>=0])/np.pi
    elif npoles != 0:
        from multipole import fit_multipole, getdata_file #, write_f_as_sum_of_poles
        print(" ### ================== ###")
        print(" ###    Multipole fit   ###")
        print(" Number of poles:", npoles)
        omegampole =  np.zeros((nkpt,nband,npoles))
        ampole =  np.zeros((nkpt,nband,npoles))
        for ik in xrange(nkpt):
            ikeff=minkpt+ik-1
            for ib in xrange(nband):
                ibeff=minband+ib-1
                print(" ik, ib", ik, ib)
                interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
                # Here we take the curve starting from eqp and then we invert it
                # so as to have it defined on the positive x axis
                # and so that the positive direction is in the 
                # increasing direction of the array index
                #if eqp[ik,ib] <= efermi:
                if eqp[ik,ib] <= 0:
                    en3 = en[en<=eqp[ik,ib]] # So as to avoid negative omegampole
                else:
                    en3 = en[en>eqp[ik,ib]] # So as to avoid negative omegampole
                #en3 = en[en<=efermi]
                im3 = abs(interpims(en3)/np.pi) # This is what should be fitted
                en3 = en3 - eqp[ik,ib]
                #en3 = en3 - efermi
                #if eqp[ik,ib] <= efermi: 
                if eqp[ik,ib] <= 0:
                    en3 = -en3[::-1] 
                    im3 = im3[::-1]
                omegai, gi, deltai = fit_multipole(en3,im3,npoles,0)
                #if np.isnan(omegai): sys.exit(1)
                #omegampole[ik,ib] = omegai + eqp[ik,ib] - efermi
                omegampole[ik,ib] = omegai 
                ampole[ik,ib] = gi/(omegampole[ik,ib])**2 
                print(" Integral test. Compare \int\Sigma and \sum_j^N\lambda_j.")
                print(" 1/pi*\int\Sigma   =", np.trapz(im3,en3))
                print(" \sum_j^N\lambda_j =", np.sum(gi))
                #plt.plot(en3,im3,"-"); plt.plot(omegai,np.pi/2*gi*omegai/deltai,"-o")
                #plt.show(); sys.exit(1)
                #e1,f1 = write_f_as_sum_of_poles(en3,omegai,gi,deltai,0)
        # Writing out a_j e omega_j
        print(" ### Writing out a_j and omega_j...")
        outname = "a_j_np"+str(npoles)+".dat"
        outfile = open(outname,'w')
        outname = "omega_j_np"+str(npoles)+".dat"
        outfile2 = open(outname,'w')
        for ipole in xrange(npoles):
            for ik in xrange(nkpt):
                for ib in xrange(nband):
                    outfile.write("%10.5f"  % (ampole[ik,ib,ipole]))
                    outfile2.write("%10.5f" % (omegampole[ik,ib,ipole]))
                outfile.write("\n")
                outfile2.write("\n")
            outfile.write("\n")
            outfile2.write("\n")
        outfile.close()
        outfile2.close()
        # Extrinsic and interference contribution
        if extinf == 1:
            extinfname = "a_wp."+str(penergy)
            amp_exinf, w_extinf = calc_extinf_corrections(origdir,extinfname,ampole,omegampole)
            print(" ### Writing out a_j_extinf...")
            outname = "a_j_np"+str(npoles)+"_extinf."+str(penergy)
            outfile = open(outname,'w')
            for ipole in xrange(npoles):
                for ik in xrange(nkpt):
                    for ib in xrange(nband):
                        outfile.write("%10.5f"  % (amp_exinf[ik,ib,ipole]))
                    outfile.write("\n")
                outfile.write("\n")
            outfile.close()
    else:
        omegampole =  np.zeros((nkpt,nband))
        ampole =  np.zeros((nkpt,nband))
    elaps2 = time.time() - elaps1 - e0
    cpu2 = time.clock() - cpu1 - c0
    #print(elaps2, cpu2)
    print(str(" Used time (elaps, cpu): %10.6e %10.6e"% (elaps2, cpu2)))
    print(" Calculating multipole exponential A...")
    dxexp=0.005 
    enexp=np.arange(enmin,enmax,dxexp)
    nenexp=np.size(enexp)
    ftot=np.zeros((nenexp))
    f=np.zeros((nkpt,nband,nenexp))
    ftot=np.zeros((np.size(enexp)),order='Fortran')
    nen = np.size(enexp)
    # With extrinsic effects
    if extinf == 1:
        from extmod_spf_mpole import f2py_calc_spf_mpole_extinf
        for ik in xrange(nkpt):
            ikeff=minkpt+ik-1
            for ib in xrange(nband):
                ibeff=minband+ib-1
                print(" ik, ib", ik, ib)
                prefac=np.exp(-np.sum(amp_exinf[ik,ib]))/np.pi*wtk[ikeff]*pdos[ib]*abs(imeqp[ik,ib])
                akb=amp_exinf[ik,ib] # This is a numpy array (slice)
                omegakb=omegampole[ik,ib] # This is a numpy array (slice)
                wkb=w_extinf[ik,ib] # This is a numpy array (slice)
                eqpkb=eqp[ik,ib]
                imkb=imeqp[ik,ib] # + w_extinf[ik,ib]/2 # extinf width added
                #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles,wkb)
                #ftot += tmpf
                tmpf = np.zeros((nen), order='Fortran')
                tmpf = f2py_calc_spf_mpole_extinf(tmpf,enexp,prefac,akb,omegakb,wkb,eqpkb,imkb) #,np.size(enexp),npoles)
                outnamekb = "spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"_extinf."+str(penergy)
                outfilekb = open(outnamekb,'w')
                for ien in xrange(nenexp):
                    outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], tmpf[ien]))
                outfilekb.close()
                ftot = ftot + tmpf
    else: # extinf == 0
        from extmod_spf_mpole import f2py_calc_spf_mpole
        for ik in xrange(nkpt):
            ikeff=minkpt+ik-1
            for ib in xrange(nband):
                ibeff=minband+ib-1
                print(" ik, ib, ikeff, ibeff", ik, ib, ikeff+1, ibeff+1)
                prefac=np.exp(-np.sum(ampole[ik,ib]))/np.pi*wtk[ikeff]*pdos[ib]*abs(imeqp[ik,ib])
                akb=ampole[ik,ib] # This is a numpy array (slice)
                omegakb=omegampole[ik,ib] # This is a numpy array (slice)
                eqpkb=eqp[ik,ib]
                imkb=imeqp[ik,ib]
                #tmpf1 = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                #print(nen, np.size(enexp))
                #tmpf = 0.0*tmpf
                if eqpkb < 0.0:
                    tmpf = np.zeros((nen), order='Fortran')
                    tmpf = f2py_calc_spf_mpole(tmpf,enexp,prefac,akb,omegakb,eqpkb,imkb) #,nen,npoles)
                    #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                else:
                    print(" This state is empty! eqpkb ik ib:",eqpkb, ikeff+1, ibeff+1)
                    #print("omegakb", omegakb)
                    omegakb=-omegakb
                    #print("-omegakb", omegakb)
                    tmpf = np.zeros((nen), order='Fortran')
                    tmpf = f2py_calc_spf_mpole(tmpf,enexp,prefac,akb,omegakb,eqpkb,imkb) #,nen,npoles)
                    #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                #if not tmpf[0]>=0: print("ik,ib,prefac,akb,omegakb,eqpkb,imkb,npoles:",ik,ib,prefac,akb,omegakb,eqpkb,imkb,npoles; sys.exit(1))
                outnamekb = "spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"."+str(penergy)
                outfilekb = open(outnamekb,'w')
                for ien in xrange(nenexp):
                    outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], tmpf[ien]))
                outfilekb.close()
                ftot = ftot + tmpf
                #print(ftot[0], tmpf[0])
    elaps2 = time.time() - elaps1 - e0
    cpu2 = time.clock() - cpu1 - c0
    #print(elaps2, cpu2)
    print(str(" Used time (elaps, cpu): %10.6e %10.6e"% (elaps2, cpu2)))
    print(" ### Writing out A(\omega)_exp...  ")
    #enexp = enexp-efermi
    if extinf == 1:
        outname = "spftot_exp"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev_np"+str(npoles)+"_extinf.dat"
        outfile = open(outname,'w')
        for i in xrange(nenexp):
            outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot[i])) # Dump string representations of arrays
        outfile.close()
    else: # extinf == 0
        outname = "spftot_exp"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev_np"+str(npoles)+".dat"
        outfile = open(outname,'w')
        for i in xrange(nenexp):
            outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot[i])) # Dump string representations of arrays
        outfile.close()
    print(" A(\omega)_exp written in", outname)
    plt.plot(enexp,ftot,label="ftot");
# Now go back to original directory
print(" Moving back to parent directory:\n", origdir)
chdir(newdir)
#title = 'Spectral function '+ 'A (' + r'$\omega $ ' + ') - '+r'$ h\nu = $'+str(penergy)+' eV'
#plt.title(title)
plt.legend(loc=2)
plt.show()
