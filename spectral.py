#!/usr/bin/env python
"""
Written by Matteo Guzzo.
This version of the script includes extrinsic and interference effects
which are taken from the file (a_wp.dat) given by Josh. 
An external fortran module is used to calculate the 
exponential spectral function. 
It has not been tested a lot but it seems to work. 
List of files needed:
- invar.in with input variables.
- _SIG for the self-energy.
- s.dat, p_even.dat, p_odd.dat, d_even.dat, etc. 
for the orbital character and symmetries.
- cs*.dat for the photon cross sections.
- hartree.dat or elda.dat and vxc.dat for the hartree energies.
- wtk.dat for the k-points weights.
- a_wp.dat for the extrinsic/interference effects and additional lifetime.
"""
from spectral_modules import *
import numpy as np;
import matplotlib.pylab as plt;
plt.figure(1)
#from scipy.interpolate import interp1d
#from scipy import optimize
import sys
from os.path import isfile, join, isdir
from os import getcwd, pardir, mkdir, chdir
#

### ============================= ###
###  ==  PROGRAM BEGINS HERE  ==  ###
### ============================= ###

# ======== READING INPUT VARIABLES ======= #
print " Reading invar file... ",
invar = {}
if len(sys.argv)>1: 
    infname = sys.argv[1]
else: 
    infname = "invar.in"
if isfile(infname):
    infile = open(infname)
    for line in infile.readlines():
        word = line.split()
        invar[word[-1]] = word[0];
#        print "invar: ", invar
    infile.close()
    if 'sigmafile' in invar: 
        sigfilename  =       invar['sigmafile'];
    else:
        sigfilename  =       "default_SIG";
    if 'minband' in invar:
        minband      =   int(invar['minband']);
    else:
        minband      =   1
    if 'maxband' in invar:
        maxband      =   int(invar['maxband']);
    else:
        maxband      =   1
    if 'minkpt' in invar:
        minkpt       =   int(invar['minkpt']);
    else:
        minkpt = 1
    if 'maxkpt' in invar:
        maxkpt       =   int(invar['maxkpt']);
    else:
        maxkpt = 1
    if 'enmin' in invar:
        enmin        = float(invar['enmin']);
    else:
        enmin = -20.0
    if 'enmax' in invar:
        enmax        = float(invar['enmax']);
    else:
        enmax = 20.0
    if 'sfactor' in invar:
        sfac         = float(invar['sfactor']);
    else:
        sfac = 1.0
    if 'pfactor' in invar:
        pfac         = float(invar['pfactor']);
    else:
        pfac = 1.0
    if 'penergy' in invar:
        penergy      =   int(invar['penergy']);
    else:
        penergy = 0
    if 'npoles' in invar:
        npoles       =   int(invar['npoles']);
    else:
        npoles = 1
    if 'calc_gw' in invar:
        flag_calc_gw =   int(invar['calc_gw']);
    else:
        flag_calc_gw = 1
    if 'calc_exp' in invar:
        flag_calc_exp =   int(invar['calc_exp']);
    else:
        flag_calc_exp = 0
    if 'extinf' in invar:
        extinf       = float(invar['extinf']);
    else:
        extinf = 0
    if 'efermi' in invar:
        efermi       = float(invar['efermi']);
    else:
        efermi = 0.0
    if 'nkpt' in invar: 
        nkpt         =   int(invar['nkpt']);
        if nkpt != maxkpt-minkpt+1: 
            print
            print " WARNING: nkpt not in accordance with minkpt and maxkpt. Please check!"
    if 'omega_p' in invar: 
        omega_p = float(invar['omega_p']);
    else:
        omega_p = 0.0
    if 'enhartree' in invar: 
        enhartree = float(invar['enhartree']);
    else:
        enhartree = None
else : 
    print "Invar file not found ('"+str(infname)+"'). Impossible to continue."
    sys.exit(1)
print "Done."
print " minband =", minband;
print " maxband =", maxband;
nband = 1 + maxband - minband;
print " nband =", nband;
print " minkpt =", minkpt;
print " maxkpt =", maxkpt;
nkpt = 1 + maxkpt - minkpt;
print " nkpt =", nkpt;
print " enmin =", enmin;
print " enmax =", enmax;
print " S prefactor:", sfac
print " P prefactor:", pfac
# Max energy in spectrum
# TODO: write a function to read this parameter from a file (Fermi energy?)
#enmax = 15. # eV
# ====== READING HARTREE ===== #
hartree = read_hartree()
#hartree = hartree # - efermi
# ======== READING WTK ======= #
wtk = read_wtk()
# ======== READING _SIG FILE ======= #
#en, res, ims = read_sigfile(nkpt,nband,sigfilename)
print " enmin, enmax"
print  enmin, enmax
enmit = enmin+efermi
enmat= enmax+efermi
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
### ===================================================== ###
print " # ------------------------------------------------ # ";
# Here we move to a subdirectory to avoid flooding-up the current directory
newdirname = "Spfunctions"
origdir = getcwd() # remember where we are
newdir = join(origdir, newdirname) # Complete path of the new directory
print " Moving into output directory:\n ", newdir
if not isdir(newdir) :
    mkdir(newdir)
chdir(newdir)
### ================================= ###
### ===== GW SPECTRAL FUNCTION ====== ###
# GW spectral function part
if flag_calc_gw == 1:
    newen, spftot = calc_spf_gw(minkpt,maxkpt,minband,maxband,wtk,pdos,en,enmin,enmax,res,ims,hartree)
    #calc_spf_gw(nkpt,nband,wtk,pdos,en,res,ims)
        ### ==== WRITING OUT GW SPECTRAL FUNCTION === ###
    #newen = newen-efermi
    print " ### Writing out A(\omega)_GW...  "
    outname = "spftot_gw"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev"+".dat"
    outfile = open(outname,'w')
    for i in xrange(np.size(newen)):
        outfile.write("%7.4f %15.10e\n"% (newen[i],spftot[i])) # Dump string representations of arrays
    outfile.close()
    print " A(\omega)_GW written in", outname
    plt.plot(newen,spftot,label="ftot_gw");
    
# ============================= ###

### ===================================== ###
### ===== EXPONENTIAL SPECTRAL FUNCTION ====== ###
if flag_calc_exp == 1:
    # Time section
    import time
    e0=time.time()
    c0=time.clock()
    elaps1=time.time() - e0
    cpu1=time.clock() - c0
    print str(" Starting time (elaps, cpu): %10.6e %10.6e"% (elaps1, cpu1))
    print " ### Calculation of exponential A...  "
    ### ==== Finding zero in res --> Eqp ===== ###
    print " Finding zeros in real parts..."
    eqp, imeqp = calc_eqp_imeqp(nkpt,nband,en,res,ims,hartree,0,minband)
    print " Test imeqp:", imeqp
    # Writing out eqp
    # Writing out imeqp
    outname = "eqp.dat"
    outfile2 = open(outname,'w')
    outname = "imeqp.dat"
    outfile3 = open(outname,'w')
    for ik in xrange(nkpt):
        for ib in xrange(nband):
            outfile2.write("%14.5f" % (eqp[ik,ib]))
            outfile3.write("%14.5f" % (imeqp[ik,ib]))
        outfile2.write("\n")
        outfile3.write("\n")
    outfile2.close()
    outfile3.close()
    if npoles==999:
        omegampole = np.ones((nkpt,nband))*omega_p
        ampole =  np.zeros((nkpt,nband))
        for ik in xrange(nkpt):
            for ib in xrange(nband):
                print " ik, ib", ik, ib
                #interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
                #if eqp[ik,ib]<=efermi:
                if eqp[ik,ib]<=0:
                    tmpen = en[ims[ik,ib]>=0]
                    tmpim = ims[ik,ib,ims[ik,ib]>=0]
                else:
                    tmpen = en[ims[ik,ib]<0]
                    tmpim = ims[ik,ib,ims[ik,ib]<0]
                ampole[ik,ib] = abs(np.trapz(tmpim,tmpen))/np.pi
                print " 1/pi*\int\Sigma   =", ampole[ik,ib]
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
        print " ### ================== ###"
        print " ###    Multipole fit   ###"
        print " Number of poles:", npoles
        omegampole =  np.zeros((nkpt,nband,npoles))
        ampole =  np.zeros((nkpt,nband,npoles))
        for ik in xrange(nkpt):
            ikeff=minkpt+ik-1
            for ib in xrange(nband):
                ibeff=minband+ib-1
                print " ik, ib", ik, ib
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
                print " Integral test. Compare \int\Sigma and \sum_j^N\lambda_j."
                print " 1/pi*\int\Sigma   =", np.trapz(im3,en3)
                print " \sum_j^N\lambda_j =", np.sum(gi)
                #plt.plot(en3,im3,"-"); plt.plot(omegai,np.pi/2*gi*omegai/deltai,"-o")
                #plt.show(); sys.exit(1)
                #e1,f1 = write_f_as_sum_of_poles(en3,omegai,gi,deltai,0)
        # Writing out a_j e omega_j
        print " ### Writing out a_j and omega_j..."
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
            #print nkpt, nband
            #print ampole.shape, ampole[:,0].size, ampole[0,:].size
            amp_exinf, w_extinf = calc_extinf_corrections(origdir,extinfname,ampole,omegampole)
            print " ### Writing out a_j_extinf..."
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
    #print elaps2, cpu2
    print str(" Used time (elaps, cpu): %10.6e %10.6e"% (elaps2, cpu2))
    print " Calculating multipole exponential A..."
    dxexp=0.005 
    enexp=np.arange(enmin,enmax,dxexp)
    nenexp=np.size(enexp)
    ftot=np.zeros((nenexp))
    fkb = np.zeros((nkpt,nband,nenexp))
   #f=np.zeros((nkpt,nband,nenexp))
    ftot=np.zeros((np.size(enexp)),order='Fortran')
    nen = np.size(enexp)
    # With extrinsic effects
    if extinf == 1:
        from extmod_spf_mpole import f2py_calc_spf_mpole_extinf
        for ik in xrange(nkpt):
            ikeff=minkpt+ik-1
            for ib in xrange(nband):
                ibeff=minband+ib-1
                print " ik, ib", ik, ib
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
                print " ik, ib, ikeff, ibeff", ik, ib, ikeff+1, ibeff+1
                prefac=np.exp(-np.sum(ampole[ik,ib]))/np.pi*wtk[ikeff]*pdos[ib]*abs(imeqp[ik,ib])
                print
                print "\n === Normalization test === " 
                print " Prefactor:", np.exp(-np.sum(ampole[ik,ib])) 
                print " Exponent:", np.sum(ampole[ik,ib]) 
                print " Exponent/npoles:", np.sum(ampole[ik,ib])/npoles
                print
                print
                akb=ampole[ik,ib] # This is a numpy array (slice)
                omegakb=omegampole[ik,ib] # This is a numpy array (slice)
                eqpkb=eqp[ik,ib]
                imkb=imeqp[ik,ib]
                #tmpf1 = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                #print nen, np.size(enexp)
                #tmpf = 0.0*tmpf
                if eqpkb < 0.0:
                    tmpf = np.zeros((nen), order='Fortran')
                    tmpf = f2py_calc_spf_mpole(tmpf,enexp,prefac,akb,omegakb,eqpkb,imkb) #,nen,npoles)
                    #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                else:
                    print " This state is empty! eqpkb ik ib:",eqpkb, ikeff+1, ibeff+1
                    #print "omegakb", omegakb
                    omegakb=-omegakb
                    #print "-omegakb", omegakb
                    tmpf = np.zeros((nen), order='Fortran')
                    tmpf = f2py_calc_spf_mpole(tmpf,enexp,prefac,akb,omegakb,eqpkb,imkb) #,nen,npoles)
                    #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                #if not tmpf[0]>=0: print "ik,ib,prefac,akb,omegakb,eqpkb,imkb,npoles:",ik,ib,prefac,akb,omegakb,eqpkb,imkb,npoles; sys.exit(1)
                outnamekb = "spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"."+str(penergy)
                outfilekb = open(outnamekb,'w')
                for ien in xrange(nenexp):
                    #outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], tmpf[ien]))
                    outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], tmpf[ien]))
                outfilekb.close()
                fkb[ik,ib] = tmpf
                ftot = ftot + tmpf  
                #print ftot[0], tmpf[0]
    print
    print(" ### Calculation of constrained retarded cumulant ### ")
   #B_crc_kb = calc_B_crc(dict_c, eqp, newen, allkb)
   #    imskb = allkb[3]
   #B_crc_kb =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,npoles))
    B_crc_kb =  np.zeros((nkpt,nband))
    from multipole import fit_multipole, getdata_file #, write_f_as_sum_of_poles
    print " ### ================== ###" 
    print " ###    Multipole fit   ###" 
    print " Number of poles:", npoles 
    print " Number of poles_crc is FIXED TO 1"
    npoles_crc = 1
   #omegampole =  np.zeros((nkpt,nband,npoles))
   #ampole =  np.zeros((nkpt,nband,npoles))
    omegampole_crc =  np.zeros((nkpt,nband,npoles_crc))
    ampole_crc     =  np.zeros((nkpt,nband,npoles_crc))
    #for ik in range(nkpt):
    #    ikeff=minkpt+ik-1
    #bdrange = vardct['bdrange']
    #kptrange = vardct['kptrange']
    #print "kptrange, bdrange ", kptrange, bdrange 
   #for ik in kptrange:
   #    for ib in bdrange:
    for ik in xrange(nkpt):
        ikeff=minkpt+ik-1
        for ib in xrange(nband):
            ibeff=minband+ib-1
   #for ik in range(imskb[:,0,0].size):
        #for ib in range(nband):
   #    for ib in range(imskb[0,:,0].size):
           #if eqp[ik,ib] < newen[npoles]:
           ##if eqp[ik,ib] > newen[-1]:
           #    omegampole[ik,ib] = omegampole[ik,ib-1]
           #    ampole[ik,ib] = ampole[ik,ib-1]
           #    print " Eqp beyond available energy range. Values from lower band are taken." 
           #    continue
           #else:
            print " ik, ib", ik, ib 
            #interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
            #print newen.shape, imskb.shape 
           #print ims[ik,ib], enexp
           #sys.exit()
            interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
            # Here we take the curve starting from efermi and then we invert it
            # so as to have it defined on the positive x axis
            # and so that the positive direction is in the 
            # increasing direction of the array index
            #if eqp[ik,ib] <= efermi:
            if eqp[ik,ib] <= 0:
                #en3 = en[en<=eqp[ik,ib]] # So as to avoid negative omegampole
                en3 = enexp[enexp >= 0] # So as to avoid negative omegampole and ambiguity btween eqp and efermi
            else:
                en3 = enexp[enexp<0] # So as to avoid negative omegampole
                #en3 = en[en>eqp[ik,ib]] # So as to avoid negative omegampole
            #en3 = en[en<=efermi]
            if en3.size == 0:
                print  
                print " WARNING: QP energy is outside of given energy range!\n"+\
                        " This state will be skipped!\n"+\
                        "You might want to modify enmin/enmax."
                print " eqp[ik,ib], newen[-1]", eqp[ik,ib] , enexp[-1] 
                continue
            im3 = abs(interpims(en3)/np.pi) # This is what should be fitted
            en3 = en3 - eqp[ik,ib]
            if eqp[ik,ib] > 0:
                en3 = -en3[::-1] 
                im3 = im3[::-1]
           #if eqp[ik,ib] <= 0:
           #    en3 = -en3[::-1] 
           #    im3 = im3[::-1]
           #### TESTING ###
           #print "ik, ib, eqp[ik,ib], en3[0], en3[-1], newen[0], newen[-1]:\n", ik, ib, eqp[ik,ib], en3[0], en3[-1], newen[0], newen[-1] 
           #import matplotlib.pylab as plt
           #plt.plot(newen, imskb[ik,ib]/np.pi,"-")
           #plt.plot(en3+eqp[ik,ib], im3,"x")
           #plt.plot(en3, im3,"o")
           #plt.show()
           #sys.exit()
           #### END TESTING ###
            omegai, lambdai, deltai = fit_multipole(en3,im3,npoles_crc)
            # HERE WE MUST CHECK THAT THE NUMBER OF POLES 
            # IS NOT BIGGER THAN THE NUMBER OF POINTS THAT HAS TO BE FITTED
            if npoles_crc > omegai.size:
                omegampole_crc[ik,ib][:omegai.size] = omegai 
                ampole_crc[ik,ib][:omegai.size] = np.true_divide(lambdai,(np.square(omegai)))
                print  
                print " WARNING: npoles used ("+str(npoles_crc)+") is larger"+\
                        " than poles x data array can give ("+str(omegai.size)+")."
               #print "WARNING: Reduce npoles. You are wasting resources!!!" 
                print " Im(Sigma) will be interpolated to obtain the desired number of poles." 
                current_size = omegai.size
                counter = 0
                while npoles_crc > current_size:
                    counter += 1
                    print  
                    print " WARNING: Arrays are too coarse." 
                    print " npoles, omegai.size:", npoles_crc, omegai.size 
                    print " Filling arrays with interpolated values..." 
                    en1 = array_doublefill(en3)
                    im1 = array_doublefill(im3)
                    en3 = en1
                    im3 = im1
                    omegai, lambdai, deltai = fit_multipole(en1,im1,npoles_crc)
                    current_size = omegai.size
                    if counter > 4:
                        print 60*"=" 
                        print " WARNING: You are trying too hard with too few points." 
                        print " The array has been interpolated more than 4 times." 
                        print " Maybe use less poles or calculate more points for Sigma?" 
                        print 60*"=" 
        #   im1 = fit_double(im3)
            else:
                omegampole_crc[ik,ib] = omegai 
                ampole_crc[ik,ib] = np.true_divide(lambdai,(np.square(omegai)))
            B_crc_kb[ik,ib] = np.sum(ampole_crc[ik,ib])
            ### Test plot for the fit ###
           #from multipole import write_f_as_sum_of_poles
          ##import matplotlib.pylab as plt
           #import pylab
           #plt.figure(2)
           #eta = 0.5
           #enlor, flor = write_f_as_sum_of_poles(en3, omegai, lambdai, deltai, eta)
           #plt.plot(enlor, flor,"-",label="sum of poles, eta: "+str(eta))
           #plt.plot(en3,im3,"-",label="ImS(e-w)")
           #plt.plot(omegai,lambdai,"go", label = "omegai, lambdai")
           #plt.plot(omegai,lambdai/deltai,"ro", label = "omegai, lambdai/deltai")
           #plt.title("ik: "+str(ik)+", ib: "+str(ib)+", npoles: "+str(npoles))
           #plt.legend()
           #pylab.savefig('imS_fit_np'+str(npoles)+'_ik'+str(ik)+'_ib'+str(ib)+'_CRC.pdf')
           #plt.close(2)
           #### END - Test plot for the fit ###
            from multipole import write_f_as_sum_of_poles
           #ampole[ik,ib] = gi
            print
            print " Integral test. Compare \int\Sigma and \sum_j^N\lambda_j." 
            print " 1/pi*\int\Sigma   =", np.trapz(im3,en3) 
            print " \sum_j^N\lambda_j =", np.sum(lambdai) 
            print " b_j:", ampole_crc[ik,ib] 
            #plt.plot(en3,im3,"-"); plt.plot(omegai,np.pi/2*gi*omegai/deltai,"-o")
            #e1,f1 = write_f_as_sum_of_poles(en3,omegai,gi,deltai,0)
    # Writing out a_j e omega_j
    print
    print " ### Writing out a_j and omega_j..." 
    outname = "a_j_np"+str(npoles_crc)+"_crc.dat"
    outfile = open(outname,'w')
    outname = "omega_j_np"+str(npoles_crc)+"_crc.dat"
    outfile2 = open(outname,'w')
    for ipole in xrange(npoles_crc):
 #      for ik in kptrange:
 #          #for ib in range(nband):
 #          for ib in bdrange:
        for ik in range(ims[:,0,0].size):
            for ib in range(ims[0,:,0].size):
                outfile.write("%10.5f"  % (ampole_crc[ik,ib,ipole]))
                outfile2.write("%10.5f" % (omegampole_crc[ik,ib,ipole]))
            outfile.write("\n")
            outfile2.write("\n")
        outfile.write("\n")
        outfile2.write("\n")
    outfile.close()
    outfile2.close()
   #return B_crc_kb
    ### Calculation of the CRC spectral function ###
   #sftot_crc, sfkb_crc = calc_sf_crc(dict_c, B_crc_kb, enexp, allkb)
   #ftot=np.zeros((nenexp))
   #f=np.zeros((nkpt,nband,nenexp))
    ftot_crc_occ=np.zeros((np.size(enexp)),order='Fortran')
    ftot2=np.zeros((np.size(enexp)),order='Fortran')
    from extmod_spf_mpole import f2py_calc_crc_mpole
    for ik in xrange(nkpt):
        ikeff=minkpt+ik-1
        for ib in xrange(nband):
            ibeff=minband+ib-1
   #for ik in kptrange:
   #    ikeff = ik + 1
   #    for ib in bdrange:
   #        ibeff = ib + 1
            print " ik, ib, ikeff, ibeff", ik, ib, ikeff + 1, ibeff + 1
            #prefac=np.exp(-np.sum(ampole[ik,ib]))/np.pi*wtk[ik]*pdos[ib]*abs(imeqp[ik,ib])
            # Experimental fix for npoles dependence
            tmp = 1/np.pi*wtk[ikeff]*pdos[ib]*abs(imeqp[ik,ib])
            exponent = - np.sum(ampole[ik,ib]) - np.sum(ampole_crc[ik,ib])
            prefac = np.exp(exponent)*tmp
            #ss1=str(np.exp( - np.sum(ampole[ik,ib])))
            #print "BBB1\t"+ss1+"\n"
            #ss2=str(np.exp( - np.sum(ampole_crc[ik,ib])))
            #print "BBB2\t"+ss2+"\n"
            #ss3=str(np.exp(exponent))
            #print "BBB3\t"+ss3+"\n"
            #prefac=np.exp(-tmp*np.trapz(imskb[ik,ib],enexp)/np.sum(omegai)*npoles)
            print
            print "\n === Normalization test === " 
            print " Prefactor*wtk*pdos*Gamma/pi:", prefac
            print " Prefactor:", prefac/tmp
            print " Exponent:", exponent 
            print " Exponent/npoles:", exponent/npoles
            print
            print
            akb=ampole[ik,ib] # This is a numpy array (slice) of npoles length
            bkb =  B_crc_kb[ik,ib]/npoles # This is a number
            omegakb=omegampole[ik,ib] # NOT THE CRC ONES! IMPORTANT!!!
            eqpkb=eqp[ik,ib]
            imkb=imeqp[ik,ib]
            #tmpf1 = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
            #print nen, np.size(enexp) 
            #tmpf = 0.0*tmpf
            if eqpkb < 0.0:
                pass
            else:
                print " This state is empty! eqpkb ik ib:",eqpkb, ikeff+1, ibeff+1 
                #print "omegakb", omegakb 
                omegakb=-omegakb
                #print "-omegakb", omegakb 
            tmpf = np.zeros((nenexp), order='Fortran')
            tmpf = f2py_calc_crc_mpole(tmpf,enexp,bkb,prefac,akb,omegakb,eqpkb,imkb) #,nen,npoles)
                #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
            #if ikeff == 24:
             # outnamekb = "before" 
             # outfilekb = open(outnamekb,'w')
             # for ien in xrange(nenexp):
             #   outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], fkb[ik,ib,ien]))
             # outfilekb.close()
            fkb[ik,ib] = fkb[ik,ib] * np.exp(-np.sum(ampole_crc[ik,ib]))
            #ss = str(np.exp(-np.sum(ampole_crc[ik,ib])))
            #print "BBB\t"+ss+"\n"
            outnamekb = "spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"_crc."+str(penergy)
            outfilekb = open(outnamekb,'w')
            for ien in xrange(nenexp):
                #outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], fkb[ik,ib,ien]))
                outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], fkb[ik,ib,ien]))
            outfilekb.close()
            outnamekb = "spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"_crc_unocc."+str(penergy)
            outfilekb = open(outnamekb,'w')
            for ien in xrange(nenexp):
                outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], tmpf[ien]))
            outfilekb.close()
            ftot2 = ftot2 + tmpf
           #sfkb_c[ik,ib] = tmpf
            #ftot_crc_occ = np.sum(np.sum(fkb,0),0)
           #ftot_crc_occ = np.sum(np.sum(fkb[ik,ib],0),0)
     #  return ftot, sfkb_c
    ftot_crc_occ = np.sum(np.sum(fkb,0),0)
    plt.plot(enexp,ftot2, label='ftot_crc_unocc')
    ftot_crc = ftot_crc_occ + ftot2
    plt.plot(enexp, ftot_crc_occ, label='ftot_crc_occ')
    plt.plot(enexp, ftot_crc, label='ftot_crc')

    elaps2 = time.time() - elaps1 - e0
    cpu2 = time.clock() - cpu1 - c0
    #print elaps2, cpu2
    print str(" Used time (elaps, cpu): %10.6e %10.6e"% (elaps2, cpu2))
    print " ### Writing out A(\omega)_exp...  "
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
        outname2 = "spftot_exp"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev_np"+str(npoles)+"_crc_occ.dat"
        outfile = open(outname2,'w')
        for i in xrange(nenexp):
            outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot_crc_occ[i])) # Dump string representations of arrays
        outfile.close()
        outname2 = "spftot_exp"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev_np"+str(npoles)+"_crc.dat"
        outfile = open(outname2,'w')
        for i in xrange(nenexp):
            outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot_crc[i])) # Dump string representations of arrays
        outfile.close()
    print
    print " A(\omega)_exp written in", outname
    print " A(\omega)_CRC written in", outname2
    plt.plot(enexp,ftot,label="ftot");
# Now go back to original directory
print " Moving back to parent directory:\n", origdir
chdir(newdir)
#title = 'Spectral function '+ 'A (' + r'$\omega $ ' + ') - '+r'$ h\nu = $'+str(penergy)+' eV'
#plt.title(title)
plt.legend(loc=2);
plt.show();
