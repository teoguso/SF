#!/usr/bin/env python
"""
### Written by Matteo Guzzo ###
### A.D. MMXIV (2014)       ###
Functions necessary to the main script for the calculation 
of spectral functions. 
"""
from __future__ import print_function
import numpy as np;
import matplotlib.pylab as plt;
from scipy.interpolate import interp1d
from scipy import optimize
import sys
from os.path import isfile, join, isdir
from os import getcwd, pardir, mkdir, chdir
#
def read_invar(infile='invar.in'):
    """
    A small function that produces a list of the input variables
    using a dictionary.
    Returns a dictionary with keys and variables from the invar file.
    List of variables with default values:
     'sigmafile': None,   # Name of file containing the self-energy
     'minband': 1, Lowest band to be used in the calculation
     'maxband': 1, Highest band to be used in the calculation
     'minkpt': 1, Lowest kpt to be used in the calculation
     'maxkpt': 1, Highest kpt to be used in the calculation
     'nkpt': 1, Number of kpt (redundant, but useful to keep in mind the actual total)
     'enmin': -20.0, Lowest energy to be used in the calculation for the sf output
     'enmax': 20.0, Highest energy to be used in the calculation for the sf output
     'sfactor': 1., Percentage of LH spectrum to be calculated [0.;1.]
     'pfactor': 1., Percentage of LV spectrum to be calculated [0.;1.]
     'penergy': 0.0, Photon energy
     'npoles': 1, Number of poles for multipole fit
     'calc_gw': 1, Enables output of GW spectral function
     'calc_exp': 0, Enables output of cumulant spectral function
     'calc_crc': 0, Enables output of constrained retarded cumulant spectral function
     'np_crc': 1, Number of poles for B coefficient in constrained retarded cumulant spectral function
     'extinf': 0, Includes extrinsic and interference effects (0, 1)
     'efermi': 0.0, Fermi energy
     'omega_p':0.0, Arbitrary plasmon frequency (used only in the case npoles=999)
     'enhartree': 0, Converts energies to eV in case they are given in Hartree (0, 1)
     'gwcode': 'abinit', Name of code
     'nspin': 0, Number of spin polarizations (0, 1, 2)
     'spin': 0, Spin calculation (0, 1)
     'is_sc': 0, Self-consistent calculation (0, 1)
     'restart': 0, Will not calculate anything, but simply use the single spf files already calculated (0, 1)
     'add_wtk': 1, Include k-point weights (0, 1)
     'plot_fit': 0, Plot the fitted ImSigma as representation of poles (0, 1)
     'coarse': 0, Use coarser grid for the spectral function (faster?) (0, 1)
    """
    var_defaults = { 
            'sigmafile': None,   
            'minband': 1, 
            'maxband': 1,
            'minkpt': 1,
            'maxkpt': 1,
            'nkpt': 1,
            'enmin': -20.0,
            'enmax': 20.0,
            'sfactor': 1.,
            'pfactor': 1.,
            'penergy': 0.0, 
            'npoles': 1,
            'calc_gw': 1,
            'calc_exp': 0,
            'calc_crc': 0,
            'np_crc': 1,
            'extinf': 0,
            'efermi': 0.0,
            'omega_p':0.0,
            'enhartree': 0,
            'gwcode': 'abinit',
            'nspin': 0, 
            'spin': 0,
            'is_sc': 0,
            'restart': 0,
            'add_wtk': 1,
            'plot_fit': 0, 
            'coarse': 0
            }
#    varlist = list((
#            'sigmafile','minband','maxnband','minkpt','maxkpt',
#            'nkpt','enmin','enmax','sfactor','pfactor','penergy',
#            'npoles','calc_gw','calc_exp','extinf','efermi',
#            'omega_p','enhartree'
#            ))
    varlist = var_defaults.keys()
    invar_dict = var_defaults.copy()
    print("read_invar :: reading file: ", infile)
    with open(infile,'r') as f:
            lines = [line.strip('\n') for line in f]
    #print(varlist)
    for var in varlist:
        for line in lines:
            last =  line.split()[-1]
            if var == last:
                #print(var, last)
                invar_dict[var] = line.split()[0]
    return invar_dict

def read_hartree():
    """
    This function takes the file 'hartree.dat'
    (or alternatively the files 'Vxc.dat' and 'Elda.dat')
    and creates a 'nkpt x nband' array containing the 
    values of the hartree energy for each state. 
    This array is returned.
    All input files are supposed to be ordered in a 
    'nkpt x nband' fashion.
    """
# TODO: write a function to read the parameters from not ad-hoc files
    import numpy as np;
    if isfile("hartree.dat"):
        print(" Reading file hartree.dat... ",end="")
        hartreefile = open("hartree.dat");
        hartree = [];
        for line in hartreefile.readlines():
            hartree.append(map(float,line.split()));
        hartreefile.close()
        print("Done.")
        hartree = np.array(hartree);
    elif isfile("E_lda.dat") and isfile("Vxc.dat"):
        print("Auxiliary file (hartree.dat) not found.")
        print("Reading files E_lda.dat and Vxc.dat... ",end="")
        Eldafile = open("E_lda.dat");
        Vxcfile = open("Vxc.dat");
        elda = [];
        vxc = [];
        for line in Eldafile.readlines():
            elda.append(map(float,line.split()));
        Eldafile.close()
        for line in Vxcfile.readlines():
            vxc.append(map(float,line.split()));
        Vxcfile.close()
        print("Done.")
        elda = np.array(elda);
        vxc = np.array(vxc);
        hartree = elda - vxc
    else : 
        print("Auxiliary file not found (hartree/E_lda/Vxc). Impossible to continue.")
        sys.exit(1)
    return hartree

def read_wtk():
    """
    This function takes the file 'wtk.dat'
    and creates an array containing the 
    values of the k-point weights for each state. 
    This array is returned.
    The input file is supposed to be a single column of
    nkpt elements. 
    """
    import numpy as np;
    if isfile("wtk.dat"):
        wtkfile = open("wtk.dat");
    else : 
        print("Auxiliary file not found (wtk.dat). Impossible to continue.")
        sys.exit(1)
    wtk = [];
    for line in wtkfile.readlines():
        wtk.append((float(line)));
    wtkfile.close()
    wtk = np.array(wtk);
    return wtk

def read_occ(maxkpt,minband,maxband):
    """
    This function takes the file 'occ.dat'
    and creates an array containing the 
    occupation number for each state. 
    This array is returned.
    All input files are supposed to be ordered in a 
    'nkpt x nband' fashion.
    """
    import numpy as np;
    if isfile("occ.dat"):
        occfile = open("occ.dat");
        occ = [];
        for line in occfile.readlines():
            occ.append(map(float,line.split()));
        occfile.close()
        occ = np.array(occ)
    else : 
        print("Auxiliary file not found (occ.dat). ")
        print("Setting all occupations to 2.0.")
        occ = 2.0*np.ones((maxkpt,maxband-minband+1))
    return occ                                 

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

def calc_nkpt_sigfile(insigfile, spin = False):
    """
    Get number of kpt from a _SIG file
    """
    print("calc_nkpt_sigfile ::",end='')
    insigfile.seek(0)
    lines = insigfile.readlines()
    nkpt = 0
    for line in lines:
        if '#' in line:
            nkpt += 1
    nkpt = nkpt/2
    if spin == 1:
        nkpt = nkpt/2
    print("nkpt= {}".format(nkpt))
    return nkpt

def read_wtk_sigfile(insigfile, spin = 0):
    """
    Get wtk (k-point weights) from a SelfXC.dat file (exciting)
    """
    print("read_wtk_sigfile ::")
    insigfile.seek(0)
    lines = insigfile.readlines()
    wtk = []
    for line in lines:
        if 'wkpt' in line:
            wtk.append(float(line.split()[-1]))
    if spin == 1:
        wtk = wtk[::2]
    print("wtk = {}".format(wtk))
    print("read_wtk_sigfile :: Done.")
    return wtk

def read_sigfile(invar_dict):
    """
    A hopefully better version.
    This has to deal with the fact that abinit does not write 
    all values of sigma for one enery on a single line, but
    instead goes on a newline when a certain limit is reached.
    also: SPIN!!! For the moment, it will skip every second kpt,
    i.e. it will read only one spin channel. 
    """
    import numpy as np;
    import glob
    from time import sleep
    print("read_sigfile :: ",end="")
    # We put the content of the file (lines) in this array
    sigfilename = invar_dict['sigmafile']
    spin = int(invar_dict['spin'])
    nspin = int(invar_dict['nspin'])
    #en=[] 
    firstbd = 0
    lastbd = 0
    nbd = 0
    #sngl_row = True # are data not split in more than 1 line?
    if sigfilename is None:
        print("File "+str(sigfilename)+" not defined.")
        sigfilename = glob.glob('*_SIG')[0]
        print("Looking automatically for a _SIG file... ",sigfilename)
    with open(sigfilename) as insigfile:
        filelines = insigfile.readlines() 
        nkpt = calc_nkpt_sigfile(insigfile,spin)
        if invar_dict['gwcode'] == 'exciting':
            invar_dict['wtk'] = read_wtk_sigfile(insigfile)
        insigfile.seek(0)
        insigfile.readline()
        line = insigfile.readline()
        firstbd = int(line.split()[-2])
        lastbd =  int(line.split()[-1])
        nbd = lastbd - firstbd + 1
        print("nbd:",nbd)
        num_cols = len(insigfile.readline().split())
        num_cols2 = len(insigfile.readline().split())
        print("numcols:",num_cols)
        print("numcols2:",num_cols2)
        if num_cols != num_cols2: 
            print()
            print(" WARNING: newlines in _SIG file.")
            print(" Reshaping _SIG file structure...")
            print(" _SIG file length (rows):", len(filelines))
            new_list = []
            nline = 0
            a = []
            b = []
            for line in filelines:
                #if line.split()[0] == "#":
                if '#' in line:
                    print(line.strip('\n'))
                    continue
                elif nline == 0: 
                    a = line.strip('\n')
                    nline += 1
                else: 
                    b = line.strip('\n')
                    new_list.append(a + " " + b)
                    nline = 0
            print("New shape for _SIG array:",np.asarray(new_list).shape)
            tmplist = []
            tmplist2 = []
            for line in new_list:
                tmplist.append(map(float,line.split())[0])
                tmplist2.append(map(float,line.split())[1:])
            for el1 in tmplist2:
                for j in el1:
                    try:
                        float(j)
                    except:
                        print(j)
            #tmplist = map(float,tmplist)
            #tmplist2 = map(float,tmplist2)
            xen = np.asarray(tmplist)
            x = np.asarray(tmplist2)
        else:
            insigfile.seek(0)
            xen = np.genfromtxt(sigfilename,usecols = 0)
            insigfile.seek(0)
            x = np.genfromtxt(sigfilename,usecols = range(1,num_cols), filling_values = 'myNaN')
    #nkpt = int(invar_dict['nkpt'])
    print("nkpt:",nkpt)
    print("spin:",spin)
    print("nspin:",nspin)
    # From a long line to a proper 2D array, then only first row
    #print(xen.shape)
    print("x.shape", x.shape)
    if spin == 1 and nspin == 0:
        nspin = 2
    else:
        nspin = 1
    print("nspin:",nspin)
    print("size(xen):",xen.size)
    print("The size of a single energy array should be",\
            float(np.size(xen))/nkpt/nspin)
    en = xen.reshape(nkpt*nspin,np.size(xen)/nkpt/nspin)[0]
    #en = xen.reshape(nkpt,np.size(xen)/nkpt)[0]
    print("New shape en:",np.shape(en))
    print("First row of x:",x[0])
    if invar_dict['gwcode'] == 'abinit':
        nb_cols = 3
    elif invar_dict['gwcode'] == 'exciting':
        nb_cols = 2
       #b = x.reshape(nkpt*nspin,np.size(x)/nkpt/nspin/nbd/3,3*nbd)
    b = x.reshape(nkpt*nspin,np.size(x)/nkpt/nspin/nbd/nb_cols,nb_cols*nbd)
    print("New shape x:",b.shape)
    y = b[0::nspin,:,0::nb_cols]
    z = b[0::nspin,:,1::nb_cols]
    res = np.rollaxis(y,-1,1)
    ims = np.rollaxis(z,-1,1)
    print("New shape res, ims:", res.shape)
    print("First and last band in _SIG file:", firstbd, lastbd)
    print(" Done.")
    return en, res, ims, (firstbd, lastbd)

def read_cross_sections(penergy):
    """
    This function should read the values for the cross sections
    given a photon energy 'penergy' from the file
    "cs'penergy'.dat" and construct an array "cs"
    containing the values.
    For now only s and p are expected, but they can be added
    seamlessly to the file. Only, the other part of the code
    using the cs array would have to be changed accordingly. 
    cs array is returned. 
    """
    import numpy as np;
    print(" ### Reading cross sections...  ")
    if int(penergy) == 0:
        cs = np.array([0.1,0.1])
    else:
        csfilename = "cs"+str(int(penergy))+".dat"
        if isfile(csfilename):
            print(" Photon energy:", penergy,"eV")
        else:
            penergy = raw_input(" File "+csfilename+" not found. Photon energy (eV): ")
            csfilename = "cs"+str(penergy)+".dat"
        cs = []
        print(" csfilename:",csfilename)
        csfile = open(csfilename,'r')
        for line in csfile.readlines():
            cs.append((float(line)));
        csfile.close()
        cs = np.array(cs)
    #print("cs:",cs.shape,cs)
    #print("cs:",np.transpose(cs),cs.shape)
    return cs

def calc_pdos(var_dct,res=None):
    """
    Calculates projected DOS for all bands contained in the
    sigma file. 
    """
    import numpy as np
    penergy = float(var_dct['penergy'])
    sfac =  float(var_dct['sfactor'])
    pfac =  float(var_dct['pfactor'])
    #nband = int(var_dct['nband'])
    nband = int(var_dct['sig_bdgw'][1]) - int(var_dct['sig_bdgw'][0])
    if penergy != 0:
        # ======== CROSS SECTIONS ======= #
        cs = read_cross_sections(penergy)
        # ====== BAND TYPE AND SYMMETRY ==== #
        sp = read_band_type_sym(sfac,pfac,nband)
        # ===== EFFECTIVE STATE-DEPENDENT PREFACTOR ==== #
        pdos = 10000.*np.dot(cs,sp)
    else:
        pdos=np.ones((nband))
    if res is not None:
        tail = np.ones(shape=(res[0,:,0].size-pdos.size))*pdos[-1] # same weight as last item of pdos
        pdos = np.concatenate((pdos,tail))
    return pdos

def read_band_type_sym(sfac,pfac,nband):
    """
    This function reads the s,p (TODO and d,f)
    composition of bands, descerning between 
    s (mirror-even) and p (mirror-odd) symmetries
    with respect to a plane normal to the surface 
    of the sample. 
    By consequence, the preparation of the input
    files s.dat, p_even.dat and p_odd.dat is crucial
    for a correct description of cross-section 
    and symmetry effects. 
    The function takes sfac and pfac as inputs, so that 
    one can simulate LH or LV light measurements. 
    A 'sp' symmetry-biased array is created, following the 
    number of band types that are considered (s and p for now). 
    sp numpy array is returned.
    """
    import numpy as np;
    print(" Reading s bands file... ",)
    if isfile("s.dat"):
        sfile =  open("s.dat",'r')
        s = map(float,sfile.read().split())
        sfile.close()
        s = np.array(s)
        print("Done.")
    else : 
        print()
        print(" WARNING: File for orbital character not found (s.dat). S character will be 1 for all bands. ")
        s = np.ones(nband)
    print(" Reading p bands file... ",)
    if isfile("p_even.dat") and isfile("p_odd.dat"):
        # This file should contain the summed contribution of all even p orbitals
        pevenfile = open("p_even.dat",'r')
        # This file should contain the summed contribution of all odd p orbitals
        poddfile =  open("p_odd.dat",'r')
        peven = map(float,pevenfile.read().split())
        podd = map(float,poddfile.read().split())
        pevenfile.close()
        poddfile.close()
        peven = np.array(peven)
        podd = np.array(podd)
        print("Done.")
    else : 
        print()
        print(" WARNING: File for orbital character not found (p_even.dat/p_odd.dat). P character will be 1 for all bands. ")
        peven = np.ones(nband)
        podd = np.ones(nband)
    # TODO: Add d bands!!!
    # Warning: in this section it is easy to confuse
    # s and p symmetry with s and p electrons! 
    # Don't worry, it's normal.
    s = sfac*s
    peven = sfac*peven
    podd = pfac*podd
    p = peven+podd
    sp = np.array([s,p])
    #print("sp:",sp)
    return sp

#def calc_spf_gw(minkpt,maxkpt,minband,maxband,wtk,pdos,en,enmin,enmax,res,ims,hartree):
def calc_sf_gw(vardct,hartree,pdos,en,res,ims):
    """
    Macro-function calling instructions necessary to calculate 
    the GW spectral function. 
    For now it writes out the single state spectral functions
    on ascii files. This should be moved to an external module
    and just return spfkb as an output variable. 
    spf (GW spectral function) is returned.
    """
    import numpy as np;
    print("calc_sf_gw :: ")
    wtk = np.array(vardct['wtk'])
    hartree = np.array(hartree)
    pdos = np.array(pdos)
   #minkpt = int(vardct['minkpt'])
   #maxkpt = int(vardct['maxkpt'])
    nkpt = res[:,0,0].size 
    #maxkpt - minkpt + 1
   #minband = int(vardct['minband'])
   #maxband = int(vardct['maxband'])
    nband =  res[0,:,0].size 
    print("nkpt, nband ", nkpt, nband)
   #bdgw = map(int, vardct['sig_bdgw'])
   #bdrange = range(minband-bdgw[0],maxband-bdgw[0]+1)
   #kptrange = range(minkpt - 1, maxkpt)
    #maxband - minband + 1
    coarse = int(vardct['coarse'])
    if coarse == 1: 
        newdx = 0.2
    else:
        newdx = 0.005
    enmin = float(vardct['enmin'])
    enmax = float(vardct['enmax'])
    if enmin < en[0] and enmax >= en[-1]:  
        newen = np.arange(en[0],en[-1],newdx)
    elif enmin < en[0]:  
        newen = np.arange(en[0],enmax,newdx)
    elif enmax >= en[-1] :  
        newen = np.arange(enmin,en[-1],newdx)
    else :  
        newen = np.arange(enmin,enmax,newdx)
    print(" ### Interpolation and calculation of A(\omega)_GW...  ")
    spftot = np.zeros((np.size(newen)));
    # Here we interpolate re and im sigma
    # for each band and k point
    spfkb = np.zeros(shape=(res[:,0,0].size,res[0,:,0].size,np.size(newen)))
    reskb = np.zeros(shape=(res[:,0,0].size,res[0,:,0].size,np.size(newen)))
    imskb = np.zeros(shape=(res[:,0,0].size,res[0,:,0].size,np.size(newen)))
    rdenkb = np.zeros(shape=(res[:,0,0].size,res[0,:,0].size,np.size(newen)))
   #spfkb = np.zeros(shape=(nkpt,nband,np.size(newen)))
   #reskb = np.zeros(shape=(nkpt,nband,np.size(newen)))
   #imskb = np.zeros(shape=(nkpt,nband,np.size(newen)))
   #rdenkb = np.zeros(shape=(nkpt,nband,np.size(newen)))
    for ik in range(nkpt):
        #ikeff = minkpt+ik-1
        print(" k point, nband = %02d %02d" % (ik,nband))
        #print(" nband = %02d " % (nband))
        for ib in range(nband):
            #ibeff = minband+ib-1
           #print(ik, ib)
           #plt.plot(en,ims[ik,ib])
            interpres = interp1d(en, res[ik,ib], kind = 'linear', axis = -1)
            interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
            tmpres = interpres(newen)
            reskb[ik,ib] = tmpres
            #redenom = newen + efermi - hartree[ik,ib] - interpres(newen)
            redenom = newen - hartree[ik,ib] - interpres(newen)
            #plt.plot(hartree[ik,ib],'o')
            rdenkb[ik,ib] = redenom
            #print("ik ib minband maxband ibeff hartree[ik,ib]", ik, ib, minband, maxband, ibeff, hartree[ik,ib])
            tmpim = interpims(newen)
            imskb[ik,ib] = tmpim
           #print("pdos ",pdos)
            spfkb_tmp = wtk[ik] * pdos[ib] * abs(tmpim)/np.pi/(redenom**2 + tmpim**2)
            #print(spfkb.shape, spfkb_tmp.shape)
            spfkb[ik,ib] = spfkb_tmp
            spftot += spfkb_tmp
    allkb = [spfkb, reskb, rdenkb, imskb]
   #plt.plot(newen,spftot)
    print("reskb.shape:",reskb.shape)
    return newen, spftot, allkb
    #return newen, spftot, spfkb,reskb, rdemkb, imskb

def write_spfkb(vardct, newen, allkb):
    """
    Does what it says. For the GW spectral function. 
    """
    import numpy as np
    print("write_spfkb :: ")
    minkpt = int(vardct['minkpt'])
    maxkpt = int(vardct['maxkpt'])
    nkpt = maxkpt - minkpt + 1
    minband = int(vardct['minband'])
    maxband = int(vardct['maxband'])
    nband = maxband - minband + 1
    bdgw = map(int, vardct['sig_bdgw'])
    bdrange = range(minband - bdgw[0], maxband - bdgw[0] + 1)
    kptrange = range(minkpt - 1, maxkpt)
    spfkb = allkb[0]
    reskb = allkb[1]
    rdenkb = allkb[2]
    imskb = allkb[3]
   #for ik in range(nkpt):
   #    #print(" k point = %02d " % (ikeff+1))
   #    ikeff = minkpt + ik 
   #    for ib in range(nband):
    for ik in kptrange:
        ikeff = ik + minkpt
        #for ib in range(nband):
        for ib in bdrange:
            ibeff = ib + minband
            outnamekb = "spf_gw-k"+str("%02d"%(ikeff))+"-b"+str("%02d"%(ibeff))+".dat"
            outfilekb = open(outnamekb,'w')
            for ien in xrange(np.size(newen)) :
                outfilekb.write("%8.4f %12.8e %12.8e %12.8e %12.8e\n" % (newen[ien], spfkb[ik,ib,ien], rdenkb[ik,ib,ien], reskb[ik,ib,ien], imskb[ik,ib,ien]))
            outfilekb.close()
    print("write_spfkb :: Done.")

def find_eqp_resigma(en, resigma, efermi):
    """
    This function is supposed to deal with the plasmaron problem 
    and calculate the quasiparticle energy once it is fed with 
    resigma = \omega - \epsilon_H - \Re\Sigma. 
    It expects an array of increasing values on the x axis 
    and it will return 
    the x value of the last resigma=0 detected. 
    It should return the value of eqp and the number of zeros
    found (useful in case there are plasmarons or for debugging). 
    If no zeros are found, it will fit resigma with a line and 
    extrapolate a value.
    """
    import numpy as np
    nzeros=0
    zeros = []
    tmpeqp = en[0]
    for i in xrange(1,np.size(resigma)):
        #print(resigma[i]*resigma[i-1] # DEBUG)
        if  resigma[i] == 0: # Yes, it can happen
            tmpeqp = en[i] 
            zeros.append(en[i])
            nzeros+=1
        elif (resigma[i]*resigma[i-1] < 0):
            tmpeqp = en[i-1] - resigma[i-1]*(en[i] - en[i-1])/(resigma[i] - resigma[i-1]) # High school formula
            zeros.append(tmpeqp)
            nzeros+=1
    if tmpeqp>efermi: 
        tmpeqp=zeros[0]
    if nzeros==0 : 
        print()
        print(" WARNING: No eqp found! ")
        def fit_func(x, a, b): 
            return a*x + b
        from scipy.optimize import curve_fit
        params = curve_fit(fit_func, en, resigma)
        [a, b] = params[0]
        if -b/a < en[-1]:
            print("WTF!!! BYE!")
            sys.exit()
        tmpeqp = -b/a
        zeros.append(tmpeqp)
    elif nzeros>1 : 
        print()
        print(" WARNING: Plasmarons! ")
    return tmpeqp, nzeros

def write_eqp_imeqp(eqp,imeqp):
    """
    Does what it says.
    """
    import numpy as np
    print("write_eqp_imeqp :: ")
    outname = "eqp.dat"
    outfile2 = open(outname,'w')
    outname = "imeqp.dat"
    outfile3 = open(outname,'w')
    for ik in xrange(np.size(eqp[:,0])):
        for ib in xrange(np.size(eqp[0,:])):
            outfile2.write("%14.5f" % (eqp[ik,ib]))
            outfile3.write("%14.5f" % (imeqp[ik,ib]))
        outfile2.write("\n")
        outfile3.write("\n")
    outfile2.close()
    outfile3.close()
    print("write_eqp_imeqp :: Done.")
    
def calc_eqp_imeqp(en,res,ims,hartree,efermi):
    """
    This function calculates qp energies and corresponding
    values of the imaginary part of sigma for a set of
    k points and bands. 
    The function find_eqp_resigma() is used here.
    eqp and imeqp are returned. 
    """
    import numpy as np
    from scipy import interp
    nkpt = np.size(res[:,0,0])
    nband = np.size(res[0,:,0])
    eqp = np.zeros((nkpt,nband))
    imeqp = np.zeros((nkpt,nband))
    hartree = np.array(hartree)
    for ik in range(nkpt):
        for ib in range(nband):
            #temparray = np.array(en - hartree[ik,ib] - res[ik,ib])
            temparray = np.array(en - hartree[ik,ib] - res[ik,ib])
            interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
            tempim = interpims(en)
            # New method to overcome plasmaron problem
            eqp[ik,ib], nzeros = find_eqp_resigma(en,temparray,efermi)
            if nzeros==0: 
                print()
                print(" WARNING: ik "+str(ik)+" ib "+str(ib)+". No eqp found!!!")
            if (eqp[ik,ib] > en[0]) and (eqp[ik,ib] < en[-1]): 
                #print(en[0], eqp[ik,ib], en[-1])
                imeqp[ik,ib] = interpims(eqp[ik,ib])
            else:
                imeqp[ik,ib] = interp(eqp[ik,ib], en, ims[ik,ib])
            ## Warning if imaginary part of sigma < 0 (Convergence problems?)
            if imeqp[ik,ib] <= 0 : 
                print()
                print(" WARNING: im(Sigma(eps_k)) <= 0 !!! ik ib eps_k im(Sigma(eps_k)) = ", ik, ib, eqp[ik,ib], imeqp[ik,ib])
    return eqp, imeqp

def calc_extinf_corrections(origdir,extinfname,ampole,omegampole):
    """
    # Here we add the extrinsic contribution. 
    # N.B.: It has to be renormalized to the number of poles!!!
    # The data are interpolated linearly with a numerical function. 
    # The fit curve of ext+inf passes by the origin. 
    The file structure is expected to be:
    #  wp, aext, ainf, aint, width
    """
    import numpy as np;
    from multipole import getdata_file #, write_f_as_sum_of_poles
    #extinfname = "a_wp.dat"
    print(" Reading extrinsic and interference contribution from file "+str(extinfname)+"...")
    en_ei, aext = getdata_file(origdir+"/"+str(extinfname))
    en_ei, ainf = getdata_file(origdir+"/"+str(extinfname),2)
    aextinf = aext+ainf
    newen_ei = []
    newen_ei.append(0.0)
    for x in en_ei.tolist():
        newen_ei.append(x)
    newen_ei = np.array(newen_ei)
    newa_ei = []
    newa_ei.append(0.0)
    for x in aextinf.tolist():
        newa_ei.append(x)
    newa_ei = np.array(newa_ei)
    # a_int from the model is in the third column
    en_ei, aint = getdata_file(origdir+"/"+str(extinfname),3)
    newa_int = []
    a_int_zero = aint[1] - en_ei[1]*(aint[0]-aint[1])/(en_ei[0]-en_ei[1])
    newa_int.append(a_int_zero)
    for x in aint.tolist():
        newa_int.append(x)
    newa_int = np.array(newa_int)
    # broadening from the model is in the fourth column
    en_ei, width = getdata_file(origdir+"/"+str(extinfname),4)
    newwmod = []
    w_zero = width[1] - en_ei[1]*(width[0]-width[1])/(en_ei[0]-en_ei[1])
    newwmod.append(w_zero)
    for x in width.tolist():
        newwmod.append(x)
    newwmod = np.array(newwmod)
    interpwidth = interp1d(newen_ei, newwmod, kind = 'linear', axis = -1)
    w_extinf = ampole.copy()
    print("omega_p, a_extinf, a_int:")
    print(newen_ei)
    print(newa_ei)
    print(newa_ei/newa_int)
    #print(en_ei, newenexin)
    #print(aextinf, newaexin)
    interpextinf = interp1d(newen_ei, newa_ei/newa_int, kind = 'linear', axis = -1)
    amp_exinf = ampole.copy()
    #print("Type(amp_exinf, ampole):", type(amp_exinf), type(ampole))
    # Mod following discussion with Josh
    amp_mean = np.mean(ampole)
    for ik in range(ampole[:,0,0].size):
        for ib in range(ampole[0,:,0].size):
            #tmpextinf = interpextinf(omegampole[ik,ib])/npoles # <-- Divided by the number of poles (normalization)!
            try: 
                w_extinf[ik,ib] = interpwidth(omegampole[ik,ib]) # Numpy array
            except ValueError:
                print()
                print(" WARNING: A value for omega_p is beyond what contained in a_wp.x.")
                print(" The last available value is taken. ")
                print("ik, ib, omegampole: ", ik, ib, omegampole[ik,ib])
                w_extinf[ik,ib] = np.interp(omegampole[ik,ib],newen_ei,newwmod)
            try: 
                tmpextinf = interpextinf(omegampole[ik,ib]) # 
            except ValueError:
                print()
                print(" WARNING: A value for omega_p is beyond what contained in a_wp.x.")
                tmpextinf = np.interp(omegampole[ik,ib],newen_ei, newa_ei/newa_int)
            # Mod following discussion with Josh
            #amp_exinf[ik,ib] += ampole[ik,ib] * tmpextinf
            amp_exinf[ik,ib] += amp_mean * tmpextinf
    return amp_exinf, w_extinf

def read_sf(vardct, pdos, approx):
    """
    This function reads the spectral functions from  the files 
    spf_gw- or spf_exp- and returns their sum, 
    according to minband, maxband and minkpt, maxkpt. 
    """
    import numpy as np
    print("read_sf :: ")
    minkpt = int(vardct['minkpt'])
    maxkpt = int(vardct['maxkpt'])
    nkpt = maxkpt - minkpt + 1
    minband = int(vardct['minband'])
    maxband = int(vardct['maxband'])
    nband = maxband - minband + 1
    extinf = int(vardct['extinf'])
    npoles = int(vardct['npoles'])
    penergy = int(vardct['penergy'])
    #wtk = np.array(vardct['wtk'])
    #hartree = np.array(hartree)
    bdrange = vardct['bdrange']
    print("bdrange", bdrange)
    kptrange = vardct['kptrange']
    print("kptrange", kptrange)
    pdos = np.array(pdos)
    if extinf == 1: 
        str_exi = "_extinf"
    else:
        str_exi = ""
    if approx == 'exp':
        end_fname = "_np"+str(npoles)+str_exi+"."+str(penergy)
    elif approx == 'gw':
        end_fname = ".dat"
    # Initialize sf
    ikeff = minkpt
    ibeff = minband
    fname = "spf_"+str(approx)+"-k"+str("%02d"%(ikeff))+"-b"+str("%02d"%(ibeff))+str_exi+end_fname
    en = np.genfromtxt(fname, usecols = 0) # sigfilename,usecols = range(1,num_cols), filling_values = 'myNaN')
    #sf = np.genfromtxt(fname, usecols = 1) # sigfilename,usecols = range(1,num_cols), filling_values = 'myNaN')
    sf = np.zeros((en.size))
    #bdgw = map(int, vardct['sig_bdgw'])
    #for ik in range(nkpt):
    for ik in kptrange:
        #print(" k point = %02d " % (ikeff+1))
        ikeff = ik + 1 
        #for ib in range(nband):
        for ib in bdrange:
            #ibeff = minband + ib
            ibeff = ib + 1 
            print("ikeff, ibeff: ",ikeff,ibeff)
            #outnamekb = "spf_gw-k"+str("%02d"%(ikeff))+"-b"+str("%02d"%(ibeff))+".dat"
            #outnamekb = "spf_"+str(approx)+"-k"+str("%02d"%(ikeff))+"-b"+str("%02d"%(ibeff))+"_np"+str(npoles)+str_exi+"."+str(penergy)
            fname = "spf_"+str(approx)+"-k"+str("%02d"%(ikeff))+"-b"+str("%02d"%(ibeff))+str_exi+end_fname
            print("fname: ",fname)
            tmp_sf = np.genfromtxt(fname, usecols = 1) # sigfilename,usecols = range(1,num_cols), filling_values = 'myNaN')
           #en, tmp_sf = np.genfromtxt(fname) # sigfilename,usecols = range(1,num_cols), filling_values = 'myNaN')
            sf += tmp_sf
           #with open(fname,'r') as ifkb:
           #    for ien in range(en.size):
           #        ifkb.write("%8.4f %12.8f\n" % (en[ien], sfkb[ik,ib,ien]))
    print("read_sf :: Done.")
    return en, sf

def calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles,wkb=None):
    """
    OBSOLETE???
    This function calculates the exponential spectral function. 
    """
    import numpy as np;
    nenexp=np.size(enexp)
    ftot = np.zeros((np.size(enexp)))
    tmpf1=np.zeros((nenexp))
    tmpf2=np.zeros((nenexp))
    tmpf3=np.zeros((nenexp))
    if wkb is None:
        wkb = np.zeros(npoles)
    for ipole in xrange(npoles):
        #print(" First order")
        tmpf2=0.
        for jpole in xrange(npoles):
            #print(" Second order")
            tmpf3=0.
            for kpole in xrange(npoles):
                ## Fourth order
                #tmpf4 = np.zeros((np.size(enexp)))
                #for ipole4 in xrange(npoles):
                    #tmpf4[:] += 1./4. * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole]  + omegampole[ik,ib,jpole]  + omegampole[ik,ib,kpole]  + omegampole[ik,ib,ipole4] )**2 + ( imeqp[ik,ib] )**2 ) )
                tmpomegap = omegakb[ipole]+omegakb[jpole]+omegakb[kpole]
                tmpgamma = (wkb[ipole]+wkb[jpole]+wkb[kpole])/2
                tmpf3+=1./3.*akb[kpole]/((enexp-eqpkb+tmpomegap)**2+(tmpgamma+imkb)**2 ) 
            tmpomegap = omegakb[ipole]+omegakb[jpole]
            tmpgamma = (wkb[ipole]+wkb[jpole])/2
            tmpf2+=1./2.*akb[jpole]*(1./((enexp - eqp[ik,ib]+tmpomegap)**2+(tmpgamma+imkb)**2 ) + tmpf3 ) 
        tmpomegap = omegakb[ipole]
        tmpgamma = (wkb[ipole])/2
        tmpf1+=1.*akb[ipole]*(1./((enexp-eqp[ik,ib]+tmpomegap)**2+(tmpgamma+imkb)**2)+tmpf2) 
    #f[ik,ib]=prefac*(1./((enexp-eqpkb)**2+(imkb)**2)+tmpf1) 
    #ftot += f[ik,ib]
    f=prefac*(1./((enexp-eqpkb)**2+(imkb)**2)+tmpf1) 
    ftot += f
    return ftot

def write_sfkb_c(vardct,en,sfkb):
    """
    Does what it says. For the cumulant spectral function. 
    """
    import numpy as np
    import matplotlib.pylab as plt
    sfkb = np.array(sfkb)
    print("write_sfkb_c :: ")
    minkpt = int(vardct['minkpt'])
    maxkpt = int(vardct['maxkpt'])
    nkpt = maxkpt - minkpt + 1
    minband = int(vardct['minband'])
    maxband = int(vardct['maxband'])
    nband = maxband - minband + 1
    bdgw = map(int, vardct['sig_bdgw'])
    bdrange = range(minband-bdgw[0],maxband-bdgw[0]+1)
    kptrange = range(minkpt - 1, maxkpt)
    extinf = int(vardct['extinf'])
    npoles = int(vardct['npoles'])
    penergy = int(vardct['penergy'])
   #plt.plot(en,sfkb[0,4])
   #plt.show()
    if extinf == 1: 
        str_exi = "_extinf"
    else:
        str_exi = ""
    for ik in kptrange:
        #print(" k point = %02d " % (ikeff+1))
       #ikeff = minkpt + ik 
        ikeff = ik + 1
        for ib in bdrange:
           #ibeff = minband + ib
            ibeff = ib + bdgw[0] 
            outnamekb = "spf_exp-k"+str("%02d"%(ikeff))+"-b"+str("%02d"%(ibeff))+"_np"+str(npoles)+str_exi+"."+str(penergy)
            print("ik,ib: ",ik,ib)
            print(" Writing on file : ", outnamekb)
            with open(outnamekb,'w') as ofkb:
                for ien in range(en.size):
                    ofkb.write("%9.4f %15.8f\n" % (en[ien], sfkb[ik,ib,ien]))
    print("write_sfkb_c :: Done.")

def calc_sf_c(vardct, hartree, pdos, eqp, imeqp, newen, allkb):
    """
    Meta-function that calls serial or para version.
    Still wishful thinking ATM, as it just calls the serial version. 
    """
    import numpy as np
    npoles = int(vardct['npoles'])
   #if npoles == 0 or npoles == 1 or npoles == 999: 
    while True:
        enexp, ftot, sfkb_c = \
                calc_sf_c_serial(\
                vardct, hartree, pdos, eqp, imeqp, newen, allkb)
        break
   #else:
   #    enexp, ftot, sfkb_c = \
   #            calc_sf_c_para(\
   #            vardct, hartree, pdos, eqp, imeqp, newen, allkb)
    return  enexp, ftot, sfkb_c

def calc_multipole(npoles, imskb, kptrange, bdrange, eqp, newen):
    """
    Function that calculates frequencies and amplitudes
    of ImSigma using the multipole model. 
    """
    import numpy as np
    from multipole import fit_multipole #, write_f_as_sum_of_poles
    print(" ### ================== ###")
    print(" ###    Multipole fit   ###")
    print(" Number of poles:", npoles)
    omegampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,npoles))
    ampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,npoles))
    for ik in kptrange:
        for ib in bdrange:
            if eqp[ik,ib] > newen[-npoles]:
                omegampole[ik,ib] = omegampole[ik,ib-1]
                ampole[ik,ib] = ampole[ik,ib-1]
                print(" Eqp beyond available energy range. Values from lower band are taken.")
                continue
            else:
                ibeff = bdrange[0] + ib - 1
                print(" ik, ib", ik, ib)
                # Here we take the curve starting from eqp and then we invert it
                # so as to have it defined on the positive x axis
                # and so that the positive direction is in the 
                # increasing direction of the array index
                if eqp[ik,ib] <= 0:
                    en3 = newen[newen<=eqp[ik,ib]] # So as to avoid negative omegampole
                    im3 = abs(imskb[ik,ib][newen<=eqp[ik,ib]]/np.pi) # This is what should be fitted
                else:
                    en3 = newen[newen>eqp[ik,ib]] # So as to avoid negative omegampole
                    im3 = abs(imskb[ik,ib][newen<=eqp[ik,ib]]/np.pi) # This is what should be fitted
                if en3.size == 0:
                    print()
                    print(" WARNING: QP energy is outside of given energy range!\n"+\
                            " This state will be skipped!\n"+\
                            "You might want to modify enmin/enmax.")
                    print(" eqp[ik,ib], newen[-1]", eqp[ik,ib] , newen[-1])
                    continue
                en3 = en3 - eqp[ik,ib]
                if eqp[ik,ib] <= 0:
                    en3 = -en3[::-1] 
                    im3 = im3[::-1]
                omegai, lambdai, deltai = fit_multipole(en3,im3,npoles)
                # HERE WE MUST CHECK THAT THE NUMBER OF POLES 
                # IS NOT BIGGER THAN THE NUMBER OF POINTS THAT HAS TO BE FITTED
                if npoles > omegai.size:
                    omegampole[ik,ib][:omegai.size] = omegai 
                    ampole[ik,ib][:omegai.size] = np.true_divide(lambdai,(np.square(omegai)))
                    print()
                    print(" WARNING: npoles used ("+str(npoles)+") is larger"+\
                            " than x data array ("+str(omegai.size)+").")
                    print(" Reduce npoles. You are wasting resources!!!")
                else:
                    omegampole[ik,ib] = omegai 
                    ampole[ik,ib] = np.true_divide(lambdai,(np.square(omegai)))
                print(" Integral test. Compare \int\Sigma and \sum_j^N\lambda_j.")
                print(" 1/pi*\int\Sigma   =", np.trapz(im3,en3))
                print(" \sum_j^N\lambda_j =", np.sum(lambdai))
    # Writing out a_j e omega_j
    print(" ### Writing out a_j and omega_j...")
    outname = "a_j_np"+str(npoles)+".dat"
    outfile = open(outname,'w')
    outname = "omega_j_np"+str(npoles)+".dat"
    outfile2 = open(outname,'w')
    for ipole in xrange(npoles):
        for ik in range(imskb[:,0,0].size):
            for ib in range(imskb[0,:,0].size):
                outfile.write("%10.5f"  % (ampole[ik,ib,ipole]))
                outfile2.write("%10.5f" % (omegampole[ik,ib,ipole]))
            outfile.write("\n")
            outfile2.write("\n")
        outfile.write("\n")
        outfile2.write("\n")
    outfile.close()
    outfile2.close()
    return omegampole, ampole

def calc_extinf(vardct, ampole, omegampole):
    """
    # Extrinsic and interference contribution
    """
    import numpy as np
    penergy = int(vardct['penergy'])
    origdir = vardct['origdir']
    extinfname = "a_wp."+str(penergy)
    amp_exinf, w_extinf = calc_extinf_corrections(origdir,extinfname,ampole,omegampole)
    print(" ### Writing out a_j_extinf...")
    outname = "a_j_np"+str(npoles)+"_extinf."+str(penergy)
    with open(outname,'w') as outfile:
        for ipole in xrange(npoles):
            for ik in range(imskb[:,0,0].size):
                for ib in range(imskb[0,:,0].size):
                    outfile.write("%10.5f"  % (amp_exinf[ik,ib,ipole]))
                outfile.write("\n")
            outfile.write("\n")
    return amp_extinf, w_extinf

def calc_sf_c_para(vardct, hartree, pdos, eqp, imeqp, newen, allkb):
    """
    STILL BROKEN!!!!
    Parallel version of calc_sf_c_serial, hopefully more modular,
    functional and polished than its predecessor. 
    The idea is to parallelize the energies over the number of 
    cores that will be detected.
    """
    print(" calc_sf_c_para :: ")
    import numpy as np;
    wtk = np.array(vardct['wtk'])
    hartree = np.array(hartree)
    pdos = np.array(pdos)
    minkpt = int(vardct['minkpt'])
    maxkpt = int(vardct['maxkpt'])
    #nkpt = maxkpt - minkpt + 1
    minband = int(vardct['minband'])
    maxband = int(vardct['maxband'])
    nband = maxband - minband + 1
    bdgw = map(int, vardct['sig_bdgw'])
    bdrange = range(minband-bdgw[0],maxband-bdgw[0]+1)
    kptrange = range(minkpt - 1, maxkpt)
    newdx = 0.005
    enmin = float(vardct['enmin'])
    enmax = float(vardct['enmax'])
    npoles = int(vardct['npoles'])
    extinf = int(vardct['extinf'])
    reskb = allkb[1]
    imskb = allkb[3]
    # Setting up multipole:
    omegampole, ampole = calc_multipole(npoles, imskb, kptrange, bdrange, eqp, newen)
    if extinf == 1:
        amp_extinf, w_extinf = calc_extinf(vardct, ampole, omegampole)
    print(" Calculating multipole exponential A...")
    dxexp=0.005 
    enexp = np.arange(enmin,enmax,dxexp)
    nenexp = np.size(enexp)
    ftot = np.zeros((np.size(enexp)),order='Fortran')
    sfkb_c = np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,nenexp))
    from extmod_spf_mpole import f2py_calc_spf_mpole, f2py_calc_spf_mpole_extinf
    for ik in kptrange:
        ikeff = ik + 1
        for ib in bdrange:
            ibeff=bdgw[0]+ib
            print(" ik, ib, ikeff, ibeff", ik, ib, ikeff, ibeff)
            tmp = 1/np.pi*wtk[ik]*pdos[ib]*abs(imeqp[ik,ib])
            prefac=np.exp(-np.sum(ampole[ik,ib]))*tmp
            akb=ampole[ik,ib] # This is a numpy array (slice)
            omegakb=omegampole[ik,ib] # This is a numpy array (slice)
            eqpkb=eqp[ik,ib]
            imkb=imeqp[ik,ib] # + w_extinf[ik,ib]/2 # extinf width added
            if eqpkb < 0.0:
                pass
            else:
                print(" This state is empty! eqpkb ik ib:",eqpkb, ikeff+1, ibeff+1)
                omegakb=-omegakb
            # The calculation can be done in both cases by the extinf version, 
            # with wkb = 0 for the intrinsic. 
            if extinf == 1: 
                akb=amp_exinf[ik,ib] # This is a numpy array (slice)
                wkb=w_extinf[ik,ib] # This is a numpy array (slice) 
            else: 
                wkb = np.zeros((akb.size))
               #tmpf = f2py_calc_spf_mpole(tmpf,enexp,prefac,akb,omegakb,eqpkb,imkb) #,nen,npoles)
            tmpf = np.zeros((nenexp), order='Fortran')
            # PARALLELISM STARTS HERE
            # PARALLELIZATION OVER ENERGIES
            print(" ==================== ")
            print(" PARALLELIZATION HERE ")
            print(" ==================== ")
            import multiprocessing as mp
            ncpu = mp.cpu_count()
           #ncpu = 3
            print("TEST cpu_count()", ncpu)
            if enexp.size%ncpu > 0: 
                bite = enexp.size/ncpu + 1
            else: 
                bite = enexp.size/ncpu
            print("bite ", bite)
            print("bite*ncpu", bite*ncpu)
            print("enexp.size", enexp.size)
            split_idx = range(bite,enexp.size,bite)
            print("split_indices:", split_idx)
            # ONLY ENERGY-DEPENDENT QUANTITIES HAVE TO BE SPLIT
            sub_enexp = np.split(enexp,split_idx)
            sub_tmpf = np.split(tmpf,split_idx)
           #sub_prefac = np.split(prefac,split_idx) 
           #sub_akb = np.split(akb,split_idx) 
           #sub_omegakb = np.split(omegakb,split_idx) 
           #sub_wkb = np.split(wkb,split_idx) 
           #sub_eqpkb = np.split(eqpkb,split_idx) 
            arglist = []
            for a,b in zip(sub_tmpf,sub_enexp):
                arglist.append((a,b,prefac,akb,omegakb,wkb,eqpkb,imkb))
           #    a = f2py_calc_spf_mpole_extinf(a,b,c,d,e,f,g,h) #,np.size(enexp),npoles)
            print("len(sub_enexp), length of chunks:", len(sub_enexp), [x.size for x in sub_enexp])
            print("len(sub_tmpf), length of chunks:", len(sub_tmpf), [x.size for x in sub_tmpf])
            # This determines the number of threads
           #pool = mp.Pool(ncpu)
           #pool.map(f2py_calc_spf_mpole_extinf,arglist)
            print(np.array(list(arglist[0])).shape)
            print(arglist[0])
            sub_tmpf[0] = f2py_calc_spf_mpole_extinf(arglist[0])
            output = mp.Queue()
            processes = [mp.Process(target = f2py_calc_spf_mpole_extinf, args = arglist[i]) for i in range(ncpu)]
            for p in processes:
                print("Starting process")
                p.start()
            for p in processes:
                print("Joining process")
                p.join()
            print("ALL GOOD SO FAR")
            results = [output.get() for p in processes]
            print(results)
            tmpf = f2py_calc_spf_mpole_extinf(tmpf,enexp,prefac,akb,omegakb,wkb,eqpkb,imkb) #,np.size(enexp),npoles)
            sfkb_c[ik,ib] = tmpf
            ftot = ftot + tmpf
    write_sftot_c(vardct, enexp, ftot)
    print(" calc_sf_c_para :: Done.")
    return enexp, ftot, sfkb_c

def write_sftot_c(vardct, enexp, ftot):
    """
    Just a small piece of code that writes two numpy arrays 
    (frequencies, spectral function) to disk.
    """
    import numpy as np
    print(" write_sf_c ::")
    minkpt = int(vardct['minkpt'])
    maxkpt = int(vardct['maxkpt'])
    #nkpt = maxkpt - minkpt + 1
    minband = int(vardct['minband'])
    maxband = int(vardct['maxband'])
    bdgw = map(int, vardct['sig_bdgw'])
    bdrange = range(minband-bdgw[0],maxband-bdgw[0]+1)
    kptrange = range(minkpt - 1, maxkpt)
    newdx = 0.005
    npoles = int(vardct['npoles'])
    extinf = int(vardct['extinf'])
    sfac = vardct['sfactor']
    pfac = vardct['pfactor']
    penergy = int(vardct['penergy'])
    if extinf == 1:
        outname = "spftot_exp"+"_kpt_"+str(minkpt)+"_"+str(maxkpt)+"_bd_"+str(minband)+"_"+str(maxband)+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev_np"+str(npoles)+"_extinf.dat"
    else: # extinf == 0
        outname = "spftot_exp"+"_kpt_"+str(minkpt)+"_"+str(maxkpt)+"_bd_"+str(minband)+"_"+str(maxband)+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev_np"+str(npoles)+".dat"
    outfile = open(outname,'w')
    with open(outname,'w') as outfile:
        outfile.write("# kpt "+str(minkpt)+" "+str(maxkpt)+"\n")
        outfile.write("# bd  "+str(minband)+" "+str(maxband)+"\n")
        for i in range(enexp.size):
            outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot[i])) # Dump string representations of arrays
    print(" write_sf_c :: Done.")

def array_doublefill(en_in):
    """
    Just doubles the array size and fills the gaps 
    with linear interpolation. 
    """
    import numpy as np
    en_out = np.zeros((2*en_in.size-1))
    j = 0
    for i in range(en_in.size - 1):
        en_out[j] = en_in[i]
        j += 1
        en_out[j] = (en_in[i+1] + en_in[i]) / 2
        j += 1
    i += 1
    en_out[j] = en_in[i]
    return en_out

def calc_B_crc(vardct, eqp, newen, allkb):
    """
    This method calculates the B coefficient for each k point
    for the CRC spectral function, using the multipole fit 
    for the greater part of the self-energy. 
    It has to include wtk to take into account for the
    degeneracy of k points. 
    """
    print(" calc_B_crc :: ")
    import numpy as np;
    wtk = np.array(vardct['wtk'])
    minkpt = int(vardct['minkpt'])
    maxkpt = int(vardct['maxkpt'])
    nkpt = maxkpt - minkpt + 1
    minband = int(vardct['minband'])
    maxband = int(vardct['maxband'])
    nband = maxband - minband + 1
    bdgw = map(int, vardct['sig_bdgw'])
    bdrange = range(minband-bdgw[0],maxband-bdgw[0]+1)
    kptrange = range(minkpt - 1, maxkpt)
    npoles = int(vardct['np_crc'])
    imskb = allkb[3]
    B_crc_kb =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,npoles))
    from multipole import fit_multipole, fit_multipole, getdata_file #, write_f_as_sum_of_poles
    print(" ### ================== ###")
    print(" ###    Multipole fit   ###")
    print(" Number of poles:", npoles)
   #omegampole =  np.zeros((nkpt,nband,npoles))
   #ampole =  np.zeros((nkpt,nband,npoles))
    omegampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,npoles))
    ampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,npoles))
    #for ik in range(nkpt):
    #    ikeff=minkpt+ik-1
    #bdrange = vardct['bdrange']
    #kptrange = vardct['kptrange']
    #print("kptrange, bdrange ", kptrange, bdrange)
    for ik in kptrange:
        for ib in bdrange:
   #for ik in range(imskb[:,0,0].size):
        #for ib in range(nband):
   #    for ib in range(imskb[0,:,0].size):
            if eqp[ik,ib] < newen[npoles]:
            #if eqp[ik,ib] > newen[-1]:
                omegampole[ik,ib] = omegampole[ik,ib-1]
                ampole[ik,ib] = ampole[ik,ib-1]
                print(" Eqp beyond available energy range. Values from lower band are taken.")
                continue
            else:
                print(" ik, ib", ik, ib)
                #interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
                #print(newen.shape, imskb.shape)
                interpims = interp1d(newen, imskb[ik,ib], kind = 'linear', axis = -1)
                # Here we take the curve starting from efermi and then we invert it
                # so as to have it defined on the positive x axis
                # and so that the positive direction is in the 
                # increasing direction of the array index
                #if eqp[ik,ib] <= efermi:
                if eqp[ik,ib] <= 0:
                    #en3 = en[en<=eqp[ik,ib]] # So as to avoid negative omegampole
                    en3 = newen[newen >= 0] # So as to avoid negative omegampole and ambiguity btween eqp and efermi
                else:
                    en3 = newen[newen<0] # So as to avoid negative omegampole
                    #en3 = en[en>eqp[ik,ib]] # So as to avoid negative omegampole
                #en3 = en[en<=efermi]
                if en3.size == 0:
                    print()
                    print(" WARNING: QP energy is outside of given energy range!\n"+\
                            " This state will be skipped!\n"+\
                            "You might want to modify enmin/enmax.")
                    print(" eqp[ik,ib], newen[-1]", eqp[ik,ib] , newen[-1])
                    continue
                im3 = abs(interpims(en3)/np.pi) # This is what should be fitted
                en3 = en3 - eqp[ik,ib]
               #if eqp[ik,ib] <= 0:
               #    en3 = -en3[::-1] 
               #    im3 = im3[::-1]
               #### TESTING ###
               #print("ik, ib, eqp[ik,ib], en3[0], en3[-1], newen[0], newen[-1]:\n", ik, ib, eqp[ik,ib], en3[0], en3[-1], newen[0], newen[-1])
               #import matplotlib.pylab as plt
               #plt.plot(newen, imskb[ik,ib]/np.pi,"-")
               #plt.plot(en3+eqp[ik,ib], im3,"x")
               #plt.plot(en3, im3,"o")
               #plt.show()
               #sys.exit()
               #### END TESTING ###
                omegai, lambdai, deltai = fit_multipole(en3,im3,npoles)
                # HERE WE MUST CHECK THAT THE NUMBER OF POLES 
                # IS NOT BIGGER THAN THE NUMBER OF POINTS THAT HAS TO BE FITTED
                if npoles > omegai.size:
                    omegampole[ik,ib][:omegai.size] = omegai 
                    ampole[ik,ib][:omegai.size] = np.true_divide(lambdai,(np.square(omegai)))
                    print()
                    print(" WARNING: npoles used ("+str(npoles)+") is larger"+\
                            " than poles x data array can give ("+str(omegai.size)+").")
                   #print("WARNING: Reduce npoles. You are wasting resources!!!")
                    print(" Im(Sigma) will be interpolated to obtain the desired number of poles.")
                    current_size = omegai.size
                    counter = 0
                    while npoles > current_size:
                        counter += 1
                        print()
                        print(" WARNING: Arrays are too coarse.")
                        print(" npoles, omegai.size:", npoles, omegai.size)
                        print(" Filling arrays with interpolated values...")
                        en1 = array_doublefill(en3)
                        im1 = array_doublefill(im3)
                        en3 = en1
                        im3 = im1
                        omegai, lambdai, deltai = fit_multipole(en1,im1,npoles)
                        current_size = omegai.size
                        if counter > 4:
                            print(60*"=")
                            print(" WARNING: You are trying too hard with too few points.")
                            print(" The array has been interpolated more than 4 times.")
                            print(" Maybe use less poles or calculate more points for Sigma?")
                            print(60*"=")
            #   im1 = fit_double(im3)
                else:
                    omegampole[ik,ib] = omegai 
                    ampole[ik,ib] = np.true_divide(lambdai,(np.square(omegai)))
                B_crc_kb[ik,ib] = ampole[ik,ib]
               #ampole[ik,ib] = gi
                print(" Integral test. Compare \int\Sigma and \sum_j^N\lambda_j.")
                print(" 1/pi*\int\Sigma   =", np.trapz(im3,en3))
                print(" \sum_j^N\lambda_j =", np.sum(lambdai))
                print("b_j:", ampole[ik,ib])
                #plt.plot(en3,im3,"-"); plt.plot(omegai,np.pi/2*gi*omegai/deltai,"-o")
                #e1,f1 = write_f_as_sum_of_poles(en3,omegai,gi,deltai,0)
    # Writing out a_j e omega_j
    print(" ### Writing out a_j and omega_j...")
    outname = "a_j_np_crc"+str(npoles)+".dat"
    outfile = open(outname,'w')
    outname = "omega_j_np_crc"+str(npoles)+".dat"
    outfile2 = open(outname,'w')
    for ipole in xrange(npoles):
 #      for ik in kptrange:
 #          #for ib in range(nband):
 #          for ib in bdrange:
        for ik in range(imskb[:,0,0].size):
            for ib in range(imskb[0,:,0].size):
                outfile.write("%10.5f"  % (ampole[ik,ib,ipole]))
                outfile2.write("%10.5f" % (omegampole[ik,ib,ipole]))
            outfile.write("\n")
            outfile2.write("\n")
        outfile.write("\n")
        outfile2.write("\n")
    outfile.close()
    outfile2.close()
    return B_crc_kb

def calc_sf_crc(dict_c, B_crc_kb, hartree, newen, allkb):
    """
    Calculation of the CRC part of the spectral function and of the
    total CRC spectral function. 
    """
    print(" calc_sf_c_serial :: ")
    import numpy as np;
    from extmod_spf_mpole import f2py_calc_spf_mpole
    wtk = np.array(vardct['wtk'])
    hartree = np.array(hartree)
    pdos = np.array(pdos)
    minkpt = int(vardct['minkpt'])
    maxkpt = int(vardct['maxkpt'])
    nkpt = maxkpt - minkpt + 1
    minband = int(vardct['minband'])
    maxband = int(vardct['maxband'])
    nband = maxband - minband + 1
    bdgw = map(int, vardct['sig_bdgw'])
    bdrange = range(minband-bdgw[0],maxband-bdgw[0]+1)
    kptrange = range(minkpt - 1, maxkpt)
    newdx = 0.005
    enmin = float(vardct['enmin'])
    enmax = float(vardct['enmax'])
    npoles = int(vardct['npoles'])
    extinf = int(vardct['extinf'])
    penergy = int(vardct['penergy'])
    for ik in kptrange:
        ikeff = ik + 1
        for ib in bdrange:
            ibeff = ib + 1
            print(" ik, ib, ikeff, ibeff", ik, ib, ikeff, ibeff)
            #prefac=np.exp(-np.sum(ampole[ik,ib]))/np.pi*wtk[ik]*pdos[ib]*abs(imeqp[ik,ib])
            # Experimental fix for npoles dependence
            tmp = 1/np.pi*wtk[ik]*pdos[ib]*abs(imeqp[ik,ib])
            prefac=np.exp(-np.sum(ampole[ik,ib]))*tmp
            #prefac=np.exp(-tmp*np.trapz(imskb[ik,ib],enexp)/np.sum(omegai)*npoles)
           #print("\n === Normalization test === ")
           #print(" Prefactor:", np.exp(-np.sum(ampole[ik,ib])))
           #print(" Exponent:", np.sum(ampole[ik,ib]))
           #print(" Exponent/npoles:", np.sum(ampole[ik,ib])/npoles,end="\n\n")
            akb=ampole[ik,ib] # This is a numpy array (slice)
            omegakb=omegampole[ik,ib] # This is a numpy array (slice)
            eqpkb=eqp[ik,ib]
            imkb=imeqp[ik,ib]
            #tmpf1 = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
            #print(nen, np.size(enexp))
            #tmpf = 0.0*tmpf
            if eqpkb < 0.0:
                pass
            else:
                print(" This state is empty! eqpkb ik ib:",eqpkb, ikeff+1, ibeff+1)
                #print("omegakb", omegakb)
                omegakb=-omegakb
                #print("-omegakb", omegakb)
            tmpf = np.zeros((nenexp), order='Fortran')
            tmpf = f2py_calc_spf_mpole(tmpf,enexp,prefac,akb,omegakb,eqpkb,imkb) #,nen,npoles)
                #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
            #outnamekb = "spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"."+str(penergy)
            #outfilekb = open(outnamekb,'w')
            #for ien in xrange(nenexp):
            #    outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], tmpf[ien]))
            #outfilekb.close()
            sfkb_c[ik,ib] = tmpf
            ftot = ftot + tmpf
        return ftot, sfkb_c

def calc_sf_c_serial(vardct, hartree, pdos, eqp, imeqp, newen, allkb):
    """
    This method takes care of the calculation of the cumulant. 
    Different values of npoles command different options:
    - npoles = 999: single-pole calculation with omega_p
    used as the plasmon frequency for every state.
    - npoles = 0: QP-only calculation. Z = 1 and satellite 
    weights are put to 0. 
    - Standard cumulant for any other value of npoles.
    """
    print(" calc_sf_c_serial :: ")
    import numpy as np;
    wtk = np.array(vardct['wtk'])
    hartree = np.array(hartree)
    pdos = np.array(pdos)
    minkpt = int(vardct['minkpt'])
    maxkpt = int(vardct['maxkpt'])
    nkpt = maxkpt - minkpt + 1
    minband = int(vardct['minband'])
    maxband = int(vardct['maxband'])
    nband = maxband - minband + 1
    bdgw = map(int, vardct['sig_bdgw'])
    bdrange = range(minband-bdgw[0],maxband-bdgw[0]+1)
    kptrange = range(minkpt - 1, maxkpt)
   #print("kptrange, bdrange ", kptrange, bdrange)
    newdx = 0.005
    enmin = float(vardct['enmin'])
    enmax = float(vardct['enmax'])
    #if enmin < en[0] and enmax >= en[-1]:  
    #    newen = np.arange(en[0],en[-1],newdx)
    #elif enmin < en[0]:  
    #    newen = np.arange(en[0],enmax,newdx)
    #elif enmax >= en[-1] :  
    #    newen = np.arange(enmin,en[-1],newdx)
    #else :  
    #    newen = np.arange(enmin,enmax,newdx)
    npoles = int(vardct['npoles'])
    extinf = int(vardct['extinf'])
    penergy = int(vardct['penergy'])
    #allkb = [spfkb,reskb, rdenkb, imskb]
    reskb = allkb[1]
    imskb = allkb[3]
    if npoles==999: # same omega_p for every state, with the intensity calculated integrating Im(Sigma)
        omega_p = float(vardct['omega_p'])
       #omegampole = np.ones((nkpt,nband))*omega_p
       #ampole =  np.zeros((nkpt,nband))
        omegampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size))*omega_p
        ampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size))
        #for ik in range(nkpt):
        #for ik in kptrange:
           #for ib in range(nband):
           #for ib in bdrange:
        for ik in range(imskb[:,0,0].size):
            for ib in range(imskb[0,:,0].size):
                print(" ik, ib", ik, ib)
                #interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
                #if eqp[ik,ib]<=efermi:
                if eqp[ik,ib]<=0:
                    tmpen = newen[imskb[ik,ib]>=0]
                    tmpim = imskb[ik,ib,imskb[ik,ib]>=0]
                else:
                    tmpen = newen[imskb[ik,ib]<0]
                    tmpim = imskb[ik,ib,ims[ik,ib]<0]
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
        from multipole import fit_multipole, fit_multipole, getdata_file #, write_f_as_sum_of_poles
        print(" ### ================== ###")
        print(" ###    Multipole fit   ###")
        print(" Number of poles:", npoles)
       #omegampole =  np.zeros((nkpt,nband,npoles))
       #ampole =  np.zeros((nkpt,nband,npoles))
        omegampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,npoles))
        ampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,npoles))
        #for ik in range(nkpt):
        #    ikeff=minkpt+ik-1
        #bdrange = vardct['bdrange']
        #kptrange = vardct['kptrange']
        #print("kptrange, bdrange ", kptrange, bdrange)
        for ik in kptrange:
            for ib in bdrange:
       #for ik in range(imskb[:,0,0].size):
            #for ib in range(nband):
       #    for ib in range(imskb[0,:,0].size):
                if eqp[ik,ib] > newen[-npoles]:
                #if eqp[ik,ib] > newen[-1]:
                    omegampole[ik,ib] = omegampole[ik,ib-1]
                    ampole[ik,ib] = ampole[ik,ib-1]
                    print(" Eqp beyond available energy range. Values from lower band are taken.")
                    continue
                else:
                    ibeff=minband+ib-1
                    print(" ik, ib", ik, ib)
                    #interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis = -1)
                    #print(newen.shape, imskb.shape)
                    interpims = interp1d(newen, imskb[ik,ib], kind = 'linear', axis = -1)
                    # Here we take the curve starting from eqp and then we invert it
                    # so as to have it defined on the positive x axis
                    # and so that the positive direction is in the 
                    # increasing direction of the array index
                    #if eqp[ik,ib] <= efermi:
                    if eqp[ik,ib] <= 0:
                        #en3 = en[en<=eqp[ik,ib]] # So as to avoid negative omegampole
                        en3 = newen[newen<=eqp[ik,ib]] # So as to avoid negative omegampole
                       #en3 = newen[newen<0.] # So as to avoid negative omegampole
                    else:
                        en3 = newen[newen>eqp[ik,ib]] # So as to avoid negative omegampole
                        #en3 = en[en>eqp[ik,ib]] # So as to avoid negative omegampole
                    #en3 = en[en<=efermi]
                    if en3.size == 0:
                        print()
                        print(" WARNING: QP energy is outside of given energy range!\n"+\
                                " This state will be skipped!\n"+\
                                "You might want to modify enmin/enmax.")
                        print(" eqp[ik,ib], newen[-1]", eqp[ik,ib] , newen[-1])
                        continue
                    im3 = abs(interpims(en3)/np.pi) # This is what should be fitted
                   #zcut = 3.0
                   #for i in range(en3.size):
                   #    if en3[i]>(eqp[ik,ib]-zcut) and en3[i]<(eqp[ik,ib]+zcut):
                   #        im3[i] = 0.
                   #import matplotlib.pylab as plt
                   #plt.plot(en3,im3,'-')
                   #plt.show()
                    en3 = en3 - eqp[ik,ib]
                    if eqp[ik,ib] <= 0:
                        en3 = -en3[::-1] 
                        im3 = im3[::-1]
                   #### TESTING ###
                   #print("ik, ib, eqp[ik,ib], en3[0], en3[-1], newen[0], newen[-1]:\n", ik, ib, eqp[ik,ib], en3[0], en3[-1], newen[0], newen[-1])
                   #import matplotlib.pylab as plt
                   #plt.plot(newen, imskb[ik,ib]/np.pi,"-")
                   #plt.plot(en3+eqp[ik,ib], im3,"x")
                   #plt.show()
                   #sys.exit()
                   #### END TESTING ###
                    omegai, lambdai, deltai = fit_multipole(en3,im3,npoles)
                    plot_fit = int(vardct['plot_fit'])
                    if plot_fit == 1:
                        from multipole import write_f_as_sum_of_poles
                        import matplotlib.pylab as plt
                        import pylab
                        plt.figure(2)
                        eta = 0.5
                        enlor, flor = write_f_as_sum_of_poles(en3, omegai, lambdai, deltai, eta)
                        plt.plot(enlor, flor,"-",label="sum of poles, eta: "+str(eta))
                        plt.plot(en3,im3,"-",label="ImS(e-w)")
                        plt.plot(omegai,lambdai,"go", label = "omegai, lambdai")
                        plt.plot(omegai,lambdai/deltai,"ro", label = "omegai, lambdai/deltai")
                        plt.title("ik: "+str(ik)+", ib: "+str(ib)+", npoles: "+str(npoles))
                        plt.legend()
                        pylab.savefig('imS_fit_np'+str(npoles)+'_ik'+str(ik)+'_ib'+str(ib)+'.pdf')
                        plt.close()
                   ## TESTING THE MULTIPOLE REPRESENTATION
                   #from multipole import write_f_as_sum_of_poles
                   #import matplotlib.pylab as plt
                   #import pylab
                   #eta = 0.01
                   #for eta in [0.1]: #, 0.1, 0.5]:
                   #    for npoles in [1,10,20,100]:
                   #        omegai, lambdai, deltai = fit_multipole_const(en3,im3,npoles)
                   #        print("ik, ib, eqp[ik,ib], en3[0], en3[-1], newen[0], newen[-1]:\n", ik, ib, eqp[ik,ib], en3[0], en3[-1], newen[0], newen[-1])
                   #        print(omegai, lambdai, deltai)
                   #        enlor, flor = write_f_as_sum_of_poles(en3, omegai, lambdai, deltai, eta)
                   #        plt.plot(enlor, flor,"-",label="sum of poles, eta: "+str(eta))
                   #        plt.plot(en3,im3,"-",label="ImS(e-w)")
                   #        plt.plot(omegai,lambdai,"go", label = "omegai, lambdai")
                   #        plt.plot(omegai,lambdai/deltai,"ro", label = "omegai, lambdai/deltai")
                   #        plt.title("ik: "+str(ik)+", ib: "+str(ib)+", npoles: "+str(npoles))
                   #        plt.legend()
                   #        pylab.savefig('imS_test_np'+str(npoles)+'_ik'+str(ik)+'_ib'+str(ib)+'_eta'+str(eta)+'.pdf')
                   #        plt.show()
                   #sys.exit()
                    # END TESTING THE MULTIPOLE REPRESENTATION 
                    # HERE WE MUST CHECK THAT THE NUMBER OF POLES 
                    # IS NOT BIGGER THAN THE NUMBER OF POINTS THAT HAS TO BE FITTED
                    if npoles > omegai.size:
                        omegampole[ik,ib][:omegai.size] = omegai 
                        ampole[ik,ib][:omegai.size] = np.true_divide(lambdai,(np.square(omegai)))
                        print()
                        print(" WARNING: npoles used ("+str(npoles)+") is larger"+\
                                " than poles x data array can give ("+str(omegai.size)+").")
                       #print("WARNING: Reduce npoles. You are wasting resources!!!")
                        print(" Im(Sigma) will be interpolated to obtain the desired number of poles.")
                        current_size = omegai.size
                        counter = 0
                        while npoles > current_size:
                            counter += 1
                            print()
                            print(" WARNING: Arrays are too coarse.")
                            print(" npoles, omegai.size:", npoles, omegai.size)
                            print(" Filling arrays with interpolated values...")
                            en1 = array_doublefill(en3)
                            im1 = array_doublefill(im3)
                            en3 = en1
                            im3 = im1
                            omegai, lambdai, deltai = fit_multipole(en1,im1,npoles)
                            current_size = omegai.size
                            if counter > 4:
                                print(60*"=")
                                print(" WARNING: You are trying too hard with too few points.")
                                print(" The array has been interpolated more than 4 times.")
                                print(" Maybe use less poles or calculate more points for Sigma?")
                                print(60*"=")
                #   im1 = fit_double(im3)
                    else:
                        omegampole[ik,ib] = omegai 
                        ampole[ik,ib] = np.true_divide(lambdai,(np.square(omegai)))
                   #ampole[ik,ib] = gi
                    print(" Integral test. Compare \int\Sigma and \sum_j^N\lambda_j.")
                    print(" 1/pi*\int\Sigma   =", np.trapz(im3,en3))
                    print(" \sum_j^N\lambda_j =", np.sum(lambdai))
                    #plt.plot(en3,im3,"-"); plt.plot(omegai,np.pi/2*gi*omegai/deltai,"-o")
                    #e1,f1 = write_f_as_sum_of_poles(en3,omegai,gi,deltai,0)
        # Writing out a_j e omega_j
        print(" ### Writing out a_j and omega_j...")
        outname = "a_j_np"+str(npoles)+".dat"
        outfile = open(outname,'w')
        outname2 = "omega_j_np"+str(npoles)+".dat"
        outfile2 = open(outname2,'w')
        for ipole in xrange(npoles):
     #      for ik in kptrange:
     #          #for ib in range(nband):
     #          for ib in bdrange:
            for ik in range(imskb[:,0,0].size):
                for ib in range(imskb[0,:,0].size):
                    outfile.write("%15.7e"  % (ampole[ik,ib,ipole]))
                    outfile2.write("%15.7e" % (omegampole[ik,ib,ipole]))
                   #outfile.write("%10.5f"  % (ampole[ik,ib,ipole]))
                   #outfile2.write("%10.5f" % (omegampole[ik,ib,ipole]))
                outfile.write("\n")
                outfile2.write("\n")
            outfile.write("\n")
            outfile2.write("\n")
        outfile.close()
        outfile2.close()
        # Extrinsic and interference contribution
        if extinf == 1:
            origdir = vardct['origdir']
            extinfname = "a_wp."+str(penergy)
            amp_exinf, w_extinf = calc_extinf_corrections(origdir,extinfname,ampole,omegampole)
            print(" ### Writing out a_j_extinf...")
            outname = "a_j_np"+str(npoles)+"_extinf."+str(penergy)
            outfile = open(outname,'w')
            for ipole in xrange(npoles):
           #    for ik in kptrange:
           #        for ib in bdrange:
                for ik in range(imskb[:,0,0].size):
                    for ib in range(imskb[0,:,0].size):
                        outfile.write("%10.5f"  % (amp_exinf[ik,ib,ipole]))
                    outfile.write("\n")
                outfile.write("\n")
            outfile.close()
    else: # npoles == 0
        omegampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size))
        ampole =  np.zeros((imskb[:,0,0].size,imskb[0,:,0].size))
       #omegampole =  np.zeros((nkpt,nband))
       #ampole =  np.zeros((nkpt,nband))
    #elaps2 = time.time() - elaps1 - e0
    #cpu2 = time.clock() - cpu1 - c0
    #print(elaps2, cpu2)
    #print(str(" Used time (elaps, cpu): %10.6e %10.6e"% (elaps2, cpu2)))
    print(" Calculating multipole exponential A...")
    dxexp=0.005 
    enexp = np.arange(enmin,enmax,dxexp)
    nenexp = np.size(enexp)
    ftot = np.zeros((np.size(enexp)),order='Fortran')
    nen = np.size(enexp)
    #sfkb_c = np.zeros((nkpt,nband,nenexp))
    sfkb_c = np.zeros((imskb[:,0,0].size,imskb[0,:,0].size,nenexp))
    ############################
    # With extrinsic effects ###
    if extinf == 1:
        from extmod_spf_mpole import f2py_calc_spf_mpole_extinf
        #for ik in range(nkpt):
        for ik in kptrange:
            ikeff = ik + 1
            #for ib in range(nband):
            for ib in bdrange:
                ibeff=bdgw[0]+ib
                print(" ik, ib, ikeff, ibeff", ik, ib, ikeff, ibeff)
               #prefac=np.exp(-np.sum(amp_exinf[ik,ib]))/np.pi*wtk[ik]*pdos[ib]*abs(imeqp[ik,ib])
                # Experimental fix for npoles dependence
                tmp = 1/np.pi*wtk[ik]*pdos[ib]*abs(imeqp[ik,ib])
                prefac=np.exp(-np.sum(ampole[ik,ib]))*tmp
               #prefac=np.exp(-tmp*np.trapz(imskb[ik,ib],enexp)/np.sum(omegai)*npoles)
                akb=amp_exinf[ik,ib] # This is a numpy array (slice)
                omegakb=omegampole[ik,ib] # This is a numpy array (slice)
                wkb=w_extinf[ik,ib] # This is a numpy array (slice)
                eqpkb=eqp[ik,ib]
                imkb=imeqp[ik,ib] # + w_extinf[ik,ib]/2 # extinf width added
                #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles,wkb)
                #ftot += tmpf
                if eqpkb < 0.0:
                    pass
                    #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                else:
                    print(" This state is empty! eqpkb ik ib:",eqpkb, ikeff+1, ibeff+1)
                    #print("omegakb", omegakb)
                    omegakb=-omegakb
                    #print("-omegakb", omegakb)
                tmpf = np.zeros((nenexp), order='Fortran')
                tmpf = f2py_calc_spf_mpole_extinf(tmpf,enexp,prefac,akb,omegakb,wkb,eqpkb,imkb) #,np.size(enexp),npoles)
               #outnamekb = "spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"_extinf."+str(penergy)
               #outfilekb = open(outnamekb,'w')
               #for ien in xrange(nenexp):
               #    outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], tmpf[ien]))
               #outfilekb.close()
                sfkb_c[ik,ib] = tmpf
                ftot = ftot + tmpf
    else: # extinf == 0
        from extmod_spf_mpole import f2py_calc_spf_mpole
        #for ik in range(nkpt):
            #for ib in range(nband):
        for ik in kptrange:
            ikeff = ik + 1
            for ib in bdrange:
                ibeff = ib + 1
                print(" ik, ib, ikeff, ibeff", ik, ib, ikeff, ibeff)
                #prefac=np.exp(-np.sum(ampole[ik,ib]))/np.pi*wtk[ik]*pdos[ib]*abs(imeqp[ik,ib])
                # Experimental fix for npoles dependence
                tmp = 1/np.pi*wtk[ik]*pdos[ib]*abs(imeqp[ik,ib])
                prefac=np.exp(-np.sum(ampole[ik,ib]))*tmp
                #prefac=np.exp(-tmp*np.trapz(imskb[ik,ib],enexp)/np.sum(omegai)*npoles)
                print("\n === Normalization test === ")
                print(" Prefactor:", np.exp(-np.sum(ampole[ik,ib])))
                print(" Exponent:", np.sum(ampole[ik,ib]))
                print(" Exponent/npoles:", np.sum(ampole[ik,ib])/npoles,end="\n\n")
                akb=ampole[ik,ib] # This is a numpy array (slice)
                omegakb=omegampole[ik,ib] # This is a numpy array (slice)
                eqpkb=eqp[ik,ib]
                imkb=imeqp[ik,ib]
                #tmpf1 = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                #print(nen, np.size(enexp))
                #tmpf = 0.0*tmpf
                if eqpkb < 0.0:
                    pass
                else:
                    print(" This state is empty! eqpkb ik ib:",eqpkb, ikeff+1, ibeff+1)
                    #print("omegakb", omegakb)
                    omegakb=-omegakb
                    #print("-omegakb", omegakb)
                tmpf = np.zeros((nenexp), order='Fortran')
                tmpf = f2py_calc_spf_mpole(tmpf,enexp,prefac,akb,omegakb,eqpkb,imkb) #,nen,npoles)
                    #tmpf = calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles)
                #outnamekb = "spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"."+str(penergy)
                #outfilekb = open(outnamekb,'w')
                #for ien in xrange(nenexp):
                #    outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], tmpf[ien]))
                #outfilekb.close()
                sfkb_c[ik,ib] = tmpf
                ftot = ftot + tmpf
                #print(ftot[0], tmpf[0])
    #elaps2 = time.time() - elaps1 - e0
    #cpu2 = time.clock() - cpu1 - c0
    #print(elaps2, cpu2)
    #print(str(" Used time (elaps, cpu): %10.6e %10.6e"% (elaps2, cpu2)))
    #print(" ### Writing out A(\omega)_exp...  ")
    #enexp = enexp-efermi
    write_sftot_c(vardct, enexp, ftot)
    print(" calc_sf_c_serial :: Done.")
    return enexp, ftot, sfkb_c

