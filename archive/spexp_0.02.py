#!/usr/bin/env python
"""
Written by Matteo Guzzo.
Last revision date: 03/12/2011
List of files needed:
- invar.in with input variables.
- _SIG for the self-energy.
- s.dat, p_even.dat, p_odd.dat, d_even.dat, etc. for the orbital character.
- cs*.dat for the photon cross sections.
- hartree.dat or elda.dat and vxc.dat for the hartree energies.
- wtk.dat for the k-points weights.
"""
import numpy as np;
import matplotlib.pylab as plt;
from scipy.interpolate import interp1d
from scipy import optimize
import sys
from os.path import isfile, join, isdir
from os import getcwd, pardir, mkdir, chdir
#
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
# Reading Hartree energies
# TODO: write a function to read the parameters from not ad-hoc files
	if isfile("hartree.dat"):
		print " Reading file hartree.dat... ",
		hartreefile = open("hartree.dat");
		hartree = [];
		for line in hartreefile.readlines():
			hartree.append(map(float,line.split()));
		hartreefile.close()
		print "Done."
		hartree = np.array(hartree);
	elif isfile("E_lda.dat") and isfile("Vxc.dat"):
		print " Auxiliary file (hartree.dat) not found."
		print " Reading files E_lda.dat and Vxc.dat... ",
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
		print "Done."
		elda = np.array(elda);
		vxc = np.array(vxc);
		hartree = elda - vxc
	else : 
		print " Auxiliary file not found (hartree/E_lda/Vxc). Impossible to continue."
		sys.exit(1)
	return hartree
hartree = read_hartree()
# End read_hartree()
#
def read_wtk():
	"""
	This function takes the file 'wtk.dat'
	and creates an array containing the 
	values of the k-point weights for each state. 
	This array is returned.
	The input file is supposed to be a single column of
	nkpt elements. 
	"""
	if isfile("wtk.dat"):
		wtkfile = open("wtk.dat");
	else : 
		print " Auxiliary file not found (wtk.dat). Impossible to continue."
		sys.exit(1)
	wtk = [];
	for line in wtkfile.readlines():
		wtk.append((float(line)));
	wtkfile.close()
	wtk = np.array(wtk);
	return wtk

def read_sigfile(nkpt,nband,sigfilename):
	if isfile(sigfilename):
		insigfile = open(sigfilename);
	else:
		print "File "+sigfilename+" not found."
		insigfile = open(raw_input("Self-energy file name (_SIG): "))
	# Indexes initialization
	ikpt = 0;
	ib = 0;
	# We put the content of the file (lines) in this array
	filelines = insigfile.readlines();
	en=[];
	istart = 0;
	# loop to prepare  the energy array
	print " Reading array of energies from first k-point in _SIG file... ",
	for line in filelines:
		if line[0:3]=='# k':
			continue
		elif line[0:3]=='# b':
			istart = 1 # Supposed to be 0 before
			continue
		else : 
			data=line.split()
			en.append(float(data[0]))
			if float(data[0]) >= enmax : break
	print "Done."
	print " Length of the energy array detected from _SIG file, first k point: "+str(len(en))
	print " len(en): "+str(len(en))
	en = np.array(en)
	print " size(en): "+str(np.size(en))
	dx = (en[-1]-en[0])/np.size(en)
	print " dx:",dx
	# Preparation of arrays
	res=np.zeros((nkpt,nband,np.size(en)));
	ims=np.zeros((nkpt,nband,np.size(en)));
	print " np.shape(res), np.shape(ims):", np.shape(res), np.shape(ims)
	print " ### ================================= #### "
	print " ### ===== Reading _SIG file ... ===== #### "
	print " ### ================================= #### "
	# Cycle over the lines of the file
	io = 0 # row/energy counter
	for line in filelines:
		icol = 0;
		ib = 0;
		if line[0:3]=='# k' :
			# Detect, in commented lines, the k point declaration
			if ikpt==nkpt :
				print " End of the reading loop: ikpt == nkpt. ikpt,nkpt: ", ikpt, nkpt; 
				print " #### ================================================ #### ";
				break;
			print " --- k point:  %02i ---" % (ikpt+1);
			ikpt = ikpt + 1;
		# Detect, in commented lines, the bands declaration
		elif line[0:3]=='# b':
			io = 0;
		# TODO: This if test is incorrect: this way it always starts from ib = 0
		elif io < np.size(en) and ib < nband :
			data=map(float,line.split());
			for col in data:
				# First element (energy) goes into en
				if icol == 0 : 
					icol = icol + 1;
					io = io + 1;
				elif (icol+2)%3 == 0 :
					res[ikpt-1,ib,io-1] = col;
					icol = icol + 1;
					continue;
				elif (icol+1)%3 == 0 :
					ims[ikpt-1,ib,io-1] = col;
					icol = icol + 1;
					continue;
				else : 
					if ib == nband-1 : break # number of bands reached
					ib = ib + 1;
					icol = icol + 1;
					continue;
			else: continue
	return en, res, ims

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
	print " ### Reading cross sections...  "
	csfilename = "cs"+str(penergy)+".dat"
	if isfile(csfilename):
		print " Photon energy:", penergy,"eV"
	else:
		penergy = raw_input(" File "+csfilename+" not found. Photon energy (eV): ")
		csfilename = "cs"+str(penergy)+".dat"
	cs = []
	print " csfilename:",csfilename
	csfile = open(csfilename,'r')
	for line in csfile.readlines():
		cs.append((float(line)));
	csfile.close()
	cs = np.array(cs)
	#print "cs:",cs.shape,cs
	#print "cs:",np.transpose(cs),cs.shape
	return cs

def read_band_type_sym(sfac,pfac):
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
	print " Reading s bands file... ",
	if isfile("s.dat"):
		sfile =  open("s.dat",'r')
		s = map(float,sfile.read().split())
		sfile.close()
		s = np.array(s)
		print "Done."
	else : 
		print
		print " WARNING: File for orbital character not found (s.dat). S character will be 1 for all bands. "
		s = np.ones(nband)
	print " Reading p bands file... ",
	if isfile("p_even.dat") and isfile("p_odd.dat"):
		# This file should contain the summed contribution of all even p orbitals
		pevenfile =  open("p_even.dat",'r')
		# This file should contain the summed contribution of all odd p orbitals
		poddfile =  open("p_odd.dat",'r')
		peven = map(float,pevenfile.read().split())
		podd = map(float,poddfile.read().split())
		pevenfile.close()
		poddfile.close()
		peven = np.array(peven)
		podd = np.array(podd)
		print "Done."
	else : 
		print
		print " WARNING: File for orbital character not found (p_even.dat/p_odd.dat). P character will be 1 for all bands. "
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
	#print "sp:",sp
	return sp

def calc_spf_gw(nkpt,nband,wtk,pdos,en,res,ims):
	"""
	Macro-function calling instructions necessary to calculate 
	the GW spectral function. 
	For now it writes out the single state spectral functions
	on ascii files. This should be moved to an external module
	and just return spfkb as an output variable. 
	spf (GW spectral function) is returned.
	"""
	# ===========INTERPOLATION ============ #
	# Interpolation section
	#print " ### ============= #### "
	#print " ### Interpolation #### "
	#print " ### ============= #### "
	# INTERPOLATION GRID DEFINED HERE #
	newdx = 0.01
	if enmin < en[0] :  
		newen = np.arange(en[0],enmax,newdx)
	else :  
		newen = np.arange(enmin,enmax,newdx)
	print " ### Interpolation and calculation of A(\omega)_GW...  "
	spftot = np.zeros((np.size(newen)));
	# Here we interpolate re and im sigma
	# for each band and k point
	for ik in xrange(np.size(ims[:,0,0])):
		print " k point = %02d " % (ik+1)
		for ib in xrange(np.size(ims[0,:,0])):
			#print " Outfile: ",outnamekb
			# Parti reali e immaginarie interpolate
			interpres = interp1d(en, res[ik,ib,:], kind = 'linear', axis =  2)
			interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
			redenom = newen - hartree[ik,ib] - interpres(newen)
	                tmpim = abs(interpims(newen))
	                spfkb = wtk[ik] * pdos[ib] * tmpim/np.pi/(redenom**2 + tmpim**2)
			spftot += spfkb 
			outnamekb = "spf_gw-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+".dat"
			outfilekb = open(outnamekb,'w')
			for ien in xrange(np.size(newen)) :
				outfilekb.write("%8.4f %12.8f %12.8f %12.8f\n" % (newen[ien], spfkb[ien], redenom[ien], tmpim[ien]))
			outfilekb.close()
	return newen, spftot

def find_eqp_resigma(en,resigma):
	"""
	This function is supposed to deal with the plasmaron problem 
	and calculate the quasiparticle energy once it is fed with 
	\omega - \epsilon_H - \Re\Sigma. 
	It expects an array of increasing values and it will return 
	the last 0 detected. 
	It should return the value of eqp. 
	"""
	import numpy as np
	import matplotlib.pylab as plt
	for i in xrange(1,np.size(resigma)):
		if (resigma[i]*resigma[i-1] < 0):
			tmpeqp = en[i-1] - resigma[i-1]*(en[i] - en[i-1])/(resigma[i] - resigma[i-1]) # High school formula
	return tmpeqp

def calc_eqp_imeqp(nkpt,nband,en,res,ims):
	"""
	This function calculates qp energies and corresponding
	values of the imaginary part of sigma for a set of
	k points and bands. 
	The function find_eqp_resigma() is used here.
	eqp and imeqp are returned. 
	"""
	eqp = np.zeros((nkpt,nband))
	imeqp = np.zeros((nkpt,nband))
	for ik in xrange(nkpt):
		for ib in xrange(nband):
			temparray = np.array(en[:] - hartree[ik,ib] - res[ik,ib,:])
			interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
			tempim = interpims(en[:])
			# New method to overcome plasmaron problem
			eqp[ik,ib] = find_eqp_resigma(en,temparray)
			imeqp[ik,ib] = interpims(eqp[ik,ib])
			## Warning if imaginary part of sigma < 0 (Convergence problems?)
			if imeqp[ik,ib] <= 0 : print " WARNING: im(Sigma(eps_k)) <= 0 !!! ik ib eps_k im(Sigma(eps_k)) = ", ik, ib, eqp[ik,ib], imeqp[ik,ib]
	return eqp, imeqp

### ============================= ###
###  ==  PROGRAM BEGINS HERE  ==  ###
### ============================= ###

# ======== READING INPUT VARIABLES ======= #
print " Reading invar file... ",
invar = {}
if isfile("invar.in"):
	infile = open("invar.in")
	for line in infile.readlines():
		word = line.split()
		invar[word[-1]] = word[0];
#		print "invar: ", invar
	infile.close()
	sigfilename  =       invar['sigmafile'];
	minband      =   int(invar['minband']);
	maxband      =   int(invar['maxband']);
	nkpt         =   int(invar['nkpt']);
	enmin        = float(invar['enmin']);
	enmax        = float(invar['enmax']);
	sfac         = float(invar['sfactor']);
	pfac         = float(invar['pfactor']);
	penergy      =   int(invar['penergy']);
	npoles       =   int(invar['npoles']);
	flag_calc_gw =   int(invar['calc_gw']);
	flag_calc_exp=   int(invar['calc_exp']);
	a_extinf     = float(invar['a_extinf']);
	extinf       =   int(invar['extinf']);
else : 
	print "Invar file not found (invar.in). Impossible to continue."
	sys.exit(1)
print "Done."
nband = 1 + maxband - minband;
print " minband =", minband;
print " maxband =", maxband;
print " nband =", nband;
print " nkpt =", nkpt;
print " enmin =", enmin;
print " enmax =", enmax;
print " S prefactor:", sfac
print " P prefactor:", pfac
# Max energy in spectrum
# TODO: write a function to read this parameter from a file (Fermi energy?)
#enmax = 15. # eV
# ======== READING WTK ======= #
wtk = read_wtk()
# ======== READING _SIG FILE ======= #
en, res, ims = read_sigfile(nkpt,nband,sigfilename)
print " ### nkpt, nband:", nkpt, nband
print " # ------------------------------------------------ # ";
# ======== CROSS SECTIONS ======= #
cs = read_cross_sections(penergy)
# ====== BAND TYPE AND SYMMETRY ==== #
sp = read_band_type_sym(sfac,pfac)
# ===== EFFECTIVE STATE-DEPENDENT PREFACTOR ==== #
pdos = 10000.*np.dot(cs,sp)
print " 10000*pdos:", pdos
print " Size(pdos):",np.size(pdos)
### ===================================================== ###
print " # ------------------------------------------------ # ";
# Here we move to a subdirectory to avoid flooding-up the current directory
newdirname = "Spfunctions"
origdir = getcwd() # remember where we are
newdir = join(origdir, newdirname) # Complete path of the new directory
print " Moving into spectral functions directory:\n ", newdir
if not isdir(newdir) :
	mkdir(newdir)
chdir(newdir)
### ================================= ###
### ===== GW SPECTRAL FUNCTION ====== ###
# GW spectral function part
if flag_calc_gw == 1:
	newen, spftot = calc_spf_gw(nkpt,nband,wtk,pdos,en,res,ims)
		### ==== WRITING OUT GW SPECTRAL FUNCTION === ###
	print " ### Writing out A(\omega)_GW...  "
	outname = "spftot_gw"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+".dat"
	outfile = open(outname,'w')
	for i in xrange(np.size(newen)):
		outfile.write("%7.4f   %15.10e\n"% (newen[i],spftot[i])) # Dump string representations of arrays
	outfile.close()
	print " A(\omega)_GW written in", outname
	plt.plot(newen,spftot,label="ftot_gw");
	
# ============================= ###

### ===================================== ###
### ===== EXPONENTIAL SPECTRAL FUNCTION ====== ###
if flag_calc_exp == 1:
	print " ### Calculation of exponential A...  "
	### ==== Finding zero in res --> Eqp ===== ###
	print " Finding zeros in real parts..."
	eqp, imeqp = calc_eqp_imeqp(nkpt,nband,en,res,ims)
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
	from multipole import fit_multipole, getdata_file #, write_f_as_sum_of_poles
	print " ### ================== ###"
	print " ###    Multipole fit   ###"
	print " Number of poles:", npoles
	omegampole =  np.zeros((nkpt,nband,npoles))
	ampole =  np.zeros((nkpt,nband,npoles))
	for ik in xrange(nkpt):
		for ib in xrange(nband):
			print " ik, ib", ik, ib
			interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
			# Here we take the curve starting from eqp and then we invert it
			# so as to have it defined on the positive x axis
			# and so that the positive direction is in the 
			# increasing direction of the array index
			en3 = en[en<=eqp[ik,ib]]
			im3 = interpims(en3)/np.pi # This is what should be fitted
			en3 = en3 - eqp[ik,ib]
			en3 = -en3[::-1] 
			im3 = im3[::-1]
			omegai, gi, deltai = fit_multipole(en3,im3,npoles,0)
			omegampole[ik,ib] = omegai
			ampole[ik,ib] = gi/(omegai**2) # The weights should not be affected by the shift
			print " Integral test. Compare \int\Sigma and \sum_j^N\lambda_j."
			print " 1/pi*\int\Sigma   =", np.trapz(im3,en3)
			print " \sum_j^N\lambda_j =", np.sum(gi)
			#e1,f1 = write_f_as_sum_of_poles(en3,omegai,gi,deltai,0)
	dxexp=0.05
	enexp=np.arange(enmin,enmax,dxexp)
	nenexp=np.size(enexp)
	ftot=np.zeros((nenexp))
	f=np.zeros((nkpt,nband,nenexp))
	print " Calculating multipole exponential A..."
	ftot=np.zeros((np.size(enexp)))
	# Extrinsic and interference contribution
	if extinf == 1:
		# Here we add the extrinsic contribution. 
		# N.B.: It has to be renormalized to the number of poles!!!
		from multipole import getdata_file #, write_f_as_sum_of_poles
		enextinf, aextinf = getdata_file(origdir+"/a_wp.dat")
		newenexin = []
		newenexin.append(0.0)
		for x in enextinf.tolist():
			newenexin.append(x)
		newenexin = np.array(newenexin)
		newaexin = []
		newaexin.append(0.0)
		for x in aextinf.tolist():
			newaexin.append(x)
		newaexin = np.array(newaexin)
		#print enextinf, newenexin
		#print aextinf, newaexin
		interpextinf = interp1d(newenexin, newaexin, kind = 'linear', axis =  2)
		amp_exinf = ampole.copy()
		#print "Type(amp_exinf, ampole):", type(amp_exinf), type(ampole)
		for ik in xrange(nkpt):
			for ib in xrange(nband):
				tmpextinf = interpextinf(omegampole[ik,ib])/npoles # <-- Divided by the number of poles!
				amp_exinf[ik,ib] += tmpextinf
	# Time section
	import time
	e0=time.time()
	c0=time.clock()
	elaps1=time.time() - e0
	cpu1=time.clock() - c0
	print " Starting time (elaps, cpu):", elaps1, cpu1
	def calc_spf_mpole(enexp,akb,omegakb,eqpkb,imkb,npoles):
		ftot = np.zeros((np.size(enexp)))
		outnamekb = "spf_exp-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+"_mpole"+str(npoles)+".dat"
		outfilekb = open(outnamekb,'w')
		tmpf1=np.zeros((nenexp))
		tmpf2=np.zeros((nenexp))
		tmpf3=np.zeros((nenexp))
		for ipole in xrange(npoles):
			#print " First order"
			tmpf2=0.
			for jpole in xrange(npoles):
				#print " Second order"
				tmpf3=0.
				for kpole in xrange(npoles):
					## Fourth order
					#tmpf4 = np.zeros((np.size(enexp)))
					#for ipole4 in xrange(npoles):
						#tmpf4[:] += 1./4. * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole]  + omegampole[ik,ib,jpole]  + omegampole[ik,ib,kpole]  + omegampole[ik,ib,ipole4] )**2 + ( imeqp[ik,ib] )**2 ) )
					tmpf3+=1./3.*akb[kpole]/((enexp-eqpkb+omegakb[ipole]+omegakb[jpole]+omegakb[kpole])**2+imkb**2 ) 
				tmpf2+=1./2.*akb[jpole]*(1./((enexp - eqp[ik,ib]+omegakb[ipole]+omegakb[jpole])**2+imkb**2 ) + tmpf3 ) 
			tmpf1+=1.*akb[ipole]*(1./((enexp-eqp[ik,ib]+omegakb[ipole])**2+imkb**2)+tmpf2) 
		#f[ik,ib]=prefac*(1./((enexp-eqpkb)**2+(imkb)**2)+tmpf1) 
		#ftot += f[ik,ib]
		f=prefac*(1./((enexp-eqpkb)**2+(imkb)**2)+tmpf1) 
		ftot += f
		for ien in xrange(nenexp) :
			outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], f[ien]))
		outfilekb.close()
		return ftot
	if extinf == 1:
		for ik in xrange(nkpt):
			for ib in xrange(nband):
				print " ik, ib", ik, ib
				prefac=np.exp(-np.sum(amp_exinf[ik,ib,:])) * wtk[ik] * pdos[ib] / np.pi * abs( imeqp[ik,ib] )
				akb=amp_exinf[ik,ib] # This is a numpy array (slice)
				omegakb=omegampole[ik,ib] # This is a numpy array (slice)
				eqpkb=eqp[ik,ib]
				imkb=imeqp[ik,ib]
				tmpf = calc_spf_mpole(enexp,akb,omegakb,eqpkb,imkb,npoles)
				ftot += tmpf
	else: # extinf == 0
		for ik in xrange(nkpt):
			for ib in xrange(nband):
				print " ik, ib", ik, ib
				prefac=np.exp(-np.sum(ampole[ik,ib,:])) * wtk[ik] * pdos[ib] / np.pi * abs( imeqp[ik,ib] )
				akb=ampole[ik,ib] # This is a numpy array (slice)
				omegakb=omegampole[ik,ib] # This is a numpy array (slice)
				eqpkb=eqp[ik,ib]
				imkb=imeqp[ik,ib]
				tmpf = calc_spf_mpole(enexp,akb,omegakb,eqpkb,imkb,npoles)
				ftot += tmpf
	elaps2 = time.time() - elaps1 - e0
	cpu2 = time.clock() - cpu1 - c0
	#print elaps2, cpu2
	print " Used time (elaps, cpu):", elaps2, cpu2
	# Writing out a_j e omega_j
	print " ### Writing out a_j and omega_j..."
	outname = "a_j_mp"+str(npoles)+".dat"
	outfile = open(outname,'w')
	outname = "omega_j_mp"+str(npoles)+".dat"
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
	print " ### Writing out A(\omega)_exp...  "
	if extinf == 1:
		outname = "spftot_exp"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev_mp"+str(npoles)+"_extinf.dat"
		outfile = open(outname,'w')
		for i in xrange(nenexp):
			outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot[i])) # Dump string representations of arrays
		outfile.close()
	else: # extinf == 0
		outname = "spftot_exp"+"_s"+str(sfac)+"_p"+str(pfac)+"_"+str(penergy)+"ev_mp"+str(npoles)+".dat"
		#outname = "spftot_exp"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+"_mpole"+str(npoles)+".dat"
		outfile = open(outname,'w')
		for i in xrange(nenexp):
			outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot[i])) # Dump string representations of arrays
		outfile.close()
	print " A(\omega)_exp written in", outname
	plt.plot(enexp,ftot,label="ftot");
# Now go back to original directory
print " Moving back to parent directory:\n", origdir
chdir(newdir)
#title = 'Spectral function '+ 'A (' + r'$\omega $ ' + ') - '+r'$ h\nu = $'+str(penergy)+' eV'
#plt.title(title)
#plt.legend(loc=2);
plt.show();
