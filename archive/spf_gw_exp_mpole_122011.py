#!/usr/bin/env python
"""
Written by Matteo Guzzo.
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

print " Reading invar file... ",
invar = {}
if isfile("invar.in"):
	infile = open("invar.in")
	for line in infile.readlines():
		word = line.split()
		invar[word[-1]] = word[0];
#		print "invar: ", invar
	infile.close()
	sigfilename =       invar['sigmafile'];
	minband     =   int(invar['minband']);
	maxband     =   int(invar['maxband']);
	nkpt        =   int(invar['nkpt']);
	enmin       = float(invar['enmin']);
	enmax       = float(invar['enmax']);
	sfac        = float(invar['sfactor']);
	pfac        = float(invar['pfactor']);
	penergy     =   int(invar['penergy']);
	npoles      =   int(invar['npoles']);
	flag_gw     =   int(invar['gw']);
	a_extinf    = float(invar['a_extinf']);
	extinf      =   int(invar['extinf']);
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
# Input file (_SIG)
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
###print a[0]
en=[];
istart = 0;
# loop to prepare  the energy array
print " Reading array of energies from first k-point in _SIG file... ",
for line in filelines:
	if line[0:3]=='# k':
#		print line,
		continue
	elif line[0:3]=='# b':
		istart = 1 # Supposed to be 0 before
#		print line,
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
#		if re.search("k =",line): 
		if ikpt==nkpt :
			print " End of the reading loop: ikpt == nkpt. ikpt,nkpt: ", ikpt, nkpt; 
			print " #### ================================================ #### ";
			break;
#		print " ### === ########### === ### ";
		print " --- k point:  %02i ---" % (ikpt+1);
#		print "line = ", line[:], 
		ikpt = ikpt + 1;
	# Detect, in commented lines, the bands declaration
	# (and read the values)
	elif line[0:3]=='# b':
#		if re.search("b =",line): 
#		print "line = ", line[:], 
# 		ib = ib + 1;
		if imposebands == 0 :
			data=line.split() 
			minband = int(data[-2]) 
			maxband = int(data[-1]) 
			nband = 1 + maxband - minband
		io = 0;
	elif io < np.size(en) and ib < nband :
#		print "io,ib,icol:",io,ib,icol
		# Split this line into array of floats
#		print "io, ib, nband, icol, data[icol]: ", io, ib, nband, icol, data[icol];
		data=map(float,line.split());
		for col in data:
		# First element (energy) goes into en
			if icol == 0 : 
#				print "icol==0,icol:",icol
				icol = icol + 1;
				io = io + 1;
	#		continue;
			elif (icol+2)%3 == 0 :
#				print "icol+2%3==0,icol:",icol
#				res[ikpt-1,ib,io-1] = en[io-1] - elda[ikpt-1,ib] + vxc[ikpt-1,ib] - col;
#				print "Test ib line 138:", ib
				res[ikpt-1,ib,io-1] = col;
#				print "icol, Reading res, res(N) = ", icol, res;
				icol = icol + 1;
				continue;
			elif (icol+1)%3 == 0 :
				ims[ikpt-1,ib,io-1] = col;
#				print "icol, Reading ims, ims(N) = ", icol, ims;
				icol = icol + 1;
				continue;
			else : 
				if ib == nband-1 : break # number of bands reached
				ib = ib + 1;
				icol = icol + 1;
				continue;
		else: continue
print " ### nkpt, nband:", nkpt, nband
print " # ------------------------------------------------ # ";

### ===================================================== ###

print " ### Cross sections...  "
csfilename = "cs"+str(penergy)+".dat"
if isfile(csfilename):
	print " Photon energy:", penergy,"eV"
else:
	penergy = raw_input(" File "+csfilename+" not found. Photon energy (eV): ")
	csfilename = "cs"+str(penergy)+".dat"
cs = []
print "csfilename:",csfilename
csfile = open(csfilename,'r')
for line in csfile.readlines():
	cs.append((float(line)));
csfile.close()
cs = np.array(cs)
#print "cs:",cs.shape,cs
#print "cs:",np.transpose(cs),cs.shape

# ======== CROSS SECTIONS ======= #
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
pdos = 10000.*np.dot(cs,sp)
print "10000*pdos:", pdos
print "Size(pdos):",np.size(pdos)
# ============= ROADWORKS ============ #

# ============= ROADWORKS ============ #
# Here we move to a subdirectory to avoid flooding-up the current directory
newdirname = "Spfunctions"
origdir = getcwd() # remember where we are
newdir = join(origdir, newdirname) # Complete path of the new directory
print " Moving into spectral functions directory:\n", newdir
if not isdir(newdir) :
	mkdir(newdir)
chdir(newdir)

newdx = 0.01
if enmin < en[0] :  
	newen = np.arange(en[0],enmax,newdx)
else :  
	newen = np.arange(enmin,enmax,newdx)

# GW spectral function part
if flag_gw == 1:
	# ===========INTERPOLATION ============ #
	# Interpolation section
	print " ### ============= #### "
	print " ### Interpolation #### "
	print " ### ============= #### "
	
	# INTERPOLATION GRID DEFINED HERE #
	
	def spf(a,b):
		spf = abs(b) / np.pi / ( a**2 + b**2 )
		return spf
	
	print " ### Interpolation and calculation of A(\omega)_GW...  "
	spftot = np.zeros(np.size(newen));
	# Here we interpolate re and im sigma
	# for each band and k point
	for ik in xrange(np.size(ims[:,0,0])):
		print " k point = %02d " % (ik+1)
	#	ikstr = 
		#spfk = np.zeros(np.size(newen))
		for ib in xrange(np.size(ims[0,:,0])):
			#print " Outfile: ",outnamekb
			# Parti reali e immaginarie interpolate
			interpres = interp1d(en, res[ik,ib,:], kind = 'linear', axis =  2)
			interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
	#		print "Test ik,ib:", ik, ib
			redenom = newen - hartree[ik,ib] - interpres(newen)
	                tmpim = interpims(newen)
	                spfkb = wtk[ik] * pdos[ib] * spf(redenom, tmpim)
			spftot += spfkb 
			outnamekb = "spf_gw-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+".dat"
			outfilekb = open(outnamekb,'w')
			for ien in xrange(np.size(newen)) :
				outfilekb.write("%8.4f %12.8f %12.8f %12.8f\n" % (newen[ien], spfkb[ien], redenom[ien], tmpim[ien]))
			outfilekb.close()
	#plt.plot(newen,spftot,'-', label='spftot_gw')
	
		### ==== WRITING OUT GW SPECTRAL FUNCTION === ###
	print " ### Writing out A(\omega)_GW...  "
	outname = "spftot_gw"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+".dat"
	outfile = open(outname,'w')
	for i in xrange(np.size(newen)):
		outfile.write("%7.4f   %15.10e\n"% (newen[i],spftot[i])) # Dump string representations of arrays
	outfile.close()
	print " A(\omega)_GW written in", outname
	
# ============================= ###

# ======== INTEGRATION ======== ###
print " ### Integration...  "
# integration
# use trapz()
def integrate_absval(f,x):
  F=np.zeros(x.size);
  for i in xrange(x.size):
   if f[i] >= 0:
    # Simpson rule seems to be problematic for f(x)~0
    #F[i]=scipy.integrate.simps(f[0:i+1],x[0:i+1]);
    # Trapezoidal rule
    #F[i]=scipy.integrate.trapz(f[0:i+1],x=x[0:i+1]);

    # we use simple summation (see FAQ 3)
    F[i]=np.sum(f[0:i+1])*(x[1]-x[0]);
  return F;

outname = "imsint.dat"
imsint = np.zeros((nkpt,nband))

outfile = open(outname,'w')
#tempim = np.zeros(np.size(en))
for ik in xrange(nkpt):
	print " k point = %02d " % (ik+1)
	for ib in xrange(nband):
		# This is a function
		interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
		# This is an array filled with the values of the interpolation function
		#tempim[:] = interpims(newen[:])
		tempim = ims[ik,ib]
		tempim[tempim<0] = 0 # negative values shall be set to 0
		# Integral of the function
		#imsint[ik,ib] = np.trapz(tempim,newen) / np.pi
		imsint[ik,ib] = np.trapz(tempim,en) / np.pi
		outfile.write("%14.5f" % (imsint[ik,ib]))
#	print
	outfile.write("\n")
#	outfile.write("%7.4f   %15.10e\n"% (newen[i],spftot[i])) # 
outfile.close()
print " IMSINT ARRAY:",
print str(imsint[0,0])+'...'

# ======== FIND MAXIMUM PART ======== ###
print " Finding maxima in imaginary parts..."
# Create array for indices of maxima of ims
imax = np.zeros((nkpt,nband), int)
for ik in xrange(nkpt):
	for ib in xrange(nband):
		imax[ik,ib] = np.argmax(ims[ik,ib])
print "Test imax(np.argmax):", imax[0,0]
print "Energy from imax(ims):",newen[imax[0,0]]

### ==== Finding zero in res --> Eqp ===== ###
print " Finding zeros in real parts..."
#from scipy.optimize import fsolve
#izerores = np.zeros((nkpt,nband), int)
eqp = np.zeros((nkpt,nband))
imeqp = np.zeros((nkpt,nband))
#tempres = np.zeros(np.size(newen))
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
	#plt.plot(en,resigma,'-')
	for i in xrange(1,np.size(resigma)):
		if (resigma[i]*resigma[i-1] < 0):
			tmpeqp = en[i-1] - resigma[i-1]*(en[i] - en[i-1])/(resigma[i] - resigma[i-1]) # High school formula
			#plt.plot(en[i-1],resigma[i-1],'o')
			#plt.plot(en[i],resigma[i],'o')
			#plt.plot(tmpeqp,0,'o')
	return tmpeqp

for ik in xrange(nkpt):
	for ib in xrange(nband):
		temparray = np.array(en[:] - hartree[ik,ib] - res[ik,ib,:])
		interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
		tempim = interpims(newen[:])
		# New method to overcome plasmaron problem
		eqp[ik,ib] = find_eqp_resigma(en,temparray)
		imeqp[ik,ib] = interpims(eqp[ik,ib])
		## Warning if imaginary part of sigma < 0 (Convergence problems?)
		if imeqp[ik,ib] <= 0 : print " WARNING: im(Sigma(eps_k)) <= 0 !!! ik ib eps_k im(Sigma(eps_k)) = ", ik, ib, eqp[ik,ib], imeqp[ik,ib]
print " Test imeqp:", imeqp
#omegaps =  np.zeros((nkpt,nband))
#lambdas =  np.zeros((nkpt,nband))
#omegaps[:,:] = eqp[:,:] - newen[imax[:,:]]
omegaps = eqp - newen[imax]
#lambdas[:,:] = imsint[:,:] / ( omegaps[:,:] )**2 
lambdas = imsint / ( omegaps )**2 
# Writing out omegaps
# Writing out lambdas
# Writing out eqp
# Writing out imeqp
outname = "omegaps.dat"
outfile = open(outname,'w')
outname = "eqp.dat"
outfile2 = open(outname,'w')
outname = "imeqp.dat"
outfile3 = open(outname,'w')
outname = "a_j.dat"
outfile4 = open(outname,'w')
for ik in xrange(nkpt):
	for ib in xrange(nband):
		outfile.write("%14.5f" % (omegaps[ik,ib]))
		outfile2.write("%14.5f" % (eqp[ik,ib]))
		outfile3.write("%14.5f" % (imeqp[ik,ib]))
		outfile4.write("%14.5f" % (lambdas[ik,ib]))
	outfile.write("\n")
	outfile2.write("\n")
	outfile3.write("\n")
	outfile4.write("\n")
outfile.close()
outfile2.close()
outfile3.close()
outfile4.close()
print " OMEGAPS ARRAY:",
print str(omegaps[0,0])+'...'
print " LAMBDAS ARRAY:",
print str(lambdas[0,0])+'...'

### ===================================== ###
# ROADWORKS #
### ===== EXPONENTIAL SPECTRAL FUNCTION ====== ###
print " ### Calculation of exponential A...  "
plt.figure(2)

from multipole import fit_multipole #, write_f_as_sum_of_poles
flag_mpole = 1
if flag_mpole == 1:
	print " ### ================== ###"
	print " ###    Multipole fit   ###"
	print " Number of poles:", npoles
	#norder = 2 # Order of the exponential development
	omegampole =  np.zeros((nkpt,nband,npoles))
	ampole =  np.zeros((nkpt,nband,npoles))
	for ik in xrange(nkpt):
		for ib in xrange(nband):
			print " ik, ib", ik, ib
			#tempim = np.zeros(np.size(en))
			interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
			# Here we take the curve starting from eqp and then we invert it
			# so as to have it defined on the positive x axis
			# and so that the positive direction is in the 
			# increasing direction of the array index
			en3 = en[en<=eqp[ik,ib]]
			im3 = interpims(en3)/np.pi # This is what should be fitted
			#plt.plot(en3,im3,'-',label="after truncation");
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
			#plt.plot(-en3,im3,'-o'); plt.plot(e1,f1,'-'); plt.plot(en2-eqp[ik,ib],tempim2,'-',label="tempim"); plt.show(); sys.exit()
			#plt.plot(en3,im3,'-o'); plt.plot(omegai,np.pi/2.*gi*omegai/deltai,'-o'); plt.show(); sys.exit()
			#sys.exit()
	dxexp = 0.05
	enexp = np.arange(enmin,enmax,dxexp)
	nenexp = np.size(enexp)
	ftot = np.zeros((nenexp))
	f =  np.zeros((nkpt,nband,nenexp))
	print " Calculating multipole exponential A..."
	ftot = np.zeros((np.size(enexp)))
	# Time section
	import time
	e0 = time.time()
	c0 = time.clock()
	elaps1 = time.time() - e0
	cpu1 = time.clock() - c0
	print " Starting time (elaps, cpu):", elaps1, cpu1
	for ik in xrange(nkpt):
		for ib in xrange(nband):
			print " ik, ib", ik, ib
			prefac = np.exp(-np.sum(ampole[ik,ib,:])) * wtk[ik] * pdos[ib] / np.pi * abs( imeqp[ik,ib] )
			akb = ampole[ik,ib] # This is a numpy array (slice)
			omegakb = omegampole[ik,ib] # This is a numpy array (slice)
			eqpkb = eqp[ik,ib]
			imkb = imeqp[ik,ib]
			outnamekb = "spf_exp-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+"_mpole"+str(npoles)+".dat"
			outfilekb = open(outnamekb,'w')
			tmpf1 = np.zeros((nenexp))
			tmpf2 = np.zeros((nenexp))
			tmpf3 = np.zeros((nenexp))
			for ipole in xrange(npoles):
				#print " First order"
				tmpf2 = 0.
				for jpole in xrange(npoles):
					#print " Second order"
					tmpf3 = 0.
					for kpole in xrange(npoles):
						## Fourth order
						#tmpf4 = np.zeros((np.size(enexp)))
						#for ipole4 in xrange(npoles):
							#tmpf4[:] += 1./4. * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole]  + omegampole[ik,ib,jpole]  + omegampole[ik,ib,kpole]  + omegampole[ik,ib,ipole4] )**2 + ( imeqp[ik,ib] )**2 ) )
						tmpf3 += 1./3.*akb[kpole] / ( ( enexp - eqpkb + omegakb[ipole]  + omegakb[jpole]  + omegakb[kpole] )**2 + ( imkb )**2 ) 
					tmpf2 += 1./2.*akb[jpole] * ( 1. / ( ( enexp - eqp[ik,ib] + omegakb[ipole]  + omegakb[jpole] )**2 + ( imkb )**2 ) + tmpf3 ) 
				tmpf1 += 1.*akb[ipole] * ( 1. / ( ( enexp - eqp[ik,ib] + omegakb[ipole] )**2 + ( imkb )**2 ) + tmpf2 ) 
			f[ik,ib] = prefac * ( 1. / ( ( enexp - eqpkb )**2 + ( imkb )**2 ) + tmpf1 ) 
			ftot += f[ik,ib]
			for ien in xrange(nenexp) :
				outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], f[ik,ib,ien]))
			outfilekb.close()
			#plt.plot(enexp,f,label="f"); plt.show()
	elaps2 = time.time() - elaps1 - e0
	cpu2 = time.clock() - cpu1 - c0
	#print elaps2, cpu2
	print " Used time (elaps, cpu):", elaps2, cpu2
	# Writing out a_j e omega_j
	print " ### Writing out a_j and omega_j..."
	outname = "ampole.dat"
	outfile = open(outname,'w')
	outname = "omegampole.dat"
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
	outname = "spftot_exp"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+"_mpole"+str(npoles)+".dat"
	outfile = open(outname,'w')
	for i in xrange(nenexp):
		outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot[i])) # Dump string representations of arrays
	outfile.close()
	print " A(\omega)_exp written in", outname
	plt.plot(enexp,ftot,label="ftot");


### ===================================== ###

# Now go back to original directory
print " Moving back to parent directory:\n", origdir
chdir(newdir)

plt.title('Spectral function '+ 'A (' + r'$\omega $ ' + ') - '+r'$ h\nu = $'+str(penergy)+' eV')
		
plt.legend(loc=2);
plt.show();
