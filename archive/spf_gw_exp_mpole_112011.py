#!/usr/bin/env python
'''
Written by Matteo Guzzo.
List of files needed:
invar.in with input variables
_SIG for the self-energy
s.dat, p_even.dat, p_odd.dat, d_even.dat, etc. for the orbital character
cs*.dat for the photon cross sections
hartree.dat or elda.dat and vxc.dat for the hartree energies
wtk.dat for the k-points weights
'''
import numpy as np;
import matplotlib.pylab as plt;
#from numpy import random;
from scipy.interpolate import interp1d
from scipy import optimize
#from scipy.signal import cspline1d, cspline1d_eval
#import re # Regexp
import sys
#from os.path import isfile, path, mkdir, chdir
from os.path import isfile, join, isdir
from os import getcwd, pardir, mkdir, chdir

#nband=input("ENTER number of bands:  ")
#nkpt=input("ENTER number of kpt:  ")
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
	plasmon2    =   int(invar['plasmon2']);
	a_extinf    = float(invar['a_extinf']);
	extinf      =   int(invar['extinf']);
	npoles      =   int(invar['npoles']);
else : 
	print "Invar file not found (invar.in). Impossible to continue."
	sys.exit(1)
print "Done."
# Input file (_SIG)
if isfile(sigfilename):
	insigfile = open(sigfilename);
else:
	print "File "+sigfilename+" not found."
	insigfile = open(raw_input("Self-energy file name (_SIG): "))
# First band to be used
#minband =   1
# Number of the upper band to be included in the spectral function (usually top valence)
#maxband =   8
# Here we decide wether to impose band parameters before or to read them from file instead
# TODO: write a function to read the parameters from files
imposebands = 1;
if imposebands == 1:
	nband = 1 + maxband - minband;
print " minband =", minband;
print " maxband =", maxband;
print " nband =", nband;
# Number of k points - this is imposed
# TODO: write a function to read the parameters from files
#nkpt = 24;
print " nkpt =", nkpt;
print " enmin =", enmin;
print " enmax =", enmax;
## Mixing parameters for s and p
#sfac = 1.0
#pfac = 1.0
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
print " np.shape(res):", np.shape(res)
ims=np.zeros((nkpt,nband,np.size(en)));
#allres=[];
#allims=[];
#exit()
#for line in Hhartreefile.readlines(): 
#	hartree.append(map(float,line.split()));
#for i in range(len(hartree[0])):
#	print hartree[:][0];
#print " hartree shape:", hartree.shape;
#print " hartree shape (len,len,):", len(hartree[:,0]),len(hartree[0,:]);
#print " hartree[:,0]:", hartree[:,0]
#print " hartree[0,:]:", hartree[0,:]
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

#spf=np.zeros((nkpt,nband,np.size(en)));
#spf[:,:,:] = ims[:,:,:]/np.pi/( res[:,:,:]**2 + ims[:,:,:]**2 )
#for ik in xrange(np.size(ims[:,0,0])):
#	for ib in xrange(np.size(ims[0,:,0])):
##		spf[:] = spf[:] + ims[ik,ib,:]/np.pi()/( res[ik,ib,:]**2 + ims[ik,ib,:]**2 )
#		spftot[:] = spftot[:] + ims[ik,ib,:]/np.pi/( res[ik,ib,:]**2 + ims[ik,ib,:]**2 )
#		continue

### ===================================================== ###

print " ### Cross sections...  "
csfilename = "cs"+str(penergy)+".dat"
if isfile(csfilename):
	print " Photon energy:", penergy,"eV"
else:
	penergy = raw_input(" File "+csfilename+" not found. Photon energy (eV): ")
	csfilename = "cs"+str(penergy)+".dat"
#penergy = 800
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
s = sfac*s
peven = sfac*peven
podd = pfac*podd
p = peven+podd
#print s
#print p
#print "p:",p
#print "s:",s
#print "p:",p
#print "cs[0]:",cs[0]
#print "cs[1]:",cs[1]
#print "s*cs[0]:",s*cs[0]
#print "s*cs[1]:",s*cs[1]
#exit()
sp = np.array([s,p])
#print "sp:",sp
pdos = 10000.*np.dot(cs,sp)
print "10000*pdos:", pdos
print "Size(pdos):",np.size(pdos)
# ============= ROADWORKS ============ #

# ============= ROADWORKS ============ #
# Interpolation section
print " ### ============= #### "
print " ### Interpolation #### "
print " ### ============= #### "

# INTERPOLATION GRID DEFINED HERE #

def spf(a,b):
	spf = abs(b) / np.pi / ( a**2 + b**2 )
	return spf
# Here we move to a subdirectory to avoid flooding-up the current directory
newdirname = "Spfunctions"
origdir = getcwd() # remember where we are
newdir = join(origdir, newdirname) # Complete path of the new directory
print " Moving into spectral functions directory:\n", newdir
if not isdir(newdir) :
	mkdir(newdir)
chdir(newdir)

print " ### Interpolation and calculation of A(\omega)_GW...  "
newdx = 0.01
#newen = np.arange(en[0],en[-1],newdx)
if enmin < en[0] :  
	newen = np.arange(en[0],enmax,newdx)
else :  
	newen = np.arange(enmin,enmax,newdx)
spfkb = np.zeros(np.size(newen));
spftot = np.zeros(np.size(newen));
# Here we interpolate re and im sigma
# for each band and k point
#interpres = np.zeros((nkpt,nband,np.size(newen)))
#interpims = np.zeros((nkpt,nband,np.size(newen)))
for ik in xrange(np.size(ims[:,0,0])):
	print " k point = %02d " % (ik+1)
#	ikstr = 
	#spfk = np.zeros(np.size(newen))
	for ib in xrange(np.size(ims[0,:,0])):
		outnamekb = "spf_gw-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+".dat"
		outfilekb = open(outnamekb,'w')
		#print " Outfile: ",outnamekb
		# Parti reali e immaginarie interpolate
		#interpres[ik,ib,:] = interp1d(en, res[ik,ib,:], kind = 'linear', axis =  2)
		#interpims[ik,ib,:] = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
#		tempres = np.array(res[ik,ib,:],float)
#		interpres[ik,ib] = 0
#		a = interp1d(en,tempres)
#		print "tempres:",tempres[0]
#		exit()
#		interpims[ik,ib,:] = interp1d(en, ims[ik,ib,:], kind = 'linear')
		interpres = interp1d(en, res[ik,ib,:], kind = 'linear', axis =  2)
		interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
#		for j in xrange(np.size(newen)):
#			spfkb =  wtk[ik] * pdos[ib] * spf(interpres(newen[j]),interpims(newen[j]))
#			spftot[j] = spftot[j] + spfkb
#		print "Test ik,ib:", ik, ib
#		redenom[:] = newen[:] - elda[ik,ib] + vxc[ik,ib] - interpres(newen[:])
		redenom = np.zeros(np.size(newen))
		im = np.zeros(np.size(newen))
		redenom[:] = newen[:] - hartree[ik,ib] - interpres(newen[:])
                im[:] = interpims(newen[:])
                spfkb[:] = wtk[ik] * pdos[ib] * spf(redenom[:], im[:])
		spftot[:] += spfkb[:] 
		for ien in xrange(np.size(newen)) :
			outfilekb.write("%8.4f %12.8f %12.8f %12.8f\n" % (newen[ien], spfkb[ien], redenom[ien], im[ien]))
		outfilekb.close()
#		spftot[:] = spftot[:] + ims[ik,ib,:]/np.pi/( res[ik,ib,:]**2 + ims[ik,ib,:]**2 )
#interpres = interp1d(en, res[0,7,:], kind = 'linear', axis =  2)
#interpims = interp1d(en, ims[0,7,:], kind = 'linear', axis =  2)

#plt.plot(newen,interpres(newen),'--');
#plt.plot(newen,interpims(newen),'--');
#plt.plot(newen,spf(interpres(newen),interpims(newen)),'-');
#plt.figure(2)
#plt.plot(newen,spftot,'-', label='spftot_gw')

#import scitools.filetable as ft
#ft.write_columns(outfile,newen,spftot)
#a = np.array([newen,spftot])
### ==== WRITING OUT GW SPECTRAL FUNCTION === ###
print " ### Writing out A(\omega)_GW...  "
outname = "spftot_gw"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+".dat"
outfile = open(outname,'w')
for i in xrange(np.size(newen)):
	outfile.write("%7.4f   %15.10e\n"% (newen[i],spftot[i])) # Dump string representations of arrays
outfile.close()
print " A(\omega)_GW wrote in", outname
#plt.plot(en,res[0,7,:],'--') #,newen,newen,'x',en,en,'o');
#plt.plot(en,ims[0,0,:],'--', label='ims[0,0]') #,newen,newen,'x',en,en,'o');
#plt.plot(en,ims[0,7,:],'--', label='ims[0,7]') #,newen,newen,'x',en,en,'o');
#plt.plot(en,ims[0,0,:],'--') #,newen,newen,'x',en,en,'o');
#plt.plot(en,ims[0,7,:],'--') #,newen,newen,'x',en,en,'o');

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
#if isfile(outname) :
#	print " imsint.dat found. Reading integral values..."
#	infile = open(outname,'r')
#	for line in hartreefile.readlines():
#		hartree.append(map(float,line.split()));
#	for ik, line in zip( xrange(nkpt), infile.readlines):
#		print " k point = %02d " % (ik+1)
#		outfile.readline("%14.5f" % (imsint[ik,ib]))
#		for ib in xrange(nband):
#	infile.read("%14.5f" % (imsint[:,:]))
#	print imsint
#	exit()

outfile = open(outname,'w')
#print "Test integral 0,0:",integrate_absval(ims[0,0,],en)
#print "Test integral 0,0:",np.trapz(abs(ims[0,0,:]),en)
tempim = np.zeros(np.size(newen))
for ik in xrange(nkpt):
	print " k point = %02d " % (ik+1)
	for ib in xrange(nband):
		# This is a function
		interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
		# This is an array filled with the values of the interpolation function
#		tempim2 = np.array(tempim)
#		tempim3 = np.array(tempim)
#		tempim2[:] = interpims2(newen[:])
#		tempim3[:] = interpims3(newen[:])
		tempim[:] = interpims(newen[:])
		tempim[tempim<0] = 0 # negative values shall be set to 0
		#plt.figure(2)
		#plt.plot(newen,tempim,label='cubic')
		#plt.plot(en,ims[ik,ib],'o')
#		plt.plot(newen,tempim2,label='quad')
#		plt.plot(newen,tempim3,label='cubic')
#		plt.legend()
#		plt.show()
#		exit()
#		imsint[ik,ib] = integrate_absval(abs(tempim),newen)
		# Integral of the function
		imsint[ik,ib] = np.trapz(tempim,newen) / np.pi
#		print imsint[ik,ib],
		outfile.write("%14.5f" % (imsint[ik,ib]))
#	print
	outfile.write("\n")
#	outfile.write("%7.4f   %15.10e\n"% (newen[i],spftot[i])) # 
outfile.close()
print " IMSINT ARRAY:",
print str(imsint[0,0])+'...'

#plt.plot(en,en[:]-elda[0,0]+vxc[0,0]-res[0,0,:],'--',label='res[0,0]') #,newen,newen,'x',en,en,'o');

# ======== FIND MAXIMUM PART ======== ###
print " Finding maxima in imaginary parts..."
# use np.max();
# Create array for indices of maxima of ims
imax = np.zeros((nkpt,nband), int)
#print imax
#imax[0:nkpt,0:nband] = np.argmax(ims[0:nkpt,0:nband,:])
for ik in xrange(nkpt):
	for ib in xrange(nband):
		interpims = interp1d(en, ims[ik,ib,:], kind = 'cubic', axis =  2)
#		print "ik, ib:",ik ,ib
		tempim[:] = interpims(newen[:])
		imax[ik,ib] = np.argmax(tempim[:])
#		plabel = 'imsmax['+str(ik)+','+str(ib)+']'
#		plt.plot(en[imax[ik,ib]],ims[ik,ib,imax[ik,ib]],'o')
print "Test imax(np.argmax):", imax[0,0]
print "Energy from imax(ims):",newen[imax[0,0]]

### ======= Finding zero in res ========= ###
print " Finding zeros in real parts..."
from scipy.optimize import fsolve
izerores = np.zeros((nkpt,nband), int)
eqp = np.zeros((nkpt,nband))
imeqp = np.zeros((nkpt,nband))
tempres = np.zeros(np.size(newen))
def find_eqp_resigma(en,resigma):
	"""
	This function is supposed to deal with the plasmaron problem 
	and calculate the quasiparticle energy once it is fed with 
	\omega - \epsilon_H - \Re\Sigma. 
	It should return the value of eqp. 
	"""
	import numpy as np
	import matplotlib.pylab as plt
	#plt.plot(en,resigma,'-')
	for i in xrange(1,np.size(resigma)):
		if (resigma[i]*resigma[i-1] < 0):
			tmpeqp = en[i-1] - resigma[i-1]*(en[i] - en[i-1])/(resigma[i] - resigma[i-1])
			#plt.plot(en[i-1],resigma[i-1],'o')
			#plt.plot(en[i],resigma[i],'o')
			#plt.plot(tmpeqp,0,'o')
	#plt.show()
	#sys.exit()
	return tmpeqp
for ik in xrange(nkpt):
	for ib in xrange(nband):
#		temparray = np.array(en[:] - elda[ik,ib] + vxc[ik,ib] - res[ik,ib,:])
		temparray = np.array(en[:] - hartree[ik,ib] - res[ik,ib,:])
		interpres = interp1d(en, temparray[:], kind = 'linear', axis =  2)
		interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
		tempres[:] = interpres(newen[:])
		tempim[:] = interpims(newen[:])
		#plt.plot(en,ims[ik,ib])
		#plt.plot(newen,tempim)
		#plt.show()
		#exit()
		#tempim[tempim<0] = 0
		#print "PLASMARONSSSSSS!!!!"
		### Here I put a temporary solution that seems to work (but I am not sure it will always do)
		#zeroguess = [en[0]+(en[-1]-en[0])/3., en[0]+(en[-1]-en[0])/2., en[0]+(en[-1]-en[0])*2./3.]
		#zeroguess = [en[0]+(en[-1]-en[0])*2./10., en[0]+(en[-1]-en[0])/2., en[0]+(en[-1]-en[0])*8./10.]
		#zeroguess = [en[0], en[-1]]
		#xzeros = fsolve(interpres,zeroguess)
		#if xzeros[0]!=xzeros[-1]:
		#	#zeroguess = [en[0]+(en[-1]-en[0])/3., en[0]+(en[-1]-en[0])/2., en[0]+(en[-1]-en[0])*2./3.]
		#	#xzeros = fsolve(interpres,zeroguess)
		#	print xzeros
		#	#zeroguess = [xzeros[0], xzeros[0]+(xzeros[-1]-xzeros[0])/3, xzeros[-1]]
		#	zeroguess = [xzeros[0], xzeros[0]+(xzeros[-1]-xzeros[0])/3, xzeros[-1]]
		#	xzeros = fsolve(interpres,zeroguess)
		#	print xzeros
		### Old zero-finding way
		izerores[ik,ib] = np.argmax(-abs(tempres[:]))
		eqp[ik,ib] = newen[izerores[ik,ib]]
		imeqp[ik,ib] = tempim[izerores[ik,ib]]
		# New method to overcome plasmaron problem
		eqp[ik,ib] = find_eqp_resigma(en,temparray)
		imeqp[ik,ib] = interpims(eqp[ik,ib])
		## New way
		#eqp[ik,ib] = interpres(xzeros[-1])
		#imeqp[ik,ib] = interpims(xzeros[-1])
		##plt.plot(newen,tempres,'-')
		##for i in xrange(len(xzeros)):
		##	plt.plot(xzeros[i],0,'o')
		##plt.show()
		##exit()
		if imeqp[ik,ib] <= 0 : print " WARNING: im(Sigma(eps_k)) <= 0 !!! ik ib eps_k im(Sigma(eps_k)) = ", ik, ib, eqp[ik,ib], imeqp[ik,ib]
print " Test res==0, i_0, en:", izerores, newen[izerores[0,0]]
print " Test imeqp:", imeqp
#exit()
#plt.plot(newen,abs(interpres(newen)),'-',label='interpres')
#plt.show()
#exit()
omegaps =  np.zeros((nkpt,nband))
lambdas =  np.zeros((nkpt,nband))
omegaps[:,:] = eqp[:,:] - newen[imax[:,:]]
lambdas[:,:] = imsint[:,:] / ( omegaps[:,:] )**2 
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
from multipole import fit_multipole, write_f_as_sum_of_poles
flag_mpole = 1
if flag_mpole == 1:
	print " ### ================== ###"
	print " ###    Multipole fit   ###"
	#mpoleparam = np.zeros((nkpt,nband,3,npoles))
	#mpoleparam = []
	#npoles =  20
	print " Number of poles:", npoles
	#norder = 2 # Order of the exponential development
	omegampole =  np.zeros((nkpt,nband,npoles))
	ampole =  np.zeros((nkpt,nband,npoles))
	for ik in xrange(nkpt):
		#mpoleparam.append([])
		for ib in xrange(nband):
			print " ik, ib", ik, ib
			#mpoleparam.append([])
			#tempim = np.zeros(np.size(newen))
			tempim = np.zeros(np.size(en))
			interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
			#tempim[:] = interpims(newen[:])
			#tempim2 = tempim[tempim>0]
			#en2 = newen[tempim>0]
			#en2 = en2-eqp[ik,ib]
			# Here we take the curve starting from eqp and then we invert it
			# so as to have it defined on the positive x axis
			# and so that the positive direction is in the 
			# increasing direction of the array index
			#en3 = en2[en2<=eqp[ik,ib]]
			#en3 = newen[newen<=eqp[ik,ib]]
			en3 = en[en<=eqp[ik,ib]]
			im3 = interpims(en3)/np.pi # This is what should be fitted
			#plt.plot(en3,im3,'-',label="after truncation");
			en3 = en3 - eqp[ik,ib]
			#en3 = -en3
			#safe_shift = 1. # this is a shift to avoid being too close to 0. It is subtracted to the values of the poles after the fit
			en3 = -en3[::-1] 
			#plt.plot(en3,im3,'-',label="after inversion and shift");
			im3 = im3[::-1]
			#omegai, gi, deltai = fit_multipole(en3-safe_shift,im3,npoles,1)
			omegai, gi, deltai = fit_multipole(en3,im3,npoles,0)
			#plt.plot(omegai,np.pi/2*gi*omegai/deltai,'-x',label="raw fit")
			#plt.plot(omegai-safe_shift,np.pi/2*gi*omegai/deltai,'-o',label="with safe_shift")
			omegampole[ik,ib] = omegai
			ampole[ik,ib] = gi/(omegai**2) # The weights should not be affected by the shift
			print " Integral test. Compare \int\Sigma and \sum_j^N\lambda_j."
			print " 1/pi*\int\Sigma   =", np.trapz(im3,en3)
			print " \sum_j^N\lambda_j =", np.sum(gi)
			#print np.shape(omegai)
			#print omegai, gi, deltai
			#plt.legend();plt.show(); sys.exit()
			#print mpoleparam
			#for i in xrange(npoles):
			#mpoleparam[ik,ib,0,:] = omegai[:]
			#mpoleparam[ik,ib,1,:] = gi[:]
			#mpoleparam[ik,ib,2,:] = deltai[:]
			#mpoleparam[ik,ib].append([omegai, gi, deltai])
			#print mpoleparam
			#print np.shape(mpoleparam)
			#e1,f1 = write_f_as_sum_of_poles(en3,omegai,gi,deltai,0)
			#plt.plot(-en3,im3,'-o'); plt.plot(e1,f1,'-'); plt.plot(en2-eqp[ik,ib],tempim2,'-',label="tempim"); plt.show(); sys.exit()
			#plt.plot(en3,im3,'-o'); plt.plot(omegai,np.pi/2.*gi*omegai/deltai,'-o'); plt.show(); sys.exit()
			#sys.exit()
	#print ampole.ravel(), omegampole.ravel()
	#for aa,bb in zip(np.ravel(omegampole), np.ravel(ampole)):print aa,bb; 
	#exit()
	dxexp = 0.05
	enexp = np.arange(enmin,enmax,dxexp)
	ftot = np.zeros((np.size(enexp)))
	f =  np.zeros((nkpt,nband,np.size(enexp)))
	print " Calculating multipole A..."
	ftot = np.zeros((np.size(enexp)))
	for ik in xrange(nkpt):
		for ib in xrange(nband):
			print " ik, ib", ik, ib
			prefac = np.exp(-np.sum(ampole[ik,ib,:])) * wtk[ik] * pdos[ib] / np.pi
			outnamekb = "spf_exp-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+"_mpole"+str(npoles)+".dat"
			outfilekb = open(outnamekb,'w')
#			import time
#			e0 = time.time()
#			c0 = time.clock()
#			for ien in enexp:
#				tmpf2 = 0
#				tmpf3 = 0
#				for ipole1 in xrange(npoles):
#					#tmpf2 = np.zeros((np.size(enexp)))
#					for ipole2 in xrange(npoles):
#						tmpf3 += np.sum(1./3.*ampole[ik,ib,:] * ( abs( imeqp[ik,ib] ) / ( ( enexp[ien] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,ipole2]  + omegampole[ik,ib,:] )**2 + ( imeqp[ik,ib] )**2 ) ) )
#					tmpf2 += np.sum( 1./2.*ampole[ik,ib,:] * ( abs( imeqp[ik,ib] ) / ( ( enexp[ien] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,:] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf3 ) )
#				tmpf1 = np.sum( 1.*ampole[ik,ib,:] * ( abs( imeqp[ik,ib] ) / ( ( enexp[ien] - eqp[ik,ib] + omegampole[ik,ib,:] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf2 ) )
#				f = prefac * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf1 ) 
#			elaps1 = time.time() - e0
#			cpu1 = time.clock() - c0
#			print elaps1, cpu1
			tmpf1 = np.zeros((np.size(enexp)))
			for ipole1 in xrange(npoles):
				#print " First order"
				tmpf2 = np.zeros((np.size(enexp)))
				for ipole2 in xrange(npoles):
					#print " Second order"
					tmpf3 = np.zeros((np.size(enexp)))
					for ipole3 in xrange(npoles):
						## Fourth order
						#tmpf4 = np.zeros((np.size(enexp)))
						#for ipole4 in xrange(npoles):
							#tmpf4[:,ipole1,ipole2,ipole3] += 1./4. * ( abs( imeqp[ik,ib] ) / ( ( en[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,ipole2]  + omegampole[ik,ib,ipole3]  + omegampole[ik,ib,ipole4] )**2 + ( imeqp[ik,ib] )**2 ) )
							#tmpf4[:] += 1./4. * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,ipole2]  + omegampole[ik,ib,ipole3]  + omegampole[ik,ib,ipole4] )**2 + ( imeqp[ik,ib] )**2 ) )
						#tmpf3[:,ipole1,ipole2] += 1./3. * ( abs( imeqp[ik,ib] ) / ( ( en[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,ipole2]  + omegampole[ik,ib,ipole3] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf4[:,ipole1,ipole2,ipole3] ) 
						tmpf3[:] += 1./3.*ampole[ik,ib,ipole3] * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,ipole2]  + omegampole[ik,ib,ipole3] )**2 + ( imeqp[ik,ib] )**2 ) ) 
					#tmpf2[:] += 1./2. * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,ipole2] )**2 + ( imeqp[ik,ib] )**2 ) ) 
					tmpf2[:] += 1./2.*ampole[ik,ib,ipole2] * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,ipole2] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf3[:] ) 
				#tmpf1[:] += 1. * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1] )**2 + ( imeqp[ik,ib] )**2 ) ) 
				tmpf1[:] += 1.*ampole[ik,ib,ipole1] * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf2[:] ) 
			f[ik,ib,:] = prefac * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf1[:] ) 
			#f[:] += prefac * ( abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf1[:] ) 
			for ien in xrange(np.size(enexp)) :
				outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], f[ik,ib,ien]))
			outfilekb.close()
			ftot[:] += f[ik,ib,:]
#			elaps2 = time.time() - elaps1 - e0
#			cpu2 = time.clock() - cpu1 - c0
#			print elaps2, cpu2
			#f[:] += prefac * ( abs( imeqp[ik,ib] ) / ( ( en[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )) 
			#plt.plot(enexp,f,label="f"); plt.show()
			#f += abs( imeqp[ik,ib] ) / ( ( en[:] - eqp[ik,ib] + omegampole[ik,ib,ipole1]  + omegampole[ik,ib,ipole2]  + omegampole[ik,ib,ipole3] )**2 + ( imeqp[ik,ib] )**2 ) + tmpf/4
			#ftot *= f
					#order0[:] = abs( imeqp[ik,ib] ) / ( ( en[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
					#order1[:] = lambdas[ik,ib] * abs( imeqp[ik,ib] ) / ( ( en[:] - eqp[ik,ib] + 1.* omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			#outfile1.write("%14.5f" % (lomegap1[ik,ib]))
			#outfile2.write("%14.5f" % (lomegap2[ik,ib]))
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
	for i in xrange(np.size(enexp)):
		outfile.write("%7.4f   %15.10e\n"% (enexp[i],ftot[i])) # Dump string representations of arrays
	outfile.close()
	print " A(\omega)_exp written in", outname
	plt.plot(enexp,ftot,label="ftot"); plt.show()
	exit()



# Lorentian fit
elif plasmon2 == 1 :
	print " ### Lorentian fit...  "
	# Fit the first set
	outname = "lomegap1.dat"
	outfile1 = open(outname,'w')
	outname = "lomegap2.dat"
	lomegap1 =  np.zeros((nkpt,nband))
	lomegap2 =  np.zeros((nkpt,nband))
	a1 =  np.zeros((nkpt,nband))
	a2 =  np.zeros((nkpt,nband))
	outfile2 = open(outname,'w')
	newen2 = np.arange(-100,en[-1],newdx)
	plt.figure(1)
	for ik in xrange(nkpt):
	#for ik in (11,):
		dummy = 0
		if nkpt%4 > 0 and nkpt > 4 :
			dummy = 1
			plt.subplot(nkpt/4,4+dummy,ik+1)
		elif nkpt <= 4 :
			dummy = 1
			plt.subplot(nkpt,dummy,ik+1)
		else :
			plt.subplot(nkpt/4,4,ik+1)
		print " k point = %02d " % (ik+1)
		for ib in xrange(nband):
		#for ib in (6,7):
			#plt.subplot(2,nband/2,ib+1) # Temporary line for testing
			print " band = %02d " % (ib+1)
			# Define anonymous functions:
			#fitfunc = lambda p, x: p[0] * p[2] / ( ( x + p[1])**2 + ( p[2] )**2 )  # Target function
			# First function is to fit small first plasmon peak
			fitfunc = lambda p, x: abs(p[0] * p[1]) / ( ( x - eqp[ik,ib] + p[2])**2 + ( p[1] )**2 )
			# Second function is to fit the second big plasmon peak, keeping fixed the first
			fitfunc2 = lambda p1, p2, x: abs(p1[0] * p1[1]) / ( ( x - eqp[ik,ib] + p1[2])**2 + ( p1[1] )**2 ) + abs(p2[0] * p2[1]) / ( ( x - eqp[ik,ib] + p2[2])**2 + ( p2[1] )**2 ) # Target function
			#fitfunc2 = lambda p, o1, x: abs(p[0] * p[1]) / ( ( x - eqp[ik,ib] + p[2])**2 + ( p[1] )**2 ) + abs(p[3] * p[4]) / ( ( x - eqp[ik,ib] + o1)**2 + ( p[4] )**2 ) # + abs(p[6] * p[8]) / ( ( x - eqp[ik,ib] + p[7])**2 + ( p[8] )**2 ) # Target function
			#fitfunc = lambda p, x: p[0] * 0.001 / ( ( x + p[1])**2 + ( 0.001 )**2 )  # Target function
			# These are the two functions to minimize
			errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
			#errfunc2 = lambda p1, p2, x, y: (fitfunc2(p1, p2, x) - y) # Distance to the target function
			errfunc2 = lambda p1, p2, c, x, y: (fitfunc2(p1, p2, x) - y)*np.exp(-(x-c)**2) # Distance to the target function
			# Parameters array: Intensity, Omega_p, InverseLifetime, ...
			#p0 = [50., omegaps[ik,ib], 0.1, 50., omegaps[ik,ib] / 4, 0.1, 50., omegaps[ik,ib] / 8, 0.1] # Initial guess for the parameters
			tempim = np.zeros(np.size(newen))
			interpims = interp1d(en, ims[ik,ib,:], kind = 'linear', axis =  2)
			tempim[:] = interpims(newen[:])
			tempim2 = tempim[tempim>0]
			en2 = newen[tempim>0]
			der = np.zeros(np.size(en2))
			#der2 = np.zeros(np.size(en2))
			der[0] = (tempim2[1]-tempim2[0])/(en2[1]-en2[0])
			der[1] = (tempim2[2]-tempim2[0])/(en2[2]-en2[0])
			der[-1] = (tempim2[-1]-tempim2[-2])/(en2[-1]-en2[-2])
			der[-2] = (tempim2[-1]-tempim2[-3])/(en2[-1]-en2[-3])
			#der2[0] = (der[1]-der[0])/(en2[1]-en2[0])
			#der2[-1] = (der[-1]-der[-2])/(en2[-1]-en2[-2])
			for i in range(2,np.size(en2)-1) :
				#der[i] = (ims[ik,ib,i+1]-ims[ik,ib,i-1])/(en[i+1]-en[i-1])
				der[i] = (tempim2[i+1]-tempim2[i-1])/(en2[i+1]-en2[i-1])
			#	der2[i-1] = (der[i]-der[i-2])/(en2[i]-en2[i-2])
			thr = eqp[ik,ib]-2*omegaps[ik,ib]/5
			#plt.plot(en,tempim,'o',label='data-orig')
			#plt.plot(en[tempim>0],tempim[tempim>0],'o',label='data')
			# Calculate parameters for small lorentzian
			p0 = [10., 0.1, omegaps[ik,ib] / 4 ] # Initial guess for the parameters
			#p1, success = optimize.leastsq(errfunc, p0, args=(en[en> thr ], ims[ik,ib,en> thr ]))
			p1, success = optimize.leastsq(errfunc, p0, args=( en2[en2> thr ], tempim2[en2> thr ]))
			#p00 = [50., 0.1, omegaps[ik,ib], p1[0], p1[1], p1[2]] # Initial guess for the parameters
			# Calculate parameters for big lorentzian, the small being fixed
			p00 = [50., 0.1, omegaps[ik,ib]] # Initial guess for the parameters
			#p11, success = optimize.leastsq(errfunc2, p00, args=(p1, en[en< thr ], ims[ik,ib,en< thr ]))
			p11, success = optimize.leastsq(errfunc2, p00, args=(p1, eqp[ik,ib] - omegaps[ik,ib], en, ims[ik,ib,:]))
			## Calculate again parameters for small lorentzian, the big being fixed
			#p00 = [50., 0.1, omegaps[ik,ib]] # Initial guess for the parameters
			#p1, success = optimize.leastsq(errfunc2, p1, args=(p11, eqp[ik,ib] - p1[2], en2[en2> thr ], tempim2[en2> thr ]))
			print " Starting parameters (Intensity, Inv_Lifetime, omega_p):"
			print p0
			print p00
			print " Optimized parameters:"
			print p1
			print p11
			lomegap1[ik,ib] = p1[-1]
			lomegap2[ik,ib] = p11[-1]
			# Calculate integral of curves
			surf_small = np.trapz(fitfunc(p1, en2[en2>eqp[ik,ib]-2*lomegap1[ik,ib]]),en2[en2>eqp[ik,ib]-2*lomegap1[ik,ib]]) / np.pi
			surf_big = np.trapz(fitfunc2(p11, p1, en),en) / np.pi
			#print " Integral of Lorentzian:", np.trapz(fitfunc2(p11, p1, en),en) / np.pi
			print " Integral of Lorentzian:", surf_big
			print " Integral of data:", np.trapz(ims[ik,ib],en) / np.pi
			a1[ik,ib] = surf_small/lomegap1[ik,ib]**2
			a2[ik,ib] = ( np.trapz(ims[ik,ib],en) / np.pi - surf_small ) / lomegap2[ik,ib]**2
			print " a1:", a1[ik,ib]
			print " a2:", a2[ik,ib]
			outfile1.write("%14.5f" % (lomegap1[ik,ib]))
			outfile2.write("%14.5f" % (lomegap2[ik,ib]))
			#plt.figure(ik)
			plt.plot(en[en<thr],ims[ik,ib,en<thr],'x',label='data')
			plt.plot(en2[en2> thr ],tempim2[en2> thr ],'o-',label='cut data')
			#plt.plot(en,fitfunc(p0,en),'-',label='starting guess')
			plt.plot(newen2,fitfunc2(p11, p1, newen2),'-',label='fit2')
			plt.plot(newen2[newen2> thr ],fitfunc(p1,newen2[newen2> thr ]),'-',label='fit1')
			plt.plot(en2[en2> thr ],der[en2> thr ],'-',label='der')
			#plt.plot(en2,der2,'-',label='der2')
			#plt.plot(en,errfunc(p1,en,tempim),'--',label='errfunc')
		outfile1.write("\n")
		outfile2.write("\n")
	outfile1.close()
	outfile2.close()
	# Writing out a1 e a2
	outname = "a1.dat"
	outfile = open(outname,'w')
	outname = "a2.dat"
	outfile2 = open(outname,'w')
	for ik in xrange(nkpt):
		for ib in xrange(nband):
			outfile.write("%14.5f"  % (a1[ik,ib]))
			outfile2.write("%14.5f" % (a2[ik,ib]))
		outfile.write("\n")
		outfile2.write("\n")
	outfile.close()
	outfile2.close()
	print " a1 array:",
	print str(a1[0])+'...'
	print " a2 array:",
	print str(a2[0])+'...'
	#plt.legend(loc=2)
# END ROADWORKS #
### ===================================== ###

### ===== EXPONENTIAL SPECTRAL FUNCTION ====== ###
print " ### Calculation of exponential A...  "
plt.figure(2)
#starten = -90.
#enden = enmax
dxexp = 0.05
enexp = np.arange(enmin,enmax,dxexp)
spfexp =  np.zeros((nkpt,nband,np.size(enexp)))
spfexptot =  np.zeros(np.size(enexp))
spfexptot0 =  np.zeros(np.size(enexp))
spfexptot1 =  np.zeros(np.size(enexp))
spfexptot2 =  np.zeros(np.size(enexp))
# 
# plasmon2 = 0 --> 1-plasmon calculation
# plasmon2 = 1 --> 2-plasmon calculation
#plasmon2 = 1
impl  =  3 # Inverse plasmon lifetime (Single plasmon-pole model)
impl1 =  3 # Inverse 1st plasmon lifetime (Double plasmon-pole model)
impl2 =  3 # Inverse 2nd plasmon lifetime (Double plasmon-pole model)
# First the case for a single plasmon (silicon)
if plasmon2 == 0 and extinf == 0 :
	print " Plasmon mode:", plasmon2
	order0 =  np.zeros(np.size(enexp))
	order1 =  np.zeros(np.size(enexp))
	order2 =  np.zeros(np.size(enexp))
	order3 =  np.zeros(np.size(enexp))
	for ik in xrange(nkpt):
		print " k point = %02d " % (ik+1)
		for ib in xrange(nband):
			outnamekb = "spf_exp-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+".dat"
			outfilekb = open(outnamekb,'w')
			prefac = np.exp(-lambdas[ik,ib]) * wtk[ik] * pdos[ib] / np.pi
			order0[:] = abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order1[:] = lambdas[ik,ib] * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 1.* omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order2[:] = lambdas[ik,ib]**2 / 2 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 2.* omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order3[:] = lambdas[ik,ib]**3 / 6 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 3.* omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			spfexp[ik,ib,:] = prefac * ( order0[:] + order1[:] + order2[:] + order3[:] )
			spfexptot0[:] += prefac * order0[:]
			spfexptot1[:] += prefac * ( order0[:] + order1[:] )
			spfexptot2[:] += prefac * ( order0[:] + order1[:] + order2[:])
			spfexptot[:] += spfexp[ik,ib,:] 
			plt.plot(eqp[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 1.*omegaps[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 2.*omegaps[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 3.*omegaps[ik,ib],0.01,'o')#label='point')
			for ien in xrange(np.size(enexp)) :
				outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], spfexp[ik,ib,ien]))
			outfilekb.close()
			continue
elif plasmon2 == 0 and extinf == 1 :
	print " Plasmon mode:", plasmon2
	print " Extrinsic and interference effect included."
	order0 =  np.zeros(np.size(enexp))
	order1 =  np.zeros(np.size(enexp))
	order2 =  np.zeros(np.size(enexp))
	order3 =  np.zeros(np.size(enexp))
	for ik in xrange(nkpt):
		print " k point = %02d " % (ik+1)
		for ib in xrange(nband):
			akb = lambdas[ik,ib] + a_extinf
			prefac = np.exp( - akb ) * wtk[ik] * pdos[ib] / np.pi
			order0[:] = abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order1[:] = akb        * abs( imeqp[ik,ib] + 1.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 1.* omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 1.* impl )
			order2[:] = akb**2 / 2 * abs( imeqp[ik,ib] + 2.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 2.* omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 2.* impl )
			order3[:] = akb**3 / 6 * abs( imeqp[ik,ib] + 3.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 3.* omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 3.* impl )
			spfexp[ik,ib,:] = prefac * ( order0[:] + order1[:] + order2[:] + order3[:] )
			spfexptot0[:] += prefac * order0[:]
			spfexptot1[:] += prefac * ( order0[:] + order1[:] )
			spfexptot2[:] += prefac * ( order0[:] + order1[:] + order2[:])
			spfexptot[:] += spfexp[ik,ib,:] 
			plt.plot(eqp[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 1.*omegaps[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 2.*omegaps[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 3.*omegaps[ik,ib],0.01,'o')#label='point')
			continue
# This is the case for 2 plasmons (graphite)
elif plasmon2 == 1 and extinf == 0 : 
	print " Plasmon mode:", plasmon2
	order0 =  np.zeros(np.size(enexp))
	order10 =  np.zeros(np.size(enexp))
	order20 =  np.zeros(np.size(enexp))
	order30 =  np.zeros(np.size(enexp))
	order40 =  np.zeros(np.size(enexp))
	order01 =  np.zeros(np.size(enexp))
	order02 =  np.zeros(np.size(enexp))
	order03 =  np.zeros(np.size(enexp))
	order04 =  np.zeros(np.size(enexp))
	order11 =  np.zeros(np.size(enexp))
	order21 =  np.zeros(np.size(enexp))
	order12 =  np.zeros(np.size(enexp))
	order13 =  np.zeros(np.size(enexp))
	order22 =  np.zeros(np.size(enexp))
	order31 =  np.zeros(np.size(enexp))
	order41 =  np.zeros(np.size(enexp))
	order32 =  np.zeros(np.size(enexp))
	order23 =  np.zeros(np.size(enexp))
	order14 =  np.zeros(np.size(enexp))
	order42 =  np.zeros(np.size(enexp))
	order33 =  np.zeros(np.size(enexp))
	order24 =  np.zeros(np.size(enexp))
	order43 =  np.zeros(np.size(enexp))
	order34 =  np.zeros(np.size(enexp))
	order44 =  np.zeros(np.size(enexp))
	spfexptot1p1 =  np.zeros(np.size(enexp))
	spfexptot1p2 =  np.zeros(np.size(enexp))
	for ik in xrange(nkpt):
		print " k point = %02d " % (ik+1)
		for ib in xrange(nband):
			outnamekb = "spf_exp-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+"-p2.dat"
			outfilekb = open(outnamekb,'w')
			prefac = np.exp( - a1[ik,ib] - a2[ik,ib] ) * wtk[ik] * pdos[ib] / np.pi
			# TEST
			#if interpims(newen[izerores[ik,ib]]) < 0 : print " WARNING: im(Sigma(eps_k)) < 0 !!! ik ib eps_k im(Sigma(eps_k)) = ", ik, ib, newen[izerores[ik,ib]], interpims(newen[izerores[ik,ib]])
			order0[:] = abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			# Diagonal terms
			order10[:] = a1[ik,ib]        * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order20[:] = a1[ik,ib]**2 /2  * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order30[:] = a1[ik,ib]**3 /6  * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order40[:] = a1[ik,ib]**4 /24 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order01[:] = a2[ik,ib]        * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order02[:] = a2[ik,ib]**2 /2  * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order03[:] = a2[ik,ib]**3 /6  * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order04[:] = a2[ik,ib]**4 /24 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			# Mixed terms now
			order11[:] = a1[ik,ib]    * a2[ik,ib]         * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order21[:] = a1[ik,ib]**2 * a2[ik,ib]    /  2 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order12[:] = a1[ik,ib]    * a2[ik,ib]**2 /  2 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order22[:] = a1[ik,ib]**2 * a2[ik,ib]**2 /  4 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order13[:] = a1[ik,ib]    * a2[ik,ib]**3 /  6 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order31[:] = a1[ik,ib]**3 * a2[ik,ib]    /  6 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order41[:] = a1[ik,ib]**4 * a2[ik,ib]    / 24 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order32[:] = a1[ik,ib]**3 * a2[ik,ib]**2 / 12 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order23[:] = a1[ik,ib]**2 * a2[ik,ib]**3 / 12 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order14[:] = a1[ik,ib]    * a2[ik,ib]**4 / 24 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order42[:] = a1[ik,ib]**4 * a2[ik,ib]**2 / 48 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order33[:] = a1[ik,ib]**3 * a2[ik,ib]**3 / 36 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order24[:] = a1[ik,ib]**2 * a2[ik,ib]**4 / 48 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order43[:] = a1[ik,ib]**4 * a2[ik,ib]**3 /144 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order34[:] = a1[ik,ib]**3 * a2[ik,ib]**4 /144 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			order44[:] = a1[ik,ib]**4 * a2[ik,ib]**4 /576 * abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			# Sum up everything
			spfexp[ik,ib,:] = prefac * ( order0[:] + order01[:] + order10[:] + order20[:] + order02[:] + order30[:] + order03[:] + order40[:] + order04[:] + order11[:] + order21[:] + order12[:] + order22[:] + order13[:] + order31[:] + order41[:] + order32[:] + order23[:] + order14[:] + order42[:] + order33[:] + order24[:] + order43[:] + order34[:] + order44[:] )
			spfexptot0[:] += prefac * order0[:]
			spfexptot1[:] += prefac * ( order0[:] + order01[:] + order10[:] )
			spfexptot2[:] += prefac * ( order0[:] + order01[:] + order10[:] + order20[:] + order02[:])
			spfexptot1p1[:] += prefac * ( order0[:] + order10[:] + order20[:] + order30[:] + order40[:] )
			spfexptot1p2[:] += prefac * ( order0[:] + order01[:] + order02[:] + order03[:] + order04[:] )
			spfexptot[:] += spfexp[ik,ib,:] 
			plt.plot(eqp[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 1.*lomegap2[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 2.*lomegap2[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 3.*lomegap2[ik,ib],0.01,'o')#label='point')
			for ien in xrange(np.size(enexp)) :
				outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], spfexp[ik,ib,ien]))
			outfilekb.close()
			continue
elif plasmon2 == 1 and extinf == 1 : 
	print " Plasmon mode:", plasmon2
	print " Extrinsic and interference effect included."
	order0 =  np.zeros(np.size(enexp))
	order10 =  np.zeros(np.size(enexp))
	order20 =  np.zeros(np.size(enexp))
	order30 =  np.zeros(np.size(enexp))
	order40 =  np.zeros(np.size(enexp))
	order01 =  np.zeros(np.size(enexp))
	order02 =  np.zeros(np.size(enexp))
	order03 =  np.zeros(np.size(enexp))
	order04 =  np.zeros(np.size(enexp))
	order11 =  np.zeros(np.size(enexp))
	order21 =  np.zeros(np.size(enexp))
	order12 =  np.zeros(np.size(enexp))
	order13 =  np.zeros(np.size(enexp))
	order22 =  np.zeros(np.size(enexp))
	order31 =  np.zeros(np.size(enexp))
	order41 =  np.zeros(np.size(enexp))
	order32 =  np.zeros(np.size(enexp))
	order23 =  np.zeros(np.size(enexp))
	order14 =  np.zeros(np.size(enexp))
	order42 =  np.zeros(np.size(enexp))
	order33 =  np.zeros(np.size(enexp))
	order24 =  np.zeros(np.size(enexp))
	order43 =  np.zeros(np.size(enexp))
	order34 =  np.zeros(np.size(enexp))
	order44 =  np.zeros(np.size(enexp))
	spfexptot1p1 =  np.zeros(np.size(enexp))
	spfexptot1p2 =  np.zeros(np.size(enexp))
	for ik in xrange(nkpt):
		print " k point = %02d " % (ik+1)
		for ib in xrange(nband):
			#outnamekb = "spf_exp-k"+str("%02d"%(ik+1))+"-b"+str("%02d"%(ib+1))+"-p2.dat"
			#outfilekb = open(outnamekb,'w')
			a1kb = a1[ik,ib] + a_extinf
			a2kb = a2[ik,ib] + a_extinf
			prefac = np.exp( - a1kb - a2kb ) * wtk[ik] * pdos[ib] / np.pi
			# TEST
			#if interpims(newen[izerores[ik,ib]]) < 0 : print " WARNING: im(Sigma(eps_k)) < 0 !!! ik ib eps_k im(Sigma(eps_k)) = ", ik, ib, newen[izerores[ik,ib]], interpims(newen[izerores[ik,ib]])
			order0[:] = abs( imeqp[ik,ib] ) / ( ( enexp[:] - eqp[ik,ib] )**2 + ( imeqp[ik,ib] )**2 )
			# Diagonal terms
			order10[:] = a1kb        * abs( imeqp[ik,ib] + 1.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] )**2 + ( imeqp[ik,ib] + 1.* impl1 )**2 )
			order20[:] = a1kb**2 /2  * abs( imeqp[ik,ib] + 2.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] )**2 + ( imeqp[ik,ib] + 2.* impl1 )**2 )
			order30[:] = a1kb**3 /6  * abs( imeqp[ik,ib] + 3.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] )**2 + ( imeqp[ik,ib] + 3.* impl1 )**2 )
			order40[:] = a1kb**4 /24 * abs( imeqp[ik,ib] + 4.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] )**2 + ( imeqp[ik,ib] + 4.* impl1 )**2 )
			order01[:] = a2kb        * abs( imeqp[ik,ib] + 1.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] + 1.* impl2 )**2 )
			order02[:] = a2kb**2 /2  * abs( imeqp[ik,ib] + 2.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] + 2.* impl2 )**2 )
			order03[:] = a2kb**3 /6  * abs( imeqp[ik,ib] + 3.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] + 3.* impl2 )**2 )
			order04[:] = a2kb**4 /24 * abs( imeqp[ik,ib] + 4.* impl ) / ( ( enexp[:] - eqp[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] + 4.* impl2 )**2 )
			# Mixed terms now
			order11[:] = a1kb    * a2kb         * abs( imeqp[ik,ib] + 1.* impl1 + 1.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 1.* impl1 + 1.* impl2 )
			order21[:] = a1kb**2 * a2kb    /  2 * abs( imeqp[ik,ib] + 2.* impl1 + 1.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 2.* impl1 + 1.* impl2 )
			order12[:] = a1kb    * a2kb**2 /  2 * abs( imeqp[ik,ib] + 1.* impl1 + 2.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 1.* impl1 + 2.* impl2 )
			order22[:] = a1kb**2 * a2kb**2 /  4 * abs( imeqp[ik,ib] + 2.* impl1 + 2.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 2.* impl1 + 2.* impl2 )
			order13[:] = a1kb    * a2kb**3 /  6 * abs( imeqp[ik,ib] + 1.* impl1 + 3.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 1.* impl1 + 3.* impl2 )
			order31[:] = a1kb**3 * a2kb    /  6 * abs( imeqp[ik,ib] + 3.* impl1 + 1.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 3.* impl1 + 1.* impl2 )
			order41[:] = a1kb**4 * a2kb    / 24 * abs( imeqp[ik,ib] + 4.* impl1 + 1.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] + 1.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 4.* impl1 + 1.* impl2 )
			order32[:] = a1kb**3 * a2kb**2 / 12 * abs( imeqp[ik,ib] + 3.* impl1 + 2.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 3.* impl1 + 2.* impl2 )
			order23[:] = a1kb**2 * a2kb**3 / 12 * abs( imeqp[ik,ib] + 2.* impl1 + 3.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 2.* impl1 + 3.* impl2 )
			order14[:] = a1kb    * a2kb**4 / 24 * abs( imeqp[ik,ib] + 1.* impl1 + 4.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 1.* lomegap1[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 1.* impl1 + 4.* impl2 )
			order42[:] = a1kb**4 * a2kb**2 / 48 * abs( imeqp[ik,ib] + 4.* impl1 + 2.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] + 2.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 4.* impl1 + 2.* impl2 )
			order33[:] = a1kb**3 * a2kb**3 / 36 * abs( imeqp[ik,ib] + 3.* impl1 + 3.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 3.* impl1 + 3.* impl2 )
			order24[:] = a1kb**2 * a2kb**4 / 48 * abs( imeqp[ik,ib] + 2.* impl1 + 4.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 2.* lomegap1[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 2.* impl1 + 4.* impl2 )
			order43[:] = a1kb**4 * a2kb**3 /144 * abs( imeqp[ik,ib] + 4.* impl1 + 3.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] + 3.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 4.* impl1 + 3.* impl2 )
			order34[:] = a1kb**3 * a2kb**4 /144 * abs( imeqp[ik,ib] + 3.* impl1 + 4.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 3.* lomegap1[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 3.* impl1 + 4.* impl2 )
			order44[:] = a1kb**4 * a2kb**4 /576 * abs( imeqp[ik,ib] + 4.* impl1 + 4.* impl2 ) / ( ( enexp[:] - eqp[ik,ib] + 4.* lomegap1[ik,ib] + 4.*  omegaps[ik,ib] )**2 + ( imeqp[ik,ib] )**2 + 4.* impl1 + 4.* impl2 )
			# Sum up everything
			spfexp[ik,ib,:] = prefac * ( order0[:] + order01[:] + order10[:] + order20[:] + order02[:] + order30[:] + order03[:] + order40[:] + order04[:] + order11[:] + order21[:] + order12[:] + order22[:] + order13[:] + order31[:] + order41[:] + order32[:] + order23[:] + order14[:] + order42[:] + order33[:] + order24[:] + order43[:] + order34[:] + order44[:] )
			spfexptot0[:] += prefac * order0[:]
			spfexptot1[:] += prefac * ( order0[:] + order01[:] + order10[:] )
			spfexptot2[:] += prefac * ( order0[:] + order01[:] + order10[:] + order20[:] + order02[:])
			spfexptot1p1[:] += prefac * ( order0[:] + order10[:] + order20[:] + order30[:] + order40[:] )
			spfexptot1p2[:] += prefac * ( order0[:] + order01[:] + order02[:] + order03[:] + order04[:] )
			spfexptot[:] += spfexp[ik,ib,:] 
			plt.plot(eqp[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 1.*lomegap2[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 2.*lomegap2[ik,ib],0.01,'o')#label='point')
			plt.plot(eqp[ik,ib] - 3.*lomegap2[ik,ib],0.01,'o')#label='point')
			#for ien in xrange(np.size(enexp)) :
			#	outfilekb.write("%8.4f %12.8f\n" % (enexp[ien], spfexp[ik,ib,ien]))
			#outfilekb.close()
			continue
print " ### Writing out A(\omega)_exp...  "
outname = "spftot_exp"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+".dat"
if plasmon2 == 1 and extinf == 0 : 
	outname1 = "spftot_exp"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+"-1p1.dat"
	outname2 = "spftot_exp"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+"-1p2.dat"
	outname  = "spftot_exp"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+"-p2.dat"
	outfile1 = open(outname1,'w')
	outfile2 = open(outname2,'w')
elif plasmon2 == 1 and extinf == 1 :
	outname = "spftot_exp"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+"-p2.extinf.dat"
elif extinf == 1 : # i.e. and plasmon2 == 0
	outname = "spftot_exp"+"_sfac"+str(sfac)+"_pfac"+str(pfac)+"_pen"+str(penergy)+".extinf.dat"
outfile = open(outname,'w')
for i in xrange(np.size(enexp)):
	outfile.write("%7.4f   %15.10e\n"% (enexp[i],spfexptot[i])) # Dump string representations of arrays
outfile.close()
print " A(\omega)_exp written in", outname
if  plasmon2 == 1 and extinf == 0 : 
	for i in xrange(np.size(enexp)):
		outfile1.write("%7.4f   %15.10e\n"% (enexp[i],spfexptot1p1[i])) # Dump string representations of arrays
		outfile2.write("%7.4f   %15.10e\n"% (enexp[i],spfexptot1p2[i])) # Dump string representations of arrays
	outfile1.close()
	outfile2.close()
	print " A(\omega)_exp (1p1) written in", outname1
	print " A(\omega)_exp (1p2) written in", outname2
# Now go back to original directory
print " Moving back to parent directory:\n", origdir
chdir(newdir)

plt.plot(newen,spftot,'-', label='spftot_gw')
plt.plot(enexp,spfexptot,'-',label='exptot')
if  plasmon2 == 1 : 
	plt.plot(enexp,spfexptot1p1,'-x',label='exptot1p1')
	plt.plot(enexp,spfexptot1p2,'-+',label='exptot1p2')
plt.plot(enexp,spfexptot2,'-',label='exptot2')
plt.plot(enexp,spfexptot1,'-',label='exptot1')
plt.plot(enexp,spfexptot0,'-',label='exptot0')
plt.title('Spectral function '+ 'A (' + r'$\omega $ ' + ') - '+r'$ h\nu = $'+str(penergy)+' eV')
		
#plt.plot(en[izerores],res[0,0,izerores],'o')
# === PLOT PARAMETERS === #
#plt.xlim(-60,10)
#plt.ylim(0,5)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);
plt.legend(loc=2);
plt.show();


