#!/usr/bin/env python
"""
Written by Matteo Guzzo 2012.
Function necessary to the main script for the calculation 
of spectral functions. 
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
# TODO: write a function to read the parameters from not ad-hoc files
	import numpy as np;
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
		print " Auxiliary file not found (wtk.dat). Impossible to continue."
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
		print " Auxiliary file not found (occ.dat). "
		print " Setting all occupations to 2.0."
		occ = 2.0*np.ones((maxkpt,maxband-minband+1))
	return occ                             	

def read_sigfile2(sigfilename,enmin,enmax,minkpt,maxkpt,minband,maxband):
	"""                                    	
	This function reads the real and imaginary parts of the self energy
	$\Sigma(\omega)$ from the file _SIG for the given bands and k points.
	It returns numpy arrays containing the energies, the real and the 
	imaginary part of the self energy.
	"""
	import numpy as np;
	if isfile(sigfilename):
		insigfile = open(sigfilename);
	else:
		print "File "+sigfilename+" not found."
		insigfile = open(raw_input("Self-energy file name (_SIG): "))
	# We put the content of the file (lines) in this array
	filelines = insigfile.readlines();
	en=[];
	# loop to prepare  the energy array
	print " Reading array of energies from first k-point in _SIG file... ",
	trigger=0
	for line in filelines:
		if line[0:3]=='# k':
			if trigger==1: break # So as to read just the first k point
			continue
		elif line[0:3]=='# b':
			firstband = int(line.split()[-2])
			if firstband>minband: 
				print "ERROR: first band available in _SIG file is higher than minband. Please check."
				sys.exit(1)
			lastband =  int(line.split()[-1])
			if lastband<maxband: 
				print "ERROR: last band available _SIG file is lower than maxband. Please check."
				sys.exit(1)
			trigger = 1
			continue
		else : 
			data=map(float,line.split())
			en.append(data[0])
			#if data[0] >= enmin and data[0] < enmax :
			#	en.append(data[0])
			#elif data[0] < enmin: continue
			#else: break
	print "Done."
	print " Length of the energy array detected from _SIG file, first k point: "+str(len(en))
	print " len(en): "+str(len(en))
	en = np.array(en)
	print " size(en): "+str(np.size(en))
	dx = (en[-1]-en[0])/np.size(en)
	print " dx:",dx
	res = []
	ims = []
	print "en[0], en[-1]"
	print en[0], en[-1]
	ik=-1
	# Loop where we actually read Sigma
	for line in filelines:
		if line[0:3] == "# k": 
			nline=0
			ik=ik+1
			print " kpt # %02d" % (ik+1)
			print line,
			res.append([])
			ims.append([])
			continue
		elif line[0:3] == "# b": print line,; continue
		#elif float(line.split()[0])>=enmin or float(line.split()[0])<enmax:
		else:
			tmplist = map(float,line.split())
			del tmplist[0]
			ib = 0
			for i in xrange(len(tmplist)):
				if i%3==0 and nline==0: 
					res[ik].append([])
					res[ik][ib].append(tmplist[i])
				elif  i%3==0:
					res[ik][ib].append(tmplist[i])
				elif (i-1)%3==0 and nline==0 : 
					ims[ik].append([])
					ims[ik][ib].append(tmplist[i])
				elif (i-1)%3==0 : 
					ims[ik][ib].append(tmplist[i])
				elif (i-2)%3==0: 
					ib = ib+1
					continue
			nline+=1
		#else:  continue
	res = np.array(res)
	ims = np.array(ims)
	ik2=0
	ib2=0
	dum1 = np.zeros((maxkpt-minkpt+1,maxband-minband+1,np.size(en)))
	dum2 = np.zeros((maxkpt-minkpt+1,maxband-minband+1,np.size(en)))
	nband = lastband-firstband+1
	for ik in xrange(minkpt-1,maxkpt):
		ib2=0
#		for ib in xrange(minband-1,maxband):
		for ib in xrange(minband-firstband,nband-(lastband-maxband)):
				#print "ik, ib, ik2, ib2, minkpt, maxkpt, minband, maxband", ik, ib, ik2, ib2, minkpt, maxkpt, minband, maxband
				dum1[ik2,ib2]=res[ik,ib]
				dum2[ik2,ib2]=ims[ik,ib]
				ib2+=1
		ik2+=1
	res = dum1
	ims = dum2
	return en, res, ims

def read_sigfile_OLD(sigfilename,enmax,minkpt,maxkpt,minband,maxband):
	"""                                    	
	This function reads the real and imaginary parts of the self energy
	$\Sigma(\omega)$ from the file _SIG for the given bands and k points.
	It returns numpy arrays containing the energies, the real and the 
	imaginary part of the self energy.
	"""
	import numpy as np;
	if isfile(sigfilename):
		insigfile = open(sigfilename);
	else:
		print "File "+sigfilename+" not found."
		insigfile = open(raw_input("Self-energy file name (_SIG): "))
	# We put the content of the file (lines) in this array
	filelines = insigfile.readlines();
	en=[];
	# loop to prepare  the energy array
	print " Reading array of energies from first k-point in _SIG file... ",
	for line in filelines:
		if line[0:3]=='# k':
			continue
		elif line[0:3]=='# b':
			continue
		else : 
			data=line.split()
			if float(data[0]) <= enmax :
				en.append(float(data[0]))
			else: break
	print "Done."
	print " Length of the energy array detected from _SIG file, first k point: "+str(len(en))
	print " len(en): "+str(len(en))
	en = np.array(en)
	print " size(en): "+str(np.size(en))
	dx = (en[-1]-en[0])/np.size(en)
	print " dx:",dx
	print " ### ===== Reading _SIG file ... ===== #### "
	# Preparation of arrays
	nkpt = maxkpt-minkpt+1
	nband = maxband-minband+1
	res=np.zeros((nkpt,nband,np.size(en)));
	ims=np.zeros((nkpt,nband,np.size(en)));
	print " np.shape(res), np.shape(ims):", np.shape(res), np.shape(ims)
	# Indexes initialization
	ikpt = 0;
	ib = 0;
	io = 0 # row/energy counter
	# Cycle over the lines of the file
	for line in filelines:
		icol = 0;
		ib = 0;
		if line[0:3]=='# k' :
			#print line,
			# Detect, in commented lines, the k point declaration
			if ikpt<minkpt-1 :
				#print " --- k point:  %02i ---" % (ikpt);
				ikpt = ikpt + 1;
				continue
			elif ikpt==maxkpt :
				print " End of the reading loop: ikpt == maxkpt. ikpt, maxkpt: ", ikpt, maxkpt; 
				print " #### ================================================ #### ";
				break;
			else:
				ikpt = ikpt + 1;
		# Detect, in commented lines, the bands declaration
		elif line[0:3]=='# b':
			io = 0;
			continue
		# TODO: This if test is incorrect: this way it always starts from ib = 0
		#elif io < np.size(en) and ib < maxband and ikpt >=minkpt-1:
		elif io < np.size(en) and ikpt >=minkpt-1:
			data=map(float,line.split());
			for col in data:
				# First element (energy) goes into en
				if icol == 0 : 
					icol = icol + 1;
					io = io + 1;
				elif (icol+2)%3 == 0 and ib>=minband-1:
					res[ikpt-minkpt,ib-minband,io-1] = col;
					icol = icol + 1;
					continue;
				elif (icol+2)%3 == 0 :
					icol = icol + 1;
					continue;
				elif (icol+1)%3 == 0 and ib>=minband-1:
					ims[ikpt-minkpt,ib-minband,io-1] = col;
					icol = icol + 1;
					continue;
				elif (icol+1)%3 == 0 :
					icol = icol + 1;
					continue;
				else : 
					if ib == maxband-1 : break # number of bands reached
					ib = ib + 1;
					icol = icol + 1;
					continue;
		else: continue
	#plt.plot(en,res[0,0]); plt.plot(en,ims[0,0]);plt.show();sys.exit(1)
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
	import numpy as np;
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
		pevenfile = open("p_even.dat",'r')
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

def calc_spf_gw(minkpt,maxkpt,minband,maxband,wtk,pdos,en,enmin,enmax,res,ims,hartree):
	"""
	Macro-function calling instructions necessary to calculate 
	the GW spectral function. 
	For now it writes out the single state spectral functions
	on ascii files. This should be moved to an external module
	and just return spfkb as an output variable. 
	spf (GW spectral function) is returned.
	"""
	import numpy as np;
	nkpt = maxkpt-minkpt+1
	nband = maxband-minband+1
	newdx = 0.01
	if enmin < en[0] and enmax >= en[-1]:  
		newen = np.arange(en[0],en[-1],newdx)
	elif enmin < en[0]:  
		newen = np.arange(en[0],enmax,newdx)
	elif enmax >= en[-1] :  
		newen = np.arange(enmin,en[-1],newdx)
	else :  
		newen = np.arange(enmin,enmax,newdx)
	print " ### Interpolation and calculation of A(\omega)_GW...  "
	spftot = np.zeros((np.size(newen)));
	# Here we interpolate re and im sigma
	# for each band and k point
	for ik in xrange(nkpt):
		ikeff = minkpt+ik-1
		print " k point = %02d " % (ikeff+1)
		for ib in xrange(0,nband):
			ibeff = minband+ib-1
			interpres = interp1d(en, res[ik,ib], kind = 'linear', axis =  2)
			interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis =  2)
	                tmpres = interpres(newen)
			#redenom = newen + efermi - hartree[ik,ib] - interpres(newen)
			redenom = newen - hartree[ik,ib] - interpres(newen)
			#print "ik ib minband maxband ibeff hartree[ik,ib]", ik, ib, minband, maxband, ibeff, hartree[ik,ib]
	                tmpim = interpims(newen)
	                spfkb = wtk[ikeff] * pdos[ib] * abs(tmpim)/np.pi/(redenom**2 + tmpim**2)
			spftot += spfkb 
			outnamekb = "spf_gw-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+".dat"
			outfilekb = open(outnamekb,'w')
			for ien in xrange(np.size(newen)) :
				outfilekb.write("%8.4f %12.8e %12.8e %12.8e %12.8e\n" % (newen[ien], spfkb[ien], redenom[ien], tmpres[ien], tmpim[ien]))
			outfilekb.close()
	return newen, spftot

def find_eqp_resigma(en,resigma):
	"""
	This function is supposed to deal with the plasmaron problem 
	and calculate the quasiparticle energy once it is fed with 
	resigma = \omega - \epsilon_H - \Re\Sigma. 
	It expects an array of increasing values on the x axis 
	and it will return 
	the x value of the last resigma=0 detected. 
	It should return the value of eqp and the number of zeros
	found (useful in case there are plasmarons or for debugging). 
	"""
	import numpy as np
	import matplotlib.pylab as plt
	nzeros=0
	tmpeqp = en[0]
	for i in xrange(1,np.size(resigma)):
		#print resigma[i]*resigma[i-1] # DEBUG
		if  resigma[i] == 0: # Yes, it can happen
			tmpeqp = en[i] 
			nzeros+=1
		elif (resigma[i]*resigma[i-1] < 0):
			tmpeqp = en[i-1] - resigma[i-1]*(en[i] - en[i-1])/(resigma[i] - resigma[i-1]) # High school formula
			nzeros+=1
	if nzeros==0 : print " WARNING: No eqp found! "
	elif nzeros>1 : print " WARNING: Plasmarons! "
	return tmpeqp, nzeros

def calc_eqp_imeqp(nkpt,nband,en,res,ims,hartree):
	"""
	This function calculates qp energies and corresponding
	values of the imaginary part of sigma for a set of
	k points and bands. 
	The function find_eqp_resigma() is used here.
	eqp and imeqp are returned. 
	"""
	import numpy as np;
	eqp = np.zeros((nkpt,nband))
	imeqp = np.zeros((nkpt,nband))
	for ik in xrange(nkpt):
		for ib in xrange(nband):
			temparray = np.array(en - hartree[ik,ib] - res[ik,ib])
			interpims = interp1d(en, ims[ik,ib], kind = 'linear', axis =  2)
			tempim = interpims(en)
			# New method to overcome plasmaron problem
			eqp[ik,ib], nzeros = find_eqp_resigma(en,temparray)
			if nzeros==0: print " ERROR: ik "+str(ik)+" ib "+str(ib)+". No eqp found!!! Bye bye!";sys.exit()
			imeqp[ik,ib] = interpims(eqp[ik,ib])
			## Warning if imaginary part of sigma < 0 (Convergence problems?)
			if imeqp[ik,ib] <= 0 : 
				print " WARNING: im(Sigma(eps_k)) <= 0 !!! ik ib eps_k im(Sigma(eps_k)) = ", ik, ib, eqp[ik,ib], imeqp[ik,ib]
#	plt.show() # DEBUG
#	sys.exit(1) # DEBUG
	return eqp, imeqp

def calc_extinf_corrections(extinfname,ampole,omegampole):
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
	print " Reading extrinsic and interference contribution from file "+str(extinfname)+"..."
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
	interpwidth = interp1d(newen_ei, newwmod, kind = 'linear', axis =  2)
	w_extinf = ampole.copy()
	print "omega_p, a_extinf, a_int:"
	print newen_ei
	print newa_ei
	print newa_ei/newa_int
	#print en_ei, newenexin
	#print aextinf, newaexin
	interpextinf = interp1d(newen_ei, newa_ei/newa_int, kind = 'linear', axis =  2)
	amp_exinf = ampole.copy()
	#print "Type(amp_exinf, ampole):", type(amp_exinf), type(ampole)
	# Mod following discussion with Josh
	amp_mean = np.mean(ampole)
	for ik in xrange(nkpt):
		for ib in xrange(nband):
			#tmpextinf = interpextinf(omegampole[ik,ib])/npoles # <-- Divided by the number of poles (normalization)!
			w_extinf[ik,ib] = interpwidth(omegampole[ik,ib]) # Numpy array
			tmpextinf = interpextinf(omegampole[ik,ib]) # 
			# Mod following discussion with Josh
			#amp_exinf[ik,ib] += ampole[ik,ib] * tmpextinf
			amp_exinf[ik,ib] += amp_mean * tmpextinf
	return amp_exinf, w_extinf

def calc_spf_mpole(enexp,prefac,akb,omegakb,eqpkb,imkb,npoles,wkb=None):
	"""
	This function calculates the exponential spectral function. 
	"""
	import numpy as np;
	ftot = np.zeros((np.size(enexp)))
	tmpf1=np.zeros((nenexp))
	tmpf2=np.zeros((nenexp))
	tmpf3=np.zeros((nenexp))
	if wkb is None:
		wkb = np.zeros(npoles)
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