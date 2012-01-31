#!/usr/bin/env python
"""
Multipole fit module. It is callable as a stand-alone script. 
"""


def getdata_file(infilename,wantedcol=1):
	"""
	This function opens a file with a given filename
	and puts (as default) the first two 
	columns into two numpy arrays
	that are returned. 
	With an optional keyword (wantedcol) one 
	can choose to wrap the n-th column of the file.
	"""
	import numpy as np
	infile = open(infilename)
	preen = []
	predata = []
	ncol = wantedcol 
	for lines in infile : 
		#print lines[0]
		if lines[0] != "#" :
			line = map(float,lines.split())
			#print line
			preen.append(line[0])
			predata.append(line[ncol])
	infile.close()
	preen = np.array(preen)
	predata = np.array(predata)
	return preen, predata

def first_inverse_moment(preen,predata):
	import numpy as np
	nbin = np.size(preen)
	fxonx = np.zeros(nbin)
	tmpen = preen
	tmpdata = predata
	dx = abs((tmpen[-1] - tmpen[0])/nbin) # or nbin + 1 ???
	for i in xrange(nbin) : 
		if tmpen[i] == 0. : 
			tmpen[i] += dx/1000
			#izero = i
			#tmpdata[i] = 0
			#print " Avoiding zero:", preen[i], predata[i] 
			break
	#tmpen[tmpen == 0] += dx/1000
	fxonx = tmpdata / tmpen
	#if izero is not None: preen[izero] = 0.
	return fxonx

def resize_en(preen, nbin) :
	"""
	This is a quite primitive solution. 
	I am not sure it will work in any situation.
	TODO: think of a smarter grid analyser.
	"""
	import numpy as np
	nbin = int(nbin)
	if np.size(preen) < float( 2 * nbin ) :
		#print " X-axis grid is too coarse for so many poles."
		print " Refining grid..."
		nx = 2*nbin+1
		print " Old dx = %g, new dx = %g." % (abs(preen[-1]-preen[0])/(np.size(preen)-1),abs(preen[-1]-preen[0])/nx)
		en = np.linspace(preen[0], preen[-1], nx)
	else :
		en = preen
	return en

def fit_multipole(preen,predata,nbin,ifilewrite=0,binmode=0):
	"""
	This function fits a curve given by some dataset (preen,predata) 
	with a given number of poles (nbin).
	preen is supposed to be all positive and with increasing
	numbers starting from its first element (not decreasing numbers). 
	It can write out a file with the calculated parameters if the 
	appropriate flag is equal to 1 (ifilewrite).
	It returns omegai, gi, deltai.
	"""
	#import matplotlib.pylab as plt
	import numpy as np
	import sys
	from scipy.interpolate import interp1d
	from scipy.integrate import fixed_quad
	nbin = int(nbin)
	eta = 0.005 # This is the Lorentzian broadening that would be used???
	safe_shift = 10. # This is a little trick to avoid x=zero which introduces errors. 
	preen = preen + safe_shift
	totalint = np.trapz(predata,preen)
	totdeltax = abs( preen[-1] - preen[0] )
	print " Totdeltax, np.size(preen), dx:", totdeltax, np.size(preen), ( preen[-1] - preen[0] ) / float( np.size(preen) - 1 )
	print " Number of poles (nbin):", nbin
	print " Total integral:", totalint 
	print " Total integral / nbin:", totalint / float(nbin)
	# This is the supposed integral within a single interval
	partint = totalint / float(nbin)
	interpdata = interp1d(preen, predata[:], kind = 'linear', axis =  2)
	# First moment
	xfx = np.zeros(np.size(preen))
	xfx = predata * preen
	interpxfx = interp1d(preen, xfx[:], kind = 'linear', axis =  2)
	# First inverse moment
	fxonx = first_inverse_moment(preen,predata)
	interpfxonx = interp1d(preen, fxonx[:], kind = 'linear', axis =  2)
	# Test plot
	#plt.plot(preen,interpdata(preen),label="data")
	#plt.plot(preen,interpxfx(preen),'-x',label="f(x)*x")
	#plt.plot(preen,(predata*preen),label="f(x)*x")
	#plt.plot(preen,interpfxonx(preen),label="f(x)/x")
	#plt.plot(preen,(predata/preen),'-x',label="f(x)/x")
	### =========================================== ###
	# Here we calculate the bins' bounds
	# First we want the x-axis grid to be finer than the density of poles
	en = resize_en(preen, nbin)
	bdensity = 1
	iboundcheck = 1
	while iboundcheck == 1:
		#print " ### ========================= ###"
		#print " ###   Calculating Delta_i     ###"
		bounds = []
		ibound = 0
		gi = []
		omegai = []
		istart = 0
		funcint = 0
		bounds.append(en[istart])
		x0 = en[istart]
		x1 = x0
		# Number of gaussians used in the integration (40 is safe, lower values are not tested too well)
		ngaussint = 20
		print " Getting poles...",
		for i in xrange(1,np.size(en)) : 
			x2 = en[i]
			x3 = ( x1 + x2 ) / 2
			(tmpint,dummy) = fixed_quad(interpdata, x0, x2,(), ngaussint)
			if tmpint > partint :
				(tmpint,dummy) = fixed_quad(interpdata, x0, x3,(), ngaussint)
				while abs( ( tmpint - partint ) / partint ) > 1E-06 :
					if tmpint < partint :
						x1 = x3
						x3 = ( x3 + x2 ) / 2
					else :
						x2 = x3
						x3 = ( x1 + x3 ) / 2
					(tmpint,dummy) = fixed_quad(interpdata, x0, x3,(), ngaussint)
				#print " Bound found. ibound, en, tmpint:", ibound, x3, tmpint
				bounds.append( x3 )
				# Formula to calculate g_i
				tmpint2, dummy = fixed_quad(interpfxonx, x0, x3,(), ngaussint)
				#gi.append(float( 2. / np.pi * tmpint2))
				tmpgi = 2. / np.pi * tmpint2 
				gi.append( tmpgi )
				# Formula to calculate omega_i
				tmpint, dummy = fixed_quad(interpxfx, x0, x3,(), ngaussint)
				tmpomegai = np.sqrt( 2/np.pi * tmpint / tmpgi ) 
				omegai.append( tmpomegai )
				#print " gi, omegai:", tmpgi, tmpomegai # Test print
				# Reset variables
				x0 = x3
				istart = i
				ibound += 1
		print "Done."
		# Add last value as the upper bound
		print " ibound       = %4i (should be %g) " % (ibound, nbin)
		print " Size(bounds) = %4i (should be %g) " % (np.size(bounds), nbin+1)
		print " Size(omegai) = %4i (should be %g) " % (np.size(omegai), nbin)
		if ibound == nbin : # i.e. if it is as it is supposed to be
			print " ibound == nbin, Fixing last value"
			bounds[-1] = en[-1]
			bounds = np.array(bounds)
			iboundcheck = 0
		# Prevent approximate integration to miss g_i and omega_i for last interval
		elif ibound == nbin - 1 :  # i.e. supposedly there is one bound missing ( ibound == nbin -1 )
			print " ibound == nbin - 1. Calculating parameters for last bin..."
			bounds.append(en[-1])
			bounds = np.array(bounds)
			x0 = bounds[-2]
			x3 = bounds[-1]
			# Formula to calculate g_i
			tmpint2, dummy = fixed_quad(interpfxonx, x0, x3,(), ngaussint)
			tmpgi = 2. / np.pi * tmpint2 
			gi.append( tmpgi )
			# Formula to calculate omega_i
			tmpint, dummy = fixed_quad(interpxfx, x0, x3,(), ngaussint)
			tmpomegai = np.sqrt( tmpint / tmpint2 ) 
			omegai.append( tmpomegai )
			print " gi, omegai:", tmpgi, tmpomegai 
			ibound += 1
			iboundcheck = 0
		else :
			print " ibound has a non-compliant value. ibound:", ibound
			iboundcheck = 1
			bdensity*=2
			print " Will retry with a higher bdensity:", bdensity
			en = resize_en(preen, bdensity*nbin)
			#sys.exit(1)
	omegai = np.array(omegai)
	gi = np.array(gi)
	# Here we assign the value as f is a sum of delta with one coefficient only (no pi/2 or else)
	gi = np.pi/2*gi*omegai
	# Here we restore the correct x axis removing safe_shift
	omegai = omegai - safe_shift
	# Uncomment to change weights to single delta function
	# gi = gi / abs(omegai)
	#print " bounds[-5:]:", bounds[-5:]
	#print " omegai[-5:]:", omegai[-5:]
	#print " Bounds:", bounds
	#print " Omegai:", omegai
	deltai = []
	sumcheck = 0
	print " Calculating deltai..."
	for i in xrange(1,np.size(bounds)) :
		deltai.append(bounds[i]-bounds[i-1])
		sumcheck += abs(bounds[i]-bounds[i-1])
	deltai = np.array(deltai)
	print " Check if sum of deltai gives the original length: ", sumcheck,
	if abs((sumcheck - abs(en[-1] - en[0])) / sumcheck) > 1E-03: 
		print
		print en[-1] - en[0]
		print "WARNING: the difference is", abs((sumcheck - abs(en[-1] - en[0])) / sumcheck)
	else: print "(OK)"
	#intcheck = np.pi/2*np.sum(gi[:]*omegai[:])
	intcheck = np.sum(gi)
	print " Check if sum of gi gives the original total integral (origint): ", intcheck, totalint
	print " ibound       = %4i (should be %g) " % (ibound, nbin)
	print " Size(bounds) = %4i (should be %g) " % (np.size(bounds), nbin+1)
	print " Size(omegai) = %4i (should be %g) " % (np.size(omegai), nbin)
	print " Size(deltai) = %4i (should be %g) " % (np.size(deltai), nbin)
	if ifilewrite == 1:
		# Print a file like Josh output
		# omega_i  gamma_i  g_i    delta_i
		#    .        .       .       .
		#    .        .       .       .
		#    .        .       .       .
		outname = "poles."+str(nbin)+".dat"
		outfile = open(outname,'w')
		# print 2-line header
		outfile.write("### number of poles: %g\n" % (nbin))
		outfile.write("### omega_i  gamma_i  g_i    delta_i \n")
		for j in xrange(np.size(omegai)) :
			outfile.write("%12.8f %12.8f %12.8f %12.8f\n" % (omegai[j], eta, gi[j], deltai[j]))
		outfile.close()
		print " Parameters written in file", outname
	return omegai, gi, deltai

def fit_multipole2(preen,predata,nbin,ifilewrite=0,binmode=0):
	"""
	Another version with constant bin separation.
	"""
	import numpy as np
	import sys
	from scipy.interpolate import interp1d
	from scipy.integrate import fixed_quad
	nbin = int(nbin)
	eta = 0.005 # This is the Lorentzian broadening that would be used???
	totalint = np.trapz(predata,preen)
	totdeltax = abs( preen[-1] - preen[0] )
	print " Totdeltax, np.size(preen), dx:", totdeltax, np.size(preen), ( preen[-1] - preen[0] ) / float( np.size(preen) - 1 )
	return omegai, gi, deltai

def write_f_as_sum_of_poles(preen,omegai,gi,deltai,ifilewrite):
	"""
	An attempt to rewrite the original function given the set of poles and weights
	"""
	import numpy as np
	def lorentz_delta(x, a, b):
		"""
		This is the lorentzian representation of a delta function. 
		x is the independent variable, a is the pole and b is the damping. 
		This is exact for b going to 0. 
		This small function can work also with numpy arrays. 
		"""
		from numpy import pi
		florentz = b / pi / ((x - a)**2 + b**2)
		return florentz
	eta = 0.001
	#for eni in preen
	#interpdata = interp1d(preen, predata[:], kind = 'linear', axis =  2)
	# Mh... Here I used a rule-of-thumb procedure.
	# TODO: think of a consistent and robust way of defining the x grid.
	eta = 0.05
	nbin = 8*np.size(omegai)/eta/100
	en = resize_en(preen, nbin)
	f = np.zeros(np.size(en))
	f1 = np.zeros(np.size(en))
	for i in xrange(np.size(omegai)):
	#for oi in omegai:
		#f[:] += 2/np.pi * gi[i] * omegai[i]**2 * lorentz(en[:]**2,omegai[i]**2,eta)
		#for j in xrange(np.size(en)):
		f1[:] += 2.* gi[i] / omegai[i] * lorentz_delta(en[:]**2, omegai[i]**2, eta)
		f[:] += gi[i] * lorentz_delta(en[:], omegai[i], eta)
	if ifilewrite == 1:
		oname = "mpfit."+str(np.size(gi))+".dat"
		ofile = open(oname,'w')
		# print 2-line header
		ofile.write("### Function reconstructed following many-pole fit. \n")
		ofile.write("### x   f(x) = \sum_i^N | \omega_i | \eta / ( ( \omega - \omega_i )^2 + \eta^2 ) \n")
		for j in xrange(np.size(en)):
			ofile.write("%15.8f %15.8f %15.8f\n" % (en[j], f[j], f1[j]))
		ofile.close()
		print " Sum of poles written in file", oname
	return en, f

if __name__ == '__main__':
	import sys
	usage = 'Usage: %s npoles infile' % sys.argv[0]
	try:
		infilename = sys.argv[2]
		nbin = sys.argv[1]
		ifilewrite = 1
	except:
		print usage
		sys.exit(1)
	preen, predata = getdata_file(infilename)
	omegai, gi, deltai = fit_multipole(preen,predata,nbin,ifilewrite)
	write_f_as_sum_of_poles(preen,omegai,gi,deltai,1)
