#!/usr/bin/env python
from __future__ import print_function
import sys
import time
import numpy as np
from numpy import random
import matplotlib.pylab as plt;
#from scipy.signal import convolve, fftconvolve

def gbroaden(x,f,sigma):
	x = np.asarray(x)
	f = np.asarray(f)
	xsize = np.size(x)
	broadf = np.zeros(f.size)
	print("gbroaden() called, x.size:", x.size)
	#dx =  x[1] - x[0] 
        dx =  ( x[-1] - x[0] ) / float(xsize - 1)
	print("gbroaden() called, gaussian window size (4 x sigma):", 4*sigma)
	for n in xrange(xsize):
		if  abs( x[n] - x[0] ) > 4*sigma : 
			nrange = n
			print("gbroaden() called, gaussian window size (4 x sigma) in x units (e.g. eV?):", x[nrange] - x[0])
			break
	if 4*sigma < dx*xsize/2 :
		# First chunk (beginning)
		for n in xrange(0,nrange):
#	#		print "Processing... "+str(int(n/xsize)*100)+"\r",
			for m in xrange(n-nrange,n+nrange):
#	#			gaub = np.exp( - ( x[n] - x[m] )**2 / 2 / sigma**2 ) / np.sqrt(2*np.pi) / sigma
				gaub = dx * np.exp( - ( dx * ( m - n ) )**2 / 2 / sigma**2 ) / np.sqrt(2*np.pi) / sigma
				if m >= 0 : 
					broadf[n] += f[m] * gaub
				else : 
					broadf[n] += f[0] * gaub
		# Last chunk (end)
		for n in xrange(xsize-nrange,xsize):
#	#		print "Processing... "+str(int(nrange+n/xsize)*100)+"\r",
			for m in xrange(n-nrange,n+nrange):
#	#			gaub = np.exp( - ( x[n] - x[m] )**2 / 2 / sigma**2 ) / np.sqrt(2*np.pi) / sigma
				gaub = dx * np.exp( - ( dx * ( m - n ) )**2 / 2 / sigma**2 ) / np.sqrt(2*np.pi) / sigma
				if m < xsize : 
					broadf[n] += f[m] * gaub
				else : 
					broadf[n] += f[-1] * gaub
		# Middle chunk (treated with the standard formula)
		for n in xrange(nrange,xsize-nrange):
#	#		print "Processing... "+str(int(n/xsize)*100)+"\r",
			for m in xrange(n-nrange,n+nrange):
				gaub = dx * np.exp( - ( dx * ( m - n ) )**2 / 2 / sigma**2 ) / np.sqrt(2*np.pi) / sigma
				#gaub = dx * np.exp( - ( x[n] - x[m] )**2 / 2 / sigma**2 ) / np.sqrt(2*np.pi) / sigma
				broadf[n] += f[m] * gaub
	else :
		# Standard formula regardless of boundaries (it should work decently in all cases)
		for n in xrange(xsize):
#	#		print "Processing... "+str(int(n/xsize)*100)+"\r",
			for m in xrange(xsize):
				gaub = dx * np.exp( - ( dx * ( m - n ) )**2 / 2 / sigma**2 ) / np.sqrt(2*np.pi) / sigma
				#gaub = dx * np.exp( - ( x[n] - x[m] )**2 / 2 / sigma**2 ) / np.sqrt(2*np.pi) / sigma
				broadf[n] += f[m] * gaub
#	af = np.trapz(f)
#	abf = np.trapz(broadf)
#	print af, abf
#	broadf = broadf * af / abf
	return broadf
def broad_file(filename):
    """
    Takes the name of a file and gives two broadened versions of it.
    """
    import numpy as np
   #infilename = sys.argv[-1]
    print("broad_file :: ")
    print("Filename: ", filename)
    infilename = filename
    file1 = open(infilename,'r')
    x = []; y = []
    for line in file1:
    	if line[0] != "#" : # Remove commented lines
    		data = line.split()
    		x.append(float(data[0]))
    		y.append(float(data[1]))
    x = np.array(x)
    y = np.array(y)
    
    yb = gbroaden(x,y,sigma)
    yb2 = gbroaden(x,y,2*sigma)
   #plt.plot(x,yb,label='broad f(x)')
   #plt.plot(x,yb2,label='2*broad f(x)')
    outfilename1 = infilename+"."+str(sigma)+".gbr"
    outfilename2 = infilename+"."+str(2*sigma)+".gbr"
    outfile1 = open(outfilename1,'w')
    outfile2 = open(outfilename2,'w')
    for (en, br1, br2) in zip(x,yb,yb2) :
    	outfile1.write("%7.4f%15.6f\n" % (en, br1))
    	outfile2.write("%7.4f%15.6f\n" % (en, br2))
    outfile1.close()
    outfile2.close()
    return 0

if __name__ == '__main__':
	# ASK FOR A BROAD VALUE IF NOT ALREADY GIVEN IN COMMAND LINE
	if  len(sys.argv) == 2:
		sigma = float(raw_input("Gaussian broadening needed: "))
		print("What I read:", sigma)
	else: sigma = float(sys.argv[1])
        for item in sys.argv[2:]:
            broad_file(item)
	
	### UNCOMMENT IF PLOTTING IS NEEDED
	#plt.legend()
	#plt.show()
