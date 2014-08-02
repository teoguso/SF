#!/usr/bin/env python
"""
Written by Matteo Guzzo. 
This script should read single ik,ib-state 
files and sum them up as needed. 
"""

def sum_kb_numpy(mink,maxk,minb,maxb,npoles=100,penergy=0,extinf=0):
	"""
	Here we actually sum up the spectra and 
	spit them out of the box (nicely).
	"""
	import numpy as np
	spftot = None
	for ik in xrange(mink,maxk+1):
		for ib in xrange(minb,maxb+1):
			#print " ik, ib:", ik,ib
			energy = []
			spfkb = []
			if extinf == 1:
				name="spf_exp-k"+str("%02d"%(ik))+"-b"+str("%02d"%(ib))+"_np"+str(npoles)+"_extinf."+str(penergy)
			else:
				name="spf_exp-k"+str("%02d"%(ik))+"-b"+str("%02d"%(ib))+"_np"+str(npoles)+".dat"
			spf_file = open(name,"r")
			for line in spf_file:
				lfloats = map(float,line.split())
				energy.append(lfloats[0])
				spfkb.append(lfloats[1])
			skb = np.array(spfkb)
			if spftot is not None:
				spftot = spftot + skb
			else:
				spftot = skb
			spfkb = None
			#del spfkb
	energy = np.array(energy)
	return energy,spftot

def sum_kb(mink,maxk,minb,maxb,npoles=100,penergy=0,extinf=0):
	"""
	Here we actually sum up the spectra and 
	spit them out of the box (nicely).
	"""
	#import numpy as np
	spftot = None
	for ik in xrange(mink,maxk+1):
		for ib in xrange(minb,maxb+1):
			#print " ik, ib:", ik,ib
			energy = []
			spfkb = []
			if extinf == 1:
				name="spf_exp-k"+str("%02d"%(ik))+"-b"+str("%02d"%(ib))+"_np"+str(npoles)+"_extinf."+str(penergy)
			else:
				name="spf_exp-k"+str("%02d"%(ik))+"-b"+str("%02d"%(ib))+"_np"+str(npoles)+"."+str(penergy)
			spf_file = open(name,"r")
			for line in spf_file:
				lfloats = map(float,line.split())
				energy.append(lfloats[0])
				spfkb.append(lfloats[1])
			#skb = np.array(spfkb)
			if spftot is not None:
				for i in xrange(len(spfkb)):
					spftot[i] = spftot[i] + spfkb[i]
			else:
				spftot = spfkb
			spfkb = None
			#del spfkb
	#energy = np.array(energy)
	return energy,spftot

if __name__ == '__main__':
	#	import sys
	#	usage = 'Usage: %s (<inputparamfile)' % sys.argv[0]
#	try:
#		continue
#		#infilename = sys.argv[2]
#		#sigma = float(sys.argv[1])
#		#ifilewrite = 1
#	except:
#		print usage 
#		sys.exit(1) 
	minkpt, maxkpt = map(int,raw_input(" min,max kpt: ").split())
	minbd, maxbd = map(int,raw_input(" min,max band: ").split())
	#minkpt = int(raw_input(" minkpt: "))
	#maxkpt = int(raw_input(" maxkpt: "))
	#minbd  = int(raw_input(" minband: "))
	#maxbd  = int(raw_input(" maxband: "))
	npoles = int(raw_input(" npoles: "))
	extinf = int(raw_input(" extinf(1/0): "))
	penergy = int(raw_input(" photon energy (eV): "))
	print "minkpt,maxkpt", minkpt,maxkpt
	print "minbd,maxbd", minbd,maxbd
	print "npoles", npoles
	print "extinf", extinf
	print "penergy", penergy
	#	name="spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+"_extinf."+str(penergy)
	#else:
	#	name="spf_exp-k"+str("%02d"%(ikeff+1))+"-b"+str("%02d"%(ibeff+1))+"_np"+str(npoles)+".dat"
	from broad import gbroaden
	sigma = 0.4
	if extinf==1:
		energy,spf = sum_kb(minkpt,maxkpt,minbd,maxbd,npoles,penergy)
		en2,spfext = sum_kb(minkpt,maxkpt,minbd,maxbd,npoles,penergy,extinf)
		brof = gbroaden(energy,spf,sigma)
		brof2 = gbroaden(en2,spfext,sigma)
	else:
		energy,spf = sum_kb(minkpt,maxkpt,minbd,maxbd,npoles,penergy)
		brof = gbroaden(energy,spf,sigma)
	#print len(energy),len(spf)
	
	#for en, sp in zip(energy,spf):
	#	print ("%10.8e %10.8e")% (en, sp)
	flag_filewrite=1
	if flag_filewrite==1:
		outname = "sumkb_k_"+str(minkpt)+"_"+str(maxkpt)+"_b_"+str(minbd)+"_"+str(maxbd)+"_np"+str(npoles)+"."+str(penergy)
		print " Output file:", outname
		outfile = open(outname,'w')
		outfile.write("### energy intrinsic extrinsic int_broad ext_broad \n")
		if extinf == 1:
			for ene,fint,fext,fib,feb in zip(energy,spf,spfext,brof,brof2):
				outfile.write("%12.9e %12.9e %12.9e %12.9e %12.9e\n" % (ene,fint,fext,fib,feb))
			outfile.close()
		else:
			for ene,fint,fib in zip(energy,spf,brof):
				outfile.write("%12.9e %12.9e %12.9e\n" % (ene,fint,fib))
			outfile.close()
	flag_plot=1
	if flag_plot==1:
		import matplotlib.pylab as plt
		plt.plot(energy,spf,label="Intrinsic")
		if extinf == 1:
			plt.plot(en2,spfext,label="Extinf")
		plt.plot(energy,brof,label="Int_broad")
		if extinf == 1:
			plt.plot(en2,brof2,label="Extinf_broad")
		plt.legend(loc="upper left")
		plt.show()
# 1 11 17 (18) 19
