#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
from __future__ import print_function
import sys,os,glob,re #,time
import numpy as np

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def outvars_parser(alist):
    """
    This function creates a dictionary using the information
    from the 'outvars' section of preprocessed variables 
    in the abinit output. 
    The input is a list of rows from the file, ideally 
    already stripped. As such, it is supposed to be
    a list of strings containing words and/or numbers,
    starting with a string, so as to fit nicely a
    dictionary. 
    All numbers are converted to float. 
    """
    value = []
    key = None
    vars_dict = {}
    init_int = ('i','j','l','m','n') # Follows fortran convention
    print("-"*3, "Parsing outvars...")
   #for x in alist:
   #    print(x)
    while(alist[0] == '' or alist[0] == ' '):
        alist.pop(0)
    if is_number(alist[0][0]):
        raise ValueError
    for row in alist:
        if row == '' or row == ' ':
            continue
        x = row.split()
        #print(x)
        # Remove single letters at beginning
        # and find keys
        if not is_number(x[0]):
            #print(key, value)
            if len(x[0])==1:
                if key: 
                    vars_dict[key]=value
                print("Anomalous line is taken care of:", x)
                x.pop(0)
                key = x.pop(0)
                value = []
            elif key: # i.e. if it's not the first time
                vars_dict[key]=value
                #print(key, value
                key = x.pop(0)
                value = []
            else:
                key = x[0]
                x.pop(0)
        #print(key
        for y in x: 
            if not is_number(y):
                value.append(y)
            elif key[0] in init_int:
                value.append(int(y))
            else:
                value.append(float(y))
    vars_dict[key]=value
    print("-"*3, "Done.")
    print(sorted(vars_dict.keys()))
    return vars_dict
    
def get_keylist(keyword,dictionary):
    """
    Returns the list of keys actually present in
    the dictionary that start with the given 
    keyword. 
    They are sorted!
    E.g.: occ--> occ2,occ3
    """
    dct = dictionary
    word = keyword
    klist = []
    for key in dct:
        #p = re.compile(word)
        p = re.compile(word+'\d+')
        if key == word:
            #print(key, dct[key]
            klist.append(key)
        elif p.match(key):
            #print(key, dct[key]
            klist.append(key)
        #print("get_keylist :: word, klist:", word, klist
        klist.sort()
    if not klist == []: return klist

def chk_gw(chunkoffile):
    """
    Checks if the given list (a dataset, 
    hopefully) contains a line 
    indicating a gw calculation.
    """
    f = chunkoffile
    response = False
    for line in f:
        if "SIGMA" in line: 
            response = True
            break
    return response

class CodeOutReader(object):
    """
    This class processes the output of an electronic
    structure code and 
    ideally puts it in a form more easily usable 
    within python, using instances and method to 
    access its data. 
    Main feature is that it provides a searchable 
    dictionary containing variables as keys and 
    actual values (as float lists or numbers) as 
    the dictionary's values.
    """
    def __init__(self, code=None, filename=None, is_sc=0):
        """
	Initialises several instances and calls 
	a bunch of test methods. 
	If no file name is given, it takes what 
	it finds in the running directory. 
	"""
#	try: 
        if filename is None:
            self.fname = raw_input("Please give a filename for the output of your GW calculation: ").strip()
            #self.fname = glob.glob('*.out')[0] # This is 'ls *.out' in current dir
        else: self.fname = filename
        self.code = code
        self.is_sc = is_sc
        # Initialize variables
        # Most of them are lists, in that they 
        # iterate over the different datasets.
	# This array should contain the starting 
	# row number for each dataset and then 
	# the number for the 'end dataset' row.
        self.ndtsets = 1
        self.version = None
        self.nversion = 0
        self.completed = False
        self.content = []
        self.var_dict = {}
        self.gw_ks_ev = []
        self.qpen = []
        self.hartree = []
# Call methods
        print(52*"=")
        print(" INITIALIZING QP CALCULATION OUTPUT PARSER... ")
        print(52*"=")
	self.get_file_content()
        print(52*"-")
        print(" QP CALCULATION OUTPUT PROCESSED. ")
        print(52*"-")

    ### METHODS HERE BELOW ###
    def __str__(self):
        """
        Outputs a decently formatted output when you try to print
        the object.
        """
        header=" + == CodeOutReader printout == + "  
        mystring = (
                header+"\n" 
                "    code       " + str(self.code)      + "\n" +
                "    version    " + str(self.version)   + "\n" + 
                "    filename   " + str(self.fname)     + "\n" +
                "    completed  " + str(self.completed) + "\n" +
                header+"\n" 
                )
        return mystring


    def get_scalar(self,label):
        """
        This method returns the scalar value(s) of 
        a requested scalar variable in the dictionary as a list. 
        If units are different than default, it 
        issues a warning. 
        """
        #ndts = self.ndtsets
        #s_e = self.dts_start_end
        dct = self.var_dict
        if dct.get(label):
            value = dct.get(label)[:]
            if type(value[0]) is not float: 
                    print("get_scalar ::","ERROR: This is no scalar.")
                    print("get_scalar ::","See?", value)
                    return None
                    #raise ValueError 
            elif not is_number(value[-1]): 
                print("get_scalar ::", "WARNING: Units used are not default, but", value[-1])
                value.pop()
                return value
            else:
                return value
        else:
            print("get_scalar ::","ERROR:", label, "not found.")
            return None
        #else: 
        #    value.pop()
        return value

    def chk_units(self):
        """
        This method simply checks whether there is a 
        string at the end of the list of values for all
        keys in the dictionary, 
        indicating that the units used for that 
        variable are not the default.
        It issues a warning if a string is found. 
        """
        print("chk_units :: ",)
        dct = self.var_dict
        for key in dct.keys():
            x = dct.get(key)
            if not is_number(x[-1]):
                print("WARNING: Units used for", key, "are not default, but", x[-1])
        print("Done.")

    def reshape_key(self,keylist,size):
        """
        This method/function reshapes a single list 
        relative to a given list of keys to a list of list of 
        a given size.
        It accepts a list of keywords, in case they are supposed 
        to be reshaped the same way. 
        WARNING: A special value of 'size' is "nband", which 
        refers to the relative number of bands for the given 
        quantity (e.g. 'occ'). 
        """
        print("reshape_key :: ")
        dct = self.var_dict
        l = list(keylist)
        if not is_number(size):
            tmp = get_keylist(size,dct)
            for key,nband in zip(l,tmp): # Should raise an error if not same size
                tmplist = [key,] # Argument must be a list
                self.reshape_key(tmplist,dct.get(nband)[0]) #recursion!!!
        else:
            klist = []
            for word in l:
                tmp = get_keylist(word,dct)
                if tmp is not None: 
                    klist.append(tmp)
                else:
                    print("WARNING:", "No key for '", word,"' found.")
            #print("klist:", klist
            for keys in klist:
                #print("key(s) to reshape, len(keys), size:", keys, len(keys), size)
                for elem in keys: 
                    bucket1 = []
                    bucket2 = []
                    value = dct.get(elem)
                    for x in value: 
                        bucket1.append(x)
                        if len(bucket1)==size: 
                            bucket2.append(bucket1) 
                            bucket1 = []
                    dct[elem] = bucket2
        print("Done.")

    def get_file_content(self):
        """
	This method keeps the file open just for the
	strictly necessary time to copy its content.  
	Lines are duly stripped of newlines. 
	"""
        with open(self.fname,'r') as f:
            self.content = [line.strip('\n') for line in f]

#######################################################################

class ExcitingOutReader(CodeOutReader):
    """
    This class processes the exciting output and 
    ideally puts it in a form more easily usable 
    within python, using instances and method to 
    access its data. 
    Main feature is that it provides a searchable 
    dictionary containing variables as keys and 
    actual values (as float lists or numbers) as 
    the dictionary's values.
    """
    def __init__(self, filename = 'EVALQP.DAT', is_sc=0):
        """
	Initialises several instances and calls 
	a bunch of test methods. 
	If no file name is given, it takes what 
	it finds in the running directory. 
	"""
        print(" ExcitingOutReader :: filename: ", filename)
        CodeOutReader.__init__(self, 'exciting', filename, is_sc)
	# Attributes
       #self.dts_start_end = []
       #self.dts_labels = []
       #self.dtsets = []
        # Methods
       #self.chk_fname()
       #self.get_version()
       #flag = self.chk_completed()
       #self.get_dtsets(flag)
       #self.get_outvars()
        self.get_var_dict()
        self.hartree = self.var_dict['hartree']
        # GW STUFF
       #self.get_gw_dts()
       #for i in self.gw_dts: 
       #self.get_gw_ks_ev()
       #self.get_gw_qpen(i,self.nversion)

    def chk_fname(self):
        """
	This method checks whether we are looking 
        a proper abinit output. 
	Also, it checks whether the file exists. 
        If not, it raises an error. 
        """
        if self.fname != "EVALQP.DAT":
	   #print(self.fname[-4:], len(self.fname)
	    print(self.fname, ":")
	    print(" The file given is not the expected exciting output file. ")
	    print(" Please check the file or simply rename it. ")
	    raise ValueError    
	elif os.path.isfile(self.fname): 
	    print(self.fname, ":")
	    print(" File name correct ")
	    print("and it also exists. ")
	else: 
	    print(self.fname, ":")
	    print(" File name correct ")
	    print("but it does not seem to exist. Bye. ")
	    raise ValueError    
	
    ### METHODS HERE BELOW ###
    def get_version(self):
        """
        Detects the version of Exciting (how?)
        This method is not working at present
        """
       #self.version = self.content[1].split()[1]
       #self.nversion = int(self.version.split('.')[0])
        self.version = 666
        self.nversion = 666
        print(" WARNING: Unable to detect code version for Exciting.")

    def get_var_dict(self):
        """
        Set up the dictionary. 
        """
       #self.var_dict = outvars_parser(self.outvars)
        self.var_dict = self.read_evalqp()
       #print(" var_dict: ", self.var_dict)
       #self.chk_units()
       #self.reshape_dict()

    def read_evalqp(self):
        """
        This method reads the EVALQP.DAT file trying to
        extract all the useful information.
        Namely it reads:
        - 
        - 
        - 
        """
        dct = {}
        nkpt = 0
        nband = 0
        kpt = []
        hartree = []
        ik = 0
        ib = 0
        print("read_evalqp :: ")
        box = self.content
        for i,line in enumerate(box):
            words = line.split()
            if 'k-point' in words:
               #print("KAPPA?")
                nkpt += 1
                kpt.append(map(float, words[-3:]))
               #print("line: ", line)
            elif 'state' in words:
                istate = 0
                hartree.append([])
            elif line == ' ':
                nband = istate
                ib = 0
                ik += 1
            else:
                istate = int(words[0])
               #print(line)
                hartree[ik].append(float(words[2])-float(words[4]))
                ib += 1
                pass
        print("nkpt: ", nkpt)
        print("len(kpt): ", len(kpt))
        print("nband: ", nband)
        print("hartree shape (ik,ib): ", len(hartree), len(hartree[0]))
        dct['nkpt'] = nkpt
        dct['kpt'] = kpt
        dct['nband'] = nband
        dct['hartree'] = hartree
        print("read_evalqp :: Done.")
        return dct

class AbinitOutReader(CodeOutReader):
    """
    This class processes the abinit output and 
    ideally puts it in a form more easily usable 
    within python, using instances and method to 
    access its data. 
    Main feature is that it provides a searchable 
    dictionary containing variables as keys and 
    actual values (as float lists or numbers) as 
    the dictionary's values.
    """
    def __init__(self, filename=None, is_sc=0):
        """
	Initialises several instances and calls 
	a bunch of test methods. 
	If no file name is given, it takes what 
	it finds in the running directory. 
	"""
        if filename is None:
            filename = glob.glob('*.out')[0] # This is 'ls *.out' in current dir
        CodeOutReader.__init__(self, 'abinit', filename, is_sc)
	# Attributes
        self.dts_start_end = []
        self.dts_labels = []
        self.dtsets = []
        # Methods
        self.chk_fname()
        self.get_version()
        flag = self.chk_completed()
	self.get_dtsets(flag)
        self.get_outvars()
        self.get_var_dict()
        # GW STUFF
        self.get_gw_dts()
        for i in self.gw_dts: 
            self.get_gw_ks_ev(i)
            self.get_gw_qpen(i,self.nversion)

    def chk_fname(self):
        """
	This method checks whether we are looking 
        a proper abinit output. 
	Also, it checks whether the file exists. 
        If not, it raises an error. 
        """
        if self.fname[-4:] != ".out":
	    print(self.fname[-4:], len(self.fname))
	    print(self.fname, ":",)
	    print(" The file given is not a standard abinit '.out' output. " )
	    print(" Please check the file or simply rename it. " )
	    raise ValueError    
	elif os.path.isfile(self.fname): 
	    print(self.fname, ":",)
	    print(" File name correct ", )
	    print("and it also exists. ")
	else: 
	    print(self.fname, ":")
	    print(" File name correct ")
	    print("but it does not seem to exist. Bye. ")
	    raise ValueError    
	
    ### METHODS HERE BELOW ###
    def get_version(self):
        """
        Detects the version of abinit. 
        """
        self.version = self.content[1].split()[1]
        self.nversion = int(self.version.split('.')[0])

    def get_var_dict(self):
        """
        Set up the dictionary. 
        """
        self.var_dict = outvars_parser(self.outvars)
        self.chk_units()
        self.reshape_dict()

    def chk_units(self):
        """
        This method simply checks whether there is a 
        string at the end of the list of values for all
        keys in the dictionary, 
        indicating that the units used for that 
        variable are not the default.
        It issues a warning if a string is found. 
        """
        print("chk_units :: ")
        dct = self.var_dict
        for key in dct.keys():
            x = dct.get(key)
            if not is_number(x[-1]):
                print("WARNING: Units used for", key, "are not default, but", x[-1])
        print("Done.")

    def reshape_dict(self):
        """
        This modules fixes a selected (user-defined) 
        set of variables into a needed shape,  
        for easier later consumption.
        - k points are reshaped from a 
        mere list of numbers to a list of triplets; 
        - bdgw* is reshaped into a list of pairs; 
        - occ* is reshaped into nkpt* lists of nband* 
        elements each; 
        """ 
        print("reshape_dict :: ")
        dct = self.var_dict
        vectors = ('acell','kpt','kptgw','rprim','symrel','xangst','xcart','xred')
        couples = ('bdgw',)
        occ = ('occ',)
        self.reshape_key(vectors,size=3)
        self.reshape_key(couples,size=2)
        self.reshape_key(occ,size='nband')
        print("Done.")

    def get_outvars(self):
        """
        This copies the selected part of the abinit output file, 
        where input variables are summarized at the beginning 
        of a run, into a list. 
	"""
	istart,iend = 0,0
        f = self.content
        for i,line in enumerate(f):
            #print(line
            if "outvars" in line and "preprocessed" in line:
                #print("OUTVARSSSS!"
                istart = i+1
            if "chkinp" in line:
                iend = i-3 # Distance checked in both abinit 5.5 and 7.10
                break
        #print(f[istart:iend]
        self.outvars = f[istart:iend]

    def get_dtsets(self,completed=True):
        """
	Collects all information on datasets in the output file, 
	such as beginning and end, and the number of them. 
	"""
        n = 0
        f = self.content
	s_e = []
        labels = []
        sets = []
        for i,line in enumerate(f):
            if "= DATASET" in line:
                if s_e: s_e[n-1].append(i-1)
		s_e.append([i])
		labels.append(line.split()[2])
                n += 1
            if "= END DATASET" in line:
		s_e[n-1].append(i)
        # If there is no end dataset:
        if completed is False: s_e[n-1].append(int(len(f)-1))
        print(s_e)
        #if n == 0: raise ValueError
        self.ndtsets = n
	self.dts_start_end = s_e
        #print(s_e)
        self.dts_labels = labels
        for i in range(len(s_e)):
            #print(f)
            sets.append(f[s_e[i][0]:s_e[i][1]])
        #print(sets[0][0])
        self.dtsets = sets
        print("get_dtsets ::", n, "dataset(s) found, start/end, labels:", s_e, labels)
    
    def chk_completed(self):
        """
	Checks whether the abinit output file 
        has been closed successfully,
	i.e. if the calculation has finished. 
	"""
        if self.content: 
	    if " Calculation completed." in self.content[-5:]: 
	        self.completed = True
                print("chk_completed :: ","Abinit calculation completed. ")
                return True
            else: 
                print("chk_completed :: ","WARNING: Calculation was not completed.")
                return False
        else: 
            print("chk_completed :: ","WARNING: abinit output has not been read yet.")

    def get_gw_dts(self):
        """
        This methods detects and stores the dataset index 
        and label (if any) where a gw 
        calculation has been performed. 
        This would be given to get_gw_ks_ev().
        """ 
        print("get_gw_dts ::")
        l = self.dtsets
        dlabels = self.dts_labels
        gwsets = []
        gwlabels = []
        for i,mylist in enumerate(l):
            if chk_gw(mylist): 
                gwsets.append(i)
                gwlabels.append(dlabels[i])
                print("GW calculation found.")
        self.gw_dts = gwsets
        self.gw_dts_labels = gwlabels
        print("Done.")

    def get_gw_ks_ev(self,iset=0):
        """
        FIX IT!!!!!!
        This methods returns a list of nkpt nband-long lists, 
        containing the ks eigenvalues (floats) printed prior to a 
        SIGMA (GW) calculation. 
        It assumes that eigenvalues are listed one k point 
        after the other (referring to 'nkpt') and they span 
        'nband' values through 10 columns per row.
        """
        print("get_gw_ks_ev ::")
        myset = self.dtsets[iset]
        mylabels = self.dts_labels
        thisnband = 'nband'+mylabels[iset]
        dct = self.var_dict
        if not dct.get(thisnband): 
            thisnband = 'nband' # Fallback to nband
        nkpt = dct.get('nkpt')[0]
        nband = dct.get(thisnband)[0]
        tmp_ks_ev = []
        if chk_gw(myset):
            tmp_kpt_ev = []
            for i,line in enumerate(myset):
                if "k  " in line and "eigenvalues" in line:
                    istart = i+1
                    break
            ikpt = 0
            ib = -1
            for line in myset[istart:]:
                words = line.split() 
                if line == '':
                    ib = -1
                elif ib == -1:
                    words.pop(0)
                    if not is_number(words[0]):
                        words.pop(0)
                    ib = 0
                if ib >= 0:
                    for word in words:
                        x = float(word)
                        if ib < (nband-1): 
                            tmp_kpt_ev.append(x)
                            ib += 1
                        else:
                            tmp_kpt_ev.append(x)
                            tmp_ks_ev.append(tmp_kpt_ev)
                            tmp_kpt_ev = [] 
                            ikpt += 1
                            ib = -1
                            break 
                if ikpt == nkpt: 
                    break 
                self.gw_ks_ev = tmp_ks_ev
            print("Done.")
        else: 
            print("Not a valid GW dataset." )
            raise ValueError 

    def get_gw_qpen(self,iset=0,version=6):
        """
        Extracts the quasiparticle energies 
        and the values for Vhartree
        (optionally of all the other quantities)
        from abinit output.
        Optional variables:
        - iset: which dataset in the dataset list
          will be considered.
        - version: which version of abinit has been used.
          From version 6 on, the qp output has an additional
          line containing the imaginary part of values.
          Accordingly, possible values are: 5, 6.
        - is_sc: Depending whether it is a self-consistent
        calculation, the number of columns in the abinit 
        output changes. 
        Expected columns in abinit 5, for G0W0:
         Band     E0      <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E
        and for a self-consistent calculation:
         Band     E_lda   <Vxclda>   E(N-1)  <Hhartree>   SigX  SigC[E(N-1)]    Z     dSigC/dE  Sig[E(N)]  DeltaE  E(N)_pert E(N)_diago
        Abinit 6 and above have not been tested. 
        """
        print("get_gw_qpen ::")
        myset = self.dtsets[iset]
       #print("iset:", iset)
       #print("dtset:", myset[:5], myset[-5:])
        mylabels = self.dts_labels
        bdgwlabel = 'bdgw'+mylabels[iset]
        #print("bdgwlabel:", bdgwlabel)
        is_sc = int(self.is_sc)
        dct = self.var_dict
        if not dct.get(bdgwlabel): 
            print(" No", bdgwlabel, "found. Falling back to bdgw.")
            bdgwlabel = 'bdgw' # Fallback to bdgw
        if 'nsppol' in dct:
            nsppol = int(dct['nsppol'][0])
        else:
            nsppol = 1
        print("nsppol: ",nsppol)
        bdgw = dct.get(bdgwlabel)
        nkpt = len(bdgw)
        nband = int(1 + bdgw[0][-1] - bdgw[0][0]) # Hopefully nband is the same for all kptgw
        #type(bdgw)
#        print("bdgwlabel:", bdgwlabel)
#        print("bdgw:",bdgw)
#        print("nkpt:",nkpt)
#        print("nband:",nband)
        #nband = 1 + dct.get(mybdgw)[-1] - dct.get(mybdgw)[0]
        #print("nband:",nband)
        istart = 0
        for i,line in enumerate(myset):
            if " k = " in line:
                istart = i
                #print("i, k point:", i, line
                break
        ik = 0
        ib = 0
        if is_sc == 0:
            elda_k = [] 
            vxc_k = []
        hartree_k = []
        qpen_k = []
        if nsppol == 2:
            nkpt = 2*nkpt
        if version >= 6.19: 
            nband = 2*nband
            ic = 0 
            qpen_im = []
            qpen_im_k = []
        for line in myset[istart:]:
           #if " k = " in line: 
           #    print(line)
           #    pass
            if ik == nkpt:
                break
            elif " k = " in line: 
               #print("k point/header:", line)
                ib = 0
                lines = []
            elif "Band" in line:
                continue
            elif line == '':
                #print("empty line", line
                continue
            elif "_gap" in line:
                #print("gap line", line
                continue
            elif "energy" in line:
                #print("energy line", line
                continue
            elif ib < nband:
#                print("ik,ib,line:", ik,ib,line
               #print(" ik, ib, nband, bdgw", ik, ib, nband, bdgw)
                words = map(float,line.split())
                lines.append(words)
                ib += 1
            if ib == nband :
                qpen = []
                if version >= 6.2: 
                    imag = lines[1::2]
                    lines = lines[0::2]
                    qpen_im = []
                if is_sc == 1: 
                    hartree = []
                    for line in lines:
                        hartree.append(line[4])
                        qpen.append(line[12])
                    hartree_k.append(hartree)
                    qpen_k.append(qpen)
                    hartree = []
                    qpen = []
                else:
                    elda = []
                    vxc = []
                    for line in lines:
                        elda.append(line[1])
                        vxc.append(line[2])
                        qpen.append(line[9])
                    elda_k.append(elda)
                    vxc_k.append(vxc)
                    qpen_k.append(qpen)
                    elda = []
                    vxc = []
                    qpen = []
                ib = 0
                ik += 1
               #print("ik += 1, ik:", ik)
               #if ik == 1:
               #    break
        if is_sc == 0:
            for a,b in zip(elda_k,vxc_k):
                tmp = []
                for c,d in zip(a,b): 
                    tmp.append(c-d) 
                hartree_k.append(tmp)
       #print(hartree_k[-1])
        if nsppol == 2:
            hartree_k = hartree_k[::2]
            qpen_k = qpen_k[::2]
       #a = np.array(hartree_k)
       #print("a.shape, nband, is_sc",  a.shape, nband, is_sc)
       #print(a[0])
        self.qpen = qpen_k
        self.hartree = hartree_k
        print("Done.")
#        print(np.array(qpen_k))
#        if version >= 6: 
#            print(np.array(qpen_im_k) )
#        print(np.array(elda_k))
#        print(np.array(vxc_k))
#        print(np.array(hartree_k))



#out0 = AbinitOutReader("sp.out")
#out1 = AbinitOutReader("tgw2_4.out")
#out0.get_gw_qpen()
#out1.get_gw_qpen()


######################################################
######################################################

############  GARBAGE FROM HERE ON ###################

######################################################
######################################################
class Garbage:
    """
    Just useless stuff. 
    """
    def get_var_old(self,label,shape='scalar'): 
        """
        This method outputs a requested variable in a 
        properly formatted way, using the variable dictionary. 
        Therefore, it puts out e.g. single value or triples
        for scalar and vectors, respectively. 
        NB: Also single values are lists, even though they are 
        comprised of a singe element. 
        The label is a string corresponding to the variable 
        in the abinit output. 
        """
        dim = 1
        if shape is 'vector': dim = 3
        ndts = self.ndtsets
        s_e = self.dts_start_end
        dct = self.var_dict
        value = []
        if ndts == 1: 
            value.append(dct[label])
        else:
            for dts in s_e[:ndts]: # Loops over dts array
                key = label+dts[1]
                if key in self.var_dict:
                    value.append(dct[key])
                elif label in dct:
                    value.append(dct[label])
                    print("get_var ::", "WARNING:", key, "not found, put to", label, ".")
                else:
                    print("get_var ::", "WARNING:", key, "not found, put to 0.")
        return value

    def get_wtk(self):
        """
	Get all values for wtk.
	TODO: finalize.
	"""
        #wtk = get_var_value("wtk",self.content)
        #wtk = get_var_value("wtk",self.outvars)
        pass

    def get_kpt(self):
        """
        BROKEN!!!
	Gets number and location of k points, 
	possibly for each dataset. 
        We need:
        ndtsets
        nkpt
        TODO: Adapt it also for kptgw?
	""" 
        s_e = self.dts_start_end
        ndts = self.ndtsets
        dct = self.var_dict
        nkpt = self.get_scalar("nkpt")
        #nkptgw = self.get_scalar("nkptgw")
        if ndts == 1: 
            nkpt.append(dct['nkpt'])
            #nkptgw.append(dct['nkptgw'])
        else: 
            for dts in s_e[:ndts]: # Loops over dts array
                key = "nkpt"+dts[1]
                if key in self.var_dict:
                    nkpt.append(dct[key])
                elif "nkpt" in dct:
                    nkpt.append(dct['nkpt'])
                    print("get_kpt ::", "WARNING:", key, "not found, put to nkpt.")
                else:
                    print("get_kpt ::", "WARNING:", key, "not found, put to 0.")
                key = "nkptgw"+dts[1]
                if key in dct:
                    nkptgw.append(dct[key])
                elif "nkptgw" in dct:
                    nkptgw.append(dct['nkptgw'])
                    print("get_kpt ::", "WARNING:", key, "not found, put to nkptgw.")
                else:
                    print("get_kpt ::", "WARNING:", key, "not found, put to 0.")

    def get_kpt_old(self):
        """
	Gets number and location of k points, 
	possibly for each dataset. 
	TODO: FIX IT
	"""
        s_e = self.dts_start_end
	f = self.content[:s_e[0]]
	n = 0
	ngw = 0
	for i in range(self.ndtsets):
	    for row in f:
	        try:
                            if ("nkpt" in row) and (row.split()[0] == 'nkpt'): 
                                n = row.split()[1]
                            if ("nkptgw" in row) and (row.split()[0] == 'nkptgw'): 
                                ngw = row.split()[1]
                except: 
                    raise ValueError
	self.nkpt = n
	self.nkptgw = ngw
	print("nkpt =", self.nkpt)
	print("nkptgw =", self.nkptgw)

    def get_kpt_manydatasets(self):
        """
	Gets number and location of k points, 
	possibly for each dataset. 
        TODO: FIX IT
	"""
        s_e = self.dts_start_end
	f = self.content[:s_e[0]]
	for i in range(self.ndtsets):
	    n = 0
	    ngw = 0
	    for row in f:
	        try:
                           if ("nkpt" in row) and (row.split()[0] == 'nkpt'): 
			       print(row)
                               n = row.split()[1]
                           if ("nkptgw" in row) and (row.split()[0] == 'nkptgw'): 
			       print(row)
                               ngw = row.split()[1]
	        except: 
	            raise ValueError
            self.nkpt = n
            self.nkptgw = ngw
	print("nkpt =", self.nkpt)
	print("nkptgw =", self.nkptgw)

def get_var_value(var,alistofrows):
    """
    BROKEN! DO NOT USE!
    This should a basic function that produces
    a list composed by a variable and its value
    from abinit output, namely from the 'outvars' 
    section at the beginning of the file. 
    The main idea is that variables are strings 
    followed by one or more numbers, until the 
    the next variable (string) comes. 
    """
    #s_e = self.dts_start_end
    f = alistofrows
    if not f: 
        print("Error: empty argument.")
        raise ValueError
    dummy = []
    print("Looking for value(s) of", var)
    for i,row in enumerate(f):
	# Double check insures we don't get stuck 
	# with empty lines. 
        if (var in row) and (row.split()[0] == var): 
            print(var, "found:", row)
            x = f[i+1].split()[0] # i.e. first string of next row
            if is_number(x) is True:
                for number in f[i+1].split(): 
                    dummy.append(float(number)) 
            else: print("BOOOOOOM")
    if not dummy: print("WARNING: variable not found.")
    return dummy

def outvars_parser_old(alist):
    """
    This function creates a dictionary using the information
    from the 'outvars' section of preprocessed variables 
    in the abinit output. 
    The input is a list of rows from the file, ideally 
    already stripped. As such, it is supposed to be
    a list of strings containing words and/or numbers,
    starting with a string, so as to fit nicely a
    dictionary. 
    All numbers are converted to float. 
    """
    value = []
    key = None
    vars_dict = {}
    init_int = ('i','j','k','l','m','n') # Follows fortran convention
    print("-"*3,"Parsing outvars...",)
    for row in alist:
        x = row.split()
        #print(x)
        for y in x:
            if not is_number(y): # A keyword is found
                if value and len(value) == 1: 
                    if key[0] in init_int:
                        vars_dict[key]=int(value[0])
                    else: 
                        vars_dict[key]=value[0] 
                    value = [] # Reset list
                elif value and key[0] in init_int:
                    vars_dict[key]=map(int,value)
                    value = []
                elif value:
                    vars_dict[key]=value
                    value = []
                else: 
                    key = y
            else: 
                value.append(float(y))
                print("KEY", key)
    vars_dict[key]=value
    print(" Done.")
    return vars_dict
    
