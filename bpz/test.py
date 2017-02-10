""" 
    Tests the program bpz.py using the HDFN redshift sample
    Usage:
    python test.py -u -b -n
    Options:
    -u skips checking the package "useful"
    -b skips checking the package "bpz_tools"
"""

import sys
import os
from Numeric import *
from MLab import *


#Define BPZPATH 
try:
    user_defined_bpzpath=os.environ['PATH']
except:
    user_defined_bpzpath=None

os.environ['BPZPATH']=os.getcwd()


usef,btools=1,1
for i in range(1,len(sys.argv)):
    if sys.argv[i][:2]=='-u': usef=0
    elif sys.argv[i][:2]=='-b': btools=0 
    else: print 'what!?'

if usef:
    #import useful
    print "TESTING MODULE useful.py"
    print ""
    print ""
    os.system('python useful.py')
    #useful.test()
    
if btools:
    #import bpz_tools
    print "TESTING MODULE bpz_tools"
    print ""
    print ""
    os.system('python bpz_tools.py')
    #bpz_tools.test()


from bpz_tools import *

print "TESTING bpz.py with the Fernandez-Soto et al. HDF-N spectroscopic sample"


os.chdir('TEST')
if os.path.exists('hdfn_z.bpz'):
    os.system('rm hdfn_z.bpz')
    
    
print 'BPZPATH=',os.environ['BPZPATH']
os.environ['BPZ']=os.environ['BPZPATH']+'/bpz.py'
print 'bpz=',os.environ['BPZ']

os.system('python '+os.environ['BPZ']+' hdfn_z.cat -CHECK yes -INTERP 2 -MAG no -PLOTS yes -INTERACTIVE yes -SPECTRA CWWSB.list -DZ 0.005')

os.chdir('..')

if user_defined_bpzpath: os.environ['BPZPATH']=user_defined_bpzpath

