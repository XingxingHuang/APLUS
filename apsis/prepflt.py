#!/usr/bin/env python
# $Id: bandcount.py,v 1.00 2012/11/02 20:15:45 Huang  $
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.00 $ '[11:-3]
__version_date__ = '$Date: 2012/11/13 $ '[7:-3]
__author__       = 'Xingxing Huang'

import pdb,os,sys,glob,popen2
from optparse import OptionParser
import xmlUtil,fUtil
import pyfits,numpy

class fltset:

    """ usage:
	   flt=fltset()
	   flt.flc2flt(dir,list)
    """

    def __init__(self,dir,list=None):
        import os, glob
        self.dir = dir
        self.list = list
	if not self.list:
	    self.list = self.flclist()

    def flclist(self):
	curdir = os.getcwd()
	os.chdir(self.dir)
	flclist = glob.glob('*_flc.fits')
	os.chdir(curdir)
        return flclist

    def flc2flt(self,outdir):
	print ">>>start to define the name"
	for files in self.list:
	    #pdb.set_trace()
	    fileslist = files.split('.')
	    if fileslist[0][-1] != "c":
		print "ERROR : check the name of "+files
	    fltname=fileslist[0][:-1]+"t.fits"
            
	    fltfile = os.path.join(outdir,fltname)
	    flcfile = os.path.join(self.dir,files)
            if os.path.isfile(fltfile):
		print fltname," exists!"
	    else:
		cmd = "ln -s "+flcfile+' '+fltfile
		sproc  = popen2.Popen3(cmd,1)
		print cmd
                print fltname," be created"
		
		

usage = 'usage: %prog sourcename [outputdir]'
version = '%prog 1.0   by xingxing'
parser = OptionParser(usage=usage, version=version)
(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_help()
    sys.exit()
dsname = args[0]
outdir=''
if len(args) == 2: 
    outdir = args[1]


flclist=[]
ingestdir = os.path.join(os.environ.get('INGEST'),dsname)
os.chdir(ingestdir)
flclist = glob.glob('*_flc.fits')           # get the list
flt=fltset(ingestdir)
if flclist:
    print ">>>*flc.fits  under $INGEST: ",len(flclist)
    if outdir:
        print ">>>Put output files into "+outdir+'\n'
	flt.flc2flt(outdir)
    else:
        flt.flc2flt(ingestdir)                  # change the files
    print "done!"
else:
    print "No *flc.fits"
