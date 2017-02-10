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
	print "Define the name"
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
