#!/usr/bin/env python
# $Id: errorout.py,v 1.00 2012/11/19 20:15:45 Huang  $
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.00 $ '[11:-3]
__version_date__ = '$Date: 2012/11/19 $ '[7:-3]
__author__       = 'Xingxing Huang'

import pdb,os,sys,string
import time

class errorout:
    """
    make a file called "err.log" to include some error messages.
    
    """
    def __init__(self,obs,cmd=None,filters=None):
    	self.obs	 = obs
    	self.root 	 = obs.root
    	self.obsName     = obs.newobs
    	self.cmd=cmd
    	self.filters=filters
    
    def printerr(self):
        j=string.find(self.root,self.obsName)
        pname=self.root[0:j] + 'err.log'
        format=str('%10s %15s %30s')
        if os.path.exists(pname):
            twf = open(pname,'a')
        else:
            twf = open(pname,'w')
            txt='This file is used for recording the band which ran fail \n'
            print >> twf, txt
            del txt
	    print >> twf, format % ("name","command","time")
        ptime=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime(time.time()))
        obs=self.obsName+'_'+self.filters
        print >> twf, format % (obs,self.cmd,ptime)
        
