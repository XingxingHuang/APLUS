#!/usr/bin/env python
# $Id: bandcount.py,v 1.00 2012/11/02 20:15:45 Huang  $
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.00 $ '[11:-3]
__version_date__ = '$Date: 2012/11/02 $ '[7:-3]
__author__       = 'Xingxing Huang'


import pdb,os
import xmlUtil,fUtil,string
import pyfits
import numpy

class bandcount:
    def __init__(self,obs,band):
	 #self.sciImageList = obs
	 self.sciImageList = obs.sciImageList
	 #self.obsFits      = '/home/zheng/aaplusa_runs/datasets/a209/a209/Images' #obs.fitsdir
	 self.obsFits      = obs.fitsdir
	 self.bandlist     = []
	 self.band         =  band
	 if self.band=='nir' or self.band=='red':
	     Detector='IR'
	 elif self.band=='opt':
	     Detector='WFC'
	 elif self.band=='uv':
            Detector='UVIS'
	 else:
            print "Incorrect band name: "+Band
            exit()		
	 self.det = Detector

    def bandcounts(self):
       """Count the filter number in band (['opt','nir','red','uv']) .
	"""
        # pdb.set_trace()

       for fitsfile in self.sciImageList:
           # Use user-defined detector for the detection image
           det = fUtil.getDetector(os.path.join(self.obsFits,fitsfile))
           if det == self.det:
               if (self.band =='red'):     
                   if string.find(fitsfile,'160') > -1 or string.find(fitsfile,'140') > -1:
			  yes = 1
			  self.bandlist.append(fitsfile)
                   else:
                       print 'Not including '+fitsfile+' in detection image.'
                       yes = 0
               else:
		      yes = 1
		      self.bandlist.append(fitsfile)
           else:
               print 'Not including '+fitsfile+' in detection image.'
               yes = 0
       print len(self.bandlist),' filters are detected under '+self.band
       print self.bandlist

