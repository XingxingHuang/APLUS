#!/usr/bin/env python
# $Id: PickCalalog.py,v 1.00 2012/12/05 20:15:45 Huang  $
#choose catalog from the context image and the catalog you define.
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.00 $ '[11:-3]
__version_date__ = '$Date: 2012/12/05 $ '[7:-3]
__author__       = 'Xingxing Huang'

import pdb,os,sys,glob,popen2
from optparse import OptionParser
import pyfits,numpy
import os, glob ,string

class Catalog:

    """ usage:

    """

    def __init__(self,catdir='./',fitsdir='./',catalogfile=None,contextfile=None,flagfile=None,contextFitsNum=2,contextHighLevel=1):
        self.catalog=catalogfile
        self.context=contextfile
        self.flag=flagfile
        self.lowvalue=contextFitsNum    #  only the pixels added by more than "lowlevel" fits files are choosed.
        self.highlevel=contextHighLevel #  only the pixels with the highest "highlevel" context value are choose.
        self.highvalue=99 
        self.catdir = catdir
        self.contextdir=fitsdir
        self.flagdir=fitsdir
        if self.catalog == None:
           self.catalog = glob.glob(os.path.join(self.catdir, '*.cat'))
        if self.context == None:
           self.context = glob.glob(os.path.join(self.contextdir, '*_context.fits'))
        if self.flag == None:
           self.flag = glob.glob(os.path.join(self.flagdir, '*_FLAG.fits'))
        self.choose = 0
        #self.highvalue = []
        
    def run(self):
    	for keys in self.catalog:
    	    fname=os.path.join(self.catdir,keys)
    	    f = open(fname,'r')
    	    print ">>> read catalog "+fname
    	    num = 1
            outall = os.path.join(self.catdir,keys+'.all')              # output catalog is put in the same directory
            outchoose = os.path.join(self.catdir,keys+'.choosed')              # output catalog is put in the same directory
            if os.path.isfile(outall):
               back = string.replace(outall,'all','allold')
               os.rename(outall,back)
               print 'old version of '+keys+'.all'+' exists! We back it up!'
            if os.path.isfile(outchoose):
               back = string.replace(outchoose,'choosed','choosedold')
               os.rename(outchoose,back)
               print 'old version of '+keys+'.choosed'+' exists! We back it up!'
            for i in f.readlines():
               #i = string.split(i,'\n')[0]
               x = float(string.split(i,' ')[0])          # Notice the catalog should be only two value in each line and the value is separated by only one blank. 
               y = float(string.split(i,' ')[1])
               print x, y
               self.choose = 1
    	       contextvalue=self.choosecontext(x,y)
    	       flagvalue=self.chooseflag(x,y)
    	              
	       logf = open(outall,'a')
               format = str('%5s %5s %10s %10s %4s' )
               print >> logf,format % (x, y, str(contextvalue),str(flagvalue),self.choose)
               logf.close() 
               
               if self.choose == 1: 
                   logf = open(outchoose,'a')
                   format = str('%5s %5s %10s %10s')
                   print >> logf,format % (x, y, str(flagvalue),str(contextvalue))
                   logf.close()   
               print '   line %3s'% num,' Context ',contextvalue,' Flag ',flagvalue,' Choosed ',self.choose
               num = num+1
    	    print "writing choosed catalog from "+keys+' complete'    
    	    f.close()   

 #   def ContextImageNum(self):
 #       
 #       for i in self.context:
 ##           # open each context file
  #          PathName=os.path.join(self.contextdir, i)
 ##           fits=pyfits.open(PathName)
 #           # define how many fits files are added
 #           l=[]
 #           for array in fits[0].data:
 #              l.extend(array)
 #           self.highvalue=DefineContextValue(max(l)).fitsnum()-self.highlevel+1
                       
    def choosecontext(self,x,y):
        ''' find out the context value from *context.fits   and determine wether the pixel is choosed or  not.'''
        contextvalue=[]
        for i in self.context:
            # open each context file
            PathName=os.path.join(self.contextdir, i)
            fits=pyfits.open(PathName)
            # define how many fits files are added
            ###l=[]
            ###for array in fits[0].data:
            ###   l.extend(array)
            ###maxvalue=max(l)
            ndrizim=fits[0].header['NDRIZIM']
            self.highvalue=ndrizim-self.highlevel+1
            # define how many fits files are added for this pixel 
            n=DefineContextValue(fits[0].data[x,y]).fitsused()
            contextvalue.append(n)
            # find out whether the pixel should be exclued 
            if (n < self.lowvalue) and (n < self.highvalue):
               self.choose = 0
               
            fits.close()
        return contextvalue

    def chooseflag(self,x,y):
        ''' find out the flag value from the *FLAG.fits '''
        flagvalue=[]
        for i in self.flag:
            PathName=os.path.join(self.flagdir, i)
            fits=pyfits.open(PathName)
            flagvalue.append(fits[0].data[x,y])
            # find out whether the pixel should be exclued 
            for i in flagvalue:
              if i == 1:
    	        self.choose = 0 
            fits.close()
        return flagvalue




class DefineContextValue:
    """Just get the value using this function"""
    def __init__(self,contextvalue):
        bin = lambda n : (n > 0) and (bin(n/2) + str(n%2)) or ''
        n=bin(contextvalue)
        self.num1=str(n).count('1')
	self.num0=str(n).count('0')
	
    def fitsnum(self):
        '''return the total number of fits files used for creating **drz.fits'''
        return self.num0 + self.num1
        
    def fitsused(self):
        '''return the number of fits files used for the pixel you measure'''
        return self.num1
 
	
        
		
		
if __name__ == '__main__':
	a=Catalog(catdir='./',fitsdir='./',contextFitsNum=2,contextHighLevel=1)
        a.run()


'''
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
       
'''
