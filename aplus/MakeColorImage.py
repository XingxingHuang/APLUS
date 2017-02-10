#!/usr/bin/env python
# $Id: MakeColorImage.py,v 1.00 2012/11/30 20:15:45 Huang  $
# ---------------------------------------------------------------------
# make color images using trilogy.py

__version__      = '$Revision: 1.00 $ '[11:-3]
__version_date__ = '$Date: 2012/11/30 $ '[7:-3]
__author__       = 'Xingxing Huang'

from trilogy import * #Trilogy
import os,sys,string,path
import pdb

class MakeColorImage:
    """
    make color images using trilogy.py
    
    """
    def __init__(self,observation,inputdir=None,outputdir=None):
    	self.inputdir 	 = inputdir
    	self.colordir 	 = outputdir
    	self.obsName     = observation
    	self.filterlist  = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w', 'f225w', 'f275w', 'f336w', 'f390w', 'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f850lp', 'f814w']
    	self.filter      = [i for i in self.filterlist]
	self.file        = ['source.in','acs.in','ir.in','uvis.in','f140.in','f160.in']
	self.irfilter    = ['f105w', 'f110w', 'f125w', 'f140w', 'f160w']
	self.acsfilter   = ['f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f850lp','f814w']
	self.uvisfilter  = ['f225w', 'f275w', 'f336w', 'f390w']
    
    def run(self):
        if self.inputdir == None:
            self.inputdir = os.path.join(path.Env().getenv('DATASETS'), self.obsName,'final','Images')
        if self.colordir == None:
            self.colordir = os.path.join(path.Env().getenv('DATASETS'), self.obsName,'final','CoImages')
        if not os.path.isdir(self.colordir):
            os.mkdir(self.colordir)
        print "Color images will be stored in "+self.colordir


        for keys in self.filterlist:
            drzfile = self.obsName+"_"+keys+"_drz.fits"
            if not os.path.isfile(os.path.join(self.inputdir,drzfile)):
                print drzfile+" not found !!"
                self.filter.remove(keys)
                if keys in self.acsfilter:
                    self.acsfilter.remove(keys)
                elif keys in self.irfilter:
                    self.irfilter.remove(keys)
                elif keys in self.uvisfilter:
                    self.uvisfilter.remove(keys)
                else:
                    print keys+' is not found in ACS IR UVIS.'
                    sys.exit()
        print "filters used:"
        print "UVIS: ",self.uvisfilter            
        print "ACS:  ",self.acsfilter            
        print "IR:   ",self.irfilter  
                  
        #pdb.set_trace()
        #source.in
        if self.acsfilter and self. irfilter and self.uvisfilter: 
            fname= os.path.join(self.colordir,'source.in')
            print "\n>>>>>>>>>make ",fname,'...'
            if os.path.isfile(fname):
                back=fname+'_old'
                os.rename(fname,back)
            f = open(fname,'w')
            print >> f,'B'
            for keys in self.uvisfilter:
                print >>f, self.obsName+'_'+keys+'_drz.fits'
            print >>f,'\n' 
            print >> f,'G'
            for keys in self.acsfilter:
                print >>f, self.obsName+'_'+keys+'_drz.fits'
            print >>f,'\n' 
            print >> f,'R'
            for keys in self.irfilter:
                print >>f, self.obsName+'_'+keys+'_drz.fits'
            print >>f,'\n' 
            
            print >> f,'indir    ',self.inputdir
            print >> f,'outname  ',self.obsName+'_source'
            print >> f,'samplesize  1000'
            print >> f,'showstamps 0'
            print >> f,'satpercent  0.001'
            print >> f,'noiselum    0.12'
            print >> f,'colorsatfac  1'
            print >> f,'deletetests  0'
            print >> f,'sampledx   330'
            f.close()
            try:
                Trilogy(infile=fname,outdir=self.colordir,**params_cl()).run()
            except Exception,err:
                print "*******ERROR: fail for "+fname
           
	#acs.in
        if len(self.acsfilter) >= 3:
            fname= os.path.join(self.colordir,'acs.in')
            print "\n>>>>>>>>>make ",fname,'...'
            if os.path.isfile(fname):
                back=fname+'_old'
                os.rename(fname,back)
            f = open(fname,'w')
            
            if len(self.acsfilter) == 3:
                print >> f,'B'
                for keys in self.acsfilter[0:1]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.acsfilter[1:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.acsfilter[-1]+'_drz.fits'
                print >>f,'\n' 
            elif len(self.acsfilter) <= 5:
                print >> f,'B'
                for keys in self.acsfilter[0:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.acsfilter[2:-1]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.acsfilter[-1]+'_drz.fits'
                print >>f,'\n' 
            else:
                print >> f,'B'
                for keys in self.acsfilter[0:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.acsfilter[2:-2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.acsfilter[-2]+'_drz.fits'
                print >>f, self.obsName+'_'+self.acsfilter[-1]+'_drz.fits'
                print >>f,'\n' 
                
            print >> f,'indir    ',self.inputdir
            print >> f,'outname  ',self.obsName+'_acs'
            print >> f,'samplesize  1000'
            print >> f,'showstamps 0'
            print >> f,'satpercent  0.001'
            print >> f,'noiselum    0.12'
            print >> f,'colorsatfac  1'
            print >> f,'deletetests  0'
            print >> f,'sampledx   330'
            f.close()
            try:
                Trilogy(infile=fname,outdir=self.colordir,**params_cl()).run()
            except Exception,err:
                print "*******ERROR: fail for "+fname
            
	#uvis.in
        if len(self.uvisfilter) >= 3:
            fname= os.path.join(self.colordir,'uvis.in')
            print "\n>>>>>>>>>make ",fname,'...'
            if os.path.isfile(fname):
                back=fname+'_old'
                os.rename(fname,back)
            f = open(fname,'w')
            
            if len(self.uvisfilter) == 3:
                print >> f,'B'
                for keys in self.uvisfilter[0:1]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.uvisfilter[1:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.uvisfilter[-1]+'_drz.fits'
                print >>f,'\n' 
            elif len(self.uvisfilter) <= 5:
                print >> f,'B'
                for keys in self.uvisfilter[0:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.uvisfilter[2:-1]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.uvisfilter[-1]+'_drz.fits'
                print >>f,'\n' 
            else:
                print >> f,'B'
                for keys in self.uvisfilter[0:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.uvisfilter[2:-2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.uvisfilter[-2]+'_drz.fits'
                print >>f, self.obsName+'_'+self.uvisfilter[-1]+'_drz.fits'
                print >>f,'\n' 
                
            print >> f,'indir    ',self.inputdir
            print >> f,'outname  ',self.obsName+'_uvis'
            print >> f,'samplesize  1000'
            print >> f,'showstamps 0'
            print >> f,'satpercent  0.001'
            print >> f,'noiselum    0.12'
            print >> f,'colorsatfac  1'
            print >> f,'deletetests  0'
            print >> f,'sampledx   330'
            f.close()
            try:
                Trilogy(infile=fname,outdir=self.colordir,**params_cl()).run()
            except Exception,err:
                print "*******ERROR: fail for "+fname
      
      	#ir.in
        if len(self.irfilter) >= 3:
            fname= os.path.join(self.colordir,'ir.in')
            print "\n>>>>>>>>>make ",fname,'...'
            if os.path.isfile(fname):
                back=fname+'_old'
                os.rename(fname,back)
            f = open(fname,'w')
            
            if len(self.irfilter) == 3:
                print >> f,'B'
                for keys in self.irfilter[0:1]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.irfilter[1:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.irfilter[-1]+'_drz.fits'
                print >>f,'\n' 
            elif len(self.irfilter) <= 5:
                print >> f,'B'
                for keys in self.irfilter[0:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.irfilter[2:-1]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.irfilter[-1]+'_drz.fits'
                print >>f,'\n' 
            else:
                print >> f,'B'
                for keys in self.irfilter[0:2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'G'
                for keys in self.irfilter[2:-2]:
                    print >>f, self.obsName+'_'+keys+'_drz.fits'
                print >>f,'\n' 
                print >> f,'R'
                print >>f, self.obsName+'_'+self.irfilter[-2]+'_drz.fits'
                print >>f, self.obsName+'_'+self.irfilter[-1]+'_drz.fits'
                print >>f,'\n' 
                
            print >> f,'indir    ',self.inputdir
            print >> f,'outname  ',self.obsName+'_ir'
            print >> f,'samplesize  1000'
            print >> f,'showstamps 0'
            print >> f,'satpercent  0.001'
            print >> f,'noiselum    0.12'
            print >> f,'colorsatfac  1'
            print >> f,'deletetests  0'
            print >> f,'sampledx   330'
            f.close()      
            try:
                Trilogy(infile=fname,outdir=self.colordir,**params_cl()).run()
            except Exception,err:
                print "*******ERROR: fail for "+fname
                
        #f140.in    
        if "f140w" in self.irfilter and 'f105w' in self.irfilter and self.acsfilter:
            fname= os.path.join(self.colordir,'f140.in')
            print "\n>>>>>>>>>make ",fname,'...'
            if os.path.isfile(fname):
                back=fname+'_old'
                os.rename(fname,back)
            f = open(fname,'w')
            print >> f,'B'
            print >>f, self.obsName+'_'+self.acsfilter[-1]+'_drz.fits'
            print >>f,'\n' 
            print >> f,'G'
            print >>f, self.obsName+'_f105w_drz.fits'
            print >>f,'\n' 
            print >> f,'R'
            print >>f, self.obsName+'_f140w_drz.fits'
            print >>f,'\n' 
            
            print >> f,'indir    ',self.inputdir
            print >> f,'outname  ',self.obsName+'_f140'
            print >> f,'samplesize  1000'
            print >> f,'showstamps 0'
            print >> f,'satpercent  0.001'
            print >> f,'noiselum    0.12'
            print >> f,'colorsatfac  1'
            print >> f,'deletetests  0'
            print >> f,'sampledx   330'
            f.close()
            try:
                Trilogy(infile=fname,outdir=self.colordir,**params_cl()).run()
            except Exception,err:
                print "*******ERROR: fail for "+fname
                
        #f160.in    
        if "f160w" in self.irfilter and 'f105w' in self.irfilter and self.acsfilter:
            fname= os.path.join(self.colordir,'160.in')
            print "\n>>>>>>>>>make ",fname,'...'
            if os.path.isfile(fname):
                back=fname+'_old'
                os.rename(fname,back)
            f = open(fname,'w')
            print >> f,'B'
            print >>f, self.obsName+'_'+self.acsfilter[-1]+'_drz.fits'
            print >>f,'\n' 
            print >> f,'G'
            print >>f, self.obsName+'_f105w_drz.fits'
            print >>f,'\n' 
            print >> f,'R'
            print >>f, self.obsName+'_f160w_drz.fits'
            print >>f,'\n' 
            
            print >> f,'indir    ',self.inputdir
            print >> f,'outname  ',self.obsName+'_f160'
            print >> f,'samplesize  1000'
            print >> f,'showstamps 0'
            print >> f,'satpercent  0.001'
            print >> f,'noiselum    0.12'
            print >> f,'colorsatfac  1'
            print >> f,'deletetests  0'
            print >> f,'sampledx   330'
            f.close()
            try:
                Trilogy(infile=fname,outdir=self.colordir,**params_cl()).run()
            except Exception,err:
                print "*******ERROR: fail for "+fname
            
        print "Done for color images!\n "

# Just make a test with the following:
if __name__ == '__main__':

   MakeColorImage('a209',inputdir='/home/zheng/aaplusa_runs/datasets/a209/final/Images',outputdir='/home/zheng/aaplusa_runs/datasets/a209/final/Images').run()  



  
