#!/usr/bin/env python

# $Id: output.py,v 1.0 2011/09/21 19:46:33 zheng Exp $
# ---------------------------------------------------------------------
# 

__version__      = '$Revision: 1.0 $ '[11:-3]
__version_date__ = '$Date: 2011/09/21 19:46:33 $ '[7:-3]
__author__       = 'W Zheng <zheng@pha.jhu.edu>'

import os, string, glob, shutil,sys
import pdb, popen2
import path
import numpy
import fUtil,xmlUtil
import pyfits
from   pUtil import logFile, ptime
# from   msg   import pMessage
from   sys   import version
from   shutil import copyfile,copy
from   xydrizzle import wcsutil,drutil #WZ
pyversion = version

class copy_drz:

    """ 
    Copy images over from acsex, wfex and uvex
    """

    def __init__(self,obs): #WZ Oct 2011
        
        os.chdir(obs.fitsdir)
        i=0
        for im in obs.sciImageList:
            substr=string.split(im,'_')
            j=len(substr)
            target=substr[j-3]
            filter=substr[j-2]
            # inpfits = pyfits.open(im)
            #detector = inpfits[0].header.get('DETECTOR')
            detector = drutil.getPrimaryKeyword(im+'[0]','DETECTOR') 
            if (detector == 'WFC'):
                indir=os.path.join('/data01/acex_runs/datasets',target,filter,'Images/')
            elif (detector == 'IR'):
                indir=os.path.join('/data01/wfex_runs/datasets',target,filter,'Images/')
            elif (detector == 'UVIS'):
                indir=os.path.join('/data02/uvex_runs/datasets',target,filter,'Images/')
            else:
                obs.logfile.write("Unlisted detector name: "+detector)
            f0 = indir +string.replace(im,'sci','drz')
            shutil.copy2(f0,im)
            f1 = indir +string.replace(obs.weightImageList[i],'sci','drz')
            shutil.copy2(f1,obs.weightImageList[i])
            f1 = indir +string.replace(obs.contextImageList[i],'sci','drz')
            shutil.copy2(f1,obs.contextImageList[i])
            f1 = indir + obs.rmsImageList[i]
            shutil.copy2(f1,obs.rmsImageList[i])
            f1 = indir +obs.flagImageList[i] 
            shutil.copy2(f1,obs.flagImageList[i])
            i=i+1
        obs.logfile.write("Drz files copied over.")            
        return


class update_tweak:

    """ 
    Read tweak.res in each working directory. Make a general "tweak.res" file
    """

    def __init__(self,obs,bandlist): #WZ Oct 2011 by xingxing
        
        os.chdir(obs.fitsdir)
        target=obs.newobs
        bname=[]
        tlist=[]
        #f1 = os.path.join('/data01/acex_runs/datasets',target,'tweak.res')
        #if os.path.isfile(f1):
        #    f = open(f1,'r')
        #    for line in f.readlines():
        #        filt = string.split(line)[0]
        #        bname.append(filt)
        #        tlist.append(line)
        #    f.close()
        #f2 = os.path.join('/data01/wfex_runs/datasets',target,'tweak.res')
        #if os.path.isfile(f2):
        #    f = open(f2,'r')
        #    for line in f.readlines():
        #        filt = string.split(line)[0]
        #        bname.append(filt)
        #        tlist.append(line)
        #    f.close()
        #f3 = os.path.join('/data01/uvex_runs/datasets',target,'tweak.res')
        #if os.path.isfile(f3):
        #    f = open(f3,'r')
        #    for line in f.readlines():
        #        filt = string.split(line)[0]
        #        bname.append(filt)
        #        tlist.append(line)
        #    f.close()
	f1 = os.path.join(os.environ['DATASETS'],target,'run_'+target+'.log')  #change later  by xingxing
#	f1 = os.path.join(os.environ['DATASETS'],target,'tweak.res')
        if os.path.isfile(f1):
            f2 = open(f1,'r')
            for line in f2.readlines():
                filt = string.split(line)[0]
                bname.append(filt)
                tlist.append(line)
            f2.close()
	else :
	    print (os.path.join(os.environ['DATASETS'],target,'run_'+target+'.log')+' does not exit !!!')
	    ##sys.exit()
	    return
        twf=open(os.path.join(obs.aligndir,'run_'+target+'.log'),'w')
        headline='     Name  \t    dX       dY    Rot (rad) ResX ResY Npt   Date'
        print >> twf, headline
        for im in obs.sciImageList:
            substr=string.split(im,'_')
            j=len(substr)
            # target=substr[j-3]
            filter=substr[j-2]
            # detector = drutil.getPrimaryKeyword(im+'[0]','DETECTOR') 
            i=0
            for line in bname:
                if string.find(line,filter)>-1:
                    index = i
                i = i + 1
            best = string.split(tlist[index])
            name=best[0]
            dx =float(best[1])
            dy =float(best[2])
            dr =float(best[3])
            rex=float(best[4])
            rey=float(best[5])
            npt=int(best[6])
            date=best[7]
            format=str('%s\t%8.2f %8.2f %10.7f %4.2f %4.2f %3d %s')
            print format % (name,dx,dy,dr,rex,rey,npt,date)
            print >> twf, format % (name,dx,dy,dr,rex,rey,npt,date)
            # pdb.set_trace()
        twf.close()
        #del f1,f2,twf,line,tlist,bname,name,dx,dy,dr,rex,rey,npt,date
        obs.logfile.write("Tweak parameters updated.")            
        return

class output:

    """ 
    Copy image and catalogs over to the Output depositiry
    """

    def __init__(self,obs,bandlist): #WZ Nov 2011 by xingxing
        
        self.outputdir = obs.prodir
        self.catdir = obs.catdir
        self.sciImageList = []
        self.rmsImageList = []
        self.whtImageList = []
        self.filterList = []
        date = string.split(ptime(),'T')[0]
        date = date.replace('-','')
        stamp= date[2:]
        pscale = numpy.rint(obs.dr_asecpix * 1000)/1000.
        ps = string.split(str(pscale),'.')[1]
        os.chdir(obs.fitsdir)
        i=0
        for im in obs.sciImageList:
            substr=string.split(im,'_')
            j=len(substr)
            target=substr[j-3]
            filter=substr[j-2]
            # inpfits = pyfits.open(im)
            #detector = inpfits[0].header.get('DETECTOR')
            detector = drutil.getPrimaryKeyword(im+'[0]','DETECTOR') 
            if (detector == 'WFC'):
                f = target+'_acs_'+filter+'_sci_'+ps+'mas_'+stamp+'.fits'
            elif (detector == 'IR'):
                f = target+'_ir_'+filter+'_sci_'+ps+'mas_'+stamp+'.fits'
            elif (detector == 'UVIS'):
                f = target+'_uvis_'+filter+'_sci_'+ps+'mas_'+stamp+'.fits'
            else:
                obs.logfile.write("Unlisted detector name: "+detector)
            outfile = os.path.join(self.outputdir,f)
            shutil.copy(im,outfile)

            im=obs.weightImageList[i]
            whtfile=string.replace(outfile,'_sci_','_wht_')
            shutil.copy(im,whtfile)

            im=obs.rmsImageList[i]
            rmsfile=string.replace(outfile,'_sci_','_rms_')
            shutil.copy(im,rmsfile)

            self.sciImageList.append(outfile)
            self.whtImageList.append(whtfile)
            self.rmsImageList.append(rmsfile)
            self.filterList.append(filter)
            i = i + 1

        for band in bandlist:
        # for band in ['nir']:
            im='detection_'+band+'.fits'
            detfile = os.path.join(self.outputdir,target+'_detection_'+band+'.fits')
            shutil.copy(im,detfile)
        obs.logfile.write("Image copying done. Start gzipping")
        os.chdir(obs.prodir)
        cmd = 'gzip *fits'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()

        for band in bandlist:                                       
            dname='detection_'+band+'.cat'
            shutil.copy(os.path.join(self.catdir, dname),target+'_'+dname)
            mname='multicolor_'+band+'.cat'
            shutil.copy(os.path.join(self.catdir, mname),target+'_'+mname)
            nname='multicolor_'+band+'.columns'
            shutil.copy(os.path.join(self.catdir, nname),target+'_'+nname)
            bname='bpz_'+band+'.cat'
            shutil.copy(os.path.join(self.catdir, bname),target+'_'+bname)
            for filt in self.filterList:
                # fname=target+'_'+filt+'_drz_'+band+'.cat'  #by xingxing
                fname=target+'_'+filt+'_'+band+'.cat'  # WZ Mar 2013
                # pdb.set_trace()
                shutil.copy(os.path.join(self.catdir, fname),fname)
        obs.logfile.write("Catalog copying done")

        cmd = 'gzip *fits'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        cmd = 'chmod 755 *'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        obs.logfile.write("Fits files compressed. Open permission set")
        del bname,dname,fname,mname,nname

        # pdb.set_trace()
        os.chdir(obs.prodir)
        os.chdir(obs.fitsdir)
        return
