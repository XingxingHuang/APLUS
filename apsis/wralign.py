#!/usr/bin/env python
import pdb
# $Id: wralign.py,v 0.1 2003/11/11 21:51:02 anderson Exp $
# ---------------------------------------------------------------------
# Wrapper module for the new align module which conforms to the old align interface

from align import *
from sys import version
from msg import pMessage
import os,pyfits,popen2
import subprocess
pyversion = version

__version__      = '$Revision: 1.34 $ '[11:-3]

class FrameSetUpClass:
    def __init__(self):
        self.outshift = None
        self.outsize = None
        self.outscale = None

    def setoutsize(self,x,y):
        self.outsize = (x,y)

    def setoutshift(self,x,y):
        self.outshift = (x,y)
      
    def setoutscale(self,outscale):
        self.outscale = outscale

class alignImage:
    """ Intended to have the same API as the old apsis AlignImage object """

    def __init__(self,obs,extRefIm=0,idcTab=None,Edgebuf=12,irafSkyToo=1,nebula=0,
                 grism=0,starclus=0,keepstep=0,notrot=None,useMinSky=0,FrameSetUp=None, 
                 sxtrthresh=None):
        self.obsAsnDict = obs.asnDict
        self.irafSkyToo = irafSkyToo
        self.deAmp      = (not keepstep)      
        self.notrot     = notrot
        self.useMinSky  = useMinSky
        self.messagedir = obs.messagedir                 # where the module message will go
        self.root       = obs.newobspath                 # root path of the observation dir
        self.obs = obs
        self.obsName    = obs.newobs
        self.allShifts  = obs.allShifts  = []            # List of shiftfile names (full path) needed by makeAsn.
        self.FrameSetUp = FrameSetUp
        # pdb.set_trace()  #class GlobalBlock Seen in align.py LINE 87
        self.GlobalBlock = GlobalBlock(obs.logfile,obs.newalign,obs.newpar,obs.newfits,obs.root,obs.pardir)
        self.logfile = self.GlobalBlock.logfile
        self.modName = self.GlobalBlock.modName
        if self.notrot:
            self.GlobalBlock.logfile.write("Align will enforce median shifts; no rotation.")
        CameraList = [WFC3UVISCamera(obs.pardir,self.modName),WFC3IRCamera(obs.pardir,self.modName),ACSWFCCamera(obs.pardir,self.modName),ACSHRCCamera(obs.pardir,self.modName)] #WZ   
        # CameraList = [ACSWFCCamera(obs.pardir,self.modName),ACSHRCCamera(obs.pardir,self.modName)]

        if self.obsAsnDict:
            self.imfiles = []
            for key in self.obsAsnDict.keys(): 
                for im in self.obsAsnDict[key]:   #class FitsFile Seen in align.py LINE 239 LINE 313
                    self.imfiles.append(FitsFile(im,key,self.GlobalBlock,useMinSky=useMinSky))  
        else:
            self.imfiles = self.obsFitsList
        
        CurExtImage = None
        if extRefIm:
            TSearchParam = ExtDefSearchParam(obs.pardir,self.modName,lowLim=0.0)
            CurExtImage = FitsFile(extRefIm,"",self.GlobalBlock,useMinSky=useMinSky,SearchParam=TSearchParam)
            CurExtImage.ChangeToRef()
            CurExtImage.ChangeToExtRef()

        self.TAlign  = Align(self.imfiles,self.allShifts,CameraList,self.GlobalBlock,SkyFile=0,idcTab=idcTab,CurExtImage=CurExtImage)
	# pdb.set_trace()
        self.TAlign.setupImages()      # class Align  seen in align.py LINE 1614
        self.outputList = self.GlobalBlock.outputList
        obs.detector = self.imfiles[0].Camera.Name
        self.outscale = self.imfiles[0].outscale
	print "Please check the [obs.detector] here  by xingxing"        

    def makeMatchCats(self,SkyFile=0):    
        self.TAlign.makeMatchCats(SkyFile,self.irafSkyToo)

    def match(self,intRef='',superalign=None,mdz_align=None,retry=0):
        TAlignParam = StandAlignParam(superalign,mdz_align)      # XXH seen in align.py LINE 1601
        if (mdz_align):
            ref=self.obs.asnDict[self.obs.asnDict.keys()[0]][0] #WZ the first image as ref          
        else:
            ref=None #WZ Nov
        self.TAlign.match(TAlignParam,self.notrot,ref,intRef,retry) # WZ
        self.obs.refcdmat = self.TAlign.refImage.refcdmat      # dictionary with final cd matrix
        self.obs.refotherkeys = self.TAlign.refImage.refotherkeys # dictionary with other keys relating to the wcs solution
        self.obs.reflogfile_add = self.TAlign.refImage.reflogfile_add

        # set shift,size frames to fit those of the external image
        if self.TAlign.refImage.extref == 1:
            shiftx = (self.TAlign.refImage.x0[0] - self.TAlign.refImage.dimx/2)*self.TAlign.refImage.outscale
            shifty = (self.TAlign.refImage.y0[0] - self.TAlign.refImage.dimy/2)*self.TAlign.refImage.outscale
            sizex = self.TAlign.refImage.dimx*self.TAlign.refImage.outscale
            sizey = self.TAlign.refImage.dimy*self.TAlign.refImage.outscale
            self.FrameSetUp.setoutshift(shiftx,shifty)
            self.FrameSetUp.setoutsize(sizex,sizey)
            self.FrameSetUp.setoutscale(self.TAlign.refImage.outscale)

    def subtractSkies(self, ForceSub=0, omitSimple=0, omitMEF=0, aveExt=0):
        self.TAlign.subtractSkies(self.useMinSky,self.deAmp,ForceSub,omitSimple,omitMEF,aveExt)

    def removeUnpicklable(self):
        self.MatchDict        = {}
        for im in self.imfiles:
            if not im.extref:          
                self.MatchDict[im.PathName] = {}
                for key in im.MatchDict.keys():
                    self.MatchDict[im.PathName][key] = im.MatchDict[key]
      
        self.TAlign.removeUnpicklable()
      
    def mkMsg(self):
        """create and write module level message for this class.
        Most of this is just compiling the info. meta is a dictionary
        of lists where each list is a list of tuples describing the
        tag lines for the particular section of the message.  This tuple 
        format conforms to that used by the xmlMessage class which is
        modeled on basic python argument passing, i.e. (key,*value,**attr).
        """
        self.meta = {}
        self.meta['module']    = []
        self.meta['meta']      = []
        self.meta['input']     = []
        self.meta['output']    = []
        self.meta['errorlist'] = []

        self.meta['module'].append(('module','name='+self.modName,'version='+__version__,'dataset='+self.obsName))
        self.meta['module'].append(('root',self.root))
        self.meta['meta'].append(('meta',))
        self.meta['meta'].append(('depend',))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','python'))
        self.meta['meta'].append(('version',pyversion.split()[0]))

        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','pyfits'))
        self.meta['meta'].append(('version',pyfits.__version__.split()[0]))


        # SExtractor info
        sub  = subprocess.Popen(['sex', '--version'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
        outp = sub.stdout.readlines()
        name = outp[0].split()[0]
        ver  = outp[0].split()[2]
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name',name))
        self.meta['meta'].append(('version',ver))
        cmdline1 = 'sex fitsfile -c self.defaultInParFile'
        self.meta['meta'].append(('commandline',cmdline1))
        del outp,sub,name,ver

        # match info
        sub  = subprocess.Popen(['match', '--version'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
        outp = sub.stdout.readlines()
        name = outp[0].split()[0][:-1]
        ver  = outp[0].split()[2]
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name',name))
        self.meta['meta'].append(('version',ver))
        # ok, here's the new, post version 0.5 command, where the part in brackets
        # [identity xsh=dx ysh=dy] is only included if the angle is expected to be small
        cmdline2 = "match file_matchin 1 2 3 self.refMatchFile 1 2 3 max_iter=10 scale=1.00 nobj=60 matchrad=str(tolmatch) trirad=0.00033 recalc medtf medsigclip=2.7 outfile=outfile [identity xsh=+str(hdr_dx)+ ysh=+str(hdr_dy)]"

        self.meta['meta'].append(('commandline',cmdline2))
        del outp,sub,name,ver

        if self.GlobalBlock.errorList:
            self.meta['errorlist'].append(('errorlist',))
            for pkg,err in self.GlobalBlock.errorList:
                self.meta['errorlist'].append(('erroritem',err,'frompkg='+pkg))

        # input section
        self.meta['input'].append(('input',))
        # 25/Jan/02 - updating meta['input'] to reflect changes in MatchDict
        for im in self.imfiles:
            if not im.extref:
                self.meta['input'].append(('file','type=image/x-fits'))
                self.meta['input'].append(('name',os.path.join("Images",im.PathName)))
                self.meta['input'].append(('attr',im.MatchDict['Nraw']          ,'name=Nraw'))
                self.meta['input'].append(('attr',im.MatchDict['Ngood']         ,'name=Ngood'))
                self.meta['input'].append(('attr',im.MatchDict['Nmatch']        ,'name=Nmatch'))
                self.meta['input'].append(('attr',im.MatchDict['xpix_shift']    ,'name=xpix_shift'))
                self.meta['input'].append(('attr',im.MatchDict['xpix_shift_err'],'name=xpix_shift_err'))
                self.meta['input'].append(('attr',im.MatchDict['ypix_shift']    ,'name=ypix_shift'))
                self.meta['input'].append(('attr',im.MatchDict['ypix_shift_err'],'name=ypix_shift_err'))
                self.meta['input'].append(('attr',im.MatchDict['xarc_shift']    ,'name=xarc_shift'))
                self.meta['input'].append(('attr',im.MatchDict['yarc_shift']    ,'name=yarc_shift'))
                self.meta['input'].append(('attr',im.MatchDict['angle']         ,'name=angle'))
                self.meta['input'].append(('attr',im.MatchDict['angle_err']     ,'name=angle_err'))
                self.meta['input'].append(('attr',im.MatchDict['skyval']        ,'name=skyval'))
                self.meta['input'].append(('attr',im.MatchDict['skyrms']        ,'name=skyrms'))

        # output section
        if self.GlobalBlock.outputList:
            self.meta['output'].append(('output',))
        for f in self.GlobalBlock.outputList.keys():
            if string.find(f,".fits") != -1:
                self.meta['output'].append(('file','type=image/x-fits'))
                self.meta['output'].append(('name',os.path.join("Images",f)))
                for pred in  self.GlobalBlock.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
            else: 
                self.meta['output'].append(('file','type=text/ascii'))
                self.meta['output'].append(('name',os.path.join("align",f)))
                for pred in  self.GlobalBlock.outputList[f]:
                    self.meta['output'].append(('predecessor',os.path.join("Images",pred)))

        # pass this dictionary to the class pMessage...
        msgFile = os.path.join(self.messagedir,self.modName+"_module.xml")
        mmsg = pMessage(self.meta)
        mmsg.writeMsg(msgFile)
        return
