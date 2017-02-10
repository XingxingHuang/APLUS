"""
     Modified W Zheng: Using header WCS if simplematch fails
     THis is an advanced version as it contains files for both ACS and WFC3
"""
# ------  new api
# 
# -- list of fits files
# -- files where the relative alignment should be good and thus can be stacked...
# -- associated with each file is parameters to find stars
# -- some outside reference image??

import pyfits,string,os,fUtil,matutil,popen2,time,wcsclass
import amputil,math,pUtil
import pdb,string,shutil #WZ
from   numpy import where,resize,less_equal,arange,zeros,abs
from   pyraf import iraf          # this is here for alternate sky values
from   xydrizzle import wcsutil,drutil #WZ
import multidrizzle

_angle_tolerance_ = 0.004     #(deg). Change from .005 on 21/Jun/02; this is
                                # now 0.15 pix at WFC edge; 0.2 pix in corner
                                # (half these for HRC); likely within error of zero.

def ReturndefShiftDict(GlobalBlock,refImage):
  defShiftDict  = {}
  defshifts = os.path.join(GlobalBlock.obsAlign,'default.shifts')
  if (not os.path.isfile(defshifts)):
    del defshifts
  else:
    GlobalBlock.logfile.write('Reading shifts from '+defshifts)
    _tmp_lines = open(defshifts).readlines()
    _tmp_refIm = ''
    for line in _tmp_lines:
      if (line[0] == '#') or (len(line.split()) < 1):
        continue
      elif ("=" in line and line.split("=")[0].split()[0].lower()[0:3] == "ref"):
        _tmp_refIm = line.split("=")[1].split()[0]
      else:
        defShiftDict[line.split()[0]] = line
    if _tmp_refIm != refImage:
      errtxt = "Error:  selected ref image differs from one specified in default shifts file."+\
               "\n\t\t\t selected: "+refImage.PathName+"   defshifts specified: "+_tmp_refIm+\
               "\n\t\t\t Use intRef flag to reselect."
      GlobalBlock.logfile.write(errtxt)
      GlobalBlock.logfile.write("-------------> Quiting match() <--------------#")
      raise Exception,errtxt
    del _tmp_lines,_tmp_refIm,defshifts
  return defShiftDict

class skyLinesClass:
    def __init__(self,SkyFile,GlobalBlock):
      self.skyLines = None
      self.GlobalBlock = GlobalBlock
      if SkyFile:
        self.skyLines = []
        GlobalBlock.logfile.write("makeMatchCats: got SkyFile flag; will look for a default.skies")
	defskies=os.path.join(self.obsFits,'default.skies')
	if (not os.path.isfile(defskies)):
          GlobalBlock.logfile.write("File not found: default.skies")
	  raise Exception, "File not found: default.skies"
	else:
	  GlobalBlock.logfile.write('Reading Sky and RMS values from '+defskies)
        _tmp_lines = open(defskies).readlines()
        for line in _tmp_lines:
          if (line[0] == '#' or (len(line.split()) < 1)):
            continue
          elif (len(line.split()) < 3):
            errtxt = "Incorrect format: default.skies must specify: Simplefits Sky RMS."
            GlobalBlock.logfile.write(errtxt)
            raise Exception, errtxt
          else:
            self.skyLines.append(line)
        del _tmp_lines,defskies
              
    def checkforfile(self,shortname):
       # Check to see if there's a sky value for the image in skyLines,
       # which is the contents of 'default.skies' if SkyFile=1 was set.
      if self.skyLines:
        for line in self.skyLines:
          if (line.split()[0].lower() == shortname.lower() or \
              line.split()[0].lower() == shortname.lower()[:-5]):
            sky=float(line.split()[1])
            rms=float(line.split()[2])
            self.GlobalBlock.logfile.write('default.skies: '+shortname+'  Sky:'+str(sky)+'  RMS:'+str(rms))
            break

class GlobalBlock:
    def __init__(self,logfile,obsAlign,obsPars,obsFits,obsRoot,parDir): #WZ
      self.logfile = logfile
      self.modName = string.split(string.split(str(self))[0],'.')[0][1:] # module name
      self.inputList = []
      self.errorList = []
      self.outputList = {}
      self.obsAlign = obsAlign      # alignment files directory
      self.obsPars = obsPars        # directory for new par sets
      self.obsFits = obsFits        # directory for fits files to align,
                                    # default.skies, 
      self.inputRoot = string.replace(obsRoot,'datasets','ingest') # Input directory #WZ
      self.parDir = parDir 

class Nothinglogfile:
  def __init__(self):
    pass
  def write(self,S):
    pass

class NothingGlobalBlock(GlobalBlock):
   def __init__(self):     
     GlobalBlock.__init__(self,Nothinglogfile(),"","","")

class KWList:
    """ KWList gives the list of keywords and values unique to a specific
        filter. """
    def __init__(self):
      self.KWList = []

    def AddKW(self,key,value):
      self.KWList.append((key,value))

    def IsEqual(self,ff):
      OK = 1
      for kwentry in self.KWList:
        (key,value) = kwentry
        cval = ff[0].header.get(key)
        if cval != value:
          OK = 0

      return OK

    def IsPresent(self,ff):
      """ Useful for testing to see if the appropriate keywords are in 
          the header at all.  Thus for error checkings. """
      OK = 1
      for kwentry in self.KWList:
        (key,value) = kwentry
        cval = ff[0].header.get(key)
        if cval == None:
          OK = 0

      return OK

    def RaiseError(self,PathName):
      Str = '%s does not have the requisite keywords, e.g., ' % (PathName)
      for kwentry in self.KWList:
        (key,value) = kwentry
        Str = Str + '%s, ' % (key)      
      raise ValueError,Str


class SearchParam:
    """ SearchParam is the class which give the parameters used to identify
    stars or other concentrated objects off of which alignment is done.        

	The Parameters specify the following limits:
	  #lowLim:     min(2.3*b_axis,fwhm) >= lowLim (e.g., 1.55  for WFC)
	  lowLim:     min(2.3*b_axis,fwhm) >= lowLim (e.g., 2.00  for UVIS)
          hiLim:      max(2.3*a_axis,fwhm) < hiLim  (e.g., 20 for WFC)
          minAxisRat: b_axis/a_axis > minAxisRat   (e.g., 0.33 for WFC)  
          EdgeBuf:    Minimum distance (pixels) from image boundary.  
          inParFile:  Name of Default inParFile
          _lower_pixthresh_: Lower Pixel Threshold """
    def __init__(self,lowLim,hiLim,minAxisRat,EdgeBuf,inParFile,
                 _lower_pixthresh_):
      self.lowLim = lowLim
      self.hiLim = hiLim
      self.minAxisRat = minAxisRat
      self.EdgeBuf = EdgeBuf
      self.inParFile = inParFile
      self._lower_pixthresh_ = _lower_pixthresh_

class Filter:
    """ idctab is the full pathname of the relevant idctab file. """
    def __init__(self,Name,idctab,KWList,SearchParam):
      self.Name = Name
      self.idctab = idctab
      self.KWList = KWList
      self.SearchParam = SearchParam

    def GetSearchParam(self):
      return self.SearchParam

    def IsEqual(self,ff):
      return self.KWList.IsEqual(ff)

    def IsPresent(self,ff):
      return self.KWList.IsPresent(ff)

    def RaiseKWError(self,PathName):
      return self.KWList.RaiseError(PathName)

class Camera:
    """ Camera is the class which stores all parameters related to a specific
        camera/instrument """
    def __init__(self,Name,KWList,updown,OrientAdj):
      self.Name = Name
      #print "self.Filterlist"
      self.FilterList = []
      self.KWList = KWList
      self.updown = updown
      self.OrientAdj = OrientAdj

    def GetName(self):
      return self.Name

    def GetSearchParam(self):
      return self.SearchParam

    def AddFilter(self,Filter):
      #print "self.Filterlist add"
      #pdb.set_trace()
      self.FilterList.append(Filter)

    def IsEqual(self,ff):
      return self.KWList.IsEqual(ff)

    def IsPresent(self,ff):
      return self.KWList.IsPresent(ff)

    def RaiseKWError(self,PathName):
      return self.KWList.RaiseError(PathName)

    def GetFilter(self,ff,PathName):
      GoodFilter = None
      KWPresent = 0
      for Filter in self.FilterList:
        if Filter.IsEqual(ff):
          GoodFilter = Filter
        if (not KWPresent) and (Filter.IsPresent(ff)):
          KWPresent = 1
      if KWPresent == 0:
        self.FilterList[0].RaiseKWError(PathName)

      if GoodFilter == None:
        raise ValueError,"Filter in file '%s' not found in %s filter list." \
              % (PathName,self.Name)

      return GoodFilter

class FitsFile:
    """ FitsFile is the module/class which does object identification /
        matching for individual fits files. 

        parameters:
          PathName = filename (path + name)
          Asn = association name
          logfile = open file to log information
          GlobalBlock = 
          SearchParam = search parameters for stars (added in EvaluateFilt)
          sci = list to accommodate scientific extensions
          NumSci = the number of scientific extensions
          ref = Is reference image
          simpleimage = Is Simple Image (not from telescope or having multiple extensions)
          outscale = 
          V2center =
          V3center =
          updown = 
          matchFile = 
    """
    def __init__(self,PathName,Asn,GlobalBlock,SearchParam = None, useMinSky = 0):
      self.PathName = PathName
      self.Asn = Asn
      self.GlobalBlock = GlobalBlock
      # recopy for logfile for convenience
      self.logfile = GlobalBlock.logfile  
      self.preMatchDict = {}
      self.ref = 0
      self.extref = 0
      self.simpleimage = 0
      self.firstext = 0      # if simpleimage == 1, then this extension contains science image
      self.useMinSky = useMinSky
      self.SearchParam = SearchParam

    def ChangeToRef(self):
      self.ref = 1

    def ChangeToExtRef(self):
      self.extref = 1

    def CheckExtRefFlags(self):
      if (self.ref == 0) or (self.extref == 0):
        raise ValueError,"ERROR -- External reference image does not have ref and extref flags set."   

    def DetermineFirstImage(self):
      ff = pyfits.open(self.PathName)
      for extnum in range(len(ff)):
        dim1 = 0
        if ff[extnum].header.get('NAXIS1'):
          try:
            dim1 = string.atof(ff[extnum].header['NAXIS1'])
          except:
            dim1 = 0
        dim2 = 0
        if ff[extnum].header.get('NAXIS2'):
          try:
            dim2 = string.atof(ff[extnum].header['NAXIS2'])
          except:
            dim2 = 0
        print dim1,dim2
        if dim1*dim2 > 10000:
          self.firstimage = extnum
          break
      ff.close()
    # pdb.set_trace()
    def EvaluateFilt(self,CameraList):
      ff = pyfits.open(self.PathName)
      KWPresent = 0
      if len(CameraList) == 0:
        raise ValueError,"ERROR -- No cameras in CameraList."

      #print "finding filter"
      #pdb.set_trace()
      if not self.SearchParam:
        self.SearchParam = None
      for Camera in CameraList:                
	if Camera.IsEqual(ff):
          self.Camera = Camera
          self.logfile.write("Found DETECTOR = "+Camera.GetName()+" in first fitsfile.")
          self.updown = Camera.updown
          try:
            Filter = Camera.GetFilter(ff,self.PathName)
          except:
            self.SearchParam = Camera.GetSearchParam()
          else:
            self.Filter = Filter
            self.SearchParam = Filter.GetSearchParam()
	if Camera.IsPresent(ff):
          KWPresent = 1

      #pdb.set_trace()
      if KWPresent == 0:
        if (not self.ref):
          CameraList[0].RaiseKWError(self.PathName)
        self.DetermineFirstImage()
        self.simpleimage = 1

      if (self.SearchParam == None) and (not self.ref):
        raise ValueError,"Detector in file '%s' not found in detector list." \
              % (self.PathName)

      ff.close()

    def DetermineExtStruct(self):
      """ Read the extension struct and populate NumSci,sci """

      if self.simpleimage == 1:
        raise Exception, "This routine (DetermineExtStruct) should not be called if simpleimage."

      ff = pyfits.open(self.PathName)
      self.sci = []                             
 # Main list to accommodate multiple sci extensions, each element is a dict
      self.NumSci=0
      for extnum in range(1,len(ff)):
      	xname = ff[extnum].header['EXTNAME'].upper()
        if xname == 'SCI':
           # start new dictionary for this sci extension
           self.sci.append({})
           self.sci[self.NumSci]['extnum'] = extnum
           self.sci[self.NumSci]['chipID'] = ff[extnum].header.get('CCDCHIP')
           self.NumSci += 1
        del extnum,xname
      if self.NumSci != len(self.sci):
        raise Exception, "BUG!"
        
      # now some chipID/SCI-exten error checking
      if self.NumSci == 0:
        raise Exception,"No science extensions found in "+self.PathName
      if self.NumSci == 1:
        if self.sci[0]['chipID'] == None:
          self.sci[0]['chipID'] = 1
      else:
        for ii in range(self.NumSci):        # multiple science exentsions
          if self.sci[ii]['chipID'] == None:
            raise ValueError,"Multi-Chip image "+firstfitsfile+ \
                    " has no CCDCHIP ID in extension "+self.sci[ii]['extnum']
          for jj in range(ii+1,self.NumSci):
            if self.sci[ii]['chipID'] == self.sci[jj]['chipID']:
              raise ValueError,"Two science extensions with same CCDCHIP in "+\
                              self.PathName

      # log what we've found out
      logtxt = 'Found '+str(self.NumSci)+' science extensions to process:'
      for ii in range(self.NumSci):
        logtxt = logtxt + '\n              '+str(self.sci[ii])
      self.logfile.write(logtxt)
      ff.close()      

    def removeUnpicklable(self):
      if self.simpleimage == 0:
        for ii in range(self.NumSci):
          del self.sci[ii]['CxMat'],self.sci[ii]['CyMat']      

    def setupIDCinfo(self,idcName):
       """ Get IDC matrix transformations for each element of self.sci[]. """
       # get the IDC filename from 1st image header, if not already specified
       # NOTE: this means that all images must use the same IDC file
       # The fUtil function getIdcTab does all the iref$ parsing and searching
       # for IDCTAB file specified in the header of the fitsfile sent to it and 
       # returns the full path name to the idc file, e.g.:
       # >>> alIm.idcName
       # '~/zheng/wfp/reffiles/wfc_idc.fits'

       if self.simpleimage == 1:
         raise Exception, "This routine (DetermineExtStruct) should not be called if simpleimage."
       self.logfile.write('Getting IDC distortion information . . .')
       if not idcName:
         firstfitsfile = self.PathName
         self.logfile.write('Taking IDCTAB specified in header of '+firstfitsfile)
         try:
           idcName = fUtil.getIdcTab(firstfitsfile)
         except Exception,err:
           self.GlobalBlock.errorList.append(('fUtil',str(err)))
           raise Exception,err
         self.logfile.write('Using ' + idcName + ' IDC table for image '+firstfitsfile)
         del firstfitsfile
            
       self.idcName = idcName
       # now that we have IDCTAB file path+name, get the IDC matrices:
       # Loop over sci extensions
       for ii in range(self.NumSci):
         chipID = self.sci[ii]['chipID']
         # pdb.set_trace()
         if self.Camera.Name == 'WFC':
           CxMatrix,CyMatrix,refdict,Norder = matutil.getACSIDCinfo(self.idcName,chipID)
         else:
           CxMatrix,CyMatrix,refdict,Norder = matutil.getIDCinfo(self.idcName,chipID)
         #CxMatrix,CyMatrix,refdict,Norder = matutil.getIDCinfo(self.idcName,chipID,filter=self.Filter.Name)
         ## for convenience later, convert 2-d 4x4 matrices to 1-d 10-element arrays
         ## and save as class attributes for later back conversion
         # NO! no need to do this anymore!
         # xcoeff,ycoeff = matutil.convertIDCtrans(CxMatrix,CyMatrix)

         # Now save important IDC info in dictionary for later use
         self.sci[ii]['CxMat']   = CxMatrix
         self.sci[ii]['CyMat']   = CyMatrix
         # self.sci[ii]['xcoeff']  = xcoeff
         # self.sci[ii]['ycoeff']  = ycoeff
         self.sci[ii]['refdict'] = refdict
         # test that pscales agree for chips
         if ii and str(refdict['PSCALE']) != str(self.sci[0]['refdict']['PSCALE']):
           raise Exception,"IDC chip output pscale's disagree."

         del CxMatrix,CyMatrix,refdict,Norder
             
       self.outscale = self.sci[0]['refdict']['PSCALE']
       #print "self.outscale should be defined! Please check it: ____________by xingxing"# by xingxing
       #print  self.outscale  
       #pdb.set_trace() #by XX
       for ii in range(1,self.NumSci):
         if self.sci[ii]['refdict']['PSCALE'] != self.outscale:
           errtxt = 'Error. Chips specify different default scales!'
           self.GlobalBlock.errorList.append(self.modName,errtxt)
           self.logfile.write(errtxt)
           raise Exception,errtxt

       self.V2center = 0.0
       self.V3center = 0.0
       for ii in range(self.NumSci):
         self.V2center += self.sci[ii]['refdict']['V2REF']
         self.V3center += self.sci[ii]['refdict']['V3REF']
       self.V2center /= self.NumSci
       self.V3center /= self.NumSci
       del ii
       return

    def recenter(self):
       if self.simpleimage == 1:
         raise Exception, "This routine (DetermineExtStruct) should not be called if simpleimage."
       ff = pyfits.open(os.path.join(self.GlobalBlock.obsFits,self.PathName), 'update')
       for ii in range(self.NumSci):
	 extnum = self.sci[ii]['extnum']
	 refwcs = wcsclass.BasicWCS(ff[extnum].header)
	 naxis1 = ff[extnum].header['NAXIS1']
	 naxis2 = ff[extnum].header['NAXIS2']
	 ff[extnum].header['CRPIX1'] = naxis1/2.
	 ff[extnum].header['CRPIX2'] = naxis2/2.
	 ra,dec = refwcs.xy2rd((naxis1/2.,naxis2/2.))
	 ff[extnum].header['CRVAL1'] = ra
	 ff[extnum].header['CRVAL2'] = dec
       ff.flush()
       ff.close()

    def returnundistortradec(self,scidex,xp,yp):
      xRef = self.sci[scidex]['refdict']['XREF']   
      yRef = self.sci[scidex]['refdict']['YREF']   
      x_asec,y_asec = \
        matutil.xyIDCtrans(xp-xRef, yp-yRef, \
                           self.sci[scidex]['CxMat'], self.sci[scidex]['CyMat'])
      x_asec = x_asec/self.outscale + self.refwcs[scidex].wcs['CRPIX1']
      y_asec = y_asec/self.outscale + self.refwcs[scidex].wcs['CRPIX2']
      (ra,dec) = self.refwcs[scidex].xy2rd((x_asec,y_asec))
      return (ra,dec)

    def setupWCS(self):
      if self.simpleimage == 0:
        raise Exception, "This routine (DetermineExtStruct) should not be called if not simpleimage."
      rf = pyfits.open(self.PathName)
      self.wcs = wcsclass.BasicWCS(rf[self.firstimage].header)
      self.NumSci = 1
      self.dimx = string.atof(rf[self.firstimage].header['NAXIS1'])
      self.dimy = string.atof(rf[self.firstimage].header['NAXIS2'])
      self.sci = {}
      self.sci[0] = {}
      self.sci[0]['extnum'] = self.firstimage
      self.outscale = pixelScale(rf[self.firstimage],self.PathName)
      self.angle = wcsangle(rf[self.firstimage],self.PathName)
      (ra,dec) = self.wcs.xy2rd((self.dimx/2,self.dimy/2))
      self.wcs.wcs['CRVAL1'] = ra
      self.wcs.wcs['CRVAL2'] = dec
      self.wcs.wcs['CRPIX1'] = self.dimx/2
      self.wcs.wcs['CRPIX2'] = self.dimy/2
      rf.close()
      del rf

    def makeMatchCat(self,skyLines,irafskytoo):
      NGood = 0
      for scidex in range(self.NumSci):
        NGood += self._makeCat(scidex, skyLines,irafskytoo, Retry=0, inParFileName=self.SearchParam.inParFile)
      Try = 0
      #if (NGood < 10):
      #  pdb.set_trace()
      while ((NGood > 4000) or (NGood < 10)) and (Try < 2):  # then do it again
        self.logfile.write("NGood = "+str(NGood)+" for "+os.path.split(self.PathName)[1]+". Adjust det pars; retry.")
        print "NGood = "+str(NGood)+" for "+os.path.split(self.PathName)[1]+"  Adjust and retry."
        # adjust params
        inpars  = pUtil.readInParFile(self.SearchParam.inParFile)
        if NGood < 10:
          if Try == 0:
            inpars["DETECT_THRESH"] = self.SearchParam._lower_pixthresh_
            inpars["ANALYSIS_THRESH"] = self.SearchParam._lower_pixthresh_
          else:
            inpars["DETECT_THRESH"] = self.SearchParam._lower_pixthresh_ / 1.3
            inpars["ANALYSIS_THRESH"] = self.SearchParam._lower_pixthresh_ / 1.3
          _cmttext = "sxt detect_thresh lowered to "+str(inpars["DETECT_THRESH"])
          print _cmttext
 
        else:
          print inpars["DETECT_THRESH"]
          if Try == 0:
            inpars["DETECT_THRESH"] = str((int(inpars["DETECT_THRESH"]) * 2))
            inpars["ANALYSIS_THRESH"] = str((int(inpars["ANALYSIS_THRESH"]) * 2))
          else:
            inpars["DETECT_THRESH"] = str((int(inpars["DETECT_THRESH"]) * 3))
            inpars["ANALYSIS_THRESH"] = str((int(inpars["ANALYSIS_THRESH"]) * 3))
          print inpars["DETECT_THRESH"]
          _cmttext = "sxt detect_thresh increased to "+str(inpars["DETECT_THRESH"])
          print _cmttext

          # pdb.set_trace()
        inParFileName   = os.path.join(self.GlobalBlock.obsPars, os.path.splitext(self.PathName)[0]+"_align.inpar")
        self._writeNewParFile(self.PathName,inpars,inParFileName)
        del _cmttext
        del inpars
                
        NGood = 0
        inParFileName=self.SearchParam.inParFile
        for scidex in range(self.NumSci):
          NGood += self._makeCat(scidex, skyLines,irafskytoo,Retry=1, inParFileName=inParFileName)
          # NGood += self._makeCat(scidex, skyLines,irafskytoo,Retry=1,inParFileName=os.path.join(self.GlobalBlock.obsPars,os.path.splitext(self.PathName)[0]+"_align.inpar"))
        Try += 1

      if NGood < 6:
        self.logfile.write("Ngood = "+str(NGood)+" still too few; will not match.")

      return

    def _writeNewParFile(self,fitsfile,inpars,inParFileName):
        """ generate a new .inpar file for align should the NgoodObjects test
	require it.
	"""
	dir,file = os.path.split(fitsfile)
	self.logfile.write('Private _writeNewParFile method called on '+ file)
	
	# Make the new inpar file named for the passed fitsfile.
	inparfile  = open(inParFileName, 'w')
	# write out the inpars.  some call to _setinpar should have made to change the self.inpars dict.
	for param, value in inpars.items():
	    inparfile.write("%-18s   %-20s\n" % (param, value))             
	inparfile.close()
	# make a new self.specialInParFile.
	self.logfile.write('Wrote new special align input parameter file: '+inParFileName)
        return

    def _makeCat(self,scidex,skyLines,irafskytoo,Retry=0,inParFileName=None):
        """ Run SExtractor on one fitsfile.  fitsfile MUST be a full pathname.
            The <ref> argument defaults to zero (false) and indicates whether this
            image is a reference image or not. Pass 1 for true.
              e.g., _makeCat('/tmp/reference.fits', ref=1) """
        fitsfile = os.path.join(self.GlobalBlock.obsFits,self.PathName)

        if os.path.isabs(fitsfile):
          pass
        else:
           self.GlobalBlock.errorList.append((self.GlobalBlock.modName,"Path is not a full path."))
           raise IOError, "Path is not a full path."
        imName  = string.split(os.path.basename(fitsfile),'.')[0]

        if not self.outscale or self.outscale < 1e-9:
	   self.GlobalBlock.errorList.append((self.GlobalBlock.modName,"_makeCat: outscale not set!"))
           self.logfile.write("Error: _makeCat: outscale not set!")
           raise Exception,"_makeCat: outscale not set!"
        
        ext = self.sci[scidex]['extnum']
        msuffix = 'matchin'
        
        if self.extref == 1:
           # set the reference image input file for matching
           matchFile    = os.path.join(self.GlobalBlock.obsAlign, imName+'_ref.'+msuffix)
           # set the reference catalog to be this file
           logtxt = 'Will run SExtractor on external Ref Image '+fitsfile+\
                 ' extension '+str(ext)+'\n \t to get match input: '+matchFile
           self.GlobalBlock.logfile.write(logtxt)
           print logtxt
           irafskytoo = 0
        else:
           # an image's input file for matching
           matchFile    = os.path.join(self.GlobalBlock.obsAlign, imName+'.'+msuffix)
           self.GlobalBlock.outputList[os.path.basename(matchFile)] = [os.path.basename(fitsfile)]

        if scidex > 0:
           matches = open(matchFile,'a')
        else:
           matches = open(matchFile,'w')
        self.matchFile = matchFile

         # get the simple fits version of the data and run sextractor on it
        simplefits,NX,NY = _getsimplefits(fitsfile, ext, self.GlobalBlock, scidex=scidex)
        #pdb.set_trace()   #by xingxing
        if Retry:
           catalog,sky,rmsnoise = _runSXtractor(simplefits,inParFileName,self.GlobalBlock,skyLines,irafskytoo)
        else:
           catalog,sky,rmsnoise = _runSXtractor(simplefits,self.SearchParam.inParFile,self.GlobalBlock,skyLines,irafskytoo)

        # make the selectSet diectionary from the header of the SExtractor catalog
        # Remember: column id is one greater than any array index you might use.
        selectSet = pUtil.makeHeaderDict(catalog)

        # get ready for the catalog-trimming loop
        if self.simpleimage == 0:
          xRef = self.sci[scidex]['refdict']['XREF']   
          yRef = self.sci[scidex]['refdict']['YREF']   

        NrawObjects = NgoodObjects = 0
        self.logfile.write('Searching ' + imName +' catalog for matchable sources...')

        # start populating the Match dictionary
        fitsKeyName = os.path.split(fitsfile)[1]  # key for the MatchDict dictionary
        # define the element to be another dictionary

        # initialize the preMatchDict[scidex] dictionary
        self.preMatchDict[scidex] = {} 
        self.preMatchDict[scidex]['fullPath'] = fitsfile
        self.preMatchDict[scidex]['simplefits'] = simplefits
        self.preMatchDict[scidex]['Raw_Xsize'] = NX
        self.preMatchDict[scidex]['Raw_Ysize'] = NY
 
        lowLim = self.SearchParam.lowLim
        hiLim = self.SearchParam.hiLim
        minAxisRat = self.SearchParam.minAxisRat
        EdgeBuf = self.SearchParam.EdgeBuf

        self.logfile.write("rej crit for "+imName+":  min: %.3f  max: %.3f  minrat: %.3f"\
                              %(lowLim,hiLim,minAxisRat))
        # pdb.set_trace()
        #i=j=0
        for line in catalog:
          if '#' not in line:
             fields = string.split(line)
             # Note: subtract 1 from the column id to get list index
             NrawObjects += 1
             ID     = int(fields[selectSet['NUMBER'] - 1])
             xp     = float(fields[selectSet['X_IMAGE']-1])
             yp     = float(fields[selectSet['Y_IMAGE']-1])
             fwhm   = float(fields[selectSet['FWHM_IMAGE']-1])
             mag    = float(fields[selectSet['MAG_BEST']-1])
             magerr = float(fields[selectSet['MAGERR_BEST']-1])
             a_axis = float(fields[selectSet['A_IMAGE']-1])
             b_axis = float(fields[selectSet['B_IMAGE']-1])
             axisRat = b_axis/a_axis
                
             if min(2.3*b_axis,fwhm) >= lowLim and max(2.3*a_axis,fwhm)< hiLim and axisRat>minAxisRat:
               #i=i+1
               if (xp > EdgeBuf and yp > EdgeBuf and \
                   xp < (NX-EdgeBuf) and yp < (NY-EdgeBuf)):
                 #j=j+1
                 # convert x,y position to arcsec using saved IDC info
 	         if self.simpleimage == 0:
                   x_asec,y_asec = matutil.xyIDCtrans(xp-xRef, yp-yRef, \
                                                      self.sci[scidex]['CxMat'], self.sci[scidex]['CyMat'])

                   if self.NumSci > 1:
                     x_asec = x_asec + (self.sci[scidex]['refdict']['V2REF'] - self.V2center)
                     y_asec = y_asec - (self.sci[scidex]['refdict']['V3REF'] - self.V3center)

                   x_opix = x_asec
                   y_opix = y_asec
                 else:
                   x_opix = (xp - self.wcs.wcs['CRPIX1']) * self.outscale
                   y_opix = (yp - self.wcs.wcs['CRPIX2']) * self.outscale

                 ###> test! 23/sep/02 test! <###
                 x_opix -= 1*self.outscale
                 y_opix -= 1*self.outscale
                 # NOTE: should actually ensure that outscale
                 # is same as drizzle outscale for this
                 ###############################

	         gdLin='%3d %10.5f %10.5f %9.4f %7.4f %8.3f %8.3f %6.3f %6.3f %5.3f %.4f\n'\
                       %(ID,x_opix,y_opix,mag,magerr,xp,yp,fwhm,a_axis,b_axis,axisRat)
	         matches.write(gdLin)
	         NgoodObjects = NgoodObjects + 1
        matches.close()
        self.preMatchDict[scidex]['skyval'] = sky
        self.preMatchDict[scidex]['skyrms'] = rmsnoise
        self.preMatchDict[scidex]['Nraw']   = NrawObjects
        self.preMatchDict[scidex]['Ngood']  = NgoodObjects
        self.logfile.write('Matchfile '+matchFile+' for '+imName+' ext '+str(ext)+' written.\n\t'+\
                            'Using '+str(NgoodObjects)+' objects of '+str(NrawObjects)+' found.')
        del catalog,sky,simplefits,NX,NY,rmsnoise,matches
        return NgoodObjects

    def delMatchDicts(self):
        """
        Delete self.MatchDict, so can remake with hdr_shifts
        """
        del self.MatchDict

    def mergeMatchDicts(self):
        """
        Merge multiple self.preMatchDict[i] dict's into 1 self.MatchDict
        """
        # just a final check of this...
        if self.NumSci != len(self.sci):
            raise ValueError,"NumSci inconsistency."

        # Get: (1) list of images; (2) list of image MatchDict Keys
        matchKeyList = self.preMatchDict[0].keys()
        
        # Loop over images 
        self.MatchDict = {}
        self.MatchDict['NumSci'] = self.NumSci
        # Loop over MatchDict keys
        for matchKey in matchKeyList:
          if matchKey == 'fullPath':
            self.MatchDict[matchKey] = self.preMatchDict[0][matchKey]
          elif matchKey == 'simplefits':
            pass
          elif matchKey[0] == 'N':
            ival = 0
            # an integer number to be summed; loop over sci extensions
            for ii in range(self.NumSci):
              ival += self.preMatchDict[ii][matchKey]
            self.MatchDict[matchKey] = ival
            del ival

          elif (self.useMinSky and (matchKey == 'skyval')):
            _tmpsky=[]
            for ii in range(self.NumSci):
              _tmpsky.append(self.preMatchDict[ii][matchKey])
            self.MatchDict[matchKey] = min(_tmpsky)
            del _tmpsky

          elif matchKey[0:3] == 'sky' or matchKey[0:3] == 'Raw' or matchKey[0:3] == 'hdr':
            fval = 0.0
            # straight average
            for ii in range(self.NumSci):
              fval += self.preMatchDict[ii][matchKey]
            self.MatchDict[matchKey] = fval/self.NumSci
            del fval

          # matchKey[0:3] == 'med' or 'fit' now only occur after this method is called
          else:
            # we don't know what this key is!
            errtxt = "Error. Unrecognized key: "+matchKey
            self.GlobalBlock.errorList.append((self.GlobalBlock.modName,errtxt))
            self.logfile.write(errtxt)
            raise Exception,errtxt
                    
        del matchKey,ii,matchKeyList
        # finished merging MatchDict's
        return

    # def _do_the_matching(self, match_cmd, matchinfile, matchf, default_scale_tol, maxfitrms):
    #     """ run match program and parse results;
    #     returns 0 if all went well, -1 if it failed.
    #     """
    #     sproc  = popen2.Popen3(match_cmd,1)
    #     output = sproc.fromchild.readlines()
    #     errs   = sproc.childerr.readlines()
    #     if errs and not (len(errs) == 1 and errs[0][0:27] == 'atFindMedtf: RMS <= 0.0, so'):
    #         self.logfile.write('match produced the following message(s) on '+matchinfile)
    #         for f in errs:
    #             self.logfile.write(string.strip(f))
    #             self.GlobalBlock.errorList.append(('match',string.strip(f)))
    # 
    #     matchline = ''
    #     for line in output:
    #         if line[0:6] == 'TRANS:':
    #             matchline = line
    #         elif line[0:6] == 'MEDTF:':
    #             pass
    #         else:
    #             self.logfile.write("Error. Unexpected output from match program:\n\t\t\t"+line)
    #             sproc.fromchild.close()
    #             sproc.childerr.close()
    #             return -1
    # 
    #     if (not matchline):
    #         self.logfile.write("Error. Match program failed on "+matchinfile)
    #         return -1
    # 
    #     # ok, it appears that match produced trans and medtf structures
    # 
    #     # process the TRANS struct first
    #     fields = string.split(matchline)
    #     if self._process_transfit(fields,default_scale_tol,maxfitrms) < 0:
    #         sproc.fromchild.close()
    #         sproc.childerr.close()
    #         return -1
    #     
    #     # for bookkeeping, it's important to write out the lines only if
    #     # all went well with the match processing
    #     matchf.write(self.PathName +' '+matchline[len('TRANS:'):])
    # 
    #     del matchline
    #     sproc.fromchild.close()
    #     sproc.childerr.close()
    #     return 0

    def _process_wcs_out(self,refImage): #WZ Dec 2010
      # feference is another similar image, to determine the RELATIVE shifts between them.
        """
        inf = open(file)
	L = inf.readlines()
        aa = string.atof(L[1].split()[1])
        bb = string.atof(L[2].split()[1])
        cc = string.atof(L[3].split()[1])
        dd = string.atof(L[4].split()[1])
        ee = string.atof(L[5].split()[1])
        ff = string.atof(L[6].split()[1])
        Sigf = 10
        NumberFitted = 10
        NumberMatches = 10
        sigx = 0.001
        sigy = 0.001
        """

        #pdb.set_trace()
        im  = refImage + '[1]'
        ref = wcsutil.WCSObject(im)
        im  = self.PathName + '[1]'
        wcs = wcsutil.WCSObject(im)
        ra = wcs.crval1
        dec = wcs.crval2
        x0, y0 = ref.rd2xy((ra,dec))
        x, y = wcs.rd2xy((ra,dec))
        dx = x - x0
        dy = y - y0
        dr = wcs.orient - ref.orient
        NX  = self.MatchDict['Raw_Xsize']
        NY  = self.MatchDict['Raw_Ysize']
        #pdb.set_trace()
        # finally, populate MatchDict fit keywords:
        self.MatchDict['Nmatch']  = -1
        self.MatchDict['Nfit']    = 0
        self.MatchDict['fit_sig'] = 0
        self.MatchDict['fit_xshift'] = -dx * self.outscale
        self.MatchDict['fit_yshift'] = -dy * self.outscale
        self.MatchDict['fit_xshift_err'] = 0.
        self.MatchDict['fit_yshift_err'] = 0.
        self.MatchDict['fit_angle']     = -dr # wcs.orient
        self.MatchDict['fit_angle_err'] = 0.
        self.MatchDict['fit_delt_angle'] = 0. # dr
        
        #self.logfile.write("Using header WCS information for alignment")
        return 0

    def _process_simplematchout(self, file):
        inf = open(file)
	L = inf.readlines()
        aa = string.atof(L[1].split()[1])
        bb = string.atof(L[2].split()[1])
        cc = string.atof(L[3].split()[1])
        dd = string.atof(L[4].split()[1])
        ee = string.atof(L[5].split()[1])
        ff = string.atof(L[6].split()[1])
        Sigf = 10
        NumberFitted = 10
        NumberMatches = 10
        sigx = 0.001
        sigy = 0.001

        # ok, now for the "derived" parameters
        fit_scale_x = math.sqrt(bb*bb + ee*ee)
        fit_scale_y = math.sqrt(cc*cc + ff*ff)
        if bb == 0.:
            angle_x  =  90.0
        else:
            angle_x  = -(180./math.pi) * math.atan2(ee,bb)
        if ff == 0.:
            angle_y  =  90.0
        else:
            angle_y  =  (180./math.pi) * math.atan2(cc,ff)

        # check to see if it's a good transformation
        NX  = self.MatchDict['Raw_Xsize']
        NY  = self.MatchDict['Raw_Ysize']
        
        errtxt = ''
        if abs(angle_x - angle_y - 360) < 1:
            angle_x -= 360
        if abs(angle_x - angle_y + 360) < 1:
            angle_x += 360

        NX  = self.MatchDict['Raw_Xsize']
        NY  = self.MatchDict['Raw_Ysize']
        r_0 = math.sqrt(1.0*NX*NY*self.NumSci) / 3    # this is roughly the mean radius for objs
        sqrtNorm  = math.sqrt(NumberMatches - 1.0)
        angle_err = (180./math.pi) * math.atan(math.sqrt(sigx**2 + sigy**2) / (r_0 * sqrtNorm))
        # NOTE: I use NumberMatches for the error estimation because the
        # errors are based on the 3-sigma clipped scatters sigx and sigy,
        # which are based on the NumberMatches objects
        
        # finally, populate MatchDict fit keywords:
        self.MatchDict['Nmatch']  = NumberMatches
        self.MatchDict['Nfit']    = NumberFitted
        self.MatchDict['fit_sig'] = Sigf
        self.MatchDict['fit_xshift'] = aa
        self.MatchDict['fit_yshift'] = dd
        self.MatchDict['fit_xshift_err'] = sigx / sqrtNorm
        self.MatchDict['fit_yshift_err'] = sigy / sqrtNorm
        self.MatchDict['fit_angle']     = (angle_x + angle_y)/2.0
        self.MatchDict['fit_angle_err'] = angle_err
        self.MatchDict['fit_delt_angle'] = (angle_x - angle_y)
        
        self.logfile.write("Success!  Matched "+str(NumberMatches)+" objects; "+"sigma_fit = "+str(Sigf))
        return 0

    # def _process_transfit(self, fields, default_scale_tol, maxfitrms):
    #     """Parse, process and write out match TRANS output line."""
    #     if len(fields) != 12 or fields[0] != 'TRANS:':
    #         errtxt = "Error in parsing match TRANS output: \n"+str(fields)
    #         self.logfile.write(errtxt)
    #         return -1
    # 
    #     # read the fit parameters into variables
    #     aa        = float(fields[1].split('=')[1])
    #     bb        = float(fields[2].split('=')[1])
    #     cc        = float(fields[3].split('=')[1])
    #     dd        = float(fields[4].split('=')[1])
    #     ee        = float(fields[5].split('=')[1])
    #     ff        = float(fields[6].split('=')[1])
    #     Sigf      = float(fields[7].split('=')[1])
    #     NumberFitted  = int(fields[8].split('=')[1])
    #     NumberMatches = int(fields[9].split('=')[1])
    #     sigx      = float(fields[10].split('=')[1])
    #     sigy      = float(fields[11].split('=')[1])
    # 
    #     # ok, now for the "derived" parameters
    #     fit_scale_x = math.sqrt(bb*bb + ee*ee)
    #     fit_scale_y = math.sqrt(cc*cc + ff*ff)
    #     if bb == 0.:
    #         angle_x  =  90.0
    #     else:
    #         angle_x  = -(180./math.pi) * math.atan2(ee,bb)
    #     if ff == 0.:
    #         angle_y  =  90.0
    #     else:
    #         angle_y  =  (180./math.pi) * math.atan2(cc,ff)
    # 
    #     # check to see if it's a good transformation
    #     NX  = self.MatchDict['Raw_Xsize']
    #     NY  = self.MatchDict['Raw_Ysize']
    #     scale_tolerance_ = default_scale_tol * 6144.0/(NX+NY)
    #     
    #     errtxt = ''
    #     if abs(angle_x - angle_y - 360) < 1:
    #         angle_x -= 360
    #     if abs(angle_x - angle_y + 360) < 1:
    #         angle_x += 360
    # 
    #     if (abs(angle_x - angle_y) > 1):
    #         errsnumpy.txt += ' x and y rotation angles disagree by > 1 deg: '+str(angle_x)+' '+str(angle_y)
    #     if (abs(1. - fit_scale_x) > scale_tolerance_) or (abs(1. - fit_scale_y) > scale_tolerance_):
    #         errtxt += ' unexpected distortion in match fit. '
    #     if max(sigx,sigy) > maxfitrms:
    #         if errtxt: errtxt += 'Also'
    #         errtxt += ' x,y rms of match fit too large:  '+str(sigx)+' '+str(sigy)+''
    #     if errtxt:
    #         errtxt = 'Error. Match: '+errtxt
    #         self.logfile.write(errtxt)
    #         self.GlobalBlock.errorList.append((self.GlobalBlock.modName,errtxt))
    #         return -1
    # 
    #     NX  = self.MatchDict['Raw_Xsize']
    #     NY  = self.MatchDict['Raw_Ysize']
    #     r_0 = math.sqrt(1.0*NX*NY*self.NumSci) / 3    # this is roughly the mean radius for objs
    #     sqrtNorm  = math.sqrt(NumberMatches - 1.0)
    #     angle_err = (180./math.pi) * math.atan(math.sqrt(sigx**2 + sigy**2) / (r_0 * sqrtNorm))
    #     # NOTE: I use NumberMatches for the error estimation because the
    #     # errors are based on the 3-sigma clipped scatters sigx and sigy,
    #     # which are based on the NumberMatches objects
    #     
    #     # finally, populate MatchDict fit keywords:
    #     self.MatchDict['Nmatch']  = NumberMatches
    #     self.MatchDict['Nfit']    = NumberFitted
    #     self.MatchDict['fit_sig'] = Sigf
    #     self.MatchDict['fit_xshift'] = aa
    #     self.MatchDict['fit_yshift'] = dd
    #     self.MatchDict['fit_xshift_err'] = sigx / sqrtNorm
    #     self.MatchDict['fit_yshift_err'] = sigy / sqrtNorm
    #     self.MatchDict['fit_angle']     = (angle_x + angle_y)/2.0
    #     self.MatchDict['fit_angle_err'] = angle_err
    #     self.MatchDict['fit_delt_angle'] = (angle_x - angle_y)
    #     
    #     self.logfile.write("Success!  Matched "+str(NumberMatches)+" objects; "+"sigma_fit = "+str(Sigf))
    #     return 0

    def _bail_out(self, matchf, quiet=0):
        """ This function is for damage control in the event of a total
        match failure, or for populating extra MatchDict elements when using
        a default.shifts and match is not being run.
        It has to do the following:
          1) write a dummy line to the matchf file
          3) populate with dummy values all the MatchDict
               parameters that would otherwise have been set by
               _process_transfit() and _process_median
        By setting NumberMatches, NumberFitted, NumberMedian all to zero, we can
        tell later whether or not we should use 'measured' or 'commanded' shifts.
        """
        # define and write out the dummy lines:
        matchline = self.PathName+'  a=0.00            b=0.0             c=0.0             d=0.0          '+\
                  '   e=0.0             f=0.0             sig=999.9      Nr=0  Nm=0   sx=999.9  sy=999.9\n'
        matchf.write(matchline)
        del matchline

        # populate MatchDict fit_ these dummy numbers if not in default shift file
        if not quiet:
            self.logfile.write('Cannot use objects in '+self.PathName+' for shifts.')
        self.MatchDict['Nmatch']  = 0
        self.MatchDict['Nfit']    = 0
        self.MatchDict['fit_sig'] = 999.9
        self.MatchDict['fit_xshift'] = 0.0
        self.MatchDict['fit_yshift'] = 0.0
        self.MatchDict['fit_xshift_err'] = 999.9
        self.MatchDict['fit_yshift_err'] = 999.9
        self.MatchDict['fit_angle']      = 0.0
        #pdb.set_trace()
        self.MatchDict['fit_angle_err']  = 180.0
        self.MatchDict['fit_delt_angle'] = 0.0

        return

    def _getDefShifts(self, shiftline):
        "Get default shifts specified in appropriate line of default.shifts file."
        
        fields = string.split(shiftline)
        errtxt = ''
        if fields[0] != self.PathName:
            errtxt += "Error: fatal progenitor inconsistency in _getDefShifts. "
        if len(fields) < 7:
            errtxt += "Error: not enough fields for "+self.PathName+" in default shifts file."
        if errtxt:
            self.logfile.write(errtxt)
            self.GlobalBlock.errorList.append((self.GlobalBlock.modName,errtxt))
            raise Exception,errtxt

        xshift     = float(fields[1])
        xshift_err = float(fields[2])
        yshift     = float(fields[3])
        yshift_err = float(fields[4])
        angle      = float(fields[5])
        angle_err  = float(fields[6])
        if len(fields) > 7:
            Nmatch     = int(fields[7])
        else:
            Nmatch     = 1

        self.logfile.write('Using '+self.PathName+' default shifts:  dx = '+\
                           str(xshift)+'  dy = '+str(yshift))
        self.MatchDict['Nmatch']  = Nmatch
        self.MatchDict['Nfit']    = Nmatch
        self.MatchDict['fit_sig'] = math.sqrt((xshift_err**2 + \
                                                             yshift_err**2
)*Nmatch)
        self.MatchDict['fit_xshift'] = xshift
        self.MatchDict['fit_yshift'] = yshift
        self.MatchDict['fit_xshift_err'] = xshift_err
        self.MatchDict['fit_yshift_err'] = yshift_err
        self.MatchDict['fit_angle']      = angle
        self.MatchDict['fit_angle_err']  = angle_err
        self.MatchDict['fit_delt_angle'] = 0.0

        return

    def _use_simplematch_for_sol(self, refMatchFile, hdr_dx, hdr_dy, hdr_rot, file, tolmatch, default_scale_tol, maxfitrms, outfile, matchf):
        # fire up simplematch ... which tries to perform a match against the
        # super star file ... which is trivial because the solution has
        # already been determined by superalign ... and simplematch merely
        # verifies it, writes the output in appropriate format for match to
        # read it, and then match is called with an input guess which
        # should be exactly right

        print "Run simplematch"
        #pdb.set_trace()
        logtxt = "hdr_xshift = %.3f, hdr_yshift = %.3f, orient_shift=%.2f" % (hdr_dx,hdr_dy,hdr_rot)
        self.logfile.write(logtxt)
        files = file + '.stars'
        filetrans = file + '.intrans'
        cmd = "simplematch %s %s %.1f %.1f %.3f %s" % (refMatchFile,files,hdr_dx,hdr_dy,hdr_rot,filetrans)
        _run_cmd('simplematch',cmd,self.GlobalBlock)

        self._process_simplematchout(filetrans)
        return 1
        
        # base = run_match + ' ' + files +' 1 2 3 '+refMatchFile +' 1 2 3 max_iter=10 scale=1.0 outfile='+outfile+\
        #        ' intrans='+filetrans+' recalc medtf '
        # self.logfile.write(base)
        # nobj=60
        # matchrad=tolmatch
        #                     
        # logtxt = "Running 1st match on "+os.path.split(file)[1]+" . . . trans guess: "
        # logtxt+= filetrans+'  nobj,mrad: '+str(nobj)+' '+str(matchrad)
        # self.logfile.write(logtxt)
        # cmd = base + 'medsigclip=2.7 nobj='+str(nobj)+' matchrad='+str(matchrad)        
        # if self._do_the_matching(cmd, file, matchf, default_scale_tol, maxfitrms) > -1:
        #     return 1
        # return -1

    # WZ Use header information for shifts
    def _use_wcs(self, refMatchFile, hdr_dx, hdr_dy, hdr_rot, file, tolmatch, default_scale_tol, maxfitrms, outfile, matchf):
        # logtxt = "hdr_xshift = %.3f, hdr_yshift = %.3f, orient_shift=%.2f" % (hdr_dx,hdr_dy,hdr_rot)
        #self.logfile.write(logtxt)
        #files = file + '.stars'
        #filetrans = file + '.intrans'
        dirname = string.split(refMatchFile,'align/')[0]+'Images/'
        refImage = dirname + string.split(string.split(refMatchFile,'align/')[1],'.ma')[0] + '.fits'
        self._process_wcs_out(refImage)
        return 1
        
        # base = run_match + ' ' + files +' 1 2 3 '+refMatchFile +' 1 2 3 max_iter=10 scale=1.0 outfile='+outfile+\
        #        ' intrans='+filetrans+' recalc medtf '
        # self.logfile.write(base)
        # nobj=60
        # matchrad=tolmatch
        #                     
        # logtxt = "Running 1st match on "+os.path.split(file)[1]+" . . . trans guess: "
        # logtxt+= filetrans+'  nobj,mrad: '+str(nobj)+' '+str(matchrad)
        # self.logfile.write(logtxt)
        # cmd = base + 'medsigclip=2.7 nobj='+str(nobj)+' matchrad='+str(matchrad)        
        # if self._do_the_matching(cmd, file, matchf, default_scale_tol, maxfitrms) > -1:
        #     return 1
        # return -1

    def DetermineShift(self,AlignParam,refim,refMatchFile,matchf,defShiftDict):
#      run_match = 'match'
      outfile   = os.path.join(self.GlobalBlock.obsAlign,'match')
      hdr_dx = self.MatchDict['hdr_xshift']
      hdr_dy = self.MatchDict['hdr_yshift']
      hdr_rot = refim.refOrient - self.MatchDict['hdr_orient']

      # pdb.set_trace()
      if self.PathName in defShiftDict.keys():
        # bail_out method here just ensures proper bookkeeping
        self._bail_out(matchf, quiet=1)
        self._getDefShifts(defShiftDict[self.PathName])
      else:
        success_match = 0
        if AlignParam.superalign != None:
          if self._use_simplematch_for_sol(refMatchFile, hdr_dx,\
                                           hdr_dy, hdr_rot, self.matchFile,\
                                           AlignParam.tolmatch,\
                                           AlignParam.default_scale_tol,\
                                           AlignParam.maxfitrms,\
                                           outfile,\
                                           matchf) > 0:
            self.logfile.write("Matching with simplematch successful.")
            success_match = 1
          else:
            errtxt = "Simplematch failure."
            self.logfile.write(errtxt)
        else: #WZ
          print "Use header WCS"
          #pdb.set_trace()
          self._use_wcs(refMatchFile, hdr_dx,hdr_dy, hdr_rot,\
                          self.matchFile,AlignParam.tolmatch,\
                          AlignParam.default_scale_tol,\
                          AlignParam.maxfitrms,outfile,matchf)
          logtxt = "Skip Superalign. Use header info"
          #GlobalBlock.logfile.write(errtxt)
          self.logfile.write(logtxt)
          # raise Exception, errtxt #WZ

#	I = 0
#        for MatchParam in AlignParam.MatchParamList: 
#          if not success_match:
#            if ((not MatchParam.identity) or (abs(hdr_rot) < 0.6)):
#              cmd = MatchParam.matchcmd(run_match,self.PathName,refMatchFile,hdr_dx,hdr_dy,outfile)
#              I = I + 1
      #         logtxt = "Running match time %i on %s" % (I, self.PathName)
      #         if self._do_the_matching(cmd, self.matchFile, matchf, MatchParam.default_scale_tol, MatchParam.maxfitrms) < 0:
      #           self.logfile.write(cmd)
      #           self.logfile.write(logtxt)
      #           self.logfile.write("Matching with guesses successful.")
      # del hdr_dx, hdr_dy, hdr_rot
      # Matching done for this file

    def setFinalShifts(self,refim,notrot):
	""" Returns a tuple of final shifts uncertainties for an image.
        The code uses the fitted shifts except in emergencies, and
        then it sets MatchDict[xpix_shift],MatchDict[ypix_shift], etc
        and returns the shifts it decides on.  It logs a warning if the
        final shifts differ significantly from the commanded ones.
        """
        # as of 6/Dec/01, shifts now in output pixels; convert for arcsec
        # xpix_shift,ypix_shift = matutil.cubinvert(xshift_asec, yshift_asec,
        #                                  self.xcoeff, self.ycoeff, verb=0)
        # however, "pixels at reference pix" no longer makes sense in the case
        # of multiple chips, which will have different input scales
        # therefore, from henceforth (30/Nov/01), these pixel shifts are
        # in units of *output pixels* !

        # the 2.5*sig tolerance is, unfortunamente, a bit arbitrary
        print "Check shifts"
        min_angle = max(_angle_tolerance_, 2.5 * self.MatchDict['fit_angle_err'])
        
        # basically, if the angle is ~0, want to go with the median shifts
        # because there should be no scale change or distortion in the
        # corrected positions used in determining the shifts.

        if self.MatchDict['Nmatch'] == 0:
            # Hopeless. Use header WCS shifts
            xpix_shift  = self.MatchDict['hdr_xshift']/self.outscale
            ypix_shift  = self.MatchDict['hdr_yshift']/self.outscale
            dxpix_shift = max(0.5/self.outscale,0.1*abs(xpix_shift))
            dypix_shift = max(0.5/self.outscale,0.1*abs(ypix_shift))
            #there's always a sci[0], and this is just header info, so
            self.MatchDict['angle']  = refim.refOrient - self.MatchDict['hdr_orient']
            self.MatchDict['angle_err'] = 9.999
	    warnline =  'WARNING: using WCS information for shifts!'
	    self.GlobalBlock.errorList.append((self.GlobalBlock.modName,warnline))
	    self.logfile.write(warnline)
            print warnline

        else:
            # use fitted shifts (or supplied default shifts if Nmed=0 and angle small)
            xpix_shift  = self.MatchDict['fit_xshift']/self.outscale
            ypix_shift  = self.MatchDict['fit_yshift']/self.outscale
            dxpix_shift = self.MatchDict['fit_xshift_err']/self.outscale
            dypix_shift = self.MatchDict['fit_yshift_err']/self.outscale
            # and leave angle nonzero
            self.MatchDict['angle']     = self.MatchDict['fit_angle'] # WZ Nov # 0.
            self.MatchDict['angle_err'] = self.MatchDict['fit_angle_err']

        # Add all the final shifts to MatchDict, then return the line
        self.MatchDict['xpix_shift']     = xpix_shift
        self.MatchDict['ypix_shift']     = ypix_shift
        self.MatchDict['xpix_shift_err'] = dxpix_shift
        self.MatchDict['ypix_shift_err'] = dypix_shift
        self.MatchDict['xarc_shift']     = xpix_shift * self.outscale 
        self.MatchDict['yarc_shift']     = ypix_shift * self.outscale 
        
        # finally, compare msd to commanded shifts:
        xptg_err = round(1e5 * (xpix_shift - self.MatchDict['hdr_xshift']/self.outscale))/1e5
        yptg_err = round(1e5 * (ypix_shift - self.MatchDict['hdr_yshift']/self.outscale))/1e5
        #pdb.set_trace()
        if max(abs(xptg_err),abs(yptg_err)) > max(0.1*max(xpix_shift,ypix_shift),10.):
            warnline = 'Warning: found x,y pointing errors of:  '+str(xptg_err)+' '+str(yptg_err)+' pix.'
            self.logfile.write(warnline)
            self.GlobalBlock.errorList.append((self.GlobalBlock.modName,warnline))


        if (self.PathName == refim.PathName): #WZ
          yes = 1
        else:
          yes = 0                      
	shifts  = (self.PathName, xpix_shift, dxpix_shift, ypix_shift, dypix_shift,\
                   self.MatchDict['angle'],self.MatchDict['angle_err'], yes) #WZ
	return shifts

    def subtractBack(self,deAmp,ForceSub=0,omitSimple=0,omitMEF=0,aveExt=0):
      if not omitSimple:
        # we just loop over and modify <simple_SCI.fits> images 
        for ii in range(self.NumSci):
          if aveExt:
            sky = self.MatchDict['skyval']
            skyrms = self.MatchDict['skyrms']
          else:
            sky = self.preMatchDict[ii]['skyval']
            skyrms = self.preMatchDict[ii]['skyrms']
            # skyrms = max(math.sqrt(sky),skyrms) # Now we get negative sky values!
            skyrms = max(math.sqrt(abs(sky)),skyrms)                    

          simplefilename = self.preMatchDict[ii]['simplefits']

          subtractextension(simplefilename,0,sky,skyrms,deAmp,'ALIGNSKY',ForceSub,self.GlobalBlock)

      if not omitMEF:
        meffilename = self.MatchDict['fullPath']
        meffits = pyfits.open(meffilename,'update')
        for ext in range(1,len(meffits)):
          # see if it's a SCI exten & hasn't had its sky subtracted yet
          xname  = meffits[ext].header['EXTNAME'].upper()
          if xname != 'SCI':
            # it's not a SCI extension - leave as unaltered
            del xname
            continue

          else:
            # need to subtract sky -- first find value to use
            if aveExt:
              sky    = self.MatchDict['skyval']
              skyrms = self.MatchDict['skyrms']
            else:
              for ii in range(self.NumSci):
                if self.sci[ii]['extnum'] == ext:
                  # this is the right extension, so set the sky
                  sky    = self.preMatchDict[ii]['skyval']
                  skyrms = self.preMatchDict[ii]['skyrms']
                  break
            ## skyrms = max(math.sqrt(sky),skyrms)
            skyrms = max(math.sqrt(abs(sky)),skyrms)
                        
            subtractextension(meffilename,ext,sky,skyrms,deAmp,'ALIGNSKY',ForceSub,self.GlobalBlock)

          del xname
          # done processing this extension
                
        # done looping over extensions -- close file
        meffits.close()
        del meffilename,meffits


    def getHeaderShift(self,refim):
         ff = pyfits.open(self.MatchDict['fullPath'])
         for ii in range(self.NumSci):
           extnum = self.sci[ii]['extnum']
	   if self.simpleimage == 0:
             ra  = ff[extnum].header['CRVAL1']
             dec = ff[extnum].header['CRVAL2']
	   else:
             ra = self.wcs.wcs['CRVAL1']
             dec = self.wcs.wcs['CRVAL2']
           
           dx_pix = 0
           dy_pix = 0
           for jj in range(refim.NumSci):
             #pdb.set_trace()
             x1,y1 = refim.refwcs[jj].rd2xy((ra,dec),hour=0)
                
             # NOTE! NOTE! NOTE! these shifts should now (19Nov2002) be
             # exact shifts in the output (rectified) frame.
             dx_pix = dx_pix + (refim.x0[jj] - x1)*refim.outscale
             dy_pix = dy_pix + (refim.y0[jj] - y1)*refim.outscale

             # NOTE! NOTE! NOTE! All shifts in the MatchDict's are numbers
             # to be *added* to the image pixel positions in order to make
             # an object's x,y coordinates the same as in the ref image;
             # thus, they are negative of dx,dy:

           dx_pix /= refim.NumSci
           dy_pix /= refim.NumSci
           self.preMatchDict[ii]['hdr_xshift'] = -dx_pix # WZ ???
           self.preMatchDict[ii]['hdr_yshift'] = -dy_pix
           if self.simpleimage == 0:
             self.preMatchDict[ii]['hdr_orient'] = ff[extnum].header['ORIENTAT'] + self.Camera.OrientAdj
           else:
             self.preMatchDict[ii]['hdr_orient'] = self.angle
         ff.close()
         del ff,dx_pix,dy_pix

    def FixWCS(self,refim):
        M_PI = 3.1415926535
        dx = (self.MatchDict['xpix_shift']*self.outscale - \
              self.MatchDict['hdr_xshift']) / refim.outscale
        dy = (self.MatchDict['ypix_shift']*self.outscale - \
              self.MatchDict['hdr_yshift']) / refim.outscale
        hdr_rot = refim.refOrient - self.MatchDict['hdr_orient']
        dtheta = hdr_rot - self.MatchDict['angle']
        if dtheta < -180:
          dtheta += 360
        if dtheta > 180:
          dtheta -= 360
	S = 'Fixing %s -- offsetting by %.2f, %.2f -- rotating by %.2f' % \
               (self.PathName,dx,dy,dtheta)
	self.logfile.write(S)
	print S

        ff = pyfits.open(self.MatchDict['fullPath'],'update')
	if ff[0].header.get('PA_V3'):
  	  ff[0].header['PA_V3'] += dtheta
	  pa_v3 = ff[0].header['PA_V3']
	if ff[0].header.get('ORIENTAT'):
  	  ff[0].header['ORIENTAT'] += dtheta

	tra = 0
	tdec = 0
	for ii in range(self.NumSci):
          extnum = self.sci[ii]['extnum']

          tra += ff[extnum].header['CRVAL1']
          tdec += ff[extnum].header['CRVAL2']

        tra /= self.NumSci
        tdec /= self.NumSci

        for ii in range(self.NumSci):
          extnum = self.sci[ii]['extnum']

	  if ff[extnum].header.get('ORIENTAT'):
  	    ff[extnum].header['ORIENTAT'] += dtheta
	  if ff[extnum].header.get('PA_V3'):
  	    ff[extnum].header['PA_V3'] += dtheta
            pa_v3 = ff[extnum].header['PA_V3']

	  cosine = math.cos(M_PI*dtheta/180.0)
	  sine = math.sin(M_PI*dtheta/180.0)

          cd11 = ff[extnum].header['CD1_1']*cosine + \
            ff[extnum].header['CD2_1']*sine
          cd21 = ff[extnum].header['CD1_1']*(-sine) + \
            ff[extnum].header['CD2_1']*cosine
          cd12 = ff[extnum].header['CD1_2']*cosine + \
            ff[extnum].header['CD2_2']*sine
          cd22 = ff[extnum].header['CD1_2']*(-sine) + \
            ff[extnum].header['CD2_2']*cosine
	  ff[extnum].header['CD1_1'] = cd11
	  ff[extnum].header['CD1_2'] = cd12
	  ff[extnum].header['CD2_1'] = cd21
	  ff[extnum].header['CD2_2'] = cd22

          tdra = 0
          tddec = 0         
          for jj in range(refim.NumSci):            
             ra1,dec1 = refim.refwcs[jj].xy2rd((dx,dy))
	     ra0,dec0 = refim.refwcs[jj].xy2rd((0,0))
                
             tdra += (ra1-ra0)
             tddec += (dec1-dec0)

	  print dx,dy,tdra,tddec
          tdra /= refim.NumSci
          tddec /= refim.NumSci

	  dra = ff[extnum].header['CRVAL1'] - tra
          ddec = ff[extnum].header['CRVAL2'] - tdec

	  cosine = math.cos(180*dtheta/M_PI) - 1

	  tdra += dra*cosine + ddec * sine
          tddec += -ddec*sine + ddec * cosine

          ff[extnum].header['CRVAL1'] += tdra
          ff[extnum].header['CRVAL2'] += tddec

        ff.close()

    def SetUpRefAngleDict(self):
        ### save these V2,V3 positions in refAngleDict
        # self.refAngleDict['V2'] = self.V2center  <-- wrong because drizzle preserves
        # self.refAngleDict['V3'] = self.V3center      1st extenstions refpix 
        self.refAngleDict = {}
        self.refAngleDict['V2'] = self.sci[0]['refdict']['V2REF']
        self.refAngleDict['V3'] = self.sci[0]['refdict']['V3REF']
        
        self.logfile.write('Taking detector center at '+\
                           str(self.V2center)+','+\
                           str(self.V3center)+' arcsec.')

	rf = pyfits.open(self.PathName)
        # save PA_V3 of refimage for later WCS construction

        try:
          self.refAngleDict['pa_v3'] = rf[0].header['PA_V3']
        except Exception,err:
          self.GlobalBlock.errorList.append((self.GlobalBlock.modName,str(err)))
          self.logfile.write("PA_V3 not found in refimage primary header; "+\
                             "got following error: \n"+str(err)+"\n")
          self.logfile.write("Looking for PA_V3 in extension header.")
          try:
            self.refAngleDict['pa_v3'] = rf[1].header['PA_V3']
          except Exception,err:
            self.GlobalBlock.errorList.append((self.GlobalBlock.modName,str(err)))
            self.logfile.write("No good: PA_V3 not in refimage ext 1 header; "+\
                               "got following error: \n"+str(err)+"\n")
            raise Exception,err
        rf.close()            
        del rf
        
    def SortOutFrame(self):
        """ Get expected shifts of all images relative to the reference image
        from the image headers. This just uses the raw WCS, not the distortion map,
        so beware."""

        if self.ref == 0:
          raise ValueError,"SortOutFrame should not be called since ref = 0."

        if self.simpleimage == 0:
          self.SetUpRefAngleDict()
          
	rf = pyfits.open(self.PathName)
        self.x0 = {}
        self.y0 = {}
        self.refwcs = {}
        self.refOrient = 0.0
        self.refotherkeys = {}
        self.reflogfile_add = []

        for ii in range(self.NumSci):
           if self.simpleimage == 0:
             extnum = self.sci[ii]['extnum']
             refwcs = wcsclass.BasicWCS(rf[extnum].header)
             self.refOrient += refwcs.wcs['ORIENTAT'] + self.Camera.OrientAdj
             self.x0[ii] = refwcs.wcs['CRPIX1']
             self.y0[ii] = refwcs.wcs['CRPIX2']
             ra0,dec0 = (refwcs.wcs['CRVAL1'],refwcs.wcs['CRVAL2'])

             _V2 = self.sci[ii]['refdict']['V2REF']
             _V3 = self.sci[ii]['refdict']['V3REF']
             PA = matutil.targroll(self.refAngleDict['pa_v3'], dec0, _V2, _V3, updown=self.updown)
	     PA_FINAL = PA
             if self.sci[ii]['refdict']['THETA']:
               PA += self.sci[ii]['refdict']['THETA']

             cdmat = matutil.makeCDdict(PA, self.outscale)
             refwcs.wcs['CD1_1'] = cdmat['CD1_1'] 
             refwcs.wcs['CD1_2'] = cdmat['CD1_2']
             refwcs.wcs['CD2_1'] = cdmat['CD2_1']
             refwcs.wcs['CD2_2'] = cdmat['CD2_2']
             self.refwcs[ii] = refwcs
             if ii == 0:
               self.refotherkeys['PA_V3'] = self.refAngleDict['pa_v3']
               self.refotherkeys['PA_FINAL'] = PA
               S = "  matutil.targroll("+str(_V3)+", "+str(dec0)+\
                           ", "+str(_V2)+", "+str(_V3)+\
                           ", updown="+str(self.updown)+")"
	       self.reflogfile_add.append(S)

               # find the actual orientation and make the CDmatrix
               self.reflogfile_add.append('    target roll angle = '+str(PA_FINAL))
             del PA,cdmat,refwcs
           else:
             self.x0[0] = self.wcs.wcs['CRPIX1']
             self.y0[0] = self.wcs.wcs['CRPIX2']
             self.refwcs[0] = self.wcs
             self.refOrient = self.angle
        
        self.refOrient /= self.NumSci
        self.refcdmat = {}
        self.refcdmat['CD1_1'] = self.refwcs[0].wcs['CD1_1']
        self.refcdmat['CD2_1'] = self.refwcs[0].wcs['CD2_1']
        self.refcdmat['CD1_2'] = self.refwcs[0].wcs['CD1_2']
        self.refcdmat['CD2_2'] = self.refwcs[0].wcs['CD2_2']
        if self.simpleimage == 1:
            self.refcdmat['EXTREF'] = 1
        else:
            self.refcdmat['EXTREF'] = 0
        rf.close()            
	del rf

class MatchParam:
    def __init__(self,identity=None,minscale=1,maxscale=1,nobj=60,trirad=None,medsigclip=2.7,default_scale_tol=0.003,matchrad=6,maxfitrms=2.0):
      self.identity = identity
      self.minscale = minscale
      self.maxscale = maxscale
      self.nobj = nobj
      self.trirad = trirad
      self.medsigclip = medsigclip
      self.default_scale_tol = default_scale_tol
      self.matchrad = matchrad
      self.maxfitrms = maxfitrms

    # def matchcmd(self,run_match,file,refMatchFile,xsh,ysh,outfile):
    #   if (self.minscale > 0.9999) and (self.maxscale < 1.0001):
    #     scalecmd = 'scale=1.00'
    #   else:
    #     scalecmd = 'min_scale=%.2f max_scale=%.2f' % (self.minscale,self.maxscale)
    #   if (self.identity):
    #     identitycmd = 'identity xsh='+str(xsh)+' ysh='+str(ysh)
    #   else:
    #     identitycmd = 'trirad=%.4f' % (self.trirad)
    # 
    #   cmd  = run_match + ' ' + file +' 1 2 3 '+ refMatchFile + \
    #          ' 1 2 3 max_iter=10 ' + scalecmd + \
    #          ' nobj=%i matchrad=%.1f' % (self.nobj,self.matchrad) + \
    #          ' %s recalc medtf ' % (identitycmd) + \
    #          'medsigclip=%.2f ' % (self.medsigclip) + \
    #          'outfile=' + outfile
    #   return cmd

class AlignParam:
    def __init__(self,superalign,mdz_align,tolmatch,default_scale_tol,maxfitrms):
      # self.MatchParamList = []
      self.superalign = superalign
      self.mdz_align = mdz_align
      self.tolmatch = tolmatch
      self.default_scale_tol = default_scale_tol
      self.maxfitrms = maxfitrms

    # def AddMatchParam(self,MatchParam):
    #   self.MatchParamList.append(MatchParam) 

class StandAlignParam(AlignParam):
    def __init__(self,superalign,mdz_align):
      AlignParam.__init__(self,superalign,mdz_align,tolmatch=6/20.0,default_scale_tol=0.003,maxfitrms=2.0/20.0)
      # self.AddMatchParam(MatchParam(identity=1,minscale=1,maxscale=1,nobj=60,medsigclip=2.7,default_scale_tol=0.003,matchrad=6/20.0))
      # self.AddMatchParam(MatchParam(identity=1,minscale=1,maxscale=1,nobj=1000,medsigclip=2.5,default_scale_tol=0.0042,matchrad=5/20.0))
      # self.AddMatchParam(MatchParam(identity=1,minscale=1,maxscale=1,nobj=1000,medsigclip=2.5,default_scale_tol=0.0042,matchrad=12/20.0))
      # self.AddMatchParam(MatchParam(identity=1,minscale=1,maxscale=1,nobj=1000,medsigclip=2.3,default_scale_tol=0.00672,matchrad=18/20.0))
      # self.AddMatchParam(MatchParam(identity=1,minscale=1,maxscale=1,nobj=1000,medsigclip=2.0,default_scale_tol=0.00672,matchrad=30/20.0))
      # self.AddMatchParam(MatchParam(identity=1,minscale=1,maxscale=1,nobj=60,trirad=0.0005,medsigclip=2.7,default_scale_tol=0.003,matchrad=6/20.0))
      # self.AddMatchParam(MatchParam(identity=1,minscale=0.99,maxscale=1.01,nobj=120,trirad=0.0005,medsigclip=2.7,default_scale_tol=0.003,matchrad=5/20.0))
      # self.AddMatchParam(MatchParam(identity=1,minscale=0.98,maxscale=1.02,nobj=200,trirad=0.001,medsigclip=2.7,default_scale_tol=0.0042,matchrad=5/20.0))


class Align:
    """ Align is the module/class which does all the work for aligning 
        a huge set of images from different instruments / filters / etc.

      Parameters:
        imfiles consists of a list of FitsFile instances
        CameraList specifies a list of different Camera objects.
        Groups contains a list of different groupings
        CurExtImage specifies an external image off of which to do
          alignments.  It should be of class FitsFile.
        GlobalBlock contains many of global parameters used for logging
          operations and specifying the overall directory structure.

    Basic usage:
       al = align.alignImage(obs)      # all the setup
       al.makeMatchCats()              # all sextractor runs; cat generation
       al.match()                      # all matching/merging/output'ing    
       al.subtractSkies()              # subtract skyval's from simple & mef images 
     """
    def __init__(self,imfiles,allShifts,CameraList,GlobalBlock,SkyFile=0,idcTab=None,CurExtImage=None,Groups=None):
      self.imfiles = imfiles
      self.allShifts = allShifts
      self.CameraList = CameraList
      self.idcTab = idcTab
      self.GlobalBlock = GlobalBlock
      self.logfile = GlobalBlock.logfile
      self.Groups = Groups
      self.refImage = None
      self.AsnList = []
      for file in self.imfiles:
        if not (file.Asn in self.AsnList):
          self.AsnList.append(file.Asn)
      
      if CurExtImage:       
        for file in self.imfiles:
          if CurExtImage.PathName == file.PathName:
            GlobalBlock.errorList.append((GlobalBlock.modName,"External Reference Image cannot be in Internal Image List!"))
            raise Exception,"External Reference Image cannot be in Internal Image List!"

        CurExtImage.CheckExtRefFlags()
        self.imfiles.append(CurExtImage)
        self.ExtImage = CurExtImage
        self.refImage = CurExtImage

    def removeUnpicklable(self):
      for fitsfile in self.imfiles:
        fitsfile.removeUnpicklable()
        
    def setupImages(self):
      for fitsfile in self.imfiles:
        if fitsfile.extref == 0:
          fitsfile.EvaluateFilt(self.CameraList)          
	else:
          fitsfile.simpleimage = 1
          fitsfile.DetermineFirstImage()
          print fitsfile.PathName
        if fitsfile.simpleimage == 0:
          fitsfile.DetermineExtStruct()
          #pdb.set_trace() #
          f = pyfits.open(fitsfile.PathName)
          self.det=f[0].header.get('DETECTOR')
          f.close()
          # if (string.find(det,'WFC')>-1: # WZ
          fitsfile.setupIDCinfo(self.idcTab)
	else:
          fitsfile.setupWCS()

    def recenterImages(self):
      for fitsfile in self.imfiles:
        if fitsfile.simpleimage == 0:
          fitsfile.recenter()

    def makeMatchCats(self,SkyFile,irafskytoo):
      skyLines = skyLinesClass(SkyFile,self.GlobalBlock)

      #pdb.set_trace()
      for fitsfile in self.imfiles:
        fitsfile.makeMatchCat(skyLines,irafskytoo)
      del skyLines

    def match(self,AlignParam,notrot,ref,intRef = None,retry=0):
	"""  Match sources in the .matchin files using match. If no refFile
        has yet been defined, it calls _chooseRefimage to get one.  
        You can specify a default shift file in the case of total match failures;
        if not an absolute path, it assumes the file is in the align directory.
        Optional parameters:
           intRef  -  specify a reference image internal to the dataset.
           retry   -  running match again without reloading or remaking MatchCats,etc.
                     Note: this only works right if not using an externally supplied
                     reference image.
        """
	self.GlobalBlock.outputList['Match.shifts'] = []
        if retry:
            self.logfile.write('\n\n#---------------> *** redoing match *** <----------------#')
        else:
            self.logfile.write('\n\n#-------------------> Begin match() <-------------------#')

        if intRef:
            errtxt=''
            if os.path.isabs(intRef):
                errtxt += "Error: internal ref image spec cannot be full path. "
            if self.refImage:
                errtxt += "Error: cannot set both external and internal ref images!"
            if errtxt:
                self.logfile.write(errtxt)
                raise Exception,errtxt
            del errtxt
        
        self._mergeMatchDicts()
        if not self.refImage:
          self.refImage = self._chooseRefimage(ref,intRef=intRef)
          if ref: #WZ
            self.refImage.PathName = ref
          self.refImage.ChangeToRef()
        self.refMatchFile = self.refImage.matchFile

        defShiftDict = ReturndefShiftDict(self.GlobalBlock,self.refImage)
        self._getHeaderShifts(self.refImage)
        # remake MatchDicts to include header shifts
        self._delMatchDicts()
        self._mergeMatchDicts()

        # pdb.set_trace()
        if AlignParam.mdz_align != None:
          self._MdzShift(AlignParam,self.refImage,self.refMatchFile,defShiftDict) #WZ                                                                                  
        else:
          if AlignParam.superalign != None:
            self._superalignprep()

          matchf = open(os.path.join(self.GlobalBlock.obsAlign,'Match.shifts'),'w')
          twf = open(os.path.join(self.GlobalBlock.obsAlign,'alignTweak.txt'),'w') #WZ Nov
          for fitsfile in self.imfiles:
            if (not self.refImage.extref) or (self.refImage.PathName != fitsfile.PathName):
              fitsfile.DetermineShift(AlignParam,self.refImage,self.refMatchFile,matchf,defShiftDict)
            ff = pyfits.open(fitsfile.PathName)
            pa= ff[0].header.get('PA_V3')
            print >> twf, fitsfile.PathName[0:9], '0.000 0.000 0.000', int(round(pa)), fitsfile.PathName[0:7]
            ff.close()
          matchf.close()
          twf.close()
          # finally, generate individual shift files for each association           
          # pdb.set_trace()                                                         
          self._finalizeShifts(notrot)

        self.logfile.write('\n#-------------------> Finished match() <-------------------#\n')

    def FixWCS(self):
      for fitsfile in self.imfiles:
        if (not self.refImage.extref) or (self.refImage.PathName != fitsfile.PathName):
          fitsfile.FixWCS(self.refImage)

    def _chooseRefimage(self,ref,intRef=None):
      """Returns a reference image file name, from the first image in the asn file #WZ
         Also sets self.refMatchFile for each ii
      """

      Maxgood=10
      # if we've been sent a ref image, use it after checking to see if it's valid
      if intRef:
        valid_intRef=0
        for image in self.imfiles:
          if image.PathName == intRef:
            valid_intRef=1
            break
        if not valid_intRef:
          errtxt = "Error: specified internal ref image not found!"
          self.GlobalBlock.logfile.write(errtxt)
          raise Exception,errtxt
        refImage = intRef
        logtxt = 'Using specified internal reference image: '+refImage.PathName+ \
                 ' with match input file(s): '
      else:
        for image in self.imfiles:
          Ngood=image.MatchDict['Ngood']
          if Ngood > Maxgood:
            newRefImage = image
            Maxgood = Ngood
        for image in self.imfiles:
          if (image.PathName == ref):
            newRefImage = image
        if not Maxgood:
          self.GlobalBlock.errorList.append((self.GlobalBlock.modName,"Cannot determine a reference image!"))
          raise ValueError,"Cannot determine a reference image!"

        refImage = newRefImage
        logtxt = 'Reference image chosen as: '+refImage.PathName+' (Ngood='+ \
                  str(Maxgood)+') with match input file: '
        
      # we have self.refImage, now get name of matchfile
      logtxt += ' '+refImage.matchFile
      self.GlobalBlock.logfile.write(logtxt)
      return refImage

    def _delMatchDicts(self):
      for fitsfile in self.imfiles:
        fitsfile.delMatchDicts()      

    def _mergeMatchDicts(self):
      for fitsfile in self.imfiles:
        fitsfile.mergeMatchDicts()

    def _getHeaderShifts(self,refim):
      refim.SortOutFrame()
      for fitsfile in self.imfiles:
        fitsfile.getHeaderShift(refim)

    def _superalignprep(self):
      # writes out 'supermatch.in' contains a list of all the input
      # star files plus initial guess for the superalign package
      supermatchfin = os.path.join(self.GlobalBlock.obsAlign,'supermatch.in')
      # supermatchfout is to contain the solutions (dx,dy,drot) for
      # all images
      supermatchfout = os.path.join(self.GlobalBlock.obsAlign,'supermatch.out')
      supermatchf = open(supermatchfin,'w')
      S = '%i %i\n' % (len(self.imfiles),self.refImage.extref)            
      supermatchf.write(S)

      # put together initial guess file
      for I in range(2):
        for fitsfile in self.imfiles:
          if ((I == 0) and (fitsfile.PathName == self.refImage.PathName)) or \
             ((I == 1) and (fitsfile.PathName != self.refImage.PathName)):
            hdr_dx = fitsfile.MatchDict['hdr_xshift']
            hdr_dy = fitsfile.MatchDict['hdr_yshift']
            hdr_rot = self.refImage.refOrient - fitsfile.MatchDict['hdr_orient']
            file = fitsfile.matchFile
            S = '%s %.3f %.3f %.3f\n' % (file,hdr_dx,hdr_dy,hdr_rot)
            supermatchf.write(S)
      supermatchf.close()
      print "Done"
          
      self.refMatchFile = self.refMatchFile + '.super'
      # self.refMatchFile is to contain the super star list from all
      # the images, after rejecting across the different images
          
      print "Run supermatch"
      #pdb.set_trace()
      superalign_cmd = "superalign %s %s %s" % \
            (supermatchfin,self.refMatchFile,supermatchfout)
      # all hard work is done by superalign / python code is just
      # a wrapper, so to speak
      _run_cmd('superalign',superalign_cmd,self.GlobalBlock)

    def _finalizeShifts(self,notrot):
        """ This function loops over each image of each association and calls
        setFinalShifts() to decide on the final best x,y shifts (median or fit)
        and rotation (only possible if fitted better than median) for each image.
        The results are then saved as yet more MatchDict entries and written out
        in the shifts_<asn>.txt files.
        """
        for Asn in self.AsnList:
          asnName = string.split(Asn,'_asn')[0]
          sftxt   = os.path.join(self.GlobalBlock.obsAlign,'shifts_'+asnName+'.txt')
  	  asnShiftFile = open(sftxt,'w')
          del sftxt
          sdict = {}
            
          # Write a shiftfile "shift.txt', break into associations later

          for fitsfile in self.imfiles:
            if fitsfile.Asn == Asn:
              if (not self.refImage.extref) or (self.refImage.PathName != fitsfile.PathName):
                #pdb.set_trace()
                shifts  = fitsfile.setFinalShifts(self.refImage,notrot)
                print "Shifts: ",shifts
                asnShiftFile.write('%s %9.4f %7.4f %9.4f %7.4f %12.6f %9.6f %d\n' %
                                   shifts)
                sdict[fitsfile.PathName] = (shifts[1],shifts[3])

          # remainder of this loop is the same
          self.allShifts.append((asnName,sdict))
  	  asnShiftFile.close()
	  del asnShiftFile,fitsfile
        del Asn
	return
                
    def subtractSkies(self, useMinSky, deAmp, ForceSub=0, omitSimple=0, omitMEF=0, aveExt=0):
        """ Subtract the sky values of all the input FITS images
        (simplefits *and* mef, unless overwritten by optional parameters).
        Update headers to reflect subtraction using ALIGNSKY keyword.
        Optional input parameters:
           ForceSub   (def 0):  if ==1, subtract sky even if ALIGNSKY keyword is 
                                   already in the header. Normally it won't subtract 
                                   the sky if this keyword is already present.
           omitSimple (def 0):  if ==1, omit all processing of simple FITS images
           omitMEF    (def 0):  if ==1, omit all processing of simple MEF images
           aveExt     (def 0):  if ==1, subtract average sky from multiple image
                                   extensions; otherwise, subtract individual sky's.
           
        There's a high danger of confusion if you change any of these options
        from their default values.
        """
        self.logfile.write('\n#-------------------> Begin subtractSkies() <-------------------#\n'+
                           'params:  ForceSub='+str(ForceSub)+' omitSimple='+str(omitSimple)+\
                           ' omitMEF='+str(omitMEF)+' aveExt='+str(aveExt)+' self.useMinSky='+str(useMinSky))
        if useMinSky and (not aveExt):
            aveExt=1
            self.logfile.write(' setting aveExt=1 for consistency with useMinSky=1')

        for fitsfile in self.imfiles:
          if (not self.refImage.extref) or (self.refImage.PathName != fitsfile.PathName):
            fitsfile.subtractBack(deAmp,ForceSub=ForceSub,omitSimple=omitSimple,omitMEF=omitMEF,aveExt=aveExt)

        self.logfile.write('\n#-------------------> Finished subtractSkies() <-------------------#\n')
        return


    def _MdzShift(self,AlignParam,refim,refMatchFile,defShiftDict):
      # WZ Run xydrizzle to determine shift/rotation Sep 2012
      # if (self.refImage.Camera.Name== 'WFC'):
      # XX
      # pdb.set_trace()
      if str.lower(self.refImage.Camera.Name)== 'wfc3/ir' :
	 print "The WFC3/IR used, using the 'match_ir.inpar.  [align.py]_________By Xingxing'"
         pname = os.path.join(self.GlobalBlock.parDir,str('match_ir.inpar'))
      else :
         print "The "+str.lower(self.refImage.Camera.Name)+" used! [align.py]__________By Xingxing"
         pname = os.path.join(self.GlobalBlock.parDir,str('match_'+str.lower(self.refImage.Camera.Name)+'.inpar'))
      
      fp = open(pname,'r')
      for line in fp.readlines():
        if (string.find(line,'#')<0):
          tmp=line.split(',')
          dd = 0;
      fp.close()
      os.chdir(self.GlobalBlock.obsFits)
      print os.getcwd()
      # pdb.set_trace()
      # f = open('tem.log','w')
      nfile=len(self.imfiles)
      for i in range(nfile):
      #  if (string.find(self.imfiles[i].PathName,self.refImage.PathName)<0):
        if (i==0):
          fname=self.imfiles[0].PathName
        else:
          tmp=fname
          fname=tmp+','+self.imfiles[i].PathName
      # fname=self.refImage.PathName+','+tmp1
      pdb.set_trace()
      # md = multidrizzle.Multidrizzle(self.refImage.Asn)
      md = multidrizzle.Multidrizzle(fname)
      # md.editpars()
      md.build() # >> f
      md.run()   # >> f
      # f.close()
      print "\nmultidrizzle done !\n"  #by xingxing

      dx = zeros(nfile,float)
      dy = zeros(nfile,float)
      dr = zeros(nfile,float)
      dn = zeros(nfile,int)
      ni = zeros(nfile,int)
      dx0=dy0=dr0=0.
      #dr_min=1e5
      # fr = open('multidrizzle.run','r')
      # sftxt   = '/Users/zheng/shifts_'+tmp2+'.txt'
      # cmd       = 'cp '+sftxt+' .'
      # subproc   = popen2.Popen4(cmd, 1)
      # errs      = subproc.fromchild.readlines()     #
      # pdb.set_trace() # set pars 
      fr = open('align.log','r')
      tmp=string.split(self.GlobalBlock.inputRoot,'ingest/')
      tmp1 = string.split(tmp[1],'/')
      tmp2 = tmp1[0]+'_'+string.upper(tmp1[1])
      sftxt   = os.path.join(self.GlobalBlock.obsAlign,'shifts_'+tmp2+'.txt')
      # sftxt   = os.path.join(self.GlobalBlock.obsAlign,'shifts.txt')   #___________by xingxing  combDither LINE 283
      fw = open(sftxt,'w')
      fmt=str('%s %9.4f  %7.4f %9.4f %7.4f %11.6f %8.5f %1d') 
      for line in fr.readlines():
        tmp=line.split(' ')
        if (len(tmp)> 3):
          fname=tmp[0]
          # fname=str.split(tmp[1],'[sci,')[0]
          for i in range(nfile):
            # pdb.set_trace() # set pars 
            if (string.find(self.imfiles[i].PathName,fname) > -1):
              dn[i] = dn[i] + 1
              # v= float(string.split(tmp[5],'=')[1])
              vx = float(tmp[1])
              if (dn[i]>1):
                self.imfiles[i].MatchDict['xpix_shift_err']= abs((vx-dx[i])/2.0)
              else:
                self.imfiles[i].MatchDict['xpix_shift_err']= 0.00067
              dx[i] = dx[i] + vx
              vy = float(tmp[2])
              if (dn[i]>1):
                self.imfiles[i].MatchDict['ypix_shift_err']= abs((vy-dy[i])/2.0)
              else:
                self.imfiles[i].MatchDict['ypix_shift_err']= 0.00067
              dy[i] = dy[i] + vy
              vr = float(tmp[3])
              # if (abs(vr) < abs(dr_min)):
              #  pdb.set_trace()
              #  dr_min=vr
              #  self.refImage.MatchDict=self.imfiles[i].MatchDict
              if (dn[i]>1):
                self.imfiles[i].MatchDict['angle_err']= abs((vr-dr[i])/2.0)
              else:
                self.imfiles[i].MatchDict['angle_err']= 0.00002
              dr[i] = dr[i] + vr
	else:
	  print "There are not enough elements in align.log!!"   #by xingxing
      for i in range(nfile):
        if (self.imfiles[i].PathName == self.refImage.PathName):
          # pdb.set_trace()
          ni[i]=1
          #dx0 = dx0 + vx
          #dy0 = dy0 + vy
          #dr0 = dr0 + vr
        # dx[i] = (dx[i] - dx0)/ dn[i]
        # dy[i] = (dy[i] - dy0)/ dn[i]
        # dr[i] = (dr[i] - dr0)/ dn[i]
        dx[i] = dx[i] / dn[i]
        dy[i] = dy[i] / dn[i]
        dr[i] = dr[i] / dn[i]
        self.imfiles[i].MatchDict['Nmatch']=10
        self.imfiles[i].MatchDict['fit_xshift']=dx[i]
        self.imfiles[i].MatchDict['fit_yshift']=dy[i]
        self.imfiles[i].MatchDict['angle']=-dr[i]
        self.imfiles[i].MatchDict['xpix_shift']=dx[i]
        self.imfiles[i].MatchDict['ypix_shift']=dy[i]
        self.imfiles[i].MatchDict['xarc_shift']=dx[i]*self.imfiles[i].outscale
        self.imfiles[i].MatchDict['yarc_shift']=dy[i]*self.imfiles[i].outscale
        print >> fw, fmt % (self.imfiles[i].PathName,dx[i],\
                            self.imfiles[i].MatchDict['xpix_shift_err'],dy[i],\
                            self.imfiles[i].MatchDict['ypix_shift_err'],dr[i],\
                            self.imfiles[i].MatchDict['angle_err'],ni[i])
        # self.imfiles[i].MatchDict['angle']=dr[i]
      fr.close()
      fw.close()
      print "parameters are writen into "+sftxt   # by xingxing
      # pdb.set_trace()
      # get the catalog into a list
      del sftxt,fname,dx,dy,dr,dn,ni,pname,fmt,nfile    #  attentions here!

    def _MdzShift(self,AlignParam,refim,refMatchFile,defShiftDict):
      # WZ Run multidrizzle to determine shift/rotation Sep 2012                    
      # if (self.refImage.Camera.Name== 'WFC'):                                     
      pname = os.path.join(self.GlobalBlock.parDir,str('match_'+str.lower(self.refImage.Camera.Name)+'.inpar'))
      fp = open(pname,'r')
      for line in fp.readlines():
        if (string.find(line,'#')<0):
          tmp=line.split(',')
          dd = 0;
      fp.close()
      os.chdir(self.GlobalBlock.obsFits)
      print os.getcwd()
      # pdb.set_trace()                                                             
      # f = open('tem.log','w')                                                     
      nfile=len(self.imfiles)
      for i in range(nfile):
      #  if (string.find(self.imfiles[i].PathName,self.refImage.PathName)<0):       
        shutil.copy(self.imfiles[i].PathName,string.replace(self.imfiles[i].PathName,'flt','bck'))
        if (i==0):
          fname=self.imfiles[0].PathName
        else:
          tmp=fname
          fname=tmp+','+self.imfiles[i].PathName
      # fname=self.refImage.PathName+','+tmp1                                       
      # pdb.set_trace()                                                             
      md = multidrizzle.Multidrizzle(fname)
      md.build()
      md.run()
      for i in range(nfile):
        os.rename(string.replace(self.imfiles[i].PathName,'flt','bck'),self.imfiles[i].PathName)
      dx = zeros(nfile,float)
      dy = zeros(nfile,float)
      dr = zeros(nfile,float)
      dn = zeros(nfile,int)
      ni = zeros(nfile,int)
      dx0=dy0=dr0=0.
      #dr_min=1e5                                                                   
      # fr = open('multidrizzle.run','r')                                           
      # sftxt   = '/Users/zheng/shifts_'+tmp2+'.txt'                                
      # cmd       = 'cp '+sftxt+' .'                                                
      # subproc   = popen2.Popen4(cmd, 1)                                           
      # errs      = subproc.fromchild.readlines()     #                             
      # pdb.set_trace() # set pars                                                  
      asnName = string.split(self.AsnList[0],'_asn')[0]
      sftxt   = os.path.join(self.GlobalBlock.obsAlign,'shifts_'+asnName+'.txt')
      altxt   = os.path.join(self.GlobalBlock.inputRoot,'align.par')
      if os.path.exists(altxt):
        shutil.copy(altxt,sftxt)
        fname=tableio.get_str(sftxt,0)
        dx,ex,dy,ey,dr,er=tableio.get_data(sftxt,(1,2,3,4,5,6))
        for i in range(len(dx)):
          for j in range(nfile):
            # pdb.set_trace() # set pars                                            
            if  (string.find(self.imfiles[j].PathName,fname[i]) > -1):
              # pdb.set_trace()                                                     
              self.imfiles[j].MatchDict['Nmatch']=10
              self.imfiles[j].MatchDict['fit_xshift']=dx[i]
              self.imfiles[j].MatchDict['fit_yshift']=dy[i]
              self.imfiles[j].MatchDict['angle']=-dr[i]
              self.imfiles[j].MatchDict['angle_err']= er[i]
              self.imfiles[j].MatchDict['xpix_shift']=dx[i]
              self.imfiles[j].MatchDict['xpix_shift_err']=ex[i]
              self.imfiles[j].MatchDict['ypix_shift']=dy[i]
              self.imfiles[j].MatchDict['ypix_shift_err']=ey[i]
              self.imfiles[j].MatchDict['xarc_shift']=dx[i]*self.imfiles[i].outscale
              self.imfiles[j].MatchDict['yarc_shift']=dy[i]*self.imfiles[i].outscale
        del ex,ey,er
      else:
        fw = open(sftxt,'w')
        fmt=str('%s %9.4f  %7.4f %9.4f %7.4f %11.6f %8.5f %1d')
        # pdb.set_trace()                                                           
        fr = open('align.log','r')
        tmp=string.split(self.GlobalBlock.inputRoot,'ingest/')
        tmp1  = string.split(tmp[1],'/')
        tmp2  = tmp1[0]+'_'+string.upper(tmp1[1])
        for line in fr.readlines():
          tmp=line.split(' ')
          if (len(tmp)> 3):
            fname=tmp[0]
            # fname=str.split(tmp[1],'[sci,')[0]                                    
            for i in range(nfile):
              # pdb.set_trace() # set pars                                          
              if  (string.find(self.imfiles[i].PathName,fname) > -1):
                dn[i] = dn[i] + 1
                # v= float(string.split(tmp[5],'=')[1])                             
                vx  = float(tmp[1])
                if (dn[i]>1):
                  self.imfiles[i].MatchDict['xpix_shift_err']= abs((vx-dx[i])/2.0)
                else:
                  self.imfiles[i].MatchDict['xpix_shift_err']= 0.00067
                dx[i] = dx[i] + vx
                vy = float(tmp[2])
                if (dn[i]>1):
                  self.imfiles[i].MatchDict['ypix_shift_err']= abs((vy-dy[i])/2.0)
                else:
                  self.imfiles[i].MatchDict['ypix_shift_err']= 0.00067
                dy[i] =dy[i] + vy
                vr = float(tmp[3])
                # if (abs(vr) < abs(dr_min)):                                       
                #  pdb.set_trace()                                                  
                #  dr_min=vr                                                        
                #  self.refImage.MatchDict=self.imfiles[i].MatchDict                
                if (dn[i]>1):
                  self.imfiles[i].MatchDict['angle_err']= abs((vr-dr[i])/2.0)
                else:
                  self.imfiles[i].MatchDict['angle_err']= 0.00002
                dr[i] = dr[i] + vr
                self.imfiles[i].MatchDict['Nmatch']=10
                # self.imfiles[i].MatchDict['xpix_shift']=dx[i]                     
                # self.imfiles[i].MatchDict['ypix_shift']=dy[i]                     
        # pdb.set_trace()                                                           
        fr.close()
        for i in range(nfile):
          if (self.imfiles[i].PathName == self.refImage.PathName):
            ni[i]=1
            #dx0 = dx0 + vx                                                         
            #dy0 = dy0 + vy                                                         
            #dr0 = dr0 + vr                                                         
            # dx[i] = (dx[i] - dx0)/ dn[i]                                          
            # dy[i] = (dy[i] - dy0)/ dn[i]                                          
            # dr[i] = (dr[i] - dr0)/ dn[i]                                          
          dx[i] = dx[i] / dn[i]
          dy[i] = dy[i] / dn[i]
          dr[i] = dr[i] / dn[i]
          self.imfiles[i].MatchDict['Nmatch']=10
          self.imfiles[i].MatchDict['fit_xshift']=dx[i]
          self.imfiles[i].MatchDict['fit_yshift']=dy[i]
          self.imfiles[i].MatchDict['angle']=-dr[i]
          self.imfiles[i].MatchDict['xpix_shift']=dx[i]
          self.imfiles[i].MatchDict['ypix_shift']=dy[i]
          self.imfiles[i].MatchDict['xarc_shift']=dx[i]*self.imfiles[i].outscale
          self.imfiles[i].MatchDict['yarc_shift']=dy[i]*self.imfiles[i].outscale
          print >> fw, fmt % (self.imfiles[i].PathName,dx[i],\
                            self.imfiles[i].MatchDict['xpix_shift_err'],dy[i],\
                            self.imfiles[i].MatchDict['ypix_shift_err'],dr[i],\
                            self.imfiles[i].MatchDict['angle_err'],ni[i])
        # self.imfiles[i].MatchDict['angle']=dr[i]                                  
        fw.close()
        del tmp,tmp1,tmp2,dn,ni,fmt
      # pdb.set_trace()                                                             
      del sftxt,fname,dx,dy,dr,pname,nfile

def _irafDitherSky(simplefits,skyest,rms,GlobalBlock,irafskytoo):
  """Uses the 'sky' task of the STSDAS Dither package to estimate a new
  sky value for the 'simplefits' image.  Returns new sky estimate."""

  # set all iraf'y setup stuff
  iraf.flpr();iraf.flpr()
  iraf.stsdas()
  iraf.analysis()
  iraf.dither()
  iraf.unlearn(iraf.dither.sky)
  iraf.reset(imtype='fits')  
  iraf.set(tmp='./')

  # set 'sky' task params
  iraf.dither.sky.input = simplefits
  iraf.dither.sky.masks = ''
  iraf.dither.sky.lower = skyest - 2.4*rms
  # iraf.dither.sky.upper = skyest + 1.5*rms  -- 12/Apr/2002,jpb
  iraf.dither.sky.upper = skyest + 1.6*rms
  iraf.dither.sky.dq = ''
  iraf.dither.sky.subsky = iraf.no         # do not subtract sky!
  # iraf.dither.sky.width = 1.9*rms  -- 12/Apr/2002,jpb
  iraf.dither.sky.width = 2.0*rms
  iraf.dither.sky.bunit = 'counts'
  iraf.dither.sky.expname = ''
  iraf.dither.sky.skyname = ''
  iraf.dither.sky.verbose = iraf.yes
  iraf.dither.sky.tempdir = './'
  iraf.dither.sky.skyvalue = 0.
  iraf.dither.sky.mode = 'h'

  # run sky using mode estimator
  iraf.dither.sky.stat = 'mode'
  iraf.dither.sky()
  sky_mode = iraf.dither.sky.skyvalue
  GlobalBlock.logfile.write(' IRAF dither.sky mode estimate: '+str(sky_mode))

  if irafskytoo == 2: 
    # run sky using mean estimator
    iraf.dither.sky.stat = 'mean'       
    iraf.dither.sky()
    sky_mean = iraf.dither.sky.skyvalue
    GlobalBlock.logfile.write(' IRAF dither.sky mean estimate: '+str(sky_mean))
    skyraf = min(sky_mean,sky_mode)
  else:
    skyraf = sky_mode
        
  return skyraf
    
def _runSXtractor(simplefits,inParFile,GlobalBlock,skyLines,irafskytoo):
  " runs sextractor on simplefits; sends back: catalog, noise, rms "
  # a quirk of SExtractor (or maybe the popen2 module factory function popen3 or
  # both) is forcing us to write the catalog to a file and then read it back in.
  # Under popen3, STDERR and STDOUT are separate descriptors but when one or the
  # the other's buffer fills it blocks the other one.  Because SExtractor writes so
  # much junk to STDERR, the buffer fills,this blocks the catalog going to STDOUT
  # and both just hang.  Though the docs say you can adjust the buffer size, this
  # seems to have no effect, nor does making things unbuffered (bufsize=0) which
  # you would think would fix things.

  # cd into the align directory to write the catalog

  curdir = os.getcwd()
  os.chdir(GlobalBlock.obsAlign)
  GlobalBlock.logfile.write('moved to align directory...')

  shortname = os.path.basename(simplefits)
  cmd       = 'sex ' + simplefits + ' -c ' + inParFile  # parfile specifies filename atemp.cat
  GlobalBlock.logfile.write('running SExtractor on: '+shortname+'  --  with command:\n\t  '+cmd)
  subproc   = popen2.Popen4(cmd, 1)
  errs      = subproc.fromchild.readlines()     # fromchild is stdout and stderr, but nothing is in stdout 
  # get the catalog into a list
  catalog = open('atemp.cat').readlines()

  sky_sxt,rms = None,None
  for f in range(len(errs)):
    if string.find(errs[f],'(M+D) Background:') != -1:
      linef = string.split(errs[f])
      GlobalBlock.logfile.write(shortname+': '+linef[1]+linef[2]+'  '+linef[3]+linef[4])
      # SExtractor seems to overestimate the sky, by (0--0.5)*sigma,
      # depending on how crowded the field is. 
      rms = float(linef[4])
      sky_sxt = float(linef[2])

  if sky_sxt == None or rms == None:
    raise Exception,"No sky or rms given by sextractor!"

  # now look for errors in the output
  # errorSearch returns a dictionary of unique error strings and the
  # number of occurrences found in the passed list (errs).

  foundErrs = pUtil.errorSearch(errs)
  if foundErrs:
    for f in foundErrs.keys():
      GlobalBlock.logfile.write("SExtractError: "+f+": occurred "+str(foundErrs[f])+" times")
      GlobalBlock.errorList.append(('SExtractor', f+": occurred "+str(foundErrs[f])+" times"))

  sxtr_dsky = rms/50.   # somewhat arbitrary lowering of sxtr sky

  # consider alternate sky's
  #pdb.set_trace()   #by xingxing
  if irafskytoo:
    sky_iraf = _irafDitherSky(simplefits,sky_sxt,rms,GlobalBlock,irafskytoo)
    sigdiff = (sky_sxt - sky_iraf)/rms
    if abs(sigdiff) > 0.85:
      skyErrTxt = 'Warning: SXT sky: %.4f  DTH sky: %.4f differ by %.2f sigma' \
                   %(sky_sxt,sky_iraf,sigdiff)
      GlobalBlock.errorList.append((GlobalBlock.modName,skyErrTxt))
      GlobalBlock.logfile.write(skyErrTxt)
      print skyErrTxt
      # in this case, go with sextractor
      # sky = sky_sxt - sxtr_dsky
      ##> jpb, 12/Apr/2002
      # no, don't: just boost up the rms
      rms = max(rms,abs(sky_sxt - sky_iraf))
    ## always go with minimum? ...
    # sky = min(sky_sxt,sky_iraf) ... No.
    # The sextractor sky may be biased high, but it
    # tends to be more consistent than iraf.dither.sky
    # so, if it's less than sky_iraf, assumes that it's
    # better; otherwise, average the two.
    if sky_sxt < sky_iraf:
      sky = sky_sxt
    else:
      # if clearly discrepant, take min value
      if abs(sigdiff) > 0.99:
        sky = min(sky_sxt,sky_iraf)
      else:
        sky = (sky_sxt + sky_iraf)/2.0
  else:
    sky = sky_sxt - sxtr_dsky

  # Check to see if there's a sky value for the image in skyLines,
  # which is the contents of 'default.skies' if SkyFile=1 was set.
  skyLines.checkforfile(shortname)

# clean up
  os.remove('atemp.cat')
  os.chdir(curdir)
  return catalog,sky,rms

def _getsimplefits(fitsfile, extnum, GlobalBlock, scidex=None):
    """
    takes:  fits filename; extension number;
              optional science index for XVER error check
    returns:  simple fits filename; NX; NY
    If sent a simple fits filename, returns name and NX,NY.
    Otherwise, decides what name the simple fits version should have
    (based on mef file name and XNAME,XVER in the given extenstion),
    and looks for that file.  If it finds it, then great; otherwise,
    it calls convertFile.  In either case, returns sf file name and NX,NY."""

    suffix = string.split(fitsfile,'.')[-1]
    base = fitsfile[:-len(suffix)-1]
    if suffix != 'fits':
        raise pyfits.FITS_SevereError,("suffix of "+fitsfile+" not 'fits'?")

    print "Making SCI and RMS images"
    fo = pyfits.open(fitsfile)
    if len(fo) == 1:
       NX = fo[0].header.get('NAXIS1')
       NY = fo[0].header.get('NAXIS2')
       fo.close()
       del fo
       return fitsfile,NX,NY
            
    else:
       xname = fo[extnum].header.get('EXTNAME')
       if xname:
         xname = xname.upper()
         xver  = fo[extnum].header.get('EXTVER')  # may not exist
         if xver:
           expectfile = base+"_"+xname+"_"+str(xver)+"."+suffix
           if scidex and (scidex+1 != xver):
                errtxt = 'WARNING! (Error?) Header EXTVER = '+str(xver)+\
                         'disagrees with sci index = '+str(scidex)+' for ext '+\
                         str(extnum)+' of image '+fitsfile
                GlobalBlock.logfile.write(errtxt)
                GlobalBlock.errorList.append((GlobalBlock.modName,errtxt))
                print errtxt
         else:
           expectfile = base+"_"+xname+"."+suffix
       else:
         expectfile = base+"_0."+suffix
       fo.close()
       shortname = os.path.basename(expectfile)
       # rmsname = string.replace(shortname,'SCI','RMS') #WZ
       # Ok, we know what it's supposed to be named - is it there?
       if not os.path.isfile(expectfile):
            # it's not there, make it ourselves
            GlobalBlock.logfile.write('Simple fits file '+shortname+' does not exist!' + \
                                      '\n    will try to break out mef file '+fitsfile)
            simplefile,warnList = fUtil.convertExtension(fitsfile,extnum+1) #WZ make RMS images
            simplefile,warnList = fUtil.convertExtension(fitsfile,extnum)
            if simplefile == None:
                raise pyfits.FITS_SevereError,("Cannot understand "+fitsfile)
            newfile = os.path.basename(simplefile)
	    GlobalBlock.outputList[newfile] = [os.path.basename(fitsfile)]
            if shortname != newfile:
                raise pyfits.FITS_SevereError,("expected: "+shortname+"  found: "+newfile)
            # if we get this far, then sf file should be there and we know name
            GlobalBlock.logfile.write('Successfully liberated '+shortname)
       else:
            GlobalBlock.logfile.write('Found simple fits file '+shortname)

       NX,NY   = fUtil.getNXNY(expectfile)
       return expectfile,NX,NY

def _run_cmd(program,cmd,GlobalBlock):
  # procedure to run superalign program and check for errors
  GlobalBlock.logfile.write(cmd)
  print cmd
  sproc  = popen2.Popen3(cmd,1)
  output = sproc.fromchild.readlines()
  errs   = sproc.childerr.readlines()
  if errs:
    GlobalBlock.logfile.write(program+' choked on '+cmd)
    for f in errs:
      GlobalBlock.logfile.write(string.strip(f))
      GlobalBlock.errorList.append((program,string.strip(f)))
  sproc.fromchild.close()
  sproc.childerr.close()
  return 0

def subtractextension(filename,whext,sky,skyrms,deAmp,skykeyword,ForceSub,GlobalBlock):
# subtract background on single extension
  simplefits = pyfits.open(filename,'update')
  colvals=None
                    
  if simplefits[whext].header.get(skykeyword) == None or ForceSub:
    # then it doesn't have sky keyword, so we'll go ahead and subtract

    if sky == 0 and (not deAmp):
      # then no need to change the data array
      pass
    else:
      if not deAmp:
        simplefits[whext].data -= sky
      else:
        step,stepsig,Nstep = amputil.getstep(simplefits,sky,skyrms,ext=whext,verb=0,dowrite=0)
        if not colvals:
          (_NY,_NX) = simplefits[whext].data.shape
          colvals   = resize(arange(_NX),(_NY,_NX))
          # colvals array handy for step subtraction
        simplefits[whext].data = \
          where(less_equal(colvals,(_NX/2-1)), \
                simplefits[whext].data+0.5*step-sky, \
                simplefits[whext].data-0.5*step-sky)
        del colvals
        colvals = None

        if simplefits[whext].header.has_key('COMMENT'):
          index = simplefits[whext].header.ascard.index_of('COMMENT')
	elif simplefits[whext].header.has_key('HISTORY'):
          index = simplefits[whext].header.ascard.index_of('HISTORY')
        else:
          index = len(simplefits[whext].header.ascard) - 1
        cmtxt = '  %s  =  %.6e  subtracted  %s' %(skykeyword, sky, \
                time.asctime(time.localtime(time.time()))[3:])
        simplefits[whext].header.ascard.insert(index,pyfits.Card('COMMENT',value=cmtxt))
        simplefits[whext].header.update(skykeyword,float(sky),before='COMMENT')

        GlobalBlock.logfile.write('Subtracting '+str(sky)+' from '+os.path.basename(filename)+\
                   ' and updating header key '+skykeyword)
        if deAmp:
          if Nstep < 1:
            GlobalBlock.logfile.write(' ... insufficient data for determining amplifier step')
          else:
            GlobalBlock.logfile.write((' ... correcting amplifier step = %8.5f +/- %.5f (Ns: %d)'\
                                      %(step,stepsig,Nstep)))
          del step,stepsig,Nstep
                                
        del index,cmtxt
                    
  else:
    # sky is already subtracted
    warnTxt = 'Warning: sky already subtracted for '+os.path.basename(filename)+\
                ' ?  subtractSky() not subtracting again.'
    GlobalBlock.logfile.write(warnTxt)
    GlobalBlock.errorList.append((GlobalBlock.modName,warnTxt))
    del warnTxt

  simplefits.close()
  del simplefits,sky,skyrms

def pixelScale(fext,fname):
  "Return the pixel scale from the WCS info of the passed fitsfile."
  pixscale = None
  # make sure we understand the units first
  ctype1 = fext.header.get('CTYPE1')
  ctype2 = fext.header.get('CTYPE2')
  if (ctype1 == None or ctype2 == None):
    raise KeyError,'Could not find CTYPE Keys in '+fname
  ctype1.strip()
  ctype2.strip()
  if (ctype1 == 'RA---TAN' and ctype2 == 'DEC--TAN'):
    cd_11 = float(fext.header.get('CD1_1'))
    cd_12 = float(fext.header.get('CD1_2'))
    cd_21 = float(fext.header.get('CD2_1'))
    cd_22 = float(fext.header.get('CD2_2'))
    pixscale = math.sqrt((cd_11**2 + cd_12**2 + cd_21**2 + cd_22**2)/2.0)
    pixscale *= 3600.
  return pixscale

def wcsangle(fext,fname):
  "Return the pixel scale from the WCS info of the passed fitsfile."
  pixscale = None
  # make sure we understand the units first
  ctype1 = fext.header.get('CTYPE1')
  ctype2 = fext.header.get('CTYPE2')
  if (ctype1 == None or ctype2 == None):
    raise KeyError,'Could not find CTYPE Keys in '+fname
  ctype1.strip()
  ctype2.strip()
  if (ctype1 == 'RA---TAN' and ctype2 == 'DEC--TAN'):
    cd_11 = float(fext.header.get('CD1_1'))
    cd_12 = float(fext.header.get('CD1_2'))
    cd_21 = float(fext.header.get('CD2_1'))
    cd_22 = float(fext.header.get('CD2_2'))
    angle = (math.atan2(cd_12,cd_22)+math.atan2(cd_21,-cd_11))/2
    angle = angle * 180/3.14159265

  return angle


class ACSDefSearchParam(SearchParam):
    def __init__(self,pardir,modName,lowLim=1.55):
      inParFileName = os.path.join(pardir,modName + '_acs.inpar')
      SearchParam.__init__(self,lowLim,20.0,0.33,12,inParFileName,3.175)
      
class ACSWFCCamera(Camera):
    def __init__(self,pardir,modName):
      self.SearchParam = ACSDefSearchParam(pardir,modName)
      TKWList = KWList()
      TKWList.AddKW("DETECTOR","WFC")
      Camera.__init__(self,"WFC",TKWList,1,2.5)
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER1","F850LP")
      TKWList.AddKW("FILTER2","CLEAR2L")
      Camera.AddFilter(self,Filter('F850LP','jref$oas1611bj_idc.fits',TKWList,ACSDefSearchParam(pardir,modName)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER1","F775W")
      TKWList.AddKW("FILTER2","CLEAR2L")
      Camera.AddFilter(self,Filter('F775W','jref$oas1611bj_idc.fits',TKWList,ACSDefSearchParam(pardir,modName)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER1","F606W")
      TKWList.AddKW("FILTER2","CLEAR2L")
      Camera.AddFilter(self,Filter('F606W','jref$oas1611bj_idc.fits',TKWList,ACSDefSearchParam(pardir,modName)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER1","CLEAR1L")
      TKWList.AddKW("FILTER2","F435W")
      Camera.AddFilter(self,Filter('F435W','jref$oas1611bj_idc.fits',TKWList,ACSDefSearchParam(pardir,modName)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER1","CLEAR1L")
      TKWList.AddKW("FILTER2","F814W")
      Camera.AddFilter(self,Filter('F814W','jref$oas1611bj_idc.fits',TKWList,ACSDefSearchParam(pardir,modName)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER1","F475W")
      TKWList.AddKW("FILTER2","CLEAR2L")
      Camera.AddFilter(self,Filter('F475W','jref$oas1611bj_idc.fits',TKWList,ACSDefSearchParam(pardir,modName)))
      del TKWList

class ACSHRCCamera(Camera):
    def __init__(self, pardir, modName):
      self.SearchParam = ACSDefSearchParam(pardir,modName)
      TKWList = KWList()
      TKWList.AddKW("DETECTOR","HRC")
      Camera.__init__(self,"ACS/HRC",TKWList,0,0.0)

class NIC3DefSearchParam(SearchParam):
    def __init__(self,pardir,modName,lowLim=0.5):
      inParFileName = os.path.join(pardir,modName + '_nic3.inpar')
      SearchParam.__init__(self,lowLim,20.0,0.15,8,inParFileName,4.0)
      
#{}class ACSDefSearchParam(SearchParam):
class WFC3DefSearchParam(SearchParam):
    #def __init__(self,pardir,modName,lowLim=1.55):
    def __init__(self,pardir,modName,type=1,lowLim=2.00): #WZ
	  #TODO: Check for file extension
      #pdb.set_trace()
      if (type ==1):
        inParFileName = os.path.join(pardir,modName + '_ir.inpar')
      else:
        inParFileName = os.path.join(pardir,modName + '_uvis.inpar')
      SearchParam.__init__(self,lowLim,20.0,0.33,12,inParFileName,3.175)
      
class ExtDefSearchParam(SearchParam):
    #def __init__(self,pardir,modName,lowLim=1.55):
    def __init__(self,pardir,modName,lowLim=2.00): #WZ
      inParFileName = os.path.join(pardir,'ext.inpar')
      SearchParam.__init__(self,lowLim,20.0,0.33,12,inParFileName,3.175)

class WFC3UVISCamera(Camera):
    def __init__(self,pardir,modName):
      #{}self.SearchParam = ACSDefSearchParam(pardir,modName)
      self.SearchParam = WFC3DefSearchParam(pardir,modName,0) #WZ
      TKWList = KWList()
      #{}TKWList.AddKW("DETECTOR","WFC")
      #{}Camera.__init__(self,"WFC",TKWList,1,2.5)
      TKWList.AddKW("DETECTOR","UVIS")
      Camera.__init__(self,"UVIS",TKWList,1,2.5)
      del TKWList

      TKWList = KWList()
      TKWList.AddKW("FILTER","F200LP")
      Camera.AddFilter(self,Filter('F200LP','iref$u1r16227i__idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,0))) # WZ TBD
      del TKWList

      TKWList = KWList()
      TKWList.AddKW("FILTER","F225W")
      Camera.AddFilter(self,Filter('F225W','iref$u1r16227i__idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,0)))
      del TKWList

      TKWList = KWList()
      TKWList.AddKW("FILTER","F275W")
      Camera.AddFilter(self,Filter('F275W','iref$u1r16227i__idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,0)))
      del TKWList

      TKWList = KWList()
      TKWList.AddKW("FILTER","F336W")
      Camera.AddFilter(self,Filter('F336W','iref$u1r16227i__idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,0)))
      del TKWList

      TKWList = KWList()
      TKWList.AddKW("FILTER","F350LP")
      Camera.AddFilter(self,Filter('F350LP','iref$u1r16227i__idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1))) # WZ TBD
      del TKWList

      TKWList = KWList()
      TKWList.AddKW("FILTER","F390W")
      Camera.AddFilter(self,Filter('F390W','iref$u1r16227i__idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,0)))
      del TKWList

      TKWList = KWList()
      TKWList.AddKW("FILTER","F606W")
      Camera.AddFilter(self,Filter('F606W','iref$u1r16227i__idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1))) # WZ TBD
      del TKWList

      TKWList = KWList()
      TKWList.AddKW("FILTER","F814W")
      Camera.AddFilter(self,Filter('F814W','iref$u1r16227i__idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1))) # WZ TBD
      del TKWList


class WFC3IRCamera(Camera):
    def __init__(self,pardir,modName):
      #{}self.SearchParam = ACSDefSearchParam(pardir,modName)
      #pdb.set_trace()
      self.SearchParam = WFC3DefSearchParam(pardir,modName)
      TKWList = KWList()
      TKWList.AddKW("DETECTOR","IR")
      # Camera.__init__(self,"WFC3/IR",TKWList,0,0.0)
      Camera.__init__(self,"IR",TKWList,0,0.0) #XXH WZ
      del TKWList

      #TKWList.AddKW("FILTER","F105W")
      #TKWList.AddKW("FILTER","F110W")

      TKWList = KWList()
      TKWList.AddKW("FILTER","F098M")
      Camera.AddFilter(self,Filter('F098M','iref$u1r16228i_idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER","F105W")
      Camera.AddFilter(self,Filter('F105W','iref$u1r16228i_idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER","F110W")
      Camera.AddFilter(self,Filter('F110W','iref$u1r16228i_idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER","F125W")
      Camera.AddFilter(self,Filter('F125W','iref$u1r16228i_idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER","F140W")
      Camera.AddFilter(self,Filter('F140W','iref$u1r16228i_idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER","F160W")
      Camera.AddFilter(self,Filter('F160W','iref$u1r16228i_idc.fits',TKWList,WFC3DefSearchParam(pardir,modName,1)))
      del TKWList
      TKWList = KWList()
      TKWList.AddKW("FILTER","F126N") # WZ Mar 2013
      Camera.AddFilter(self,Filter('F126N','iref$w3m18525i_idc.fits',TKWList,WFC3DefSearchParam(pardir,modName)))
      del TKWList


def DoAlign(imfiles,CameraList,ExtImage,logfile,newalign,newpar,newfits,pardir,idcTab=None):
  GB = GlobalBlock(logfile,newalign,newpar,newfits)
  fitsfiles = []
  for im in imfiles:
    fitsfiles.append(FitsFile(im,'A',GB,0))
    
  CurExtImage = None
  if ExtImage:
    TSearchParam = ExtDefSearchParam(pardir,GB.modName,lowLim=0.0)
    CurExtImage = FitsFile(ExtImage,"",GB,SearchParam=TSearchParam)
    CurExtImage.ChangeToRef()
    CurExtImage.ChangeToExtRef()

  allShifts = []
  TAlign  = Align(fitsfiles,allShifts,CameraList,GB,CurExtImage=CurExtImage,idcTab=idcTab)
  TAlign.setupImages()
  TAlign.recenterImages()
  TAlign.makeMatchCats(0,0)
  TAlignParam = StandAlignParam(1)
  TAlign.match(TAlignParam,0,'',0)
  TAlign.FixWCS()




