#!/usr/bin/env python

#by xingxing  2012/10/12
#from the ingest.py




import os,string,glob
import pdb
import path
import fUtil,xmlUtil
import pyfits
import re           #AKS
import convertCts   #AKS
from   pUtil import logFile
from   msg   import pMessage
from   sys   import version
from  shutil import copyfile,copy

class DataSet:

    def __init__(self,newobs,ownIraf=0):
        self.modName    = string.split(string.split(str(self))[0],'.')[0][1:] # this is the module name
        self.newobs     = newobs                                              # new observation's name
        self.base       = path.Env()
        self.configdir  = self.base.getenv('CONFIGS')
        self.pardir     = self.base.getenv('PARS')
        self.ingestdir  = os.path.join(self.base.getenv('INGEST'), newobs)   # where it is in the INGEST dir
        self.newobspath = os.path.join(self.base.getenv('DATASETS'), newobs, 'final') # new observation's path under $PIPELINE/DATASETS
        self.newpardir     = os.path.join(self.newobspath, 'par')               # parfile dir
        self.fitsdir    = os.path.join(self.newobspath, 'Images')            # the dir where all data is sent
        self.catdir    = os.path.join(self.newobspath, 'Catalogs')
        self.newpar     = os.path.join(self.newobspath, 'par')               # parfile dir
        self.newfits    = os.path.join(self.newobspath, 'Images')            # the dir where all data is sent
        self.newcats    = os.path.join(self.newobspath, 'Catalogs')   # WARNING: notice the above three elements! __________by xingxing
        self.newalign   = os.path.join(self.newobspath, 'align')
        self.prodir    = os.path.join(self.newobspath, 'Output') #WZ
        self.reddir    = os.path.join(self.newobspath, 'Red') #WZ
        self.aligndir   = os.path.join(self.newobspath, 'align')
        self.messagedir = os.path.join(self.newobspath, 'Messages')
        self.root       = self.newobspath
        self.ownIraf    = ownIraf
        self.errorList  = []
        self.inputList  = []
        self.outputList = {}                                                  # outputList is now a dictionary to accomodate 
                                                                              # predecessor images.
        # add by xingxing.  files like these:
	##obs.sciImageList=['a209_f110w_drz.fits','a209_f390w_drz.fits','a209_f814w_drz.fits']
	##obs.contextImageList=['a209_f110w_drz_context.fits','a209_f390w_drz_context.fits','a209_f814w_drz_context.fits']
	##obs.weightImageList=['a209_f110w_drz_weight.fits','a209_f390w_drz_weight.fits','a209_f814w_drz_weight.fits']
	##obs.flagImageList=['a209_f110w_FLAG.fits','a209_f390w_FLAG.fits','a209_f814w_FLAG.fits']
	##obs.rmsImageList=['a209_f110w_RMS.fits','a209_f390w_RMS.fits','a209_f814w_RMS.fits']
	self.datasetsdir=self.base.getenv('DATASETS')  #path
        self.sciImageList=[]
	self.contextImageList=[]
	self.weightImageList=[]
	self.flagImageList=[]
	self.rmsImageList=[]

        # this instance variable is initialised here
        # and again in the buildObs method or not depending on the
        # value of the method call counter.  self.fitslist must be
        # reset if buildObs() is called multiple times though if that
        # is being done, the caller is doing something extraordinary.
        self.fitslist   = []

        # this counter will track the number of calls of the buildObs method.
        self.buildObsCallCounter = 0



    def getObsPath(self):
        """ return the path of the DataSet's root directory. """
        return self.newobspath

    def getFitsPath(self):
        """ return the path of the DataSet's Images directory. """
        return self.fitsdir

    def getCatPath(self):
        """ return the path of the DataSet's Catalogs directory. """
        return self.catdir

    def getProdPath(self):
        """ return the path of the DataSet's products directory. """
        return self.prodir #WZ

    def getProdPath(self):
        """ return the path of the DataSet's products directory. """
        return self.reddir #WZ

    def getParPath(self):
        """ return the path of the DataSet's par directory. """
        return self.newpardir

    def buildObs(self):
        """ set up the new observation directory structure. 
        """
        datasets_dir = self.base.getenv('DATASETS')
        if not os.path.isdir(datasets_dir):
            print 'Cannot find the $DATASETS directory....You should run acex programe first.' #xingxing
	    sys.exit()

        if not os.path.isdir(self.newobspath):
            os.mkdir(self.newobspath)
            os.mkdir(self.fitsdir)
            self.logfile = logFile(self.newobspath)           # Initiate the logfile
            #self.asnDict  = fUtil.makeAsnDict(self.ingestdir)  #Returns a dictionary of asn table names with an embedded list of associated files.

            # buildObsCallCounter tracks the number of times this method is called.
            # Previously every call of this method would reset the fitslist attr.
            # An effort to allow users to use their own fitslists of images
            # for wfp processing of non-wfp filter images.  However, fitslist
            # will still get zeroed if this method is called more than once.    
            #self.buildObsCallCounter += 1
            #if self.buildObsCallCounter == 1:
            #    pass
            #else :
            #    self.fitslist = []
            self.fitslist = []
                
            #for key in self.asnDict.keys():    # delete by xingxing.  will be combined with copy_drz.py later ____________
            #    self.fitslist.append(key)
            #    for file in self.asnDict[key]:
            #        self.fitslist.append(file)
            #for i in self.fitslist:
            #    try:
            #        copyfile(os.path.join(self.ingestdir,i)f,os.path.join(self.fitsdir,i))
            #    except IOError,err:
            #        self.logfile.write("An IOError has occurred in the copyfile call")
            #        self.logfile.write("IOError:"+str(err))
            #        raise IOError,err

            self.inputList = self.fitslist             # inputList is for the mkMsg() method.  
            self.logfile.write('Ingest Data moved to observation FITS dir complete.')
            os.mkdir(self.newpardir)            
            os.mkdir(self.catdir) 
            os.mkdir(self.prodir) #WZ 
            os.mkdir(self.reddir) #WZ 
            os.mkdir(self.aligndir)

           # # get any default.* optional input files           # delete by xingxing.
           # self.defaultlist = glob.glob(self.ingestdir+'/default*')
           # for deffile in self.defaultlist:
           #     copy(deffile,self.fitsdir)

            # read the $PIPELINE/configs/login.cl file, adjust the home and userid settings 
            # and write it to the images dir and make a uparm dir.  See Bugzilla bug # 2077 
            # as to why this is being done.
            if not self.ownIraf:
                irafLoginLines = open(os.path.join(self.configdir,"login.cl")).readlines()
                newLoginFile   = open(os.path.join(self.fitsdir,"login.cl"),"w")

                for line in irafLoginLines:
                    if string.find(line,"set\thome\t\t=") != -1:
                        newLoginFile.write("set\thome\t\t= \""+self.fitsdir+"/\"\n")
                    elif string.find(line,"set\tuserid\t\t=") != -1:
                        newLoginFile.write("set\tuserid\t\t= \""+os.environ["USER"]+"\"\n")
                    else:
                        newLoginFile.write(line)
                newLoginFile.close()
                os.mkdir(os.path.join(self.fitsdir,"uparm"))

                # set the env var MY_IRAF_HOME to be the images dir of the dataset
                os.environ["MY_IRAF_HOME"] = self.fitsdir
            #
            #pdb.set_trace()
            self._setIraf()
            if os.path.isfile(os.path.join(self.fitsdir,"default.shifts")):
                os.rename(os.path.join(self.fitsdir,"default.shifts"),os.path.join(self.aligndir,"default.shifts"))
                self.logfile.write("Warning: buildObs method found a default.shifts file. Moved to align.")
            os.mkdir(self.messagedir)
            self.logfile.write('Directory build complete.')

           # # AKS - Include convert cts/s -> cts script here   # XX  needn't
           # self.logfile.write('Converting counts/sec to counts.')
           # for asn in self.asnDict.keys():
           #   for file in self.asnDict[asn]:
           #     convertCts.convert_fits(self, os.path.join(self.fitsdir, file))

            #for i in self.fitslist:
            #    # Skip the association file marked by '_asn.fits'
            #    if not re.search('_asn', i, re.IGNORECASE):
            #        convertCts.convert_fits(self, os.path.join(self.fitsdir, i))

            #pdb.set_trace()
            #self._rotateHdr() #AKS

        else:
	    print "============================================================================by xingxing" #WZ Apr 2013
	    print "   Named Datasets already exists in $DATASETS. Are you sure to continue ?!!!"
	    print "    Skip the step of copy files, please delete other files except the original fits"
	    print "    Skip the *Detection Image :      >>>>detection_nir.fits  and other files"
	    # print "    Skip the *Detection Catalog  :   >>>>detection_nir.cat" 
	    # print "    Skip the *Filter Catalogs :      >>>>a209_f110w_drz.cat"
	    print "============================================================================by xingxing"
	    self.logfile = logFile(self.newobspath)
	    #######pdb.set_trace()

	mainList=[]      # ['f110w','f390w',...]
	asnDict = fUtil.makeAsnDict(os.path.join(os.environ['INGEST'], self.newobs)) #WZ Nov
	for key in asnDict.keys():           # fitslist
  	  filter=string.split(key,'_')[1]
  	  if str.lower(filter[0])=='f':
  	      mainList.append(filter)
 	  else:
 	      print "Mistake: ",filter
 	mainList.sort() #WZ
    	# print filter
	print str(len(mainList))+' filters are found:'+ str(mainList)

	self.datasetsdir2=os.path.join(self.datasetsdir,self.newobs)
	#tmpaplusfiles=os.listdir(self.datasetsdir2)  #improved  by xingxing
	#self.logfile.write("choose directories under "+self.datasetsdir2+"( only f***w !!) ====IMPORTANT====")
	#for file in tmpaplusfiles:
	#    if (len(file)!=5) or (str.lower(file[0]) !='f') or (str.lower(file[-1])!='w'):
	#        print "__"+file+"__ deleted from the band list."
	#    else:
	#	self.aplusfiles.append(file)
	#	print "__"+file+"__ included ."
	#print self.aplusfiles
 
        #for ii in self.aplusfiles :
        for ii in mainList :
	    ii=str.lower(ii)
            self.sciImageList.append(self.newobs+'_'+ii+'_'+'drz.fits')
	    self.contextImageList.append(self.newobs+'_'+ii+'_'+'drz_context.fits')
	    self.weightImageList.append(self.newobs+'_'+ii+'_'+'drz_weight.fits')
	    self.flagImageList.append(self.newobs+'_'+ii+'_'+'FLAG.fits')
	    self.rmsImageList.append(self.newobs+'_'+ii+'_'+'RMS.fits')
        #print "You can check obs.sciImageList here:"
        #print self.sciImageList
        
	self.skip = 1 	
	for file in self.sciImageList:
	    if not os.path.isfile(os.path.join(self.fitsdir,file)):
		self.skip=0
                print ""+file+" has not been copied! "
	    else:
                print ""+file+" already exists! "
        if self.skip==1:
            print "All image already exist! We will skip all the copy step! ====IMPORTANT===="
        else :
            self.logfile.write("Copy files over: *RMS* *FLAG* and *drz*fits. This would take a while ")
            for i in range(0,len(mainList)) :
	        copyfile(os.path.join(self.datasetsdir2,mainList[i],'Images',self.sciImageList[i]),os.path.join(self.fitsdir,self.sciImageList[i]))
	        copyfile(os.path.join(self.datasetsdir2,mainList[i],'Images',self.contextImageList[i]),os.path.join(self.fitsdir,self.contextImageList[i]))
	        copyfile(os.path.join(self.datasetsdir2,mainList[i],'Images',self.weightImageList[i]),os.path.join(self.fitsdir,self.weightImageList[i]))
	        copyfile(os.path.join(self.datasetsdir2,mainList[i],'Images',self.flagImageList[i]),os.path.join(self.fitsdir,self.flagImageList[i]))
	        copyfile(os.path.join(self.datasetsdir2,mainList[i],'Images',self.rmsImageList[i]),os.path.join(self.fitsdir,self.rmsImageList[i]))
	        self.logfile.write(mainList[i]+" band images be copied!")
	self.logfile.write("Image copy complete")
        return


    def convertData(self):
        """Convert an observation's multi-extension fits data to simple fits format.
        Also calls the xml markup routines to markup the images.
        """
        self.logfile.write('converting data to simple fits format...')
        for asn in self.asnDict.keys():
            for file in self.asnDict[asn]:
                pred_images = []
                pred_images.append(file)
                ffile = os.path.join(self.fitsdir,file)
                self._keycheck(ffile)
                try:
                    newfiles,warnings = fUtil.convertFile(ffile)   # convertFile returns a list of the new files it made.
                    for i in newfiles:                     # This attaches the predecessor image to each breakout file
                        self.outputList[os.path.basename(i)] = pred_images
                    if warnings:
                        for line in warnings:
                            self.errorList.append(('fUtil',line))
                            self.logfile.write(line)
                except IOError,err:
                    print IOError,err
                    self.errorList.append(('fUtil',str(err)))
                    continue
                except Exception,err:
                    print Exception,err
                    self.errorList.append(('fUtil',str(err)))
                    continue
        self.logfile.write('Data conversion complete.')
        return

    def mkMsg(self):
        """create and write module level message for this class.
        Most of this is just compiling the info. meta is a dictionary
        of lists where each list is a list of tuples describing the
        tag lines for the particular section of the message.  This tuple 
        format conforms to that used by the xmlMessage class which is
        modeled on basic python argument passing, i.e. (key,*value,**attr).
        """
        self.meta = {}
        self.meta['module']= []
        self.meta['meta']  = []
        self.meta['input'] = []
        self.meta['output']= []
        self.meta['errorlist'] = []

        self.meta['module'].append(('module','name='+self.modName,'version='+__version__,'dataset='+self.newobs))
        #instname = string.split(string.split(str(self))[0],'.')[1]
        self.meta['module'].append(('root',self.root))
        #self.meta['module'].append(('instance', instname))
        self.meta['meta'].append(('meta',))
        self.meta['meta'].append(('depend',))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','python'))
        self.meta['meta'].append(('version',pyversion.split()[0]))
        self.meta['meta'].append(('pkg',))
        self.meta['meta'].append(('name','pyfits'))
        self.meta['meta'].append(('version',pyfits.__version__.split()[0]))

        if self.errorList:
            self.meta['errorlist'].append(('errorlist',))
            for pkg,err in self.errorList:
                self.meta['errorlist'].append(('erroritem',err,'frompkg='+pkg))

        # input section
        if self.inputList:
            self.meta['input'].append(('input',))
            for f in self.inputList:
                if string.find(f,"_asn") != -1:
                    self.meta['input'].append(('file','type=image/x-fits'))
                    self.meta['input'].append(('name',os.path.join("Images",f)))
                else:
                    self.meta['input'].append(('file','type=image/x-fits'))
                    self.meta['input'].append(('name',os.path.join("Images",f)))

        # output section
        if self.outputList:
            self.meta['output'].append(('output',))
            for f in self.outputList.keys():
                if string.find(f,".xml") != -1:
                    self.meta['output'].append(('file','type=text/xml'))
                    self.meta['output'].append(('name',os.path.join("Images",f)))
                    for pred in  self.outputList[f]:
                        self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
                elif string.find(f,".fits") != -1:
                    self.meta['output'].append(('file','type=image/x-fits'))
                    self.meta['output'].append(('name',os.path.join("Images",f)))
                    for pred in  self.outputList[f]:
                        self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
                else:
                    self.meta['output'].append(('file','type=text/ascii'))
                    self.meta['output'].append(('name',os.path.join("Images",f)))
                    for pred in  self.outputList[f]:
                        self.meta['output'].append(('predecessor',os.path.join("Images",pred)))
        

        # pass this dictionary to the class pMessage...
        msgFile = os.path.join(self.messagedir,self.modName+"_module.xml")
        mmsg = pMessage(self.meta)
        mmsg.writeMsg(msgFile)
        return

#########################################################################################
############################# "private" helper functions ################################
#########################################################################################

    def _setIraf(self):
        """This does nothing but import drizzle in order to grab the login.cl for
        the first time, before anything else can.  Reasons?  believe me, you don't 
        want to know.  It looks like an egregious hack...and maybe it is, but
        we're dealing with iraf here....
        """
        os.chdir(self.fitsdir)
        import xydrizzle
        return

    def _keycheck(self,file):
        """helper function to check keyword conditions in a fits header.
        The argument, file, is a path to a fits file.
        """
        pri_keys = ['SIMPLE', 
                'BITPIX' ,
                'NAXIS',
                'EXTEND' ,
                # 'NEXTEND', 
                'EQUINOX', 
                'EXPTIME' , 
                'FILENAME',
                'TELESCOP',
                'INSTRUME',
                'DETECTOR',
                'FILTER',
                'IDCTAB'  ,
                'DATE-OBS',
                ]
        ext_keys = ['BITPIX',
                'NAXIS',
                'NAXIS1', 
                'NAXIS2',
                'CRVAL1',
                'ORIENTAT', 
                'CRVAL2',
                'CRPIX1',
                'CRPIX2',
                'CD1_1' ,
                'CD1_2' ,
                'CD2_1' ,
                'CD2_2' ,
                'CTYPE1',
                'CTYPE2',
                'EXTNAME',
                'EXTVER'
                ]
        # These keys will only be found in the science extension header
        sciext_keys = ['PHOTZPT' ,
                   'PHOTFLAM',
                   'PHOTPLAM',
                   ]

        filename = os.path.basename(file)
        fobj = pyfits.open(file)
        for i in range(len(fobj)):
            head = fobj[i].header
            # checking the primary header keywords
            if i == 0:
                for key in pri_keys:
                    val = self._lookup(i,head,key,filename)     # eventually this value will be checked.
            # extension keywords. The SCI extension needs to be check for an extra
            # set of keywords (sciext_keys) which other extensions won't have.
            else:
                if head['EXTNAME'] == 'SCI':
                    for key in sciext_keys:
                        val = self._lookup(i,head,key,filename) # eventually this value will be checked.
                for key in ext_keys:
                    val = self._lookup(i,head,key,filename)     # eventually this value will be checked.
        return


    def _lookup(self,i,head,key,filename):
        """try and except clause on a keyword for a header."""
        try:
            return head[key]   
        except KeyError:
            warntxt = "WARNING: "+key+" keyword missing in extension "+str(i)+": "+filename
            self.logfile.write(warntxt)
            self.errorList.append((self.modName,warntxt))
            print pyfits.FITS_Warning,warntxt


    # AKS, WZ
    def _rotateHdr(self):
        """rotate the idcv2ref/idcv3ref data in each header of UVIS image."""

        self.logfile.write('Rotating PA_V3 header data')
        for asn in self.asnDict.keys():
          for file in self.asnDict[asn]:
            fitsfile = pyfits.open(file, mode='update')
            instr = fitsfile[0].header.get('DETECTOR')
            if instr == 'UVIS': #WZ
                for ext in fitsfile:
                    pa = ext.header.get('PA_V3')
                    if pa != None:
                        ext.header.update('PA_V3', pa + 180.0)
                        print 'Updating ', file, '[', ext.name, ']'
                        print "New PA_V3 value - ", ext.header.get('PA_V3')
            fitsfile.close()

#        deg = -45.0
#        angle = (deg) * 3.141519 / 180.0

#              theta = ext.header.get('IDCTHETA')
#              if theta != None:
#                ext.header.update('IDCTHETA', 0.0)
#
#              (v2, v3) = (ext.header.get('IDCV2REF'), ext.header.get('IDCV3REF'))
#              if v2 != None and v3 != None:
#                (v2, v3) = (v2 * math.cos(angle) + v3 * math.sin(angle), 
#                           -v2 * math.sin(angle) + v3 * math.cos(angle))
#              ext.header.update('IDCV2REF', v2)
#              ext.header.update('IDCV3REF', v3)
#
#              print 'pa_v3' , ext.header.get('PA_V3'), \
#               'idctheta', ext.header.get('IDCTHETA'), \
#               'idcv2ref', ext.header.get('IDCV2REF'), \
#               'idcv3ref', ext.header.get('IDCV3REF')
              






















