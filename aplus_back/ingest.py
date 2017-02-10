#!/usr/bin/env python

# $Id: ingest.py,v 1.17 2005/10/27 19:46:33 anderson Exp $
# ---------------------------------------------------------------------
# the buildObs method now does a copytree from the ingest area, since
# align is now munging  the original fits images (it adds a keyword and
# subtracts the sky levels).  Once this is done, we can't then use those
# hacked files to re-run any pipeline tests.  So we just do a copy now and
# leave the originals in the INGEST area.  This, of course, takes a  lot
# longer.  But, this needed to happen anyway, since you cannot rename across
# partitions and it is quite possible that INGEST and DATASETS could be 
# setup on different partitions.
# K Anderson
# 11-02-02

__version__      = '$Revision: 1.0 $ '[11:-3]
__version_date__ = '$Date: 2010/8/20 19:46:33 $ '[7:-3]
__author__       = 'W Zheng <zheng@pha.jhu.edu>'

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
pyversion = version

class DataSet:

    """ the class defines a new observation directory structure.
    The constructor builds the paths for the new observation
    (newobs).  This tells anything receiving one of these objects
    where all the data is. The buildObs method creates the directory
    structure.

    """

    def __init__(self,newobs,ownIraf=0):
        self.modName    = string.split(string.split(str(self))[0],'.')[0][1:] # this is the module name
        self.newobs     = newobs                                              # new observation's name
        self.base       = path.Env()
        self.configdir  = self.base.getenv('CONFIGS')
        self.pardir     = self.base.getenv('PARS')
        self.ingestdir  = os.path.join(self.base.getenv('INGEST'), newobs)   # where it is in the INGEST dir
        self.newobspath = os.path.join(self.base.getenv('DATASETS'), newobs) # new observation's path under $PIPELINE/DATASETS
        self.newpardir     = os.path.join(self.newobspath, 'par')               # parfile dir
        self.fitsdir    = os.path.join(self.newobspath, 'Images')            # the dir where all data is sent
        self.catdir    = os.path.join(self.newobspath, 'Catalogs'
)
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
        """ set up the new observation directory structure. Copy the fits data from
        the Ingest area to the new observations dir under DATASETS. Copied files are
        determined from the asnDict dictionary.  Also, this method sets up an IRAF 
        environment for the pipeline, removing dependency on the user's environment,
        hopefully.
        """
        datasets_dir = self.base.getenv('DATASETS')
        if not os.path.isdir(datasets_dir):
            print 'Cannot find the $DATASETS directory....WFP will try to create it.'
            try:
               os.makedirs(datasets_dir)
               print 'Created $DATASETS directory %s' % datasets_dir
            except OSError, error:
               print error
               sys.exit()

        if not os.path.isdir(self.newobspath):
            os.mkdir(self.newobspath)
            os.mkdir(self.fitsdir)
            self.logfile = logFile(self.newobspath)           # Initiate the logfile
            self.asnDict  = fUtil.makeAsnDict(self.ingestdir)

            # buildObsCallCounter tracks the number of times this method is called.
            # Previously every call of this method would reset the fitslist attr.
            # An effort to allow users to use their own fitslists of images
            # for wfp processing of non-wfp filter images.  However, fitslist
            # will still get zeroed if this method is called more than once.
        
            self.buildObsCallCounter += 1
            if self.buildObsCallCounter == 1:
                pass
            else :
                self.fitslist = []
                
            for key in self.asnDict.keys():
                self.fitslist.append(key)
                for file in self.asnDict[key]:
                    self.fitslist.append(file)
            for i in self.fitslist:
                try:
                    copyfile(os.path.join(self.ingestdir,i),os.path.join(self.fitsdir,i))
                except IOError,err:
                    self.logfile.write("An IOError has occurred in the copyfile call")
                    self.logfile.write("IOError:"+str(err))
                    raise IOError,err

            self.inputList = self.fitslist             # inputList is for the mkMsg() method.  
            self.logfile.write('Ingest Data moved to observation FITS dir complete.')
            os.mkdir(self.newpardir)            
            os.mkdir(self.catdir) 
            os.mkdir(self.prodir) #WZ 
            os.mkdir(self.reddir) #WZ 
            os.mkdir(self.aligndir)

            # get any default.* optional input files
            self.defaultlist = glob.glob(self.ingestdir+'/default*')
            for deffile in self.defaultlist:
                copy(deffile,self.fitsdir)

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

            # AKS - Include convert cts/s -> cts script here
            self.logfile.write('Converting counts/sec to counts.')
            for asn in self.asnDict.keys():
              for file in self.asnDict[asn]:
                convertCts.convert_fits(self, os.path.join(self.fitsdir, file))

            #for i in self.fitslist:
            #    # Skip the association file marked by '_asn.fits'
            #    if not re.search('_asn', i, re.IGNORECASE):
            #        convertCts.convert_fits(self, os.path.join(self.fitsdir, i))

            #pdb.set_trace()
            self._rotateHdr() #AKS

        else:
            raise NameError,('Named observation dir already exists.\n')
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
        import drizzle
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
              


