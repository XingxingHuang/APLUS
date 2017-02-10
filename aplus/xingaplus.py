#!/usr/bin/env python

# $Id: APLUS,v 3.00 2012/10/15 as/XX  $
# ---------------------------------------------------------------------

__version__      = '$Revision: 3.00 $ '[11:-3]
__version_date__ = '$Date: 2012/10/15 $ '[7:-3]
__author__       = "Xingxing Huang, <xhuang25@pha.jhu.edu>"


#-------------------------------------------------------------------------------------------#
#                               Handle Cl Options                  
#-------------------------------------------------------------------------------------------#


try:
    import drexVersion
    drex_release = drexVersion.version
    drex_ver_date = drexVersion.version_date
    print 'APLUS v' + drex_release + ' (' + drex_ver_date + ')'
except ImportError:
    drex_release = __version__
    print "This is not a released version of APLUS."
    print "This is revision " +  __version__ + " of APLUS.  (by XX)"

import sys,os,string,getopt,pdb

mod  = os.path.basename(sys.argv[0])

usage= '\n\tUsage: '+ mod + ' [options] obs\n\n'\
       '\t--crlower\tUse dangerously low driz_cr rejection thresholds.\n'\
       '\t--intRef=<im>\tUse fits file <im> in the dataset as reference image.\n'\
       '\t--extRef=<im>\tUse external fits file <im> as reference image.\n'\
       '\t--idctab=<file>\tUse fits IDC table <file> instead of IDCTAB from header.\n'\
       '\t--justgrism\tThis dataset aligned as a pure grism observation.\n'\
       '\t--starcluster\tUses the starcluster parameter set for alignment.\n'\
       '\t--alpixthr=<sig>\tUse this sxtr detection thresh for align.\n'\
       '\t--noAve\t\tExtension averaging will *not* be done during sky subtraction.\n'\
       '\t--nosub\t\tNo sky subtraction will be done.\n'\
       '\t--minSky\tUse minimum sky value from multi extension images.\n'\
       '\t--keepstep\tDo not try to remove amplifier discontinuity for WFC3/UVIS.\n'\
       '\t--skyfile\tLook for a default.skies listing sky values.\n'\
       '\t--noclean\tIntermediate drizzle products will not be deleted.\n'\
       '\t--noCRreject\tDo not reject cosmic rays (1-pass drizzling).\n'\
       '\t--richfield\tUse rich (cluster) field deblend parameters.\n'\
       '\t--justDriz\tWill quit immediately after drizzling.\n'\
       '\t--smartstack\tUse optimal median stacking during drizzle.\n'\
       '\t--noContext\tWill not produce context images in final drizzling.\n'\
       '\t--keepbits=<number>\tSum of pixel data quality flags to keep during drizzle.\n'\
       '\t--noGSC \tWill not attempt to correct image WCS by webquery to GSC2.\n'\
       '\t--nocats\tNo catalogs should be produced for this pipeline run.\n'\
       '\t--notrot \tDo not use rotation in transform found by \'match\' (use med shifts).\n'\
       '\t--notrim\tDo not trim detection weight where Nim=1.\n'\
       '\t--dfilt=<filt> \tUse the distortion model in the IDCTAB for this filter. (Currently disabled)\n'\
       '\t--excludefilt=<filtname> \tCSV list of filters which will not be included in .\n'\
       '\t                         \tthe detection image.\n'\
       '\t--padfac=<f> \tPad linear output size by this factor.\n'\
       '\t--outsize=<nx>,<ny>\tUse this output N_x,N_y *pix* (comma separated) image size.\n'\
       '\t--outshift=<dxA>,<dyA>\tApply these overall x,y *ASEC* shifts (csv) in output frame.\n'\
       '\t--dr_asecpix=<outscale>\tPixel scale of final drizzled image (asec/pix).\n'\
       '\t--dr_pixfrac=<pixfrac>\tDriz pixfrac ("dropsize") param used for final drizzling.\n'\
       '\t--dr_kernel=<kernel>\tDriz kernel (square|point|turbo|tophat|lanczos2|lanczos3|gaussian).\n'\
       '\t--mask=<mask_file>\tExplicitly mask out specific regions in certain input files.\n'\
       '\t--noXml \tDo not do XML markup of data products.\n'\
       '\t--OwnIraf \tUse the login.cl file from the users environment.\n'\
       '\t--testonly \tTest configuration; exit before directory build.\n'\
       '\t--debug \tPrint full traceback when an exception occurs.\n'\
       '\t--help\t\tThis message.\n\n'

long_options = ['help',
                'debug',
                'testonly',
                'noAve',
                'minSky',
                'nocats',
                'noGSC',
                'noContext',
                'justDriz',
                'smartstack',
                'excludefilt=',
                'nosub',
                'notrim',
                'noXml',
                'starcluster',
                'justgrism',
                'keepstep',
                'notrot',
                'richfield',
                'noclean',
                'crlower',
                'intRef=',
                'extRef=',
                'dr_asecpix=',
                'dr_pixfrac=',
                'dr_kernel=',
                'keepbits=',
                'idctab=',
                'skyfile',
                'noCRreject',
                'ownIraf',
                'padfac=',
                # 'dfilt=',
                'outsize=',
                'outshift=',
                'noContext',
                'mask='
                ]

try:
    opts, arg = getopt.getopt(sys.argv[1:],'',long_options)
except getopt.GetoptError:
    # print help information and exit:
    sys.exit(usage)

if not arg:
    sys.exit(usage)

# Only ONE observation (argument) can be specified
if len(arg) != 1:
    sys.exit(usage)

observation = arg[0]
#if not os.path.isdir(os.path.join(os.environ["DATASETS"],observation)):  #by xingxing
if not string.find(os.environ["DATASETS"],observation):
    sys.exit("No observation, "+observation+", found in DATASETS.")    

# initialise the switches.  Default state is off (None).

debug   = 1
testonly= None
noAve   = None
minSky  = None
nocats  = None
noGSC   = 1 # WZ 
justDriz= None
smartstack= None
nosub   = 1 # WZ  No sky subtraction will be done.
notrim  = None
padfac  = None
outsize = None
outshift= None
noXml   = 1
ownIraf = None
starclus= None
alpixthr= None
grism   = None
keepstep= None
notrot  = None
noclean = None
crlower = None
intRef  = None
extRef  = None
idctab  = None
dfilt   = None
skyfile = None
richfield= None
noContext= None
excludefilt= None
dr_asecpix = 0.065 # None
dr_pixfrac = None
dr_kernel  = 'square'
keepbits   = None
noRej      = 1
superalign = 0   #Will not attempt to reject CR's
noContext = None
MaskFile = None

dr_kernlist = ['square','point','turbo','tophat','lanczos2','lanczos3','gaussian']

cl_switches = []
if opts:
    for o, a in opts:
        if a:
            cl_switches.append(o+"="+a)
        else:
            cl_switches.append(o)
            
        if o in ("--help",):
            sys.exit(usage)
        print '\n',
        
        if o in ("--noAve",):
            print "Running subtractSkies without extension averaging."
            noAve = 1
            continue
        
        if o in ("--debug",):
            print "Using debugging mode."
            debug = 1
            continue
        
        if o in ("--testonly",):
            print "Will test configuration and exit before directory build."
            testonly = 1
            continue
        
        if o in ("--minSky",):
            print "Will use minimum sky value for mef data."
            minSky = 1
            continue
        
        if o in ("--nocats",):
            print "No catalogs will be written for this run."
            nocats = 1
            continue
        
        if o in ("--noGSC",):
            print "Will not attempt to correct image WCS by webquery to GSC2."
            noGSC = 1
            continue
        
        if o in ("--noContext",):
            print "Will not produce context images during final drizzling."
            noContext = 1
            continue
        
        if o in ("--justDriz",):
            print "Will quit immediately after drizzling."
            justDriz = 1
            continue
        
        if o in ("--richfield",):
            print "Using rich field deblend parameter set for source detection."
            richfield = 1
            continue	
        
        if o in ("--nosub",):
            nosub = 1
            print "Sky subtraction will not be done."
            continue
        
        if o in ("--ownIraf",):
            ownIraf = 1
            print "DANGER! You've set secret 'ownIraf' flag to use your own iraf setup...."
            continue
        
        if o in ("--starcluster",):
            starclus = 1
            print "Treating dataset as a starcluster"
            continue
        
        if o in ("--justgrism",):
            grism = 1
            print "Treating dataset as a grism observation for align."
            continue
        
        if o in ("--excludefilt",):
            excludefilt=string.split(a,',')
            print "Excluding filter(s) ",excludefilt,"from detection image."
            continue
        
        if o in ("--notrot",):
            notrot = 1
            print "Will not allow match transforms with rotation (does not affect header trans)."
            continue
        
        if o in ("--keepstep",):
            keepstep = 1
            print "Will not try to remove residual amp step in UVIS data."
            continue
        
        if o in ("--noclean",):
            noclean = 1
            print "No deletion intermediate drizzle products."
            continue
        
        if o in ("--notrim",):
            notrim = 1
            print "Will not trim detection weight image where Nim=1 (detects edge CRs)."
            continue
	    
        if o in ("--smartstack",):
            print "Using optimal median stacking during drizzle."
            smartstack = 1
            notrim = 1
            print "Turning on notrim option for use with smartstack. Will not trim detection weight image."
            continue
        
        if o in ("--padfac",):
            padfac = float(a)
            print "Using",padfac,"to rescale output linear image size."
            continue
        
        if o in ("--alpixthr",):
            alpixthr = float(a)
            print "Using",alpixthr,"for SXTR pixthresh in align."
            continue
        
        if o in ("--outsize",):
            if ',' not in a:
                sys.exit(usage)
            __nx,__ny = (int(a.split(',')[0]),int(a.split(',')[1]))
            if min(__nx,__ny) < 2:
                sys.ext('Outsize dimensions must be > 1')
            outsize = (__nx,__ny)
            del __nx,__ny
            print "Using Nx,Ny = %d,%d for output image size." %(outsize[0],outsize[1])
            continue
        
        if o in ("--outshift",):
            if ',' not in a:
                sys.exit(usage)
            __dx,__dy = (float(a.split(',')[0]),float(a.split(',')[1]))
            outshift = (__dx,__dy)
            del __dx,__dy
            print "Using dx,dy = %G,%G arcsec for overall output image shift." %(outshift[0],outshift[1])
            continue
        
        if o in ("--noXml",):
            noXml = 1
            print "Will not do XML markup."
            continue
        
        if o in ("--crlower",):
            crlower = 1
            print "Using lower CR rejection thresholds."
            continue
        
        if o in ("--noCRreject",):
            noRej = 1
            print "Will not attempt to reject CR's."
            continue
        
        if o in ("--intRef",):
            intRef = a
            print "Using",intRef,"as internal match reference image."
            continue
        
        if o in ("--extRef",):
            extRef = a
            print "Using",extRef,"as external match reference image."
            continue
        
        if o in ("--dr_asecpix",):
            dr_asecpix = float(a)
            print "Using",dr_asecpix,"for output pixel scale in final images."
            continue
        
        if o in ("--dr_pixfrac",):
            dr_pixfrac = float(a)
            print "Using",dr_pixfrac,"for final drizzle pixfrac (dropsize) param."
            continue
        
        if o in ("--dr_kernel",):
            dr_kernel = a
            if a not in dr_kernlist:
                sys.exit('drizzle kernel must be one of:\n  '+str(dr_kernlist))
            else:
                print "Using",dr_kernel,"drizzle kernel."
            continue
        
        if o in ("--keepbits",):
            keepbits = int(a)
            if keepbits < 0 or keepbits > 16384:
                sys.exit('bit total must be > 0 and < 16384')
            else:
                print "Will keep flagged bits indicated by "+str(keepbits)
            continue
        
        if o in ("--idctab",):
            idctab = a
            print "Using",idctab,"as IDCTAB reference file."
            continue
        
        if o in ("--noContext",):
            noContext = 1
            print "Will not generate context images."
            continue
	    
        if o in ("--mask",):
            MaskFile = a
            print "Using",MaskFile,"as mask file."
            continue
        
        if o in ("--skyfile",):
            skyfile = 1
            print "Will look for a default.skies file."
            continue


if padfac and outsize:
    sys.exit("\nAhem. 'padfac' and 'outsize' switches are mutually exclusive.")

#-------------------------------------------------------------------------------------------#
#                         End command line options handler    
#-------------------------------------------------------------------------------------------#

# Some imports 

print "\nLoading modules..." 

import time, string
import glob
from   shutil import rmtree
import ingest, msg, filters
from   pUtil import ptime,jarIt,rmFiles
import buildAplusDir # by xingxing
#from multiprocessing import Process # by xingxing  

#if dfilt:
#    WFC3filters = filters.WFC3FilterInfo()
#    if dfilt not in WFC3filters.all:
#        sys.exit('Sorry, filter '+dfilt+' not a recognized WFC3 filter.\n')
#    wheel = WFC3filters.getFilterAttr(dfilt,'wheel')
#    if wheel == '1':
#        dfilts = (dfilt,'CLEAR')
#    elif wheel == '2':
#        dfilts = ('CLEAR',dfilt)
#    else:
#        sys.exit('Sorry, '+dfilt+' appears to be an SBC filter.\n')
#    print "Will use IDCTAB row for filter combination: "+str(dfilts)
#    del wheel,WFC3filters
#else:
#    dfilts = (None,None)

# For now use the default coefficients until filter dependent distortion coefficients
# are implenented in the new align module
dfilts = (None,None)

#-------------------------------------------------------------------------------------------#
def jar(obdict,obsdir):
    """ufunc for the pipeline which will write out the pipeline objects which
    exist when called. obdict is a dictionary of objects:names which is built
    as the pipeline progresses.  obsdir is the observation's directory, i.e. obs.root 
    """
    pdir = os.path.join(obsdir,'Picklejar')

    # call the jarIt function which gets an object to pickle and a path to put them in, 
    # This is just so i can recreate a pipeline run if needed.
    # Of course, for this to be useful, the actual data will need to be online.

    print "pickling pipeline objects..."
    try:
	os.mkdir(pdir)
    except OSError:             # OSError will likely be from the pdir already existing
	pass
    for obj in obdict.keys():
	file = obdict[obj]
	jarIt(os.path.join(pdir,file),obj)  # ufunc jarIt from pUtil
    return
#-------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------#

# function to clear out the dataset's Images directory of the temp iraf junk
# See Bugzilla bug # 2077 as to why this is present.

def irafCleanup(obs):
    print "Cleaning up iraf environment..."
    obs.logfile.write("Cleaning up iraf environment...")
    curdir = os.getcwd()
    os.chdir(obs.fitsdir)
    try:
        rmtree("pyraf")              # rmtree func from shutil module
        obs.logfile.write("Images/pyraf dir removed.")
    except OSError,err:
        print "removing pyraf directory failed:"
        print err
    try:
        rmtree("uparm")
        obs.logfile.write("Images/uparm dir removed.")
    except OSError,err:
        print "removing uparm directory failed:"
        print err
    try:
        os.remove("login.cl")
        obs.logfile.write("Images/login.cl file removed.")
    except OSError,err:
        print "failed to remove login.cl file:"
        print err
    obs.logfile.write("Done.")
    return

#-------------------------------------------------------------------------------------------#

# function for printing the full traceback to stdout
import traceback

def show_tb():
    """Prints full traceback"""
    err_head = "%s ERROR %s" % ("-"*30, "-"*30)
    print err_head
    traceback.print_exc(file=sys.stdout)
    print err_head

#-------------------------------------------------------------------------------------------#



#-------------------------------------------------------------------------------------------#
#                                     PIPELINE BEGIN                                          
#-------------------------------------------------------------------------------------------#

print "\nSetting new observation's directory structure..."
print "initialising various lists for the pipeline run."
obdict      = {}
errorList   = []
obs         = buildAplusDir.DataSet(observation,ownIraf=ownIraf)         #by xingxing
obdict[obs] = "obs"

#pdb.set_trace()
print "Starting directory build and copy files from each filter later."
try:
    obs.buildObs()  #by xingxing  only build the direction under $DATASETS.
except NameError,err:
    sys.exit("Named Dataset already exists in $DATASETS")
except IOError,err:
    sys.exit("buildObs threw an IOError exception:\n"+str(err))
except Exception,err:
    warntxt = "buildObs method failed."
    errorList.append(warntxt)
    errorList.append(str(err))
    print warntxt
    if debug:
        show_tb()
    obs.mkMsg()
    msg.runMessage(cl_switches,errorList,obs)
    jar(obdict,obs.root)
    sys.exit(str(err))

print os.getenv('DATASETS')+observation," directory built."

# We actually do not need the files generated so far. Copy the original from working directories
#print "Copy files over: *RMS* *FLAG* and *sci*fits. This would take a while \n"
#product.copy_drz(obs)          # this step has already been included in buildObs()  by xingxing 
#os.chdir(obs.root)  #change later

#-------------------------------------------------------------------------------------------#
#				Define obs
#-------------------------------------------------------------------------------------------#
print "You can check the obs here: by xingxing"
# dir(obs)    check by xingxing
#['__doc__', '__init__', '__module__', '_keycheck', '_lookup', '_rotateHdr', '_setIraf', 'aligndir', 'base', 'buildObs', 'buildObsCallCounter', 'catdir', 'configdir', 'convertData', 'errorList', 'fitsdir', 'fitslist', 'getCatPath', 'getFitsPath', 'getObsPath', 'getParPath', 'getProdPath', 'ingestdir', 'inputList', 'messagedir', 'mkMsg', 'modName', 'newobs', 'newobspath', 'newpardir', 'outputList', 'ownIraf', 'pardir', 'prodir', 'reddir', 'root']

#pdb.set_trace()
#obs.newobspath   	= '/home/zheng/aaplus_runs/datasets/a209/a209'
#obs.newobs		='a209'
#obs.newpardir          ='/home/zheng/aaplus_runs/datasets/a209/a209/par'
#obs.fitsdir            ='/home/zheng/aaplus_runs/datasets/a209/a209/Images'  # path to fits files
#obs.messagedir		='/home/zheng/aaplus_runs/datasets/a209/a209/Messages'
#obs.configdir		='/home/zheng/aaplus/reffiles/configs'
#obs.logfile
#obs.pardir		='/home/zheng/aaplus/reffiles/pars'

##obs.sciImageList=['a209_f110w_drz.fits','a209_f390w_drz.fits','a209_f814w_drz.fits']
##obs.contextImageList=['a209_f110w_drz_context.fits','a209_f390w_drz_context.fits','a209_f814w_drz_context.fits']
##obs.weightImageList=['a209_f110w_drz_weight.fits','a209_f390w_drz_weight.fits','a209_f814w_drz_weight.fits']
##obs.flagImageList=['a209_f110w_FLAG.fits','a209_f390w_FLAG.fits','a209_f814w_FLAG.fits']
##obs.rmsImageList=['a209_f110w_RMS.fits','a209_f390w_RMS.fits','a209_f814w_RMS.fits']

#-------------------------------------------------------------------------------------------#
# Now import the rest of the modules.

import combFilter          #wralign, combDither, combFilter
import filterCatalog, detectionCatalog, colorCatalog, photoz
import time
from  astrometer import WebQueryError

#-------------------------------------------------------------------------------------------#
#                               Detection Image (combFilter.py)
#-------------------------------------------------------------------------------------------#

# The rest calls will go through two rounds: optical and near-infrared: 
#for band in ['opt','nir','red','uv']:                 ## CHANGE HERE TO CHOOSE DIFFERENT BANDS!
#for band in ['opt','nir','uv']: 
bandlist=['nir']
for band in bandlist:

    from bandcount import bandcount # by xingxing  count the filter number in the band. 
    tmp=bandcount(obs,band)
    tmp.bandcounts()
    if len(tmp.bandlist) < 1 :
	print "\nNone filter is detected under "+band+"\n"
	pass
    else :
	print "\n",len(tmp.bandlist)," filters under "+band+" are detected:"
	print tmp.bandlist,"\n"

    TimeDetectionImage=time.time()  # XX
    print "==================================================================="
    print "Now making detection image in ", band," band"
    print "==================================================================="
    try:
        detImage = combFilter.detectionImage(obs,Band=band,excludeFilters=excludefilt,noContext=noContext)
        obdict[detImage] = "detection_"+band
    except combFilter.detectionImageError,err:
        print err.value
        obs.logfile.write("detectionImage constructor threw a detectionImageError exception")
        obs.logfile.write("APLUS, v"+drex_release+" exiting at  "+ptime()+"\n")
        errorList.append(err.value)
        errorList.append(warntxt)
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage)    
        sys.exit("Exit spawned by detectionImageError Exception. \nNo non-grism images in dataset.")
    except Exception,err:
        warntxt = "Error encountered making detectionImage object..."+band
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage)
        jar(obdict,obs.root)
        sys.exit(str(err))

    detImage.setpars()
    try:
        detImage.getStats(band)
    except Exception,err:
        warntxt = "Error encountered in getStats method..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        detImage.mkMsg()           # write the module message
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage)
        jar(obdict,obs.root)
        sys.exit(str(err))

    # detectionImage class now has 4 methods available to make the detection Image.
    # This latest uses an inverse variance algorithm, pixel by pixel.  

    try:
	if os.path.isfile(os.path.join(obs.fitsdir,"detection_"+band+".fits")):
	    obs.logfile.write("detection_"+band+".fits exits!!")
	    obs.logfile.write("skip detImage.mkInvarIm() for saving time  by xingxing ======IMPORTANT=======")
	else:
            detImage.mkInvarIm()   #XX making the detection Image               
    except Exception,err:
        warntxt = "Error encountered making detection image..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        detImage.mkMsg()           # write the module message
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage)
        jar(obdict,obs.root)
        sys.exit(str(err))

    if not noXml:
        try:
            detImage.writeXml()
        except Exception,err:
            warntxt = "Error encountered marking up products..."
            errorList.append(warntxt)
            errorList.append(str(err))
            print warntxt
            if debug:
                show_tb()
            detImage.mkMsg()           # write the module message
            msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage)
            jar(obdict,obs.root)
            sys.exit(str(err))
    
    try:
        detImage.mkMsg()           # write the module message
    except Exception,err:
        warntxt = "Error encountered marking up products..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage)
        jar(obdict,obs.root)
        sys.exit(str(err))

    jar(obdict,obs.root)
    print "Detection Images made."

    # n=irDetection.prepIrList(obs)

    ## build a remove list for cleaning out the input images to save space.
    # rmList = []
    # for key in obs.asnDict.keys():
    #    rmList += obs.asnDict[key]
    #pdb.set_trace() #WZ move down
    
    # calling irafCleanup fn to clean up the iraf junk...

    # irafCleanup(obs)

    ## command line switch --nocats is checked here.  If true, bug out.
    if nocats:
        # call the new pUtil function rmFiles to remove the input images.
        # This is a space saving action only.
        print "Calling rmFiles to remove input from the disk."
        obs.logfile.write("Calling rmFiles to remove input from the disk.")
        nfiles = rmFiles(obs.fitsdir,rmList)
        obs.logfile.write(str(nfiles)+" files removed successfully.")
        print nfiles,"files removed successfully."
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage)
        obs.logfile.write("End of pipeline run as requested. No catalogs written.")
        donestr = "APLUS, v" + drex_release + ", finished processing " + obs.newobs + " at " + ptime()
        obs.logfile.write(donestr + "\n")
        print donestr
        sys.exit("End of pipeline run as requested. No catalogs written.")

        #-------------------------------------------------------------------------------------------#
        #                                     End Of Imaging 
        #-------------------------------------------------------------------------------------------#


        #-------------------------------------------------------------------------------------------#
        #                                     Catalogs Begin
        #-------------------------------------------------------------------------------------------#


            #---------------------------------------------------------------------------------#
            #                      Detection Catalog (detectionCatalog.py)
            #---------------------------------------------------------------------------------#
    TimeDetectionCat=time.time()  # XX
    print "Making source catalog for the detection image..."
    obs.logfile.write("Begin cataloging...")
    # Pass the richfield argument to the detection catalog constructor.
    # richfield indicates the rich field (cluster) par set is to be used.

    #    print "Default catalog parameter set to be used."
    #    print "Internal module determination of set UVIS, IR or WFC"
    #    obs.logfile.write("default parameter set will be used") 
    obs.band=band
    if (band =='opt'):            
        obs.detector='WFC'
    else:
        if (band =='uv'):
            obs.detector='UVIS'
        else:
            obs.detector='IR'
    try:
        dcat = detectionCatalog.detectionCatalog(obs,richfield)
        obdict[dcat] = "detectionCatalog"
    except Exception,err:
        warntxt = "Error encountered making detectionCatalog object " + band
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        dcat.mkMsg()           # write the module message
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage)
        jar(obdict,obs.root)
        sys.exit(str(err))

    ### The responsibility for running SExtractor using the rms and weight images
    ### falls on the pipeline....

    wgtfile = string.replace(dcat.detectionImage,'.fits','_wht.fits')
    dcat.setinpar('WEIGHT_TYPE','MAP_WEIGHT')
    dcat.setinpar('WEIGHT_IMAGE',wgtfile)
    dcat.setinpar('INTERP_TYPE','NONE')
    dcat.setpars()

    try:
	if os.path.isfile(os.path.join(obs.catdir,"detection_"+band+".cat")):
	    obs.logfile.write("detection_"+band+".cat exits!!")
	    obs.logfile.write("skip dcat.run() for saving time by xingxing ======IMPORTANT=======")
	else:
            dcat.run()
    except Exception,err:
        warntxt = "Error encountered in run method..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        dcat.mkMsg()           # write the module message
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat)
        jar(obdict,obs.root)
        sys.exit(str(err))

    if not noXml:
        try:
            dcat.writeXml()
        except Exception,err:
            warntxt =  "Error encountered in writeXml method..."
            errorList.append(warntxt)
            errorList.append(str(err))
            print warntxt
            if debug:
                show_tb()
            dcat.mkMsg()           # write the module message
            msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat)
            jar(obdict,obs.root)
            sys.exit(str(err))

    try:
        dcat.mkMsg()           # write the module message
    except Exception,err:
        warntxt = "Error encountered writing module message..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat)
        jar(obdict,obs.root)
        sys.exit(str(err))

    jar(obdict,obs.root)
    print "Detection catalogs and messages written."

#-------------------------------------------------------------------------------------------#
#                            Filter Catalogs (filterCatalog.py)
#-------------------------------------------------------------------------------------------#
# now the filter catalogs.  The run2() method should always be used in the
# creation of filter catalogs when a detectionImage catalog is available. run()
# can be called if this is not the case.
#-------------------------------------------------------------------------------------------#

    TimeFilterCat=time.time()  # XX
    print "Making source catalogs for the filter images."
    print "This will take a bit..."
    # pdb.set_trace()
    try:
        fcat = filterCatalog.filterCatalog(obs,richfield)
        obdict[fcat]= "filterCatalogs"
    except Exception,err:
        warntxt = "Error encountered making filterCatalog object..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        fcat.mkMsg()           # write the module message 
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat)
        jar(obdict,obs.root)
        sys.exit(str(err))

### The responsibility for running SExtractor using the rms and weight images
### falls on the pipeline....
### Here we set the parameters to use the weight map of the detection image
### detIm.detWgtName and the RMS map of the filters which came out of 
### combDither (see obs.rmsImageList).  Each filter has an associated
### rms map.  We have to set each filter's parameter file individually.
### This is using the reconfigured setpars method which now takes a fitsfile
### from the caller and writes the appropriate par file.

    fcat.setinpar('WEIGHT_TYPE','MAP_WEIGHT,MAP_RMS')
    fcat.setinpar('INTERP_TYPE','NONE')
    # This is where one can manipolate the parameters
    for i in range(len(fcat.sciImageList)):
        im      = fcat.sciImageList[i]
        rms_im  = fcat.rmsImageList[i]
        rmsfile = rms_im
        fcat.setinpar('WEIGHT_IMAGE',wgtfile+','+rmsfile)   # uses the same weight map, wgtfile, as defined above.
        substring=string.split(fcat.inpars['WEIGHT_IMAGE'],'_wht')
        wname = string.split(substring[0],'detection_')[0] + 'detection_'+band+'_wht' + substring[1]
        fcat.setinpar('WEIGHT_IMAGE',wname)
        # pname=string.replace(fcat.inpars['PARAMETERS_NAME'],'_sci.','_opt.')
        # fcat.setinpar('PARAMETERS_NAME',pname)
        # cname=fcat.obsCats+'/'+string.replace(im,'_sci.fits','_opt.cat')
        # fcat.setinpar('CATALOG_NAME',cname)
        #fcat.setinpar('DETECT_THRESH',obs.dethresh)  # Override the fcan values, as values are different for NIR and OPT  WZ
        #fcat.setinpar('ANALYSIS_THRESH',obs.anathresh)
        #fcat.setpars(im,band) # Forget a file
        curdir = os.getcwd()  # chdir by xingxing
        os.chdir(obs.newfits)
        fcat.setpars(im,band)
	print "fcat.setpars "+im+" with "+band+"done!!"
        os.chdir(curdir)
    del wname
    try:
	#print "skip fcat.run2(band) for saving time_________________by xingxing"  ##
	#can not skip totally for fcat.catalogList will be used in colorCatalog.py by XX
	#I try to skip in colorCatalog.py, you can check it!!!
        fcat.run2(band)
    except Exception,err:
        warntxt = "Error encountered in run2 method..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
            fcat.mkMsg()           # write the module message 
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat)
        jar(obdict,obs.root)
        sys.exit(str(err))

    # if (band != 'opt'): # Run SExtractor over optical detection image
    #    fcat.setinpar('WEIGHT_IMAGE',wgtfile+','+rmsfile) 
    #    fcat.setpars(im,band) # Forge a file
    #    fcat.run2(band)
        
    if not noXml:
        try:
            fcat.writeXml()
        except Exception,err:
            warntxt =  "Error encountered in writeXml method..."
            errorList.append(warntxt)
            errorList.append(str(err))
            print warntxt
            if debug:
                show_tb()
                fcat.mkMsg()           # write the module message 
            msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat)
            jar(obdict,obs.root)
            sys.exit(str(err))

    try:
        fcat.mkMsg()           # write the module message 
    except Exception,err:
        warntxt = "Error encountered writing module message..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        print err
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat)
        jar(obdict,obs.root)
        sys.exit(str(err))

    jar(obdict,obs.root)
    print "Filter Catalogs written."

#-------------------------------------------------------------------------------------------#
#                            MultiColor Catalog (colorCatalog.py)
#-------------------------------------------------------------------------------------------#

    TimeMultiColorCat=time.time()  # XX
    print "Instantiating a colorCatalog object for the multicolor catalog..."
    try:
        ccat = colorCatalog.colorCatalog(obs,fcat.catalogList,band) #WZ
        obdict[ccat]= "colorCatalog"
    except Exception,err:
        warntxt = "Error encountered making colorCatalog object..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat)
        jar(obdict,obs.root)
        sys.exit(str(err))

    try:
        ccat.run(band)
    except Exception,err:
        warntxt = "Error encountered in colorCatalog run method..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
            ccat.mkMsg()           # write the module message 
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat)
        jar(obdict,obs.root)
        sys.exit(str(err))

    if not noXml:
        try:
            ccat.writeXml()
        except Exception,err:
            print "Error encountered in writeXml method..."
            warntxt =  "Error encountered in writeXml method..."
            errorList.append(warntxt)
            errorList.append(str(err))
            print warntxt
            if debug:
                show_tb()
            ccat.mkMsg()           # write the module message 
            msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat)
            jar(obdict,obs.root)
            sys.exit(str(err))

    try:
        ccat.mkMsg()           # write the module message 
    except Exception,err:
        warntxt = "Error encountered writing module message..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        if debug:
            show_tb()
        msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat)
        jar(obdict,obs.root)
        sys.exit(str(err))

    jar(obdict,obs.root)
    # """

    TimeCatEnd=time.time()  # XX
    print '\n     Step                    Elapsed time'
    print '1  Detection Image:      '+str(TimeDetectionCat-TimeDetectionImage)+' sec'
    print '2  Detection Catalog:    '+str(TimeFilterCat-TimeDetectionCat)+' sec'
    print '3  Filter Catalog:       '+str(TimeMultiColorCat-TimeFilterCat)+' sec'
    print '4  MultiColor Catalog:   '+str(TimeCatEnd-TimeMultiColorCat)+' sec \n\n'  #XX


import photoz
#-------------------------------------------------------------------------------------------#
#                          Photometric Redshift (aka BPZ) (photoz.py)
#-------------------------------------------------------------------------------------------#
print "==================================================================="
print "Instantiating a photoZ object for the bpz catalog..."
print "==================================================================="
try:
        bpzcat = photoz.photoZ(obs,band)
        obdict[bpzcat] = "bpzCatalog"
except Exception,err:
        warntxt = "Error encountered making bpzCatalog object..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        msg.runMessage(drex_release,cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat)
        jar(obdict,obs.root)
        sys.exit(str(err))

print "Running BPZ..."
try:
        bpzcat.runBpz(band)
except photoz.FiltersetError,err:
        warntxt = "Filterset Error encountered in call to runBpz."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        print "Note: No BPZ catalog was written."
        bpzcat.mkMsg()           # write the module message 
        msgfile = msg.runMessage(drex_release,cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat,bpzcat)
        if archive:
            print "Sending run message as an archive request."
            msg.writeArchMsg(msgfile)    
            jar(obdict,obs.root)
            obs.logfile.write("APLUS, release "+drex_release+", finished processing "+obs.newobs+" at "+ptime()+"\n")
            print "APLUS, v",drex_release,", finished processing",obs.newobs,"at",ptime()
            sys.exit()
except RuntimeError,err:
        warntxt = "Error encountered running bpz..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        print str(err)
        bpzcat.mkMsg()           # write the module message 
        msg.runMessage(drex_release,cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat,bpzcat)
        jar(obdict,obs.root)
        sys.exit(str(err))
except Exception,err:
        warntxt = "Error encountered running bpz..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        print str(err)
        bpzcat.mkMsg()           # write the module message 
        msg.runMessage(drex_release,cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat,bpzcat)
        jar(obdict,obs.root)
        sys.exit(str(err))

if not noXml:
        try:
            bpzcat.writeXml()
        except Exception,err:
            warntxt = "Error encountered in writeXml method..."
            errorList.append(warntxt)
            errorList.append(str(err))
            print warntxt
            bpzcat.mkMsg()           # write the module message 
            msg.runMessage(drex_release,cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat,bpzcat)
            jar(obdict,obs.root)
            sys.exit(str(err))

try:
        bpzcat.mkMsg()           # write the module message 
except Exception,err:
        warntxt = "Error encountered writing module message..."
        errorList.append(warntxt)
        errorList.append(str(err))
        print warntxt
        msg.runMessage(drex_release,cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat,bpzcat)
        jar(obdict,obs.root)
        sys.exit(str(err))

jar(obdict,obs.root)
print "BPZ catalog written."
print "===================================================================\n"
    # """
# This is the end of the band loop: 'opt','nir', and 'red'

#-------------------------------------------------------------------------------------------#
#                                   Cataloging End
#-------------------------------------------------------------------------------------------#

obs.logfile.write("\nMaking output files")
import product
#obs.ref.pscale = dr_asecpix #    Set up wcs for the output image, using the field of shift list new reference: center position
obs.dr_asecpix = dr_asecpix  # by xingxing
output=product.output(obs,bandlist)   #    Copy image and catalogs over to the Output depositiry
# product.update_tweak(obs,bandlist)    #    Read tweak.res in each working directory. Make a general "tweak.res" file

# processing complete. Now, clean up the input files...do not need them hanging around.
# call the new pUtil function rmFiles to remove the input images. This is a disk space saving 
# action only.

pdb.set_trace()
## build a remove list for cleaning out the input images to save space.
#rmList = []
#for key in obs.asnDict.keys():
#    rmList += obs.asnDict[key]
#pdb.set_trace()
#for fname in os.listdir(os.getcwd()):
#    print fname
#    if fname[0]=='_':
#        rmList.append(fname)
#    if string.find(fname,"temp") !=-1:
#        rmList.append(fname)
#    if string.find(fname,"SCI_") !=-1:
#        rmList.append(fname)
#    if string.find(fname,"BACK") !=-1:
#        rmList.append(fname)
#    if string.find(fname,"asn") !=-1:
#        rmList.append(fname)

#pdb.set_trace()
    
# calling irafCleanup fn to clean up the iraf junk...

# irafCleanup(obs)

#print "Calling rmFiles to remove input from the disk."
#obs.logfile.write("Calling rmFiles to remove input from the disk.")
#nfiles = rmFiles(obs.fitsdir,rmList)
#obs.logfile.write(str(nfiles)+" files removed successfully.")
#print nfiles,"files removed successfully."
os.chdir(obs.prodir)

#---------------------------------------------------------------------------------------#
#
#                                   Run Message
#
#---------------------------------------------------------------------------------------#


print "Writing run-level message..."

# the runMessage function takes any number of objects as arguments to make the run message.

msgfile = msg.runMessage(cl_switches,errorList,obs,alImage,drzImage,detImage,dcat,fcat,ccat)

obs.logfile.write("runMessage.xml file written for "+obs.newobs)
donestr = "APLUS, v" + drex_release + ", finished processing " + obs.newobs + " at " + ptime()
obs.logfile.write(donestr + "\n")
print "run message written."
sys.exit(donestr)


#---------------------------------------------------------------------------------------#
#
#                                  Run Message End
#
#---------------------------------------------------------------------------------------#





