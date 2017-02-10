#!/usr/bin/env python

# $Id: msg.py,v 1.11 2004/10/18 22:32:44 anderson Exp $
# ---------------------------------------------------------------------

__version__      = '$Revision: 1.11 $ '[11:-3]
__version_date__ = '$Date: 2004/10/18 22:32:44 $ '[7:-3]
__author__       = "Ken Anderson, anderson@pha.jhu.edu"


try:
    import apsisVersion
    apsis_release = apsisVersion.version
except ImportError:
    apsis_release = __version__

import os, os.path, pdb
import string
import path
from   types   import *
from   shutil  import copy
from   time    import sleep
from   xmlUtil import xmlMessage
from   pUtil   import ptime

########################################################################################
################################## pMessage class ######################################
########################################################################################

class pMessage:
    """ this class provides for a pMessage object from which pipeline
    xml messages can be made.  The constructor receives a dictionary
    of of lists which are built (usually) by the mkMsg methods of the
    various pipeline modules.
    """
    def __init__(self,metaInfo):
        """constructor for pMessage.  metadict is a dictionary of lists
        provided by the caller.  
        format of the dictionary:

        {'module':[(key1,*val1,**attr),(key2,*val2,**attr),...]
         'meta'  :[(key1,*val1,**attr),(key2,*val2,**attr),...]
         'input' :[(key1,*val1,**attr),(key2,*val2,**attr),...]
         'output':[(key1,*val1,**attr),(key2,*val2,**attr),...]
         ...
        }
        The tuples are intended to be python-like argument lists to accomodate
        an optional value and optional keyword parameters which the 
        xmlMessage class uses. This class will parse this dictionary and then
        call for an xmlMessage object to be written from the above information.
        The tuples in the above lists will be able to be passed directly to the
        various methods of the xmlMessage object.
        """
        self.metaInfo = metaInfo

        # we make an ordered list of the keys in the metaInfo dictoinary

        self.keys = ['module','meta','input','output','errorlist']
        self.metaElements = ['configuration','depend']
    
    def writeMsg(self,outfile):
        """ Make and write the pipeline message from the passed information.
        """
    
        # call for an xmlMessage object to build the message and then write
        # the message
    
        self.xmlMsg = xmlMessage()
        self.xmlMsg.docTypeTag("Module")
        self.xmlMsg.startTag('modulemessage',date=ptime())
        endTagList = []
        endTagList.append('modulemessage')
        for key in self.keys:

            # Don't write an end tag for the module level tags
            # this needs to wrap the whole message
            while len(endTagList) > 2:
                self.xmlMsg.endTag(endTagList.pop())
            tagList = self.metaInfo[key]
            configTag = None
            for tag in tagList:
                name = tag[0]
                if name == self.metaElements[0]:
                    configTag = 1
                if name == self.metaElements[1] and configTag:
                    self.xmlMsg.endTag(endTagList.pop())
                    configTag = None
                val = ''
                attr = ''
                kdic = {}             # this will be a dictionary of tag attributes for use in apply()

                # A problem with this method of writing the xml just based solely on tag structure is that
                # one cannot determine that some element is not inclusive of another.  This, of course, is
                # what a DTD does for you.  So, in this case, an explicit element <errs> is specified because 
                # there is no a priori way to determine that it is not an element of the <depend> element. 
 
                if endTagList:
                    if name == endTagList[-1]:
                        self.xmlMsg.endTag(endTagList.pop())

                for item in tag[1:]:
                    if type(item) == StringType:
                        if " " in item:            # if there are spaces in item, it MUST be a value string.
                            val = item
                            continue
                        elif "=" in item:
                            part1,part2 = string.split(item,"=")  # attributes in the tuple
                            attr += part1+'='+part2+' '
                        else:
                            val = item
                    else:
                        val = str(item)
                # OK, now we have a tag name (name), maybe a value (val) and maybe some attributes (attr)
                # pass these to the xmlMessage object.  

                if attr:
                    pairs = string.split(attr)
                    for i in range(len(pairs)):
                        atnam,atval = string.split(pairs[i],"=")
                        kdic[atnam]=atval
                    if val:
                        apply(self.xmlMsg.startTag,(name,val),kdic)
                    else:
                        apply(self.xmlMsg.startTag,(name,),kdic)
                        endTagList.append(name)
                else:
                    if val:
                        self.xmlMsg.startTag(name,val)
                    else:
                        self.xmlMsg.startTag(name)
                        endTagList.append(name)
                        
            if key != 'module':
                self.xmlMsg.endTag(endTagList.pop())
        while endTagList:
            self.xmlMsg.endTag(endTagList.pop())

        self.xmlMsg.write(outfile)
        return

##################################End pMessage class################################################


################################## runMessage function##############################################
# I need to refactor this into the pMessage class somehow.  
####################################################################################################

def runMessage(switches,errorlist=None,*pObjs):
    """make a run-level message for a run of the pipeline.  The caller
    (pipeline) sends this function a couple lists (command line switches, errorList)
    and a bunch of pipeline module objects which have everything needed to build the 
    run message.
    """
    # this figures out what the dataset name is for the dataset attribute
    # on the pipelinemessage tag
    #pdb.set_trace()
    if pObjs:
        for obj in pObjs:
            root    = obj.root
            logfile = os.path.basename(obj.logfile.logFileName)
            if obj.modName == "ingest" or 1:    # by xingxing
                obsname = obj.newobs
                break
            else:
                obsname = obj.obsName
                break
    else:
        obsname = "N/A"
    runMsg = xmlMessage()
    runMsg.docTypeTag("Runmessage")
    runMsg.startTag("pipelinemessage",version=apsis_release,date=ptime(),dataset=obsname)
    runMsg.startTag("meta")
    runMsg.startTag("user",os.environ['USER'])
    runMsg.startTag("host",os.uname()[1])
    runMsg.startTag("description","Basic Imaging")
    if switches:
        for sw in switches:
            runMsg.startTag("switch",sw)
    runMsg.startTag("root",root)
    runMsg.startTag("modulelist")
    for ob in pObjs:
        runMsg.startTag("module",ob.modName,message=os.path.join(os.path.basename(ob.messagedir),ob.modName+"_module.xml"))
    runMsg.endTag("modulelist")
    runMsg.endTag("meta")
    runMsg.startTag("input")
    runMsg.startTag("dataset",name=obsname)
    for f in obj.fitslist:
        if string.find(f,"_asn") != -1:
            runMsg.startTag("file","Images/"+f,type="image/x-fits")
        else:
            runMsg.startTag("file","Images/"+f,type="image/x-fits")
    runMsg.endTag("dataset")
    runMsg.endTag("input")
    runMsg.startTag("output")
    #
    # OK, now we search through the outputList first picking out the stuff 
    # we don't want in the archive and signal that.  Need a better way to do
    # this to allow for a selection.
    # RMS and detectionWeight images now appear commented.  This will cause
    # these images to be marked archive="yes" in the runMessage. See Bugzilla
    # bug #1245. The FLAG images are now being marked as archive="no" . See
    # same bug.
    # K Anderson 10-may-2002

    for obj in pObjs:
        for f in obj.outputList.keys():
            if string.find(f,"_shifts_asn_fits.xml") != -1:
                runMsg.startTag("file","Images/"+f, type="text/xml",     archive="no",module=obj.modName)
                continue
            elif string.find(f,"_shifts_asn.fits") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"_cr.fits") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"_SCI_") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"_ERR_") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"_DQ_") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"_FLAG.fits") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"_FLAG_fits.xml") != -1:
                runMsg.startTag("file","Images/"+f, type="text/xml",archive="no",module=obj.modName)
                continue
            elif string.find(f,"medriz_") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"_drz_sci_context.fits.CON") != -1:
                runMsg.startTag("file","Images/"+f, type="text/xml",archive="no",module=obj.modName)
                continue
            elif string.find(f,"Edgemask.fits") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"detectionImage_APER.fits") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="no",module=obj.modName)
                continue
            elif string.find(f,"detectionImage_APER_fits.xml") != -1:
                runMsg.startTag("file","Images/"+f, type="text/xml",    archive="no",module=obj.modName)
                continue
            elif string.find(f,"fits.xml") != -1:
                runMsg.startTag("file","Images/"+f, type="text/xml",    archive="yes",module=obj.modName)
                continue
            elif string.find(f,"fits.hdr") != -1:
                runMsg.startTag("file","Images/"+f, type="text/xml",    archive="no",module=obj.modName)
                continue
            elif string.find(f,".fits") != -1:
                runMsg.startTag("file","Images/"+f, type="image/x-fits",archive="yes",module=obj.modName)
                continue
            elif string.find(f,".xml") != -1:
                runMsg.startTag("file","Catalogs/"+f,type="text/xml",   archive="yes",module=obj.modName)
                continue
            elif string.find(f,".cat.old") != -1:
                runMsg.startTag("file","Catalogs/"+f,type="text/ascii" ,archive="no",module=obj.modName)
                continue            
            elif string.find(f,".cat") != -1:
                runMsg.startTag("file","Catalogs/"+f,type="text/ascii" ,archive="yes",module=obj.modName)
                continue
            elif string.find(f,".matchin") != -1:
                runMsg.startTag("file","align/"+f,  type="text/ascii",  archive="no",module=obj.modName)
                continue
            elif string.find(f,"Match") != -1:
                runMsg.startTag("file","align/"+f,  type="text/ascii",  archive="no",module=obj.modName)
                continue
            else:
                runMsg.startTag("file",f,type="text/ascii",archive="no",module=obj.modName)
    
    runMsg.startTag("file",logfile,type="text/ascii",archive="yes")
    runMsg.endTag("output")
    if errorlist:
        runMsg.startTag("errorlist")
        for i in errorlist:
            runMsg.startTag("erroritem",i)
        runMsg.endTag("errorlist")
    runMsg.endTag("pipelinemessage")
    msgfile = os.path.join(root,obsname+'_runMessage.xml')
    runMsg.write(msgfile)
    return msgfile

#############################################################################################################

def msgLookUp(obs):
    """ poll a directory for new messages from the
    archive. when found, take some action."""

    root = path.Env()
    arch_messages = root.getenv('ARCHMSG')
    processed_messages = os.path.join(arch_messages,'processed')

    while 1:
        msgs = []   
        allfiles  = os.listdir(arch_messages)

        for afile in allfiles:
            if os.stat(os.path.join(arch_messages,afile))[4] == 6011:
                msgs.append(afile)
        if msgs:
            for arch_reply in msgs:
                newMessage = open(os.path.join(arch_messages,arch_reply))
                action = newMessage.readlines()
                print "=========="
                print "Message from acsybase. Request for the following action:\n"
                for i in range(len(action)):
                    print string.rstrip(action[i])
                os.rename(os.path.join(arch_messages,arch_reply), os.path.join(processed_messages,arch_reply))
                print "\nArchive message processed and filed.\n=========="
            break
        else:
            sleep(10)
    
    return

def writeArchMsg(file):
    """ copy (uses shutil module) a file to the ARCHIVE directory."""
    root = path.Env()
    arch_req = root.getenv('ARCH_REQ')
    copy(file,arch_req)
    
    return
