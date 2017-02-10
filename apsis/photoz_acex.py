#!/usr/local/bin/python -u

# $Id: photoz.py,v 1.11 2004/05/05 19:16:43 anderson Exp $
# ---------------------------------------------------------------------
# $Log: photoz.py,v $
# Revision 1.11  2004/05/05 19:16:43  anderson
# fixed typo "sdterr" to "stderr" line 133.
#
# Revision 1.10  2003/10/29 20:31:18  anderson
# addition of better exception handling if bpz throws one.  Usually, all bpz
# does is send some error to stderr, but this module was not acting upon that.
#
# Revision 1.9  2003/10/02 20:58:03  anderson
# revision addresses Bugzilla issues
# #2726 - Implement new BPZ templates
# #2727 - Implement new filter transmission curves
# New templates and curves were delivered and are now incorporated in
# bpz.1.94e.  This revision of photoz adjusts the cli cmd to use these new
# files.  Also, a new helper method, _getBpzVersion, to reap the version
# of bpz, which only appears in the directory name.
#
# Revision 1.8  2002/09/23 19:49:36  anderson
# Reveision to address Bugzilla bug #1795
# - photoz failed to raise a FiltersetError on some single filter WFC datasets
# The _checkfile method now checks to the "M_0" line in the columns file
# which it did not do before.
#
# Revision 1.7  2002/08/20 20:16:53  anderson
# Revision to use a new command line passed to bpz
# From Benitez, Aug 2002.
#
# Revision 1.6  2002/06/27 20:33:07  anderson
# update to mkMsg method for bpz versioning
#
# Revision 1.5  2002/06/20 21:07:44  anderson
# Revision to provide consistent filenames in messages with relative
# pathnames from the path described in the <root> element of each
# module message.  This was an update to this module's mkMsg method.
#
# Revision 1.4  2002/06/20 16:46:47  anderson
# Re: Bugzilla bug #1396:
# The photoz module now has a FiltersetError exception
# class which the private method _checkfile will raise
# if only one filter is present in the dataset.
#
# Revision 1.3  2002/05/16 18:46:13  anderson
# Bringing all apsis modules to a common rev level.
# K Anderson
#
# Revision 1.2  2002/05/15 19:29:40  anderson
# fixed a problem in the __version__=$Revision: 1.11 $ cvs keyword
# interpolation line.  For some reason this line actually read
# __version__ = $$ ,which did not populate the attribute
# with anything.
#
# Revision 1.1.1.1  2002/05/03 17:39:34  anderson
# initial revision and release of the ACS Science pipeline, apsis
#
# ---------------------------------------------------------------------
#
# added the mkMsg method to write out this module's metadata message.
#

__version__      = '$Revision: 1.11 $ '[11:-3]
__version_date__ = '$Date: 2004/05/05 19:16:43 $ '[7:-3]
__author__       = "Ken Anderson, anderson@pha.jhu.edu"


import os, string, popen2
import xmlUtil
from   pUtil import fillHeader,ptime
from   msg   import pMessage
from   sys   import version
pyversion = version
import pdb

class photoZ:
    """
    This module also provides the photoZ class which has only the run method to 
    create a photometric red-shift catalog.  Constructor is passed an ACS pipeline 
    DataSet object and looks for a multicolor catalog and a .columns
    file throws an exception if these file do not exist.  This multicolor
    catalog should be created by a call to the colorCatalog class (see colorcatalog.py).

    usage:

    zcat = photoz.photoZ(obs)
    zcat.run()
    zcat.writeXml()
    zcat.mkMsg()

    Look for the bpz.cat catalog in the particular DataSet's catalog directory.
    """
    def __init__(self, obs):
	self.modName    = string.split(string.split(str(self))[0],'.')[0][1:]   # this is the module name
	self.root       = obs.newobspath                                        # root path of the observation dir
	self.obsName    = obs.newobs
	self.obsCats    = obs.newcats
	self.messagedir = obs.messagedir                                        # where the module message will go.
	self.colorCat   = os.path.join(self.obsCats,'multicolor.cat')
	self.logfile    = obs.logfile                                           # append log entries to this file
	self.bpzCat     = os.path.join(self.obsCats,'bpz.cat')
	self.errorList  = []
	self.inputList  = []
	self.outputList = {}

	self.inputList.append('multicolor.cat')
	self.inputList.append('multicolor.columns')

    def runBpz(self):
	"""Runs BPZ on the multicolor catalog file using the .columns file
	it finds in obs' Catalogs directory.
	"""
        # Need to test that the multicolor catalog actually is multicolor
	# i.e. it has more than one filter's worth of data.  Some bpz utilities
	# will not handle a data file with only one set of filter data.
	# So, check the multicolor.cat file.
	file = os.path.join(self.obsCats,'multicolor.columns') 
	self._checkfile(file) 

        self.logfile.write('Starting photometric redshift determination...')
        print 'Starting photometric redshift determination...'
	bpz = os.path.join(os.environ['P_BPZPATH'],'bpz.py ')
	cmd = 'python ' + bpz + self.colorCat + ' -VERBOSE no -INTERP 2 -DZ 0.005 -SPECTRA CWWSB_Benitez2003.list'
        print cmd
        # pdb.set_trace()
	try:
	    subproc = popen2.Popen3(cmd, 1)
	    stderr_lines = subproc.childerr.read()  
	    if stderr_lines:
                catname = os.path.basename(self.colorCat)
                warntxt = "BPZ produced the following message(s) on "+catname+":"
		self.logfile.write(warntxt)
		print warntxt
		self.errorList.append((self.modName,stderr_lines))
		self.logfile.write(stderr_lines)
                print stderr_lines
                raise RuntimeError,stderr_lines
	    else:
		self.logfile.write('Photometric redshift catalog written')		
		print 'Photometric redshift catalog written'	       
        except Exception, err:
	    self.errorList.append((self.modName,str(err)))
	    raise Exception, err
	else:
	    base = string.split(self.colorCat,'.')[0]
	    os.rename(base+'.bpz',self.bpzCat)                        # mgmt request the bpz filename changed to bpz.cat
	    self.outputList[os.path.basename(self.bpzCat)] = self.inputList
	return

    def writeXml(self):

        """
        writeXmlCatalog marks source catalogs with the pipeline
        protocol markup. Writes the .xml file into the catalogs dir
        of the observation object.

        """
	self.logfile.write('Generating the XML catalog for BPZ output.')

	base = string.split(self.bpzCat,'.')[0]
	xmlCatalogName  =  base +'.xml'
	xmlUtil.xmlStartCat(xmlCatalogName,self.obsName) 
	pardict = self._bpzpars(self.bpzCat)
        xmlUtil.xmlStartConfig(xmlCatalogName)
	try:
	    xmlUtil.xmlPars(pardict,xmlCatalogName) 
	    self.logfile.write('Input Paramaters prepended to xml catalog file.')
	except Exception,err:
	    self.logfile.write('function xmlPars encountered a problem.')
	    self.errorList.append((self.modName,str(err))) 
	    raise Exception,err
        xmlUtil.xmlEndConfig(xmlCatalogName)
        
	try:
	    xmlUtil.xmlData(self.bpzCat,xmlCatalogName)
	    self.logfile.write('XML catalog written for ' + self.bpzCat) 
	except Exception,err:
	    self.logfile.write('function xmlify encountered a problem.')
	    self.errorList.append((self.modName,str(err)))  
	    raise Exception,err
	xmlUtil.xmlEndCat(xmlCatalogName)
	self.outputList[os.path.basename(xmlCatalogName)] = [os.path.basename(self.bpzCat)]
	return


    def mkMsg(self):
	"""create and write module level message for this class.
	Most of this is just compiling the info. meta in a dictionary
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

	self.meta['module'].append(('module','name='+self.modName,'version='+__version__,'dataset='+self.obsName)) 
	#instname = string.split(string.split(str(self))[0],'.')[1]
	self.meta['module'].append(('root',self.root))
	#self.meta['module'].append(('instance', instname))
	self.meta['meta'].append(('meta',))
	self.meta['meta'].append(('depend',))
	self.meta['meta'].append(('pkg',)) 
	self.meta['meta'].append(('name','python')) 
	self.meta['meta'].append(('version',pyversion.split()[0]))
	self.meta['meta'].append(('pkg',)) 
	self.meta['meta'].append(('name','bpz'))
        bpzVersion = self._getBpzVersion()
	self.meta['meta'].append(('version',bpzVersion))


	if self.errorList:
	    self.meta['errorlist'].append(('errorlist',))
	    for pkg,err in self.errorList:
		self.meta['errorlist'].append(('erroritem',err,'frompkg='+pkg))

	# input section
	self.meta['input'].append(('input',))
	for f in self.inputList:
	    if string.find(f,".xml") != -1:
		self.meta['input'].append(('file','type=text/xml'))
		self.meta['input'].append(('name',os.path.join("Catalogs",f)))
	    else:
		self.meta['input'].append(('file','type=text/ascii'))
		self.meta['input'].append(('name',os.path.join("Catalogs",f)))

	# output section
	if self.outputList:
	    self.meta['output'].append(('output',))
	for f in self.outputList.keys():
	    if string.find(f,".xml") != -1:
		self.meta['output'].append(('file','type=text/xml'))
		self.meta['output'].append(('name',os.path.join("Catalogs",f)))
		for pred in  self.outputList[f]:
		    self.meta['output'].append(('predecessor',os.path.join("Catalogs",pred)))
	    else:
		self.meta['output'].append(('file','type=text/ascii'))
		self.meta['output'].append(('name',os.path.join("Catalogs",f)))
		for pred in  self.outputList[f]:
		    self.meta['output'].append(('predecessor',os.path.join("Catalogs",pred)))
		

	# pass this dictionary to the class pMessage...
	msgFile = os.path.join(self.messagedir,self.modName+"_module.xml")
	mmsg = pMessage(self.meta)
	mmsg.writeMsg(msgFile)
	return



    def splice(self):
	"""Method splices the photo-z catalog (.bpz) file and the multicolor.cat
	   file to produce a final photometric redshift catalog.  Raises an exception
	   if one or both files cannot be found.  All these files will be in dir defined 
	   by self.obsCats.
        """

	self.bpzCat = os.path.join(self.obsCats,'bpz.cat')
	if not os.path.exists(self.colorCat):
	    raise IOError, "Multicolor catalog file not found."
	elif not os.path.exists(self.bpzCat):
	    raise IOError, "BPZ catalog file not found."
	
	# Use the fillHeader function to get a list of header lines from each catalog.

	bpzHeaders   = fillHeader(self.bpzCat)
	colorHeaders = fillHeader(self.colorCat)
	allH = bpzHeaders + colorHeaders
	
	# delete the extra element containing the 'NUMBER' column.

	for i in range(len(allH)):
	    col,name = allH[i]
	    if name  == 'NUMBER':
		del allH[i]
		break

	# Renumber the columns via a counter

	for i in range(len(allH)):
	    col,name = allH[i]
	    allH[i]  = (i+1,name)

	# open the new catalog file and write these headers

	newCat = open(os.path.join(self.obsCats,'final_photometry.cat'),'w')
	newCat.write('## Photometry Catalog for Observation: ' + self.obsName +'\n')
	newCat.write('## Generated by the ACS Pipeline, ' + ptime() + '\n')
	newCat.write('##\n')

	f1 = open(self.bpzCat)
	while 1:
	    line = f1.readline()
	    fields = string.split(line)
	    if fields[0] == '##':
		newCat.write(line)
	    else:
		break

	f1.close()
	del f1

	for col,name in allH:
	    newCat.write('# ' + str(col) + '\t' + name + '\n')

	# slurp up the data from each catalog.

	cat1 = open(self.bpzCat).readlines()
	cat2 = open(self.colorCat).readlines()
	
	# grab just the data lines
	cat1Data = []
	cat2Data = []

	for line in cat1:
	    if '#' in line:
		pass
	    else:
		cat1Data.append(string.rstrip(line))

	# Delete the extra field 'NUMBER' from colorCat data as was done (above) for the header.

	for line in cat2:
	    if '#' in line:
		pass
	    else:
		fields = string.split(string.rstrip(line))
		del fields[0]
		newline = string.joinfields(fields)
		cat2Data.append(newline)

	# Write the concatenated line to the new catalog

	if len(cat1Data) != len(cat2Data):
	    raise IndexError,("Catalog length mismatch.")

	for i in range(len(cat1Data)):
		newCat.write(cat1Data[i] + '  ' + cat2Data[i] + '\n')
	newCat.close()
	return


##########################################################################################################
#################################### "private" helper methods ############################################
##########################################################################################################


    def _bpzpars(self,cat):
	"""returns a parameter dictionary (pardict) which is made from a bpz catalog 
	header on the ## lines in the catalog.  This paramter dictionary is to be sent
	to the xml functions to markup the parameter sets for the catalogs.
	"""
	pardict = {}
	catalog = open(cat)
	while 1:
	    newline = string.split(catalog.readline())
	    if '###' in newline:
		pass
	    elif '##' in newline:
		pardict[newline[1]] = newline[3] 
	    else: break
	return pardict

    def _checkfile(self,file):
	"""check that input data contains more than one set of filter data."""
	num_line = 0
	for line in open(os.path.join(self.obsCats,'multicolor.columns')).readlines():
	    fields = string.split(line)
	    # continue looking for an actual line of data. Ignore the M_0 line if it is there.
	    if len(fields) == 0    : continue
	    if fields[0]   == '##' : continue
	    if fields[0]   == '# ' : continue
	    if fields[0]   == 'M_0': continue
	    num_line += 1
	if num_line <= 1:
	    errtxt = "BPZ requires more than one set of filter data."
	    self.errorList.append((self.modName,"ERROR: "+errtxt))
	    self.logfile.write("ERROR: "+errtxt)
	    raise FiltersetError,errtxt
	return

    def _getBpzVersion(self):
        os.chdir(os.environ["P_BPZPATH"])
        bpzPath = os.getcwd()
        bpzId   = os.path.split(bpzPath)[1]
        ver     = bpzId.split("z")[1][1:]
        return ver

class FiltersetError(IOError):
    """Raise this error when there is only one filter in the dataset."""
