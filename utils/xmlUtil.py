#!/usr/bin/env python

# $Id: xmlUtil.py,v 1.16 2005/10/27 20:05:41 anderson Exp $
# ---------------------------------------------------------------------
# ACS pipeline xml writing tools


__version__      = '$Revision: 1.16 $ '[11:-3]
__version_date__ = '$Date: 2005/10/27 20:05:41 $ '[7:-3]
__author__       = "Ken Anderson, anderson@pha.jhu.edu"


import os,string
import pyfits                          
import filters,fUtil

from   pUtil import fillHeader,ptime
from   types import *

def markupImage(fitsfile,dataset=None,datalist=None):
    """mark up a fits file header with xml tags."""

    if not os.path.exists(fitsfile):
        raise IOError,"Can't find file"

    dir,file = os.path.split(fitsfile)
    parts    = string.split(file,'.')
    newfile  = parts[0]+'_'+parts[1]+'.xml'

    fo = pyfits.open(fitsfile)
    try:
        del fo[0].header.ascard["EXTEND"]
    except KeyError:
        pass
    
    # OK, we now have a pyfits object which we can work on
    # open an xmlMessage object and start the tagging.

    xmlheader = xmlMessage()
    xmlheader.docTypeTag("Image")
    if dataset == None:
        xmlheader.startTag('image',name=file)
    else:
        xmlheader.startTag('image',name=file,dataset=dataset)

    # A new requirement on the markup is that the header tag must now
    # use attributes arraytype and arrayversion to indicate whether
    # a data array (if any) is a SCI, ERR or DQ data array as specified
    # by the original extension keywords EXTNAME and EXTVER.  Now
    # this must be handled differently for mutli-extension fits files
    # which will have these keywords and simple fits files which will
    # not.  The info sought must be gleened from the file name in the case
    # of a simple fits file, should it be there.  arraytype and arrayversion
    # will default to SCI and 1 should these values not be found in the
    # file name.

    # So, handle a multi-extension file:
    if len(fo) > 1:
        for f in range(len(fo)):
            if f == 0:
                xmlheader.startTag('header',arraytype="None",arrayversion="None")
            else:
                xmlheader.startTag('header',arraytype=fo[f].header['EXTNAME'],arrayversion=fo[f].header['EXTVER'])
                xmlheader.startTag('ascards')
            ccount = 0
            hcount = 0
            for card in fo[f].header.ascardlist():
                key     = card.key
                value   = card.value
                comment = card.comment
                if key == 'HISTORY':
                    hcount += 1
                    xmlheader.startTag(string.lower(key),value.rstrip(),linenum=hcount)
                    continue
                if key == 'COMMENT':
                    ccount += 1
                    xmlheader.startTag(string.lower(key),value.rstrip(),linenum=ccount)
                    continue
                if type(value) == StringType:
                    realvalue = string.strip(value)
                    if not key and not realvalue:
                        xmlheader.emptyTag('blankcard')
                        continue
                    if not key and realvalue:
                        xmlheader.startTag('sectioncard',realvalue)
                        continue
                if comment:
                    xmlheader.startTag('keyword',value,comment=comment.strip(),name=key)
                else:
                    xmlheader.startTag('keyword',value,name=key)
            xmlheader.endTag('ascards')
            xmlheader.endTag('header')

    # handle a simple fits file (len == 0) and search for arraytype values in the name
    elif len(fo) == 1:
        fileNameVal = fo[0].header['FILENAME']
        flist = string.split(fileNameVal,"_")
        if len(flist) != 4:                           # if len of flist not 4 set the array values to default.
            xmlheader.startTag('header',arraytype="SCI",arrayversion="1")
        else:
            if flist[2] == 'SCI' or flist[2] == 'ERR' or flist[2] == 'DQ':
                ver = flist[3][:1]
                xmlheader.startTag('header',arraytype=flist[2],arrayversion=ver)
            else:
                # use the default, since we don't know what it is now
                xmlheader.startTag('header',arraytype='SCI',arrayversion='1')
        xmlheader.startTag('ascards')
        ccount = 0
        hcount = 0
        for card in fo[0].header.ascardlist():
            key     = card.key
            value   = card.value
            comment = card.comment
            if key == 'HISTORY':
                hcount += 1
                xmlheader.startTag(string.lower(key),value.rstrip(),linenum=hcount)
                continue
            if key == 'COMMENT':
                ccount += 1
                xmlheader.startTag(string.lower(key),value.rstrip(),linenum=ccount)
                continue
            if type(value) == StringType:
                realvalue = string.strip(value)
            if not key and not realvalue:
                xmlheader.emptyTag('blankcard')
                continue
            if not key and realvalue:
                xmlheader.startTag('sectioncard',realvalue)
                continue
            if comment:
                xmlheader.startTag('keyword',value,comment=comment.strip(),name=key)
            else:
                xmlheader.startTag('keyword',value,name=key)
        xmlheader.endTag('ascards')
        xmlheader.endTag('header')

    else:
        raise IOError,"Invalid length on this fits file: "+fitsfile

    
    if datalist:
        for path in datalist:
            datafile = os.path.basename(path)
            xmlheader.emptyTag('data',href=datafile,type='image/x-fits')
    else:
        xmlheader.emptyTag('data',href=file,type='image/x-fits')
        xmlheader.endTag('image')
        xmlheader.write(os.path.join(dir,newfile))
    del fo
    return newfile

def markupTable(fitstable,dataset=None):
    """receives a fits binary table and converts it to xml markup."""
    if not os.path.exists(fitstable):
        raise IOError,"Can't find file"

    dir,file = os.path.split(fitstable)
    parts    = string.split(file,'.')
    newfile  = parts[0]+'_'+parts[1]+'.xml'
    fo = pyfits.open(fitstable)
    
    if not fo[0].header['FILETYPE'] == 'ASN_TABLE':
        raise IOError,"file is not a fits table"

    # OK, we now have a pyfits object which we can work on
    # open an xmlMessage object and start the tagging. 

    xmlheader = xmlMessage()
    if dataset == None:
        xmlheader.startTag('asntable',name=file)
    else:
        xmlheader.startTag('asntable',name=file,dataset=dataset)
    
    # Markup the headers
    
    for f in range(len(fo)):
        xmlheader.startTag('header',arrayversion=f,arraytype=fo[f].name)
        xmlheader.startTag('ascards')
        ccount = 0
        hcount = 0
        for card in fo[f].header.ascardlist():
            key     = card.key
            value   = card.value
            comment = card.comment
            if key == 'HISTORY':
                hcount += 1
                xmlheader.startTag(string.lower(key),value.rstrip(),linenum=hcount)
                continue
            if key == 'COMMENT':
                ccount += 1
                xmlheader.startTag(string.lower(key),value.rstrip(),linenum=ccount)
                continue
            if type(value) == StringType:
                realvalue = string.strip(value)
                if not key and not realvalue:
                    xmlheader.emptyTag('blankcard')
                    continue
                if not key and realvalue:
                    xmlheader.startTag('sectioncard',realvalue)
                    continue
            if comment:
                xmlheader.startTag('keyword',value,comment=comment.strip(),name=key)
            else:
                xmlheader.startTag('keyword',value,name=key)
        xmlheader.endTag('ascards')
        xmlheader.endTag('header')

        # Now mark up the data table

        if not fo[f].data:
            continue
        xmlheader.startTag('data')
        xmlheader.startTag('table',id=fo[0].header['ASN_TAB'])

        # It is quite possible that some of these table data attributes (eg, _names)
        # may evaporate in future versions of pyfits.
 
        if fo[f].data:
            tabdata = fo[f].data
            field_names = tabdata._names
            for i in range(len(tabdata)):
                xmlheader.startTag('tr')
                for j in field_names:
                    xmlheader.startTag('td',tabdata[i].field(j))
                xmlheader.endTag('tr')    
            xmlheader.endTag('table')
            xmlheader.endTag('data')
    xmlheader.endTag('asntable')
    xmlheader.write(os.path.join(dir,newfile))
    del fo
    return newfile

############################## Utility classes for messaging ################################

class xmlMessage:
    """ provide for an xmlMessage object with methods to build your own
    xml message.  Call the StartTag method to add a start tag with optional
    value and attributes.  EndTag adds an end tag.  This is a dumb class, you 
    have to build a well-formed xml doc.  
    Eg.,
    archmessage = xmlMessage()   # get an xmlMessage object
    archmsg = xmlMessage()
    archmsg.startTag('ArchiveMessage')
    archmsg.startTag('run',pipeType='Imaging',version='dPipe')
    archmsg.startTag('date',ptime())
    archmsg.endTag('run')
    archmsg.endTag('ArchiveMessage')
    file = '/tmp/egmessage.xml'
    archmsg.write(file)

    This will write a file named egmessage.xml in /tmp and will look like this:

    <?xml version="1.0"?>
    <ArchiveMessage>
    <run version="dPipe" pipeType="Imaging">
    <date>2001-08-24T12:19:49EDT</date>
    </run>
    </ArchiveMessage>
    """

    def __init__(self):
        """Initiate the xmlMessage object."""
        self.msg = []

    def docTypeTag(self,dtd_type):
        """adds a doctype tag of the form <!DOCTYPE docRoot SYSTEM ...
        to the xml document. The passed dtd specifies which dtd type is to
        referenced in the DOCTYPE tag.  This dtd type must be one of 

        Catalog
        Colorcatalog
        Image
        Module
        Runmessage

        This should only be called once. This tag will be inserted into the
        beginning of the document no matter when it is called but for
        clarity, the caller should call this method once and then only
        at the beginning of the markup process.
        """
        valid_dtds = ["Catalog","Colorcatalog","Image","Module","Runmessage"]

        if dtd_type not in valid_dtds:
            raise TypeError,"Specified DTD not a valid type."

        if   dtd_type == "Catalog" or dtd_type == "Colorcatalog": docRoot = "catalog"
        elif dtd_type == "Image"      : docRoot = "image"
        elif dtd_type == "Module"     : docRoot = "modulemessage"
        elif dtd_type == "Runmessage" : docRoot = "pipelinemessage"

        self.msg.insert(0,('!DOCTYPE '+docRoot+' SYSTEM \n\t\"http://acs.pha.jhu.edu/science/pipeline/DTD/'+dtd_type+'.dtd\"',))
    
    def startTag(self,tag,value='',**attr):
        """adds a tag, a tag-value pair or a tag-value pair w/ attributes to 
        the msg list.  Only add the attributes if present. Value is converted
        to StringType object.
        """
        if attr:
            self.msg.append((tag,str(value),attr))
        else:
            self.msg.append((tag,str(value)))
        return 

    def endTag(self,tag):
        """add an end tag, i.e </tag>, to the msg list."""
        self.msg.append(('/'+tag,))
        return

    def emptyTag(self,tag,**attr):
        """add an empty tag, i.e. <tag/>, to the msg list."""
        if attr:
            self.msg.append((tag,attr,'/'))
        else:
            self.msg.append((tag+'/',))
        return

    def write(self,file):
        """Write the file as an xml message, where each item of the msg list
        is a line in the xml doc. Interpreted thusly:
        a single tag (no value nor attributes): <tag>
               tag-value pair (no attributes) : <tag> value </tag>
                      tag-value w/ attributes : <tag attr='blah'> value </tag>
        end tag (added via the endTag method) : </tag>
                                    empty tag : <tag/>
        """
        xfile = open(file,'w')
        xfile.write('<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n')
        spc = ''
        for tag in self.msg:
            if len(tag) == 3 and tag[-1] == '/':                        # this handles an emptyTag w/ attributes
                astring = ''
                for att in tag[1].items():
                    astring = astring+' '+str(att[0])+'="'+str(att[1])+'"'
                xfile.write(spc+'<'+tag[0]+astring+tag[2]+'>\n')
            elif len(tag) == 3:                                         # this handles a startTag, key-value w/ attributes
                astring = ''
                for att in tag[2].items():
                    astring = astring+' '+str(att[0])+'="'+str(att[1])+'"'
                if not tag[1]:
                    xfile.write(spc+'<'+tag[0]+astring+'>\n')
                    spc += '  '
                else:
                    xfile.write(spc+'<'+tag[0]+astring+'>'+tag[1]+'</'+tag[0]+'>\n')
            elif len(tag) == 2:                                         # this handles a startTag, key-value 
                if tag[1]:
                    xfile.write(spc+'<'+tag[0]+'>'+tag[1]+'</'+tag[0]+'>\n')
                else: 
                    xfile.write(spc+'<'+tag[0]+'>\n')
                    spc += '  '
            else:                                                       # everything else
                if '/' in tag[0][0]:
                    slist = list(spc)
                    spc = ''
                    slist.pop()
                    slist.pop()
                    for i in slist:
                        spc += i
                    xfile.write(spc+'<'+tag[0]+'>\n')
                else:
                    xfile.write(spc+'<'+tag[0]+'>\n')
        return
    

############################## XML catalog functions ##############################################
#
# The follwing functions are designed to allow a catalog-type file (in the ACS pipeline context)
# to be encoded in xml mark up.  The xmlStartCat function should be called first to initialise
# the xml file.  The xmlPars and xmlData functions can be optional, but you'll wind up with an
# xml file with nothing in it. The xmlPars function will add the parameter set passed to it via the
# pardict dictionary.
#
###################################################################################################

def xmlStartCat(outfile,name,imgfile=None):
    """Initialise the xml header of a catalog file. imgfile is a passed fits file.
    """
    xmlcatalog = open(outfile,'w+')
    xmlcatalog.write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n")
    xmlcatalog.write('<!DOCTYPE catalog SYSTEM ' +\
            '\n\t"http://acs.pha.jhu.edu/science/pipeline/DTD/Catalog.dtd">\n')
    if imgfile == None:
        xmlcatalog.write("<catalog type=\"xml/text\" dataset=\""+name+"\" date=\""+ptime()+"\">\n")
    else:
        img     = os.path.basename(imgfile)
        tel     = fUtil.getTel(imgfile)
        inst    = fUtil.getInstr(imgfile)
        det     = fUtil.getDetector(imgfile)
        filter1 = fUtil.getFilter1(imgfile)
        filter2 = fUtil.getFilter2(imgfile)
        acs     = filters.ACSFilterInfo()
        if filter1 and filter2:
            if filter1 not in acs.clear and filter2 not in acs.clear:
                xmlcatalog.write("<catalog type=\"xml/text\" imgsrc=\"" +img+ \
                         "\" dataset=\""+name+"\" telescope=\"" +tel+  \
                         "\" instrument=\""+inst+"\" detector=\"" +det+\
                         "\" filtername=\""+filter1+"/"+filter2+"\" date=\""+ptime()+"\">\n")
            elif filter1 in acs.clear:
                xmlcatalog.write("<catalog type=\"xml/text\" imgsrc=\"" +img+ \
                         "\" dataset=\""+name+"\" telescope=\"" +tel+  \
                         "\" instrument=\""+inst+"\" detector=\"" +det+\
                         "\" filtername=\""+filter2+"\" date=\""+ptime()+"\">\n")
            elif filter2 in acs.clear:
                xmlcatalog.write("<catalog type=\"xml/text\" imgsrc=\"" +img+ \
                         "\" dataset=\""+name+"\" telescope=\"" +tel+  \
                         "\" instrument=\""+inst+"\" detector=\"" +det+\
                         "\" filtername=\""+filter1+"\" date=\""+ptime()+"\">\n")
        else:
            print "No filters found in image:",img,"for markup."
            xmlcatalog.write("<catalog type=\"xml/text\" dataset=\""+name+"\" imgsrc=\"" \
                     +img+ "\" date=\""+ptime()+"\">\n")
    xmlcatalog.close()
    return

def xmlEndCat(outfile):
    """close off an xml catalog file."""
    xmlcatalog = open(outfile,'a')
    xmlcatalog.write("</catalog>")
    xmlcatalog.close()
    return

def xmlStartConfig(outfile):
    """Write the start tag for the configuration element.  This allows
    multiple sets of parameters to be sent to the markup.  This should only
    be called once. """
    
    xmlcatalog = open(outfile,'a')
    xmlcatalog.write("<configuration>\n")
    xmlcatalog.close()
    return

def xmlPars(pardict,outfile):
    """This function prepends xml markup of the input paramaters of a catalog-type file.
    function receives a dictionary of keyword-value pairs and an output filename. This will be
    the same filename the xmlify function receives and appends the catalog markup."""

    filenames = ['CATALOG_NAME','PARAMETERS_NAME','FLAG_IMAGE','COLUMNS','INPUT','OUTPUT']

    format2='<parameter name=\"%s\">%s</parameter>\n'
    #format3='<parameter name=\"%s\" catagory=\"%s\">%s</parameter>\n'

    xmlcatalog = open(outfile,'a')
    for name,value in pardict.items():
        if name in filenames:
            value = os.path.basename(value)
            xmlcatalog.write(format2 % (name,value))

    else:
        xmlcatalog.write(format2 % (name,value))
    xmlcatalog.close()           
    return

def xmlEndConfig(outfile):
    """Write the end tag for the configuration element. This should only
    be called once. """
    
    xmlcatalog = open(outfile,'a')
    xmlcatalog.write("</configuration>\n")
    xmlcatalog.close()
    return

    
def xmlData(infile,outfile,excludeList=None):
    """
    xmlData takes a SExtractor catalog file and marks it with the pipeline/archive
    protocol markup. Must use at least the xmlStartCat function above to initialise
    the xml file and use the xmlEndCat function to close things up nicely.
    Caller can pass a list of fields to exclude from the xml markup, if that is
    desired, if not, the list is None.

    Usage:
      xmlData('catalog','outputfile'[,excludeList=['blah']])

      where catalogType is one of SExtractor, BPZ, etc.
      (maybe more later), catalog is the catalog's
      filename, and outputfile is just that.

      
    """

    if excludeList:
        for field in excludeList:
            print "Note:",field,"being excluded from xml markup of",os.path.basename(infile)

    catalog = open(infile)
    rows    = catalog.readlines()
    catalog.close()
    
    headerList = []
    headerList = fillHeader(infile)

    # Number of sources is the first field of the last row of the rows list.
    lastrow = rows[len(rows)-1]
    # This is just getting rid of trailing zeros
    det_sources= string.split(string.split(lastrow)[0],'.')[0]         
    xmlcatalog = open(outfile,'a')
    xmlcatalog.write(" <data Nobjects=\"%s\">\n" % det_sources)

    # the function fillHeader returns a list of tuples  from the header of
    # the Sextractor catalog.  These tuples may be of length 2 or 3, depending
    # on whether a particular parameter has multiple occurrences (eg. FLUX_APER).
    # The formats will be used to handle those cases for writing the xml catalog.
    
    format2='\t<field name=\"%s\">%s</field>\n'
    format3='\t<field name=\"%s\" number=\"%s\">%s</field>\n'

    for f in range(len(rows)):
        if '#' in rows[f]: continue
        cols = string.split(rows[f])
        source = cols[0]
        xmlcatalog.write("  <object number=\"" + source + "\">\n")
        for j in range(1,len(headerList)):
            param = cols[j]
            if param.lower() == "nan":
                param = "-9999"
                print "Found a nan in catalog ",os.path.basename(infile)
                print "set to -9999 in the markup."
            if len(headerList[j]) == 2:
                column_no,name = headerList[j]
                if excludeList and name in excludeList:
                    continue
                xmlcatalog.write(format2 % (name,param))

            elif len(headerList[j]) == 3:
                column_no,name,number = headerList[j]
                if excludeList and name in excludeList:
                    continue
                xmlcatalog.write(format3 % (name,str(number),param))
            else:
                raise TypeError, "Incompatible tuple size"

        # In the case where some parameter has multiple columns  of output and is
        # the last parameter listed, this for loop attachs that last parameter name
        # to all remaining values. This behaviour is seen in the case where the parameter
        # VIGNET is the last parameter listed and all subsequent columns in the Sextractor
        # catalog are elements of the VIGNET array.
        
        for j in range(len(headerList),len(cols)):
            param = cols[j]
            xmlcatalog.write(format2 % (name,param))
        xmlcatalog.write("  </object>\n")
        #break
    xmlcatalog.write(" </data>\n")
    xmlcatalog.close()
    return 


def xmlColorData(infile,outfile,excludeList=None):
    """
    xmlColorData takes a color catalog file produced by the colorCatalog module and marks it with the 
    pipeline/archive protocol markup. Must use at least the xmlStartCat function above to initialise
    the xml file and use the xmlEndCat function to close things up nicely.
    This function is a variation of the xmlData function above but designed specifically for 
    the multicolor catalog produced by the colorCatalog module.  This was necessitated by the
    non-generalised nature of the multicolor catalog, where column headers were dependant on 
    the nature of the data (i.e. different filter names).

    Usage:
      xmlColorData('catalogfile','outputfile')

    K Anderson 02-01-02
      
    """

    if excludeList:
        for field in excludeList:
            print "Note:",field,"being excluded from xml markup",os.path.basename(infile)
   
    catalog = open(infile)
    rows    = catalog.readlines()
    catalog.close()
    
    headerList = []
    headerList = fillHeader(infile)

    # Number of sources is the first field of the last row of the rows list.
    lastrow = rows[len(rows)-1]
    # This is just getting rid of trailing zeros
    det_sources= string.split(string.split(lastrow)[0],'.')[0]         
    xmlcatalog = open(outfile,'a')
    xmlcatalog.write(" <data Nobjects=\"%s\">\n" % det_sources)

    # the function fillHeader returns a list of tuples  from the header of
    # the color catalog.  The fillHeader function _may_ return tuples of length
    # 2 or 3, but in  the case of the color catalog, length 2 is expected only.
    # The crux of the next block is to find all occurances of different filter names
    # which might be present in a color catalog header and keep a count of those
    # filter names.  it does this by appending new filter names to the list fList.
    # This count will be an attribute of the xml markup.  This is done for the general 
    # case where MAG and MAGERR columns are not necessarily paired.

    format2='\t<field name=\"%s\">%s</field>\n'
    
    for f in range(len(rows)):
        if '#' in rows[f]: continue
        cols = string.split(rows[f])
        source = cols[0]
        xmlcatalog.write("  <object number=\"" + source + "\">\n")
        for j in range(1,len(headerList)):
            param = cols[j]
            colId = headerList[j]
            if colId[1] in excludeList:
                continue
            if string.find(colId[1],"_BPZ") != -1 or string.find(colId[1],"_APER_") != -1:
                # the length of the filter lines will be based on this format:
                # single filter,  len == 6: HST_ACS_HRC_F250W_APER_CORR
                # single filter,  len == 6: HST_ACS_HRC_F250W_MAG_BPZ
                # crossed filters len == 7: HST_ACS_HRC_POL120UV_F435W_MAGERR_BPZ
                # raises and exception if len is not one of these values 'cause if'n
                # it ain't, sumpin's horrible wrong.
                filterId   = string.split(colId[1],"_")
                if len(filterId) == 6:
                    filterName = filterId[0]+"_"+filterId[1]+"_"+filterId[2]+"_"+filterId[3]
                    fieldName  = filterId[4]+"_"+filterId[5]
                    xmlcatalog.write(format2 % (fieldName,param))
                elif len(filterId) == 7:
                    filterName = filterId[0]+"_"+filterId[1]+"_"+filterId[2]+"_"+filterId[3]+"_"+filterId[4]
                    fieldName  = filterId[5]+"_"+filterId[6]
                    xmlcatalog.write(format2 % (fieldName,param))
                else:
                    raise SyntaxError,"Error: unexpected filter syntax in multicolor catalog header."
            else:
                column_no,name = colId
                xmlcatalog.write(format2 % (name,param))                    
        xmlcatalog.write("  </object>\n")
    xmlcatalog.write(" </data>\n")
    xmlcatalog.close()
    return 

############################## End XML catalog functions ##############################################
