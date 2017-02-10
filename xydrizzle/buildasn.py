""" buildAsn.py - Module to generate an ACS-style association
table from all matching files from the current working directory.

This module will also recognize shiftfiles should they be specified
and add the offsets to the ASN table.

W.J. Hack 9 April 2001
    *** INITIAL VERSION

W.J Hack 24 April 2001
    Version 0.1b - Modified to support FITS version 0.4.2

W.J Hack, 1 May 2001
    Version 0.2  - Computes delta from given shifts file and header info
    
WJH, 20 Nov 2001:
    Version 0.3 - Added ROTATION column to default ASN table created.
WJH, 4 Dec 2001:
    Added help() function and printing of Version number.
WJH, 21 Jan 2002:
    Fixed bugs in writing out ASN table with default columns.
    Also added check to insure that files are found to build ASN table.    

WJH, 8 Mar 2002:
    Added MEMPRSNT column and converted to use PyFITS v0.6.2 and numpy.

WJH, 12 June 2002 (0.5a):
    Changed 'pyfits.version' to 'pyfits.__version__' to keep up with PyFITS.
"""

import os, sre,string
import time
import pyfits
import numpy
from numpy import char as chararray


# List of supported default file types
# It will look for these file types by default
# when trying to recognize input rootnames.
EXTLIST =  ['_crj.fits','_flt.fits','_sfl.fits']

__version__ = '0.5a (12-June-2002)'

_prihdr = pyfits.Header([pyfits.Card('SIMPLE', pyfits.TRUE,'Fits standard'),
                pyfits.Card('BITPIX  ',                    16 ,' Bits per pixel'),
                pyfits.Card('NAXIS   ',                     0 ,' Number of axes'),
                pyfits.Card('ORIGIN  ',  'NOAO-IRAF FITS Image Kernel July 1999' ,'FITS file originator'),
                pyfits.Card('IRAF-TLM',  '18:26:13 (27/03/2000)' ,' Time of last modification'),
                pyfits.Card('EXTEND  ',pyfits.TRUE ,' File may contain standard extensions'),
                pyfits.Card('NEXTEND ',                     1 ,' Number of standard extensions'),
                pyfits.Card('DATE    ',  '2001-02-14T20:07:57',' date this file was written (yyyy-mm-dd)'),
                pyfits.Card('FILENAME',  'hr_box_asn.fits'            ,' name of file'),
                pyfits.Card('FILETYPE',  'ASN_TABLE'          ,' type of data found in data file'),
                pyfits.Card('TELESCOP',  'HST'                ,' telescope used to acquire data'),
                pyfits.Card('INSTRUME',  'ACS   '             ,' identifier for instrument used to acquire data'),
                pyfits.Card('EQUINOX ',                2000.0 ,' equinox of celestial coord. system'),
                pyfits.Card('ROOTNAME',  'hr_box  '              ,' rootname of the observation set'),
                pyfits.Card('PRIMESI ',  'ACS   '             ,' instrument designated as prime'),
                pyfits.Card('TARGNAME',  'SIM-DITHER'                     ,'proposer\'s target name'),
                pyfits.Card('RA_TARG ',                    0. ,' right ascension of the target (deg) (J2000)'),
                pyfits.Card('DEC_TARG',                    0. ,' declination of the target (deg) (J2000)'),
                pyfits.Card('DETECTOR',  'HRC     '           ,' detector in use: WFC, HRC, or SBC'),
                pyfits.Card('ASN_ID  ',  'hr_box  '           ,' unique identifier assigned to association'),
                pyfits.Card('ASN_TAB ',  'hr_box_asn.fits'         ,' name of the association table')])

def buildAsnTable(asnroot,suffix=None,shiftfile=None):
    """ Create an association table from all suitable files in
    the current directory selected by their suffix.  By default,
    it will include all 'crj.fits' files.
    
    It would also add the columns XOFFSET, YOFFSET and ROTATION 
    and populate them with values from a shiftfile if one was specified 
    and can be read.
    
    This function returns the name of the DITHER PRODUCT that was
    created, then inserted into the ASN table. It will be 
    the rootname of asnroot with '_dth' appended, following the
    general ACS pipeline filename conventions.
     
    """
    # print out version number for reference
    print 'BuildAsn Version '+__version__
    
    # Set default suffix
    if suffix == None:
        suffix = 'crj.fits'
    
    # Find all files in current directory with given suffix
    flist = _findFiles(suffix)

    if len(flist) == 0: 
        print 'NOTE: Please check suffix used to build list of files.'
        _errmsg = 'Could NOT find any files to build into an association!'
        raise Exception,_errmsg
        
    # Read in values from shiftfile, if present
    shiftdict = {}
    if shiftfile != None:
        # Build dictionary of shift values
        # Set to 'None' upon problems with file
        shiftdict = _readShiftFile(shiftfile,shiftdict)
        # Find file with minimum shift, assume it is reference file
        ref_file = _findRefFile(shiftdict)
        # Compute deltas by subtracting commanded offsets from headers
        # from total offsets in shiftdict
        _computeDeltaShifts(shiftdict,ref_file)
        
    _output = asnroot+'_asn.fits'
    
    # Open output ASN table file
    fasn = pyfits.HDUList()

    # Build ASN file primary header
    _buildAsnPrimary(fasn,_output,flist[0])
        
    # Make table records 
    outTable = _buildTableRows(flist,shiftdict, asnroot, suffix)
    
    # Create an extension HDU which contains a table
    exthdu = _makeTableHDU(outTable)

    fasn.append(exthdu)
    
    # Close ASN table file
    fasn.writeto(_output)
    fasn.close()
    del fasn

    # Generate output product name
    dthprod = _buildNewRootname(asnroot,extn='_dth',suffix=suffix)

    return dthprod

def _findFiles(suffix):
    """ Builds list of all files which contain suffix.
        Suffix should not contain wildcards; instead, they
        should be like 'crj.fits' or '.dat' or 'coeffs'.
    """
    _ldir = os.listdir('.')
    
    # format suffix for use in reg expression
    _indx = string.find(suffix,'.')
    if _indx > 0:
        extn = suffix[:_indx]+'[.]'+suffix[_indx+1:]
    else:
        extn = suffix
    regpatt = '.*'+extn+'.*$'
    
    # compile regular expression
    _reg = sre.compile(regpatt)
    
    # build file list
    flist = []
    for file in _ldir:
        if _reg.match(file):
            flist.append(file)
    
    return flist
    
def _updateKeyword(key,inhdr,outhdr,default='UNKNOWN'):
    """ Safely updates keyword key in outhdr from value in inhdr.
        Uses value given by 'default' if keyword is not found in 
        input header.
    """
    try:
        _keyw = inhdr[key]
    except KeyError:
        _keyw = default
    outhdr[key] = _keyw  
    

def _readShiftFile(filename,sdict):
    """ 
    Function which will read in either a 'shiftfind'
    output file or 'avshift' output file and extract the best XSH and YSH
    for each file. It builds a dictionary based on the filenames.
    sdict = {'image1':(xsh,ysh),...}.  These values are absolute
    shifts which need to have the intended offset subtracted off before
    being written into the 'XOFFSET'/'YOFFSET' columns of the output ASN table.
    """
    lines = []
    
    fshift = open(filename,'r')
    # Read first line
    line = fshift.readline()   
    ntokens = len(string.split(line))
    # Determine what kind of shift file we are working with
    if string.find(line,'Assume') > -1:
        # We are working with average shift file
        xoffindx = 3
        yoffindx = 4
    elif  ntokens == 5:
        # We are working with a shifts file from 'shiftfind'
        xoffindx = 1
        yoffindx = 3
    elif ntokens == 6:
        # We are working with a shifts file from 'shiftfind'
        # with a group number column as column 2
        xoffindx = 2
        yoffindx = 4
    else:
        # We are not working with a file we understand
        # Simply return 'None' and continue
        xoffindx = None
        yoffindx = None
        fshift.close()
        return None

    # Build dictionary now...    
    # Search all lines which have text
    while line != '':
        # Ignore commentary lines that begin with '#'
        if string.find(line,"#") < 0:
            # Any that don't begin with #, append to list
            lines.append(line)
        line = fshift.readline()
    # Done reading in file, so close it now.
    fshift.close()

    # Parse the data lines and populate the dictionary
    for line in lines:
        lsplit = string.split(line)
        sdict[lsplit[0]] = (float(lsplit[xoffindx]),float(lsplit[yoffindx]))        
    return sdict

def _findRefFile(sdict):
    """ Function which identifies the image from an association 
        which serves as the reference, zero-offset position.
        It simply returns the filename for the image with the 
        smallest total shift.
    """
    min_shift = 99999.
    reffile = None
    for img in sdict.keys():
        offx = sdict[img][0]
        offy = sdict[img][1]
        shift = numpy.sqrt(pow(offx,2) + pow(offy,2))
        
        if shift < min_shift: 
            min_shift = shift
            reffile = img
    
    return reffile


def _computeDeltaShifts(sdict,ref_file):
    """ Function which subtracts commanded offsets from headers
        from total offsets found in sdict.
        Commanded offsets will always be relative to ref_file.
    """ 
    # Determine pointing for reference image
    # Open reference image
    fref = pyfits.open(ref_file)
    # Get pointing position from SCI extension header
    ref_hdr = fref['SCI'].header
    ref_pos = (ref_hdr['CRVAL1'],ref_hdr['CRVAL2'])
    ref_cd = (ref_hdr['CD1_1'],ref_hdr['CD1_2'])
    # Close and delete reference image object
    fref.close()
    del fref
    pscale = pow((pow(float(ref_cd[0]),2)+pow(float(ref_cd[1]),2)),0.5) * 3600.

    for img in sdict.keys():
        # Open science image
        fimg = pyfits.open(img)
        scihdr = fimg['SCI'].header
        # Extract commanded position from image header 
        img_pos = (scihdr['CRVAL1'],scihdr['CRVAL2'])
        # Done with image: close and delete FITS object
        fimg.close()
        del fimg

        # Compute commanded shift here: image minus reference
        # This has to be in units of pixels
        pos = (Numeric.array(img_pos) - Numeric.array(ref_pos))/pscale
        pos = pos * 3600. * Numeric.cos(ref_pos[1]) * Numeric.array([-1.0,1.0])

        # Compute delta offset: total minus commanded
        delta_pos = Numeric.array(sdict[img]) - pos
        
        # Replace total shift with delta in shift dictionary 
        sdict[img] = tuple(delta_pos)

# Build table rows from input data
def _buildTableRows(filelist, shiftdict, dthprod,suffix=None):
    """
        Create dictionary with each entry corresponding to a list of
        values for each column.  
        
    """

    xsh,ysh,rot,name,mtype,mprsnt = [],[],[],[],[],[]
    ot = ord('T')
    of = ord('F')
    
    for file in filelist:
        try:
            xshift = float(shiftdict[file][0])
        except KeyError:
            xshift = 0.0
        xsh.append(xshift)
        
        try:
            yshift = float(shiftdict[file][1])
        except KeyError:    
            yshift = 0.0
        ysh.append(yshift)
        
        try:
            rotation = float(shiftdict[file][2])
        except KeyError:
            rotation = 0.0
        rot.append(rotation)
        
        # Build rootname for file
        name.append(_buildNewRootname(file,suffix=suffix))
        mprsnt.append(ot)
        mtype.append("EXP-DTH")
        
    # Now add product row, based on output
    name.append(dthprod)
    xsh.append(0.0)
    ysh.append(0.0)
    rot.append(0.0)
    mtype.append("PROD-DTH")
    mprsnt.append(of)
    table = {'xsh':xsh,'ysh':ysh,'rot':rot,'name':name,'mtype':mtype,'mprsnt':mprsnt}

    return table

def _buildNewRootname(filename,extn=None,suffix=None):
    """ Build rootname for a new file.
        Use 'extn' for new filename if given, otherwise
        does NOT append a suffix at all.
        Search for suffix if given.
                
        Does NOT check to see if it exists already.  
        Will ALWAYS return a new filename.
    """
    # Search known suffixes to replace ('_crj.fits',...)
    extlist = EXTLIST
    # Also, add a default where '_dth.fits' replaces
    # whatever extension was there ('.fits','.c1h',...)
    extlist.append('.')
    if suffix != None:
        extlist.append(suffix)
        
    for suffix in extlist:
        _indx = string.find(filename,suffix)
        if _indx > 0: break

    if _indx < 0:
         # default to entire rootname
        _indx = len(filename)
    
    if extn != None:
        output = filename[:_indx]+extn
    else:
        output = filename[:_indx]
        
    return output
   
# Make a Table extension HDU
def _makeTableHDU(data):
    
    # Syntax used with PyFITS 0.4.2:
    #_format = "s24,s14,s6,f32,f32,f32"
    #recData = pyfits.record.record(data, _format)
    # Do this when the Table classes are updated
    #hdu = pyfits.Table(recData, name = "ASN")
    
    col1 = pyfits.Column(name='MEMNAME',format='24a',array=chararray.array(data['name']))
    col2 = pyfits.Column(name='MEMTYPE',format='14a',array=chararray.array(data['mtype']))
    col3 = pyfits.Column(name='MEMPRSNT',format='1b',array=numpy.array(data['mprsnt']).astype(numpy.uint8))
    col4 = pyfits.Column(name='XOFFSET',format='1r',unit='arcsec',array=numpy.array(data['xsh']))
    col5 = pyfits.Column(name='YOFFSET',format='1r',unit='arcsec',array=numpy.array(data['ysh']))
    col6 = pyfits.Column(name='ROTATION',format='1r',unit='degrees',array=numpy.array(data['rot']))
    
    hdu = pyfits.new_table([col1,col2,col3,col4,col5,col6],nrows=len(data['name']))

    return hdu
    
def _buildAsnPrimary(fasn,output,img1):
    """ Creates complete Primary header for ASN table based on template.
    'fasn' is the file handle of the newly opened ASN table file.
    Uses keyword values from 'img1' to update ASN header keywords.
    """
    origin_str = 'PyFITS Version '+pyfits.__version__
    
    # Format time values for keywords IRAF-TLM, and DATE 
    _ltime = time.localtime(time.time())
    tlm_str = time.strftime('%H:%M:%S (%d/%m/%Y)',_ltime)
    date_str = time.strftime('%Y-%m-%dT%H:%M:%S',_ltime)
    
    
    # Build PRIMARY HDU
    _hdu = pyfits.PrimaryHDU(header=_prihdr)
    fasn.append(_hdu)
    newhdr = fasn['PRIMARY'].header
    # Verify that input image is a FITS file...
    try:
        # Open img1 to obtain keyword values for updating template
        fimg1 = pyfits.open(img1)
        prihdr = fimg1['PRIMARY'].header
        _updateKeyword('INSTRUME',prihdr,newhdr)
        _updateKeyword('PRIMESI',prihdr,newhdr)
        _updateKeyword('TARGNAME',prihdr,newhdr)
        _updateKeyword('DETECTOR',prihdr,newhdr)
        _updateKeyword('RA_TARG',prihdr,newhdr)
        _updateKeyword('DEC_TARG',prihdr,newhdr)
        # All done with input image. Close it now.
        fimg1.close()
        del fimg1
    except:
        pass
        
    # Update Primary header with values from input image
    newhdr['IRAF-TLM']=tlm_str
    newhdr['DATE'] = date_str
    newhdr['ORIGIN'] = origin_str
    _indx = string.find(output,'.')
    if _indx < 1:
        _indx = len(output)
    newhdr['ROOTNAME'] = output[:_indx]
    newhdr['FILENAME'] = output
    newhdr['ASN_ID'] = output[:_indx]
    newhdr['ASN_TAB'] = output

def help():
    _str = """ Create an association table from all suitable files in
    the current directory selected by their suffix. It returns the
    filename of the final product for the table, created by appending
    '_drz.fits' to the given rootname.
        
    The columns XOFFSET, YOFFSET and ROTATION are added to this table
    and populated with values from a shiftfile, if one was specified 
    and can be read. It should be in the format used by 'dither.shiftfind'
    or 'dither.avshift' from the STSDAS DITHER package.  These values should
    be DELTA offsets in units of arcseconds in RA and Dec.
    
    SYNTAX:
        dthprod = buildAsn.buildAsnTable(rootname,suffix=None,shiftfile=None)
      where
        rootname  - user-specified string for creating the output table name
        suffix    - look for files with these user-specified suffixes
                    as inputs for table (DEFAULT: crj,sfl,and flt.fits)
        shiftfile - name of file containing shifts to be used by PyDrizzle
        dthprod   - full filename of final product for association table
    EXAMPLE: 
      Create ASN table 'mymosaic_asn.fits' from all dth.fits images.
        >>> import buildAsn
        >>> dthprod = buildAsn.buildAsnTable('mymosaic',suffix='flt.fits')  
    """
    
    print _str
    print 'BuildAsn Version '+__version__
