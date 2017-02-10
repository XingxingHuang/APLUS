""" mkasntable.py - Module to generate an ACS-style association
table from all matching files from the current working directory.


DKM 30 Apr 2002:
    Hacked from buildasn.py by W. Hack
    Modified buildAsnTable to accept a list of images and take an
    association type argument
    Dumped the shift stuff
    Renamed mkasntable.py
"""

import os, re,string, pdb
import time
import pyfits
import numpy
from numpy import char as chararray


# List of supported default file types
# It will look for these file types by default
# when trying to recognize input rootnames.
EXTLIST =  ['_crj.fits','_flt.fits','_flc.fits','_sfl.fits']

__version__ = '1.0 (15-Jun-2012)'

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

def buildAsnTable(asnroot,suffix=None,flist=None,asntype=None):
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
    #print 'BuildAsn Version '+__version__
    
    # Set default suffix
    if suffix == None:
        suffix = 'crj.fits'
    
    # Find all files in current directory with given suffix
    if flist == None:
	    flist = _findFiles(suffix)

    if len(flist) == 0: 
        print 'NOTE: Please check suffix used to build list of files.'
        _errmsg = 'Could NOT find any files to build into an association!'
        raise Exception,_errmsg
    
    atypes = ['DRIZZLE', 'CRJFLT','CALACS']
    
    if asntype not in atypes:
        print 'ERROR: Association type not found'
        _errmsg = 'Undefined association type: '+asnytpe
        raise Exception, _errmsg
        
    _output = asnroot+'_asn.fits'
    
    # Open output ASN table file
    fasn = pyfits.HDUList()

    # Build ASN file primary header
    _buildAsnPrimary(fasn,_output,flist[0])
        
    # Make table records 
    outTable = _buildTableRows(flist, asnroot, asntype, suffix)
    
    # Create an extension HDU which contains a table
    exthdu = _makeTableHDU(outTable)

    fasn.append(exthdu)
    
    # Close ASN table file
    fasn.writeto(_output)
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
    _reg = re.compile(regpatt)
    
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
    

# Build table rows from input data
def _buildTableRows(filelist, dthprod, asntype, suffix=None):
    """
        Create dictionary with each entry corresponding to a list of
        values for each column.  
        
    """

    xsh,ysh,rot,name,mtype,mprsnt = [],[],[],[],[],[]
    ot = ord('T')
    of = ord('F')
    
    for file in filelist:
        xshift = 0.0
        xsh.append(xshift)
        
        yshift = 0.0
        ysh.append(yshift)
        
        rotation = 0.0
        rot.append(rotation)
        
        # Build rootname for file
        name.append(_buildNewRootname(file,suffix=suffix))
        # DKM 2002-04-30: added support for flt and raw fits files
        fsuffix = file[-8:]
        if fsuffix == "crj.fits":
            fslist = [f for f in filelist if f.find("crj.fits") != -1]
            memtype = "PROD-CR"+str(fslist.index(file)+1)
            mprsnt.append(ot)
            mtype.append(memtype)
        elif fsuffix == "flt.fits":
            fslist = [f for f in filelist if f.find("flt.fits") != -1]
            memtype = "EXP-DTH"+str(fslist.index(file)+1)
            mprsnt.append(ot)
            mtype.append(memtype)
        elif fsuffix == "flc.fits": #WZ
            fslist = [f for f in filelist if f.find("flc.fits") != -1]
            memtype = "EXP-DTH"+str(fslist.index(file)+1)
            mprsnt.append(ot)
            mtype.append(memtype)
        else:
            if len(filelist) == 1:
                mprsnt.append(ot)
                mtype.append("EXP-DTH")
            else:
                fslist = [f for f in filelist if f.find("raw.fits") != -1]
                memtype = "EXP-CRJ"+str(fslist.index(file)+1)
                mprsnt.append(ot)
                mtype.append(memtype)
        
    # Now add product row, based on output
    name.append(dthprod)
    xsh.append(0.0)
    ysh.append(0.0)
    rot.append(0.0)
    if asntype == 'CALACS':
        mtype.append("PROD-CRJ")
        mprsnt.append(ot)
    else:
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
    
    col1 = pyfits.Column(name='MEMNAME',format='35A',array=chararray.array(data['name']))
    col2 = pyfits.Column(name='MEMTYPE',format='14A',array=chararray.array(data['mtype']))
    #pdb.set_trace()
    col3 = pyfits.Column(name='MEMPRSNT',format='1B',array=numpy.array(data['mprsnt']).astype(numpy.int8))
    col4 = pyfits.Column(name='XOFFSET',format='E',unit='arcsec',array=numpy.array(data['xsh']))
    col5 = pyfits.Column(name='YOFFSET',format='E',unit='arcsec',array=numpy.array(data['ysh']))
    col6 = pyfits.Column(name='ROTATION',format='E',unit='degrees',array=numpy.array(data['rot']))
    
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
