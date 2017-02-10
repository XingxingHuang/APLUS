"""
General file functions designed for use with PyDrizzle.
These functions only rely on IRAF booleans 'yes' and 'no' and PyFITS.
"""
import pyfits
import string,os, pdb
import numpy as N
from pyraf import iraf

yes = iraf.yes
no = iraf.no

# Required keywords to build an extension with WCS info
DTH_KEYWORDS=['CD1_1','CD1_2', 'CD2_1', 'CD2_2', 'CRPIX1', 
'CRPIX2','CRVAL1', 'CRVAL2', 'CTYPE1', 'CTYPE2']

#################
#
#
#		Generic Functions
#
#
#################		
def DEGTORAD(deg):
	return (deg * N.pi / 180.)

def RADTODEG(rad):
	return (rad * 180. / N.pi)

def DIVMOD(num,val):
    if isinstance(num,N.ndarray):
        # Treat number as numpy object
        _num = N.remainder(num,val)
    else:
        _num = divmod(num,val)[1]
    return _num 


#
#Revision History:
#    Nov 2001: findFile upgraded to accept full filenames with paths, 
#               instead of working only on files from current directory. WJH
#
# Base function for  
#   with optional path.

def findFile(input):

    """ Search a directory for full filename with optional path. """
    
    _fdir,_fname = os.path.split(input)

    if _fdir == '':
        _fdir = os.curdir
        
    flist = os.listdir(_fdir)

    found = no
    for name in flist:
        if not string.find(name,_fname):
            found = yes
            continue
    return found
    

def checkFileExists(filename,directory=None):
    """ Checks to see if file specified exists in current
        or specified directory. Default is current directory.
        Returns 1 if it exists, 0 if not found.
    """
    if directory == None: directory = '.'
    _ldir = os.listdir(directory)

    _exist = 0
    # for each file in directory...
    for file in _ldir:
        # compare filename with file
        if string.find(file,filename) > -1:
            _exist = 1
            break
            
    return _exist


def copyFile(input,output,replace=None):

    """ Copy a file whole from input to output. """

    _found = findFile(output)
    if not _found or (_found and replace ):
        finput = open(input,'r')
        fout = open(output,'w')
        fout.writelines(finput.readlines())
        fout.close()
        finput.close()
        del finput
        del fout


def findKeywordExtn(ft,keyword,value):

    """
    This function will return the index of the extension in 
    a multi-extension FITS file which contains the desired keyword 
    with the given value.
    """

    i = 0
    extnum = 0
    # Search through all the extensions in the FITS object
    for chip in ft:
        hdr = chip.header
        # Check to make sure the extension has the given keyword
        if hdr.has_key(keyword):
            # If it does, then does the value match the desired value
            # MUST use 'string.strip' to match against any input string!
            if string.strip(hdr[keyword]) == value:
                extnum = i
                break
        i = i + 1
    # Return the index of the extension which contained the 
    # desired EXTNAME value.
    return extnum

    
def findExtname(fimg,extname,extver=None):
    """ This function returns the list number of the extension
	    corresponding to EXTNAME given.  
    """
    i = 0
    extnum = 0
    for chip in fimg:
        hdr = chip.header
        if hdr.has_key('EXTNAME'):
            if string.strip(hdr['EXTNAME']) == string.upper(extname):
                if extver == None or hdr['EXTVER'] == extver:
                    extnum = i
                    break
        i = i + 1
    return extnum

def rAsciiLine(ifile):

    """ Returns the next non-blank line in an ASCII file. """
    
    _line = string.strip(ifile.readline())
    while len(_line) == 0:
        _line = string.strip(ifile.readline())
    return _line

# This function read the IDC table and generates the two matrices with
# the geometric correction coefficients.  
#
# 	INPUT: FITS object of open IDC table
# 	OUTPUT: coefficient matrices for Fx and Fy
#
# def readIDCtab (tabname, chip=1, direction='forward'):
# get this function from matutil now...
def OBSOLETE_readIDCtab (tabname, chip=1, direction='forward'):

    try:
        ftab = pyfits.open(tabname)
    except:
        raise IOError,"IDC table '%s' not valid as specified!" % tabname

    #First thing we need, is to read in the coefficients from the IDC
    # table and populate the Fx and Fy matrices.

    # Read FITS header to determine order of fit, i.e. k
    order = ftab['PRIMARY'].header['NORDER']

    fx = N.zeros(shape=(order+1,order+1),dtype=N.float32)
    fy = N.zeros(shape=(order+1,order+1),dtype=N.float32)

    #Determine row from which to get the coefficients.
    # How many rows do we have in the table...
    fshape = ftab[1].data.shape
    colnames = ftab[1].data._names
    row = -1

    # Loop over all the rows looking for the one which corresponds
    # to the value of CCDCHIP we are working on...
    for i in xrange(fshape[0]):
        if 'DETCHIP' in colnames:
            detchip = ftab[1].data.field('DETCHIP')[i]
        else:
            detchip = 1

	# Below is the pydirzzle 3.0 solution to the problem of SBC idctabs
	# having a detchip of -999 in SBC idctabs.  Better would be to explicity
	# set the detchip to -999 if the detector is 'SBC'.
	# Below is the old, legacy pydrizzle 2.6 code

	if 'DIRECTION' in colnames:			
            direct = string.lower(ftab[1].data.field('DIRECTION')[i])
        else:
            direct = 'forward'

        if string.strip(direct) == string.strip(direction):
            if int(detchip) == int(chip) or int(detchip) == -999:
                row = i
		break

#         if 'DIRECTION' in colnames:			
#             direct = string.lower(ftab[1].data.field('DIRECTION')[i])
#         else:
#             raise LookupError,'fileutil: readIDCtab did not find valid DIRECTION'

#         if (string.strip(direct) == string.strip(direction)) and (int(detchip) == int(chip)):
#                 row = i
# 		print '_readIDCtab breaking on row ',i, 'detchip,direction: ',detchip,direct
# 		break
	
    if row < 0:
        print 'Row corresponding to DETCHIP of ',chip,' was not found!'
        raise LookupError

    refpix = {}
    refpix['XREF'] = ftab[1].data.field('XREF')[row]
    refpix['YREF'] = ftab[1].data.field('YREF')[row]
    refpix['XSIZE'] = ftab[1].data.field('XSIZE')[row]
    refpix['YSIZE'] = ftab[1].data.field('YSIZE')[row]
    refpix['PSCALE'] = ftab[1].data.field('SCALE')[row]
    refpix['V2REF'] = ftab[1].data.field('V2REF')[row]
    refpix['V3REF'] = ftab[1].data.field('V3REF')[row]
    refpix['XDELTA'] = 0.0
    refpix['YDELTA'] = 0.0
    refpix['centered'] = no

    # Now that we know which row to look at, read coefficients into the
    #	numeric arrays we have set up... 
    # Setup which column name convention the IDCTAB follows
    # either: A,B or CX,CY
    if 'CX10' in ftab[1].data._names:
        cxstr = 'CX'
        cystr = 'CY'
    else:
        cxstr = 'A'
        cystr = 'B'
    
    for i in xrange(order+1):
        if i > 0:
            for j in xrange(i+1):
                xcname = cxstr+str(i)+str(j)
                ycname = cystr+str(i)+str(j)
                fx[i,j] = ftab[1].data.field(xcname)[row] 
                fy[i,j] = ftab[1].data.field(ycname)[row] 

    ftab.close()			
    del ftab
    
    # Return arrays and polynomial order read in from table.
    # NOTE: XREF and YREF are stored in Fx,Fy arrays respectively.
    return fx,fy,refpix,order

def readAsnTable(fname,output,prodonly=yes):
    """
     This function reads the filenames/rootnames and shifts from a FITS
     ASN table. 

     Column names expected are:
       MEMNAME     - rootname of each member
       MEMTYPE     - type of member in association(PROD-* or EXP-*)
       XOFFSET     - shift in X for this observation
       YOFFSET     - shift in Y for this observation
       ROTATION    - additional rotation to be applied 
       SCALE	   - scale image by this value

     This will return a nested dictionary corresponding to this
     association table.
     Observation dictionary: {'xsh':0.,'ysh':0.,'rot':0.,'scale':1.}
     Product dictionary: {'output':type, 'memname1':dict, 'memname2':dict,...}
	       where dict: Observation dictionary
     The dictionary using EXP* will be:
       p = {'output':'dthname','members':{
		    'name1':{'xsh':0.,...}, 'name2':{'xsh':0.,...},...
		       }
		    }
     You get a list of input names using 'p.keys()'.

     Parameters:
	       output:     output - desired name of output image (None,'EXP', or user-filename)
	       prodonly:   yes - only select MEMTYPE=PROD* as input observations
				       no - use MEMTYPE=EXP* as input observations

     Output: If none is specified by user, the first 'PROD-DTH' filename found in the 
 		    ASN table will be used.  If there is no 'PROD-DTH' entry, the first 'PROD' 
		    entry will be used instead.  Finally, if 'output' = 'EXP', the value 'EXP'
		    will be passed along and interpreted as a switch to use the input filename
		    as the output resulting in every input creating a separate output.
    """
    # Initialize this dictionary for output
    asndict = {'output':None,'members':{}}

    #print "Read ASN table"
    #pdb.set_trace()
    # Open the table...
    try:
        ftab = pyfits.open(fname)
    except:
        raise IOError,"Association table '%s' not valid as specified!" % fname

    tablen = ftab[1].data.shape[0]
    colnames = ftab[1].data._names
    # Set a flag to specify whether ASN has a PROD-DTH member
    dthprod = no
    # Now, put it together with rootname and type...
    for row in xrange(tablen):
        # Read in required columns for each row
        if 'MEMNAME' in colnames and 'MEMTYPE' in colnames:
            # We need to make sure no EOS characters remain part of
            # the strings we read out...
            mname = string.split(ftab[1].data.field('MEMNAME')[row],'\0',1)[0]
            mtype = string.split(ftab[1].data.field('MEMTYPE')[row],'\0',1)[0]
            memname = string.strip(mname)
            memtype = string.strip(mtype)
            memrow = row
        else:
            print 'Association table incomplete: required column(s) MEMNAME/MEMTYPE NOT found!'
            raise LookupError

        # Do we care about this entry?
        # Entries that should be used to build DTH product are:
        #  PROD-RPT, PROD-CRJ, EXP-DTH
        if string.find(memtype,'PROD') < 0 and string.find(memtype,'EXP-DTH') < 0:
            if prodonly == yes:
                # We are looking at an EXP* entry we don't want...
                continue

        memdict = {}
        # Keep track of which order they were read in by their row number
        memdict['row'] = memrow
        memdict['xoff'] = 0.
        memdict['yoff'] = 0.
        memdict['rot'] = 0.
        # Read in optional data from columns
        # X offset
        if 'XOFFSET' in colnames:
            memdict['delta_x'] = ftab[1].data.field('XOFFSET')[row]
        else:
            memdict['delta_x'] = 0.

        # Y offset
        if 'YOFFSET' in colnames:			
            memdict['delta_y'] = ftab[1].data.field('YOFFSET')[row]
        else:
            memdict['delta_y'] = 0.

        # Rotation angle
        if 'ROTATION' in colnames:			
            memdict['delta_rot'] = ftab[1].data.field('ROTATION')[row]
        else:
            memdict['delta_rot'] = 0.

        # Scale: output pixel size 
        if 'SCALE' in colnames:			
            memdict['scale'] = ftab[1].data.field('SCALE')[row]
        else:
            memdict['scale'] = 1.


        # Build the shifts dictionary now...
        if string.find(memtype,'PROD') < 0 and prodonly == no:
            # We want to use this EXP* entry.
            asndict['members'][memname] = memdict
        elif memtype == 'PROD-DTH':
            # We have found a dither product specified in ASN
            # Not to be used for input, but
            # has one already been specified as the final output?
            if dthprod == no:
                if output == None:
                    # Use default output name
                    asndict['output'] = memname
                else:
                    # Use user-specified output name here
                    asndict['output'] = output
                dthprod = yes
        else:
            # We are working with a PROD* entry...
            if prodonly == yes:
                asndict['members'][memname] = memdict

        # Set up a default output filename 
        # This will be overwritten by a different output
        # name if a PROD-DTH entry is found in the ASN table
        # and 'output' was not specified by the user.
        # Useful for CR-SPLIT/REPEAT-OBS ASN tables.
        if asndict['output'] == None:
            if output == None:
                asndict['output'] = memname
            else:
                asndict['output'] = output

    # Finished reading all relevant rows from table
    ftab.close()
    del ftab

    return asndict

def buildDthProduct (pardict, output, extlist=None, outlist=None, wcs=None):
    """ 
    Parameters:
        pardict - a dictionary containing (at least):
            'data','outdata','outweight','outcontext'
            where 'data' serves as the name of the backup template image,
                and the others are filenames of the drizzle products.
        output - filename of the final combined output image
        extlist - list of EXTNAME's to be searched in the template image
        outlist - list of EXTNAME's to be used for naming the output extensions
                  
    This function will package the two or three output files from
    'drizzle' or 'tdrizzle' into a single output multi-extension 
    FITS file.
    It uses a pre-existing multi-extension FITS file as a template
    for the output product. 
    The list 'extlist' contains the names of the extensions from the
    template image which will be used to build the output file. 
    
    A number of keywords are also updated based on values from the
    drizzle products, including the CD matrix, CRPIX[1,2], CRVAL[1,2].
    In addition, the template images will be of different size than
    the drizzle products, yet this is handled automatically by PyFITS.
    
    NOTE:
    The output file will ALWAYS contain 'SCI','WHT', and 'CTX' 
    extensions.
    """
    # Set up default extlist list
    if extlist == None:
        extlist = ('SCI','ERR','DQ')
    if outlist == None:
        outlist = ('SCI','WHT','CTX')

    # Get default headers from multi-extension FITS file
    # If input data is not in MEF FITS format, it will return 'None'
    # and those headers will have to be generated from drizzle output
    # file FITS headers.
    # NOTE: These are HEADER objects, not HDUs
    prihdr,scihdr,errhdr,dqhdr = getTemplates(output,pardict['data'],extlist)
    
    if prihdr == None:
        # Open FITS image and use its Primary header
        fpri = pyfits.open(pardict['outdata'])
        prihdr = pyfits.Header(cards=fpri[0].header.ascard.copy())
        fpri.close()
        del fpri

    # Setup primary header as an HDU ready for appending to output FITS file
    prihdu = pyfits.PrimaryHDU(header=prihdr,data=None)
    
    # Start by updating PRIMARY header keywords...
    prihdu.header.update('EXTEND',pyfits.TRUE)
    prihdu.header.update('NEXTEND',3)
    prihdu.header.update('FILENAME', output)
    
    # Open the dither output SCI image	
    fsci = pyfits.open(pardict['outdata'])

    # Get the total exposure time for the image
    # If not calculated by PyDrizzle and passed through
    # the pardict, then pulled from the template image.
    inhdr = fsci[0].header
    
    if pardict.has_key('texptime'):
        _exptime = pardict['texptime']
        _expstart = pardict['expstart']
        _expend = pardict['expend']
    else:
        _exptime = inhdr['EXPTIME']
        _expstart = inhdr['EXPSTART']
        _expend = inhdr['EXPEND']
    
    prihdu.header.update('EXPTIME', _exptime)
    prihdu.header.update('TEXPTIME',_exptime)
    prihdu.header.update('EXPSTART', _expstart)
    prihdu.header.update('EXPEND', _expend)
    
    # Update DITHCORR calibration keyword if present
    # Remove when we can modify FITS headers in place...
    if prihdu.header.has_key('DITHCORR') > 0:
        prihdu.header['DITHCORR'] = 'COMPLETE'

    # Now, build the output file
    fo = pyfits.open(output,'append')
    # Add primary header to output file...
    fo.append(prihdu)    

    # Now, build SCI extension HDU
    if scihdr == None:
        scihdr = inhdr

    if wcs != None:
        # Update ORIENTAT based on PyDrizzle product's value
        # since 'drizzle' itself doesn't update that keyword.
        scihdr.update('ORIENTAT',wcs.orient)
    # Does this need to be explicitly created or can the pre-existing
    # one simply be appended?
    #       
    hdu = pyfits.ImageHDU(data=fsci[0].data,header=scihdr,name=outlist[0])
    updateDTHKeywords(hdu.header,fsci[0].header,output)

    fo.append(hdu)
    fo.flush()
    
    fsci.close()
    del fsci[0].data
    del fsci
    del hdu
    
    # Write out the WEIGHT image in the ERR array extension
    fweight = pyfits.open(pardict['outweight'])
    if errhdr == None:
        errhdr = fweight[0].header
    hdu = pyfits.ImageHDU(data=fweight[0].data,header=errhdr,name=outlist[1])
    updateDTHKeywords(hdu.header,fweight[0].header,output)

    fo.append(hdu)
    fo.flush()

    fweight.close()
    del fweight[0].data
    del fweight
    del hdu
    
    # Write out the Context image (if any was created)
    cfile = pardict['outcontext']
    _ctx = yes
    if findFile(cfile):
        fctx = pyfits.open(pardict['outcontext'])
        if dqhdr == None:
            dqhdr = fctx[0].header
        hdu = pyfits.ImageHDU(data=fctx[0].data,header=dqhdr,name=outlist[2])
        updateDTHKeywords(hdu.header,fctx[0].header,output)

    else:
        _ctx = no
        
        # Use the SCI HDU for the shape and build a default array
        imarr = N.ones(shape=fo[1].data.shape,dtype=N.int16)
        if dqhdr == None:
            dqhdr = scihdr
            dqhdr.update('EXTVER',scihdr['EXTVER'])
            
        hdu = pyfits.ImageHDU(data=imarr,header=dqhdr,name=outlist[2])
        print 'Dither Product: writing out empty context extension.'

    fo.append(hdu)
    # Close the output and template file
    print 'Finished creating FINAL dither product ',output	

    fo.close()
    del fo[1].data
    del fo
    if _ctx:
        fctx.close()
        del fctx[0].data
        del fctx

def getTemplates(oname,tname,extlist):
    # Obtain default headers for output file
    # If the output file already exists, use it 
    # If not, use an input file for this information.
    #
    # NOTE: Returns 'pyfits.Header' objects, not HDU objects!
    #
    fname = None
    if checkFileExists(oname):
        fname = oname
    else:
        fname = tname

    if fname != None and string.find(fname,'.fits') > 0:
        # Open an calibrated ACS image as a template
        _indx = string.find(tname,'[')
        if _indx > 0:
            template = tname[:_indx]
        else:
            template = tname
        ftemplate = pyfits.open(template)
        
        # Setup which keyword we will use to select each 
        # extension...
        _extkey = 'EXTNAME'

        #
        # Now, extract the headers necessary for output (as copies)
        # 1. Find the SCI extension in the template image
        # 2. Make a COPY of the extension header for use in new output file
        prihdr = pyfits.Header(cards=ftemplate['PRIMARY'].header.ascard.copy())
        extnum = findKeywordExtn(ftemplate,_extkey,extlist[0])
        scihdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())
        
        extnum = findKeywordExtn(ftemplate,_extkey,extlist[1])	
        errhdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())
        extnum = findKeywordExtn(ftemplate,_extkey,extlist[2])	
        dqhdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())

        ftemplate.close()
        del ftemplate

    else:
        # Create default headers from scratch
        prihdr = None
        scihdr = None
        errhdr = None
        dqhdr = None
    
    # Now, safeguard against having BSCALE and BZERO
    try:
        del scihdr['bscale']
        del scihdr['bzero']
        del errhdr['bscale']
        del errhdr['bzero']
        del dqhdr['bscale']
        del dqhdr['bzero']
    except:
        # If these don't work, they didn't exist to start with...
        pass
            

    # At this point, check errhdr and dqhdr to make sure they 
    # have all the requisite keywords (as listed in updateDTHKeywords).
    # Simply copy them from scihdr if they don't exist...
    if errhdr != None and dqhdr != None:
        for keyword in DTH_KEYWORDS:
            if not errhdr.has_key(keyword):
                errhdr.update(keyword,scihdr[keyword])
            if not dqhdr.has_key(keyword):
                dqhdr.update(keyword,scihdr[keyword])
            
    return prihdr,scihdr,errhdr,dqhdr    

def updateDTHKeywords(hdr,extnhdr,filename):
    """ Update header keywords in output to reflect values from
     the dither product.
    """
    for keyword in DTH_KEYWORDS:
        hdr.update(keyword,extnhdr[keyword])
    hdr.update('EXTVER', 1)
    
