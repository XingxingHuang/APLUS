#!/usr/bin/env python
# W Zheng  
import glob,math
from optparse import OptionParser
import pyfits
import mkasntable as buildasn
import os, sys, string, pdb
import shutil
from sets import Set
# import multidrizzle 
from   xydrizzle import drutil
from   pytools import wcsutil
import numpy
def setRange(list): #WZ, AKS
        """ Determine the shape of output image, and update ref parameters
        """
        ra, dec = [], []
	pscale=0.065
        for name in list:
		im = name + '[1]'
		wcs = wcsutil.WCSObject(im)
                xy = [0., 0.]
                r, d = wcs.xy2rd(xy)
                ra.append(r)
                dec.append(d)

                xy = [0., wcs.naxis2 + 1.]
                r, d = wcs.xy2rd(xy)
                ra.append(r)
                dec.append(d)

                xy = [wcs.naxis1 + 1., 0.]
                r, d = wcs.xy2rd(xy)
                ra.append(r)
                dec.append(d)

                xy = [wcs.naxis1 + 1., wcs.naxis2 + 1.]
                r, d = wcs.xy2rd(xy)
                ra.append(r)
                dec.append(d)        
                rmax, dmax = max(ra), max(dec)
                rmin, dmin = min(ra), min(dec)
                #pdb.set_trace()

        #Find the maximun and minimum ra/dec values for the 
        #current set of images
        	
        rmax, dmax = max(ra), max(dec)
        rmin, dmin = min(ra), min(dec)
        #Get image center and use the maximum
        #distance from it to the edge to calculate final image's dimensions
	ref = wcs.copy()
	ref.orient = 0.
	ref.crval1=(rmax+rmin)/2.
	ref.crval2=(dmax+dmin)/2.
	#pdb.set_trace()
        dx = abs(rmax - rmin)*3600/pscale
        dy = abs(dmax - dmin)*3600/pscale
        X = int(math.ceil(dx / 100.) * 100)
        Y = int(math.ceil(dy / 100.) * 100)        
        ref.crpix1 = X/2.
        ref.crpix2 = Y/2.
        ref.naxis1 = X
        ref.naxis2 = Y
	ref.cd11=-pscale/3600.
	ref.cd12=0.
	ref.cd21=0.
	ref.cd22=pscale/3600.

        return ref # [X, Y]

usage = 'usage: %prog dataset_name input_dir'
version = '%prog 0.0'
parser = OptionParser(usage=usage, version=version)
(options, args) = parser.parse_args()
if len(args) != 2:
    parser.print_help()
    sys.exit()
dsname = args[0]
indir = args[1]
if not os.path.exists(args[1]):
    parser.error('Directory %s not found!' % (args[1]))
    sys.exit()
fimgs = glob.glob(os.path.join(indir, '*flt.fits'))
acfilts = [] 
uvfilts = [] 
irfilts = []
acimgs = [] 
uvimgs = [] 
irimgs = []
# get all filters
rootpath=os.getcwd()
for f in fimgs:
    fobj = pyfits.open(f,'update')
    hdr = fobj[0].header
    det = hdr['detector']
    # pa = hdr['PA']
    if 'WFC' in det:
        flt = hdr['filter1']
        if 'CLEAR' in flt:
            flt = hdr['filter2']
	fobj[0].header.update('FILTER',flt)
        acimgs.append(f)
        acfilts.append(flt)
    else:
        flt = hdr['filter']        
	#pdb.set_trace()
	# fobj[0].header.update('FILTER1',flt)
	# fobj[0].header.update('FILTER2','None')
	# del rmsfitsobj[0].header.ascard["EXTEND"]
        if 'IR' in det:
            irimgs.append(f)
            irfilts.append(flt)
        else:
            if 'UVIS' in det:
                uvimgs.append(f)
                uvfilts.append(flt)
    fobj.close()
print len(fimgs),' input images'
ref=setRange(fimgs)
parfile = open("ref.par","w") 
#format=str('%4d %4d %10.7f %10.7f')
#print >> parfile, format % (ref.naxis1,ref.naxis2,ref.crval1,ref.crval2)
str="NAXIS1 "+str(ref.naxis1)
parfile.write(str)
#print "NAXIS1 ", ref.naxis1 >> parfile
#print "NAXIS2 ", ref.naxis2 >> parfile
#print "CRPIX1 ", ref.crpix1 >> parfile
#print "CRPIX2 ", ref.crpix2 >> parfile
#print "CRVAL1 ", ref.crval1 >> parfile
#print "CRVAL22 ", ref.crval2 >> parfile
#% (ref.naxis1,ref.naxis2,ref.crval1,ref.crval2)
parfile.close()
#pdb.set_trace()
#refimg='null.fits'
#fitsobj = pyfits.HDUList()
#fitsobj.append(pyfits.PrimaryHDU())
#temfits=pyfits.open(f)
#temfits.close()
#fits[0].header.update('BITPIX',16)
#flgfits[0].header.update('NAXIS1',ref.naxis1)
#flgfits[0].header.update('NAXIS2',ref.naxis2)
#fitsobj[0].header = temfits[1].header
#data = numpy.zeros(shape=(ref.naxis1,ref.naxis2))
#fitsobj[0].data   = data
#if os.path.isfile(refimg):
#    os.remove(refimg)
#fitsobj.writeto(refimg)
if (len(acfilts)>0):
    print "Working on ACS images"
uflts = Set(acfilts)
ff = zip(acimgs, acfilts)
# build asn for each filter
for uf in uflts:
    asns = []
    name = '%s_%s' % (dsname, uf)
    fname=name+'_asn.fits'
    if os.path.exists(fname):
        os.remove(fname)
    #pdb.set_trace()
    subdir=os.path.join(os.getcwd(),uf.lower())
    print subdir
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for f in ff:
        fn, flt = f
        if flt == uf:
            shutil.copy(fn,subdir)
            asns.append(os.path.basename(fn))
    buildasn.buildAsnTable(name, flist=asns, asntype='DRIZZLE')
    shutil.move(fname,subdir)
if (len(irfilts)>0) :
    print "Working on IR images"
uflts = Set(irfilts)
ff = zip(irimgs, irfilts)
# build asn for each filter
for uf in uflts:
    asns = []
    name = '%s_%s' % (dsname, uf)
    fname=name+'_asn.fits'
    if os.path.exists(fname):
        os.remove(fname)
    subdir=os.path.join(os.getcwd(),uf.lower())
    print subdir
    #pdb.set_trace()
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for f in ff:
        fn, flt = f
        if flt == uf:
            shutil.copy(fn,subdir)
            asns.append(os.path.basename(fn))
    buildasn.buildAsnTable(name, flist=asns, asntype='DRIZZLE')
    shutil.move(fname,subdir)
if (len(uvfilts)>0):
    print "Working on UVIS images"
uflts = Set(uvfilts)
ff = zip(uvimgs, uvfilts)
# build asn for each filter
for uf in uflts:
    asns = []
    name = '%s_%s' % (dsname, uf)
    fname=name+'_asn.fits'
    if os.path.exists(fname):
        os.remove(fname)
    subdir=os.path.join(os.getcwd(),uf.lower())
    print subdir
    #pdb.set_trace()
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for f in ff:
        fn, flt = f
        if flt == uf:
            #shutil.move(fn,subdir)
            shutil.copy(fn,subdir)
            asns.append(os.path.basename(fn))
    buildasn.buildAsnTable(name, flist=asns, asntype='DRIZZLE')
    shutil.move(fname,subdir)

print "Done"

