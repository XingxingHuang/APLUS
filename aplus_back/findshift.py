#!/usr/bin/env python
# $Id: combDither.py,v 2.1.0 2011/02/04 06:50:33 wz Exp $
# Accommodate different input names
# ---------------------------------------------------------------------

__version__      = '$Revision: 2.1 $ '[11:-3]
__version_date__ = '$Date: 2011/02/04  $ '[7:-3]
__author__       = 'Wei Zheng <zheng@pha.jhu.edu>'
# change to "drizzle"

import os,sys,string,glob,math
import numpy
import pyfits,matutil
import drizzle
import xmlUtil,fUtil,pyblot,augmask,tableio #WZ
import astrometer
import pdb, popen2 #WZ
from   drizzle import wcsutil,drutil #WZ
from   pUtil import ptime
from   msg   import pMessage
from   sys   import version
pyversion = version
from   pyraf import iraf

def prepMatch(reffile,imlist): #WZ
        """ Calculate shift/rotation 
            new reference: center position
        """
        de=100. # distance from edge
        circle = 500. # matching range in pixel
        rmin = 4. # FWHM minimum
        rmax = 15. # 10. # 6.
        # clx = 0.8 # must be a star
        pardir = '/home/zheng/drex/reffiles/pars'

        fitsfile = pyfits.open(reffile)
        detector = fitsfile[0].header.get('DETECTOR')
        fitsfile.close()
        #pdb.set_trace()
        weightfile = reffile.split(".fits")[0]+'_weight.fits'
        cmd = 'cp '+weightfile+' temp_weight.fits'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        if (detector =='WFC'):
            cmd = 'sex '+reffile+' -c '+pardir+'/tweak_wfc.inpar'
        else:
            if (detector =='UVIS'):
                cmd = 'sex '+reffile+' -c '+pardir+'/tweak_uvis.inpar'
            else:
                cmd = 'sex '+reffile+' -c '+pardir+'/tweak_ir.inpar'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        cmd = 'cp temp.cat temp0.cat'
        sproc  = popen2.Popen3(cmd,1)
        output = sproc.fromchild.readlines()
        errs   = sproc.childerr.readlines()
        x0,y0,mag0,err0,fwhm0,ra0,dec0,cls0=tableio.get_data('temp.cat',(1,2,3,4,5,6,7,8))
        xmin0=min(x0) + de
        xmax0=max(x0) - de
        ymin0=min(y0) + de
        ymax0=max(y0) - de
        #twf = open("tweak.dat","w")
        format=str('%4d %8.2f %8.2f %6.2f %6.2f %4.1f')
        #format2=str('%s %6.3f %5.3f %10.7f')
        for im in imlist:
            #outfile = im.split("_sci")[0]+'_ini.fits'
            fitsfile = pyfits.open(im)
            detector = fitsfile[0].header.get('DETECTOR')
            fitsfile.close()
            weightfile = im.split(".fits")[0]+'_weight.fits'
            cmd = 'cp '+weightfile+' temp_weight.fits'
            sproc  = popen2.Popen3(cmd,1)
            output = sproc.fromchild.readlines()
            errs   = sproc.childerr.readlines()
            #pdb.set_trace()
            if (detector =='WFC'):
                cmd = 'sex '+im+' -c '+pardir+'/tweak_wfc.inpar'
            else:
                if (detector =='UVIS'):
                    cmd = 'sex '+im+' -c '+pardir+'/tweak_uvis.inpar'
                else:
                    cmd = 'sex '+im+' -c '+pardir+'/tweak_ir.inpar'
            sproc  = popen2.Popen3(cmd,1)
            output = sproc.fromchild.readlines()
            errs   = sproc.childerr.readlines()
            if errs:
                print 'Program choked on '+cmd
            else:
                n,x,y,mag,err,fwhm,ra,dec,cls=tableio.get_data('temp.cat',(0,1,2,3,4,5,6,7,8))
                # pdb.set_trace()
                flag0 = numpy.zeros(len(x0),dtype=int)
                flag = numpy.zeros(len(n),dtype=int)
                index= numpy.zeros(len(n),dtype=int)
                xmin=min(x) + de
                if (xmin < xmin0):
                    xmin = xmin0
                xmax=max(x) - de
                if (xmax > xmax0):
                    xmax = xmax0
                ymin=min(y) + de
                if (ymin < ymin0):
                    ymin = ymin0
                ymax=max(y) - de
                if (ymax > ymax0):
                    ymax = ymax0
                print 'Min/max: ', xmin, xmax, ymin, ymax
                fname0 = im.split("_drz")[0]+'_1.cat'
                fname1 = im.split("_drz")[0]+'_2.cat'
                f0=open(fname0,"w")
                f=open(fname1,"w")
                k=0
                for i in range(len(x)):
                    # if (mag[i] < 90 and fwhm[i] > rmin and fwhm[i] < rmax):
                    if (mag[i] < 90):
                        for j in range(len(x0)):
                            if (abs(x[i]-x0[j]) < circle and abs(y[i]-y0[j]) < circle and mag0[j] < 90 and fwhm0[j] > rmin and fwhm0[j] < rmax):
                            #if (abs(x[i]-x0[j]) < circle and abs(y[i]-y0[j]) < circle and mag0[j] < 90):
                                index[i]=j
                                flag[i] = flag[i] + 1
                                flag0[j] = flag0[j] + 1
                                k= k +1
                print k, " sources found in ", im
                k=m=0
                for i in range(len(x)):
                    if (flag[i] >= 1):
                        #pdb.set_trace()
                        j=index[i]
                        if (flag0[j] >= 1 and fwhm0[j] < rmax and mag[i] < 90):
                        #if (flag0[j] >= 1 and mag[i] < 90):
                            print >> f, format % (k,x[i],y[i],mag[i],err[i],fwhm[i])
                            print >> f0, format % (k,x0[j],y0[j],mag0[j],err0[j],fwhm0[j])
                            k=k+1
                    else:
                        m= m+1
                print k, " sources used in ", im #,", ",m," missed"
            #pdb.set_trace()
            f0.close()
            f.close()
        return

os.chdir('/home/zheng/test')
format=str('%s %6.2f %6.2f %10.7f')
twf = open("tweak.par","w")
circle=25 # 15
match= 0
clx = 0.8
inp_list=tableio.get_str('drz.lis',(0))
twf = open("tweak.par","w")
format=str('%s %6.2f %5.2f %10.7f')

refim= inp_list[0]
pdb.set_trace()

prepMatch(refim, inp_list)

for i in range(len(inp_list)): #WZ
    im = inp_list[i]
    if string.find(refim,im) > -1: 
        dx = dy = dr = 0.0 # This is the file self
    else:
        fname0 = im.split("_drz")[0]+'_1.cat'
        fname1 = im.split("_drz")[0]+'_2.cat'
        cmd = 'match '+fname1+' 1 2 3 '+fname0+' 1 2 3 scale=1.0 matchrad='+str(circle)
        pdb.set_trace()
        sproc  = popen2.Popen3(cmd,1)
        matchout = str(sproc.fromchild.readlines())
        errs   = sproc.childerr.readlines()
        if errs:
            print 'Program choked on '+cmd
            dx=dy=dr = 0.
            match=-1
        else:
            #cmd = 'tweak '+' '+fname1+' '+fname0+' '+str(x)+' '+str(y)
            cmd = 'tweak '+' '+fname1+' '+fname0+' 750 750'
            print "Calculating tweaking parameters for "+im+" ..."
            sproc  = popen2.Popen3(cmd,1)
            matchout = str(sproc.fromchild.readlines())
            errs   = sproc.childerr.readlines()
            if errs:
                print 'Program choked on '+cmd
                dx=dy=dr= 0.
                match=-1
            else:
                dr=float(matchout[2:].split()[0])
                dx=float(matchout[2:].split()[1])
                dy=float(matchout[2:].split()[2])
                rms=matchout.split('rms:')[1]
            #self.logfile.write('Tweak: '+inpname+' dX= '+str(dx)+' dY= '+str(dy)+' dR= '+str(dr))
            print >> twf, format % ('****',dx,dy,dr) # image_name, dX, dY, rot      
twf.close()
print "Tweak.par generated"
#pdb.set_trace()

