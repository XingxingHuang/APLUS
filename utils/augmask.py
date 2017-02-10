#!/usr/bin/env python

# $Id: augmask.py,v 1.3 2003/01/31 20:50:25 anderson Exp $
# ---------------------------------------------------------------------
# augmask.py, created 29-Jan-03
# make a mask to cover areas where only a single image (with Ncombine=1)
# is available.
"""
Mental power tended to corrupt, and absolute intelligence tended to
corrupt absolutely, until the victim eschewed violence entirely in
favor of smart solutions to stupid problems.    -- Piers Anthony
"""
__version__      = '$Revision: 1.3 $ '[11:-3]
__version_date__ = '$Date: 2003/01/31 20:50:25 $ '[7:-3]
__author__       = "J. Blakeslee, jpb@pha.jhu.edu"
import os,string,numpy
import pyfits,fUtil

class augmask:
    """ Make a mask to cover wherever there's only a single non-CR-split image.
    Note, should be run from the images directory.
          augie = augmask.augmask(drob.parlists, logfile=logfile)
          augie.make()
    """
    def __init__(self, augimdict, ubername, logfile=None, clean_up=1):
        # self.augimdict[tab]['ctxlist'] = []
        
        self.modName    = string.split(string.split(str(self))[0],'.')[0][1:]
        self.logfile    = logfile
        self.augimdict  = augimdict
        self.removeList = []
        self.clean_up   = clean_up
        self.ubername   = ubername

        print self.modName,'version',__version__
        if self.logfile:
            self.logfile.write('Instantiating '+self.modName+' version '+__version__)
        
    def make(self):
        """Make the augmasks and uber-augmask (edgemask). """
        augmasks = []

        for asn in self.augimdict.keys():
            maskname = ''
            for ctx in self.augimdict[asn]['ctxlist']:
                ff = pyfits.open(ctx)
                if ff[0].header.get('NCOMBINE') == 1: 
                #if ff[0].header.get('NCOMBINE') <= 2: #WZ
                    if not maskname:
                        _NX = ff[0].header['NAXIS1']
                        _NY = ff[0].header['NAXIS2']
                        maskname = self.augimdict[asn]['maskname']
                        augmasks.append(maskname)
                        ### open the new mask file
                        if self.logfile:
                            self.logfile.write(' making '+maskname)
                        gg = fUtil.createfits(maskname,_NX,_NY,16)
                        self.removeList.append(maskname)
                        ### add contents of new mask
                    else:
                        ### sanity check
                        if _NX != ff[0].header['NAXIS1'] or _NY != ff[0].header['NAXIS2']:
                            raise Exception,'Individual Context ims of diff sizes?!'
                    gg[0].data += ff[0].data
                ff.close()
                del ff

            if maskname:
                gg[0].data = numpy.where(numpy.equal(gg[0].data,1),0,1).astype(numpy.int16)
                gg.close()
                del gg,maskname


        ### now make final asn-merged mask, if possible
        edgemask = ''
        if not augmasks:
            self.logfile.write("No Images suitable for edge masking.")
            self.logfile.write("No aug or edge masks made.")
            return

        for mask in augmasks:
            ff = pyfits.open(mask)
            if _NX != ff[0].header['NAXIS1'] or _NY != ff[0].header['NAXIS2']:
                raise Exception,'Individual augmasks of diff sizes?!'
            
            if not edgemask:
                edgemask = self.ubername
                if self.logfile:  self.logfile.write('Making '+edgemask)
                ## open new mask file
                mm = fUtil.createfits(edgemask,_NX,_NY,16)
                mm[0].data = ff[0].data

            else:
                # note, unary version (*=) doesn't work
                mm[0].data = mm[0].data * ff[0].data
            ff.close()
            
        mm[0].data = numpy.where(numpy.equal(mm[0].data,1),1,0).astype(numpy.int16)
        mm.close()
        del mm

        # final cleanup
        if self.clean_up:
            if self.logfile:
                self.logfile.write('Removing separate asn augmasks.')
            for file in self.removeList:
                try: os.remove(file)
                except: pass

        return

