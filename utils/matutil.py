#!/usr/bin/env python

# $Id: matutil.py,v 1.17 2010/08/29 20:19:42 Zheng Exp $
# ---------------------------------------------------------------------

"""
    matutil.py: utility functions for manipulating IDC distortion transformations.
    Run a separate task for non-ACS images
    #############################################################################
    AVAILABLE FUNCTIONS:
    
    getIDCinfo(IDCname, chip, direction='forward'):
        Returns: Cx,Cy,refdict,Norder.
    
    readIDCtxt(txtfile):
        Read the IDC coefficients from a text file with format:
           Cx[0,0]  Cy[0,0]
           Cx[0,1]  Cy[0,1]
             ...      ...
        Returns a tuple consisting of just the Cx,Cy matrices.
    
    getIDCxy_abstrans(x, y, IDCname, chip, verb=1,direction='forward'):
        Returns:  transformed <x,y> position.
    
    xyIDCtrans(x, y, Cx, Cy, verb=0):
        Returns:  transformed (distorted or undistorted) <x,y> position.
    
    rotateIDC(Cx,Cy,thetDeg):
        Rotate the Cx,Cy matrices by specified angle thetdeg (degrees).
        Remember that the V2,V3 coordinate system is left-handed, whereas
        the distortion-corrected x,y coordinate system is right-handed.
        If you want to rotate in the V2,V3 system, you need to use the
        negative of the angle.  So, know what you're really asking for.
        Returns rotated Cx,Cy tuple.
    
    convertIDCtrans(Cx,Cy,docheck=0,verb=0,direction='forward'):
        Returns: (xcoeffs,ycoeffs) 1-d transformation arrays.
    
    xyIDCinvert(xout, yout, Cx, Cy, verb=1, Err=InversionTol):
        Returns:  x, y  -  the distorted grid <x,y> position [relative to
                     XREF,YREF] which, when undistorted, ends up back at the
                     specified (xout, yout), to within <Err> tolerance.
    
    printIDC1d(Ctuple):
        Prints elements of each of a *tuple* of IDCTAB matrices.
    
    printIDCmat(Cmat):
        Prints an IDCTAB matrix is a nicer form than a simple 'print'.
        
    troll(paV3,dec,v2,v3):
        Input:  PA_V3(OTA), declination, V2, V3
        Returns:  PA_V3(target)

    #############################################################################
    EXAMPLES:
    
    --> from matutil import *
    --> idc3 = "wfc_idc3ord.fits"
    --> idc4 = "wfc4ord_17may02_idc.fits"
    
    ### Get the matrices and other info with getIDCinfo (calls _readIDCtab)
    
    --> c3x,c3y,refdict3,Norder = getIDCinfo(idc3,2)
    --> c4x,c4y,refdict4,Norder = getIDCinfo(idc4,2)
    --> print c3x.shape,c3y.shape
    (4, 4) (4, 4)
    --> print c4x.shape,c4y.shape
    (5, 5) (5, 5)
    
    ### Some ways of printing out the matrices in a nicer form than a
    ### simple "print" statement:
    
    --> printIDC1d((c3x,c3y))
     C[0,0]      0.0000000e+00     0.0000000e+00  
     C[1,0]      1.7180000e-03     5.0678998e-02  
     C[1,1]      4.9779002e-02     1.5829999e-03  
     C[2,0]      9.9285998e-08    -3.6764999e-07  
     C[2,1]     -2.5413999e-07     3.0679001e-07  
     C[2,2]      4.3161000e-07    -8.1236998e-08  
     C[3,0]     -1.0546000e-13    -1.9603000e-11  
     C[3,1]     -2.6058000e-11    -4.7987998e-12  
     C[3,2]     -2.6021999e-12    -2.7361001e-11  
     C[3,3]     -2.2858999e-11     3.7143999e-12  
    
    --> printIDCmat(c3y)    
       0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00  
       5.067900e-02    1.583000e-03    0.000000e+00    0.000000e+00  
      -3.676500e-07    3.067900e-07   -8.123700e-08    0.000000e+00  
      -1.960300e-11   -4.798800e-12   -2.736100e-11    3.714400e-12  
    
    --> printIDCmat(c4x) ; printIDCmat(c4y)    
       0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00  
       1.975460e-03    4.981722e-02    0.000000e+00    0.000000e+00    0.000000e+00  
       9.361248e-08   -2.480072e-07    4.220050e-07    0.000000e+00    0.000000e+00  
      -4.206189e-13   -2.571897e-11   -2.950768e-12   -2.325217e-11    0.000000e+00  
      -2.884212e-16   -6.280946e-16    2.225005e-15    6.075070e-16    1.173640e-15  
    
       0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00  
       5.027357e-02    1.486294e-03    0.000000e+00    0.000000e+00    0.000000e+00  
      -3.599211e-07    3.047960e-07   -7.519793e-08    0.000000e+00    0.000000e+00  
      -2.079035e-11   -4.851303e-12   -2.567549e-11    3.867417e-12    0.000000e+00  
      -1.328133e-15   -1.363596e-15   -1.475043e-15   -3.620283e-16   -8.563381e-16  
    
    ### Rotate the IDC transformation by 30 degrees - BEWARE - see source code!
    
    --> cxrot,cyrot = rotateIDC(c4x,c4y,30.)
    --> printIDCmat(cxrot) ; printIDCmat(cyrot)
       0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00  
       2.684758e-02    4.388613e-02    0.000000e+00    0.000000e+00    0.000000e+00  
      -9.888975e-08   -6.238253e-08    3.278681e-07    0.000000e+00    0.000000e+00  
      -1.075944e-11   -2.469893e-11   -1.539319e-11   -1.820327e-11    0.000000e+00  
      -9.138464e-16   -1.225744e-15    1.189390e-15    3.451024e-16    5.882329e-16  
     
       0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00  
       4.255046e-02   -2.362144e-02    0.000000e+00    0.000000e+00    0.000000e+00  
      -3.585070e-07    3.879646e-07   -2.761258e-07    0.000000e+00    0.000000e+00  
      -1.779466e-11    8.658133e-12   -2.076024e-11    1.497537e-11    0.000000e+00  
      -1.005986e-15   -8.668613e-16   -2.389927e-15   -6.172792e-16   -1.328431e-15  
    
    ### Transform an x,y *absolute* position (takes IDC name, not matrices):
    
    --> xn,yn= mm.getIDCxy_abstrans(3448,1724,idc3,2,verb=0)
    --> print xn,yn
    71.454543587 37.6154171062
    
    ### Transform an x,y pixel position **>relative to the reference pixel<**
    ### (this should probably be the more usual way of doing it):
    
    --> xp = 3448 - refdict4['XREF']
    --> yp = 1724 - refdict4['YREF']
    --> print xp,yp
    1400.0 700.0
    --> xyIDCtrans(1400,700,c4x,c4y)
    (71.678688505609841, 37.205813087226296)
    
    ### Convert position in arcsec back to pix (both relative to refpix)
    
    --> xp,yp=mm.xyIDCinvert(100,100,c3x,c3y)
    1942.87336095 1912.51672202
    1929.02918517 1930.21165942
    1929.02588383 1930.21620649
    1929.02588383 1930.21620649
    Linear extrapolation was off by 0.0967%
    
    ### Check if the inverse is correct by using xyIDCinvert to
    ### inverse transform the output from xyIDCtrans:
    
    --> xp1,yp1 = 2000.,1000.
    --> xn,yn = xyIDCtrans(xp1,yp1,c4x,c4y)
    --> print xn,yn
    102.676838398 53.0662295451
    
    --> xp2,yp2 = xyIDCinvert(xn,yn,c4x,c4y,verb=1)
    2021.58417063 995.782873999
    2000.00343411 999.998480465
    2000.0 1000.0
    2000.0 1000.0
    Linear extrapolation was off by 0.7808%

    ### convert roll angles from V1 axis position to aperture posn
    --> v2w1=258.023; v3w1 = 196.40408
    --> v2w2=254.255; v3w2 = 301.1133
    
    --> print targroll(roll,dec,v2w1,v3w1)
    30.1547847837
    --> print targroll(roll,dec,v2w2,v3w2)       
    30.178557734
    --> print targroll(roll,dec,v2w2,v3w2,updown=1)
    210.178557734
    
    #############################################################################
"""

__version__      = '$Revision: 1.17 $ '[11:-3]
__version_date__ = '$Date: 2006/07/21 20:19:42 $ '[7:-3]
__author__       = "J Blakeslee, jpb@pha.jhu.edu"

import string,math
import pdb
import pyfits,numpy
InversionTol = 1e-8
InverMaxIter = 500

def getACSIDCinfo(IDCname, chip, direction='forward',verbose=1,
               filter1=None,filter2=None):
    """    This function takes IDC file name, chip ID, and optional direction.
    Then it just calls _readIDCtab() and does some error checks.
    Returns: Cx,Cy,refdict,Norder
    e.g.,
    >> Cx,Cy,refdict,Norder = matutil.getIDCinfo(idcfile,2)
    """
    # pdb.set_trace()
    Cx,Cy,refdict,Norder = _readACSIDCtab(IDCname,chip,direction,verbose=verbose,
                                       filter1=filter1,filter2=filter2)
    if(Cx.shape != Cy.shape or Cx.shape[0] != Cx.shape[1]):
        raise ValueError, "Severe x,y IDC Table shape inconsistency!"
    if(Norder != (Cx.shape[0] - 1)):
        raise ValueError, "Norder inconsistent with matrix sizes!"
    return Cx,Cy,refdict,Norder

def getIDCinfo(IDCname, chip, direction='forward',verbose=1,
               filter=None):
#              filter1=None,filter2=None):
    """    This function takes IDC file name, chip ID, and optional direction.
    Then it just calls _readIDCtab() and does some error checks.
    Returns: Cx,Cy,refdict,Norder
    e.g.,
    >> Cx,Cy,refdict,Norder = matutil.getIDCinfo(idcfile,2)
    """
    Cx,Cy,refdict,Norder = _readIDCtab(IDCname,chip,direction,verbose=verbose,
                                       filter=filter)

    if(Cx.shape != Cy.shape or Cx.shape[0] != Cx.shape[1]):
        raise ValueError, "Severe x,y IDC Table shape inconsistency!"
    if(Norder != (Cx.shape[0] - 1)):
        raise ValueError, "Norder inconsistent with matrix sizes!"
    return Cx,Cy,refdict,Norder


def readIDCtxt(txtfile):
    """  Read the IDC coefficients from a text file with format:
       Cx[0,0]  Cy[0,0]
       Cx[0,1]  Cy[0,1]
       Cx[1,1]  Cy[1,1]
       Cx[0,2]  Cy[0,2]
       Cx[1,2]  Cy[1,2]
       Cx[2,2]  Cy[2,2]
       Cx[0,3]  Cy[0,3]
         ...      ...
    Comment and blank lines are ignored.
    Returns a tuple consisting of just the Cx,Cy matrices.
    """
    xco = []
    yco = []
    _rawlist = open(txtfile).readlines()
    for line in _rawlist:
        if (len(line.split()) < 1 or line.split()[0][0] == '#'):
            continue
        if len(line.split()) != 2:
            raise Exception,"Only one Cx[i,j] Cy[i,j] pair per line!"
        xco.append(float(line.split()[0]))
        yco.append(float(line.split()[1]))
    del _rawlist,line

    Ncoeff = len(xco)
    _tmp = (math.sqrt(9+8*(Ncoeff-1)) - 3)/2.0
    Norder = int(round(_tmp))
    if abs(Norder - _tmp) > 1e-9:
        raise Exception,"Bad number of coefficients in "+txtfile
    del _tmp
    Cx = numpy.zeros((Norder+1,Norder+1),numpy.float64)
    Cy = numpy.zeros((Norder+1,Norder+1),numpy.float64)

    k=0
    for i in range(Norder+1):
        for j in range(0,i+1):
            Cx[i,j] = xco[k]
            Cy[i,j] = yco[k]
            k += 1
            if k > Ncoeff:
                raise ValueError,"k>Ncoeff? That can't be!"

    return Cx,Cy

            
def getIDCxy_abstrans(x, y, IDCname, chip, verb=1,direction='forward'):
    """  Read IDC table, translate x,y absolute position.
    Takes:  <x,y> position (*absolute* - don't subtract off XREF,YREF);
            <IDCname> IDC FITS table file name;
            chip ID; optional direction (defaults to 'forward')
    Returns:  transformed <x,y> position.
    """
    # get the IDC table coefficients and check for weirdness
    Cx,Cy,refdict,Norder = getIDCinfo(IDCname,chip,direction)
    if(Cx.shape != Cy.shape or Cx.shape[0] != Cx.shape[1]):
        raise ValueError, "Severe x,y IDC Table shape inconsistency!"
    if(Norder != (Cx.shape[0] - 1)):
        raise ValueError, "Norder inconsistent with matrix sizes!"
    
    xref = refdict['XREF']
    yref = refdict['YREF']

    if verb:
        print 'Norder =',Norder, 'XREF,YREF=',xref,yref
        print Cx ; print Cy
    # subtract off the x,y zeropoints XREF,YREF (same as CRPIX1,2?)
    # note, this only affects the local copy of x,y
    x = float(x - xref)
    y = float(y - yref)
    xnew = 0.
    ynew = 0.

    # Note, the Cx[0,0] and Cy[0,0] terms should be zero for forward
    # direction, according to the IDCTab ISR.
    if ((Cx[0,0] != Cy[0,0] or Cx[0,0] != 0.) and direction=='forward'):
        raise ValueError, "Zeroth order IDCtab terms should be 0."
    
    if direction == 'forward':
        for i in range(1,Norder+1):
            for j in range(0,i+1):
                # Note: it's x^j and y^(i-j)
                imj = i-j
                xterm = Cx[i,j] * x**j * y**imj
                yterm = Cy[i,j] * x**j * y**imj
                if verb: print i,j,xterm,yterm
                xnew = xnew + xterm
                ynew = ynew + yterm
    else:
        raise ValueError, "Don\'t know how to do \'%s\' transforms."%direction
    return xnew,ynew

def xyIDCtrans(x, y, Cx, Cy, verb=0):
    """  Function to transform an x,y position using IDCtab distortion solution.
    Takes:  <x,y> posn; <Cx,Cy> IDC-like matrices; <Norder> polynomial order
    Returns:  transformed (distorted or undistorted) <x,y> position.
    **************************************************************************
    ***************************>    NOTE WELL!    <***************************
    ***  Assumes that input <x,y> are relative to the "reference position" ***
    ***  THEREFORE, to transform raw (x,y) pixel positions, YOU, dear user,***
    ***  MUST do the following:                                            ***
    ***            x = x - refdict['XREF']                                 ***
    ***            y = y - refdict['YREF']                                 ***
    ***            xc,yc = xyIDCtrans(x, y, Cx, Cy)                        ***
    ***  where <refdict> is the python dictionary returned by the          ***
    ***  drutil.readIDCtab() method or, equivalently, getIDCinfo() above   ***
    ***  (this dictionary is populated from the IDC FITS Table header).    ***
    **************************************************************************
    **************************************************************************"""
    if(Cx.shape != Cy.shape or Cx.shape[0] != Cx.shape[1]):
        raise ValueError, "Severe x,y IDC Table shape inconsistency!"
    x = float(x)
    y = float(y)
    xnew = ynew = 0.

    if (Cx[0,0] != Cy[0,0] or Cx[0,0] != 0.):
        print "WARNING: Zeroth order IDCtab terms not 0."

    for i in range(Cx.shape[0]):
        for j in range(0,i+1):
            # Note: it's x^j and y^(i-j)
            imj = i-j
            xterm = Cx[i,j] * x**j * y**imj
            yterm = Cy[i,j] * x**j * y**imj
            xnew += xterm
            ynew += yterm
            if verb: print '%d %d %14.9f %13.9f %14.9f %14.9f' \
               %(i,j,xterm,yterm,xnew,ynew)

    return xnew,ynew

def rotateIDC(Cx,Cy,thetDeg):
    """ Rotate the Cx,Cy matrices by specified angle thetdeg (degrees).
    Seems simple, but remember that the distortion solution is ultimately
    tied to the V2,V3 coordinate system, which is left-handed, whereas
    the distortion-corrected x,y coordinate system is right-handed.
    If you want to rotate in the V2,V3 system, you need to use the
    negative of the angle.  So, know what you're really asking for.
    """
    if(Cx.shape != Cy.shape or Cx.shape[0] != Cx.shape[1]):
        raise ValueError, "Severe x,y IDC Table shape inconsistency!"
    
    thetRad = thetDeg * math.pi/180.
    cosAng = math.cos(thetRad)
    sinAng = math.sin(thetRad)
    
    rotCx =  cosAng * Cx +  sinAng * Cy
    rotCy = -sinAng * Cx +  cosAng * Cy
    return rotCx,rotCy


def convertIDCtrans(Cx,Cy,docheck=0,verb=0,direction='forward'):
    """  Convert IDCtab-style Cx,Cy matrices to a 1-d  transformation
    coefficient array, the form used in standard drizzle implementation.
    The returned arrays could also be sent to cubinvert().
    ***> Input:   IDC FITS table file name; chip ID; check indices?
    ***> Returns: <xcoeffs[0:],ycoeffs[0:]> 1-d transformation arrays.
    e.g., the mapping for an Norder=3 IDCTAB (a 4x4 matrix) will be:
          xcoeff[0] = Cx[0,0]     ycoeff[0] = Cy[0,0]
          xcoeff[1] = Cx[1,1]     ycoeff[1] = Cy[1,1]
          xcoeff[2] = Cx[1,0]     ycoeff[2] = Cy[1,0]
          xcoeff[3] = Cx[2,2]     ycoeff[3] = Cy[2,2]
          xcoeff[4] = Cx[2,1]     ycoeff[4] = Cy[2,1]
          xcoeff[5] = Cx[2,0]     ycoeff[5] = Cy[2,0]
          xcoeff[6] = Cx[3,3]     ycoeff[6] = Cy[3,3]
          xcoeff[7] = Cx[3,2]     ycoeff[7] = Cy[3,2]
          xcoeff[8] = Cx[3,1]     ycoeff[8] = Cy[3,1]
          xcoeff[9] = Cx[3,0]     ycoeff[9] = Cy[3,0]
    """
    if(Cx.shape != Cy.shape or Cx.shape[0] != Cx.shape[1]):
        raise ValueError, "Severe x,y IDC Table shape inconsistency!"
    Norder = Cx.shape[0] - 1
    Ncoeff = (Norder+1)*(Norder+2)/2
    if docheck or verb:
        print 'Converting Norder =',Norder,'matrix to Ncoeff =',Ncoeff,'array.'

    xcoeff = numpy.zeros((Ncoeff,),numpy.float64)
    ycoeff = numpy.zeros((Ncoeff,),numpy.float64)
    if docheck:
        karray = numpy.zeros((Ncoeff,3),numpy.int16)

    isum=0
    for i in range(Norder+1):
        isum = isum+i
        for j in range(0,i+1):
            k = isum + (i-j)
            xcoeff[k] = Cx[i,j]
            ycoeff[k] = Cy[i,j]
            if docheck: karray[k] = (k,i,j)
    if docheck:
        for k in range(Ncoeff):
            print ('coeff[%2d]  =  C[ %d , %d ]' \
                   %(karray[k,0],karray[k,1],karray[k,2])),
            print '  =  %16.7e  %16.7e'%(xcoeff[k],ycoeff[k])

    return xcoeff, ycoeff

def xyIDCinvert(xout, yout, Cx, Cy, verb=1, Err=InversionTol):
    """  Function to transform back a dedistorted <x_c, y_c> position using the
    inverse of the 'forward' IDCtab distortion solution. The inverse is solved
    for iteratively using the Newton-Raphson method (successive approximation).

    Takes:
       xout,yout - (real) output (dedistorted) x,y from the distortion model.
        Cx, Cy   - (real) 2-d IDC arrays
         verb    - (int) 1=> print x,y at each iteration; 0=> don't print
         Err     - (real) accuracy required of the inversion (in X and Y)
    Returns:
         x, y  -  distorted grid <x,y> position [relative to XREF,YREF] which,
                 when undistorted, ends up back at the specified (xout, yout),
                 to within <Err> tolerance.
    *************************************************************************
    ************************>   NOTE THIS WELL!   <**************************
    ***  The returned x,y "untransformed" or "redistorted" pixel positons ***
    ***  are relative to the "reference position".  THEREFORE to get the  ***
    ***  absolute (x,y) pixel positions on chip, do the following:        ***
    ***            xp,yp = xyIDCinvert(x_arcsec,y_arcsec,Cx,Cy)           ***
    ***            xp += refdict['XREF']                                  ***
    ***            yp += refdict['YREF']                                  ***
    ***  where <refdict> is returned with Cx,Cy,Norder by getIDCinfo().   ***
    ***  NOTE ALSO:  this routine now no longer converts matrices to 1-d  ***
    ***    but just works directly with 2-d Cx,Cy.                        ***
    *************************************************************************
    *************************************************************************"""
    if(Cx.shape != Cy.shape or Cx.shape[0] != Cx.shape[1]):
        raise ValueError, "Cx,Cy IDC Matrices should be NxN and of same size."

    # First check to see if the transformation is pathological;
    # if not, use the linear order as the initial approximation.
    if(Cx[1,1] == 0. or Cx[1,0]*Cy[1,1] == Cx[1,1]*Cy[1,0]):
        _xtmp,_ytmp = xyIDCtrans(xout,yout,Cx,Cy)
        x = 2.0*xout - _xtmp
        y = 2.0*yout - _ytmp
        del _xtmp,_ytmp
    else:
        y0numer = Cy[1,1]*(xout - Cx[0,0]) - Cx[1,1]*(yout - Cy[0,0])
        y0denom = Cx[1,0]*Cy[1,1] - Cx[1,1]*Cy[1,0]
        y = y0numer/y0denom
        x = (xout - Cx[0,0] - Cx[1,0]*y)/Cx[1,1]
    if verb:
        xinit = x; yinit=y
        print x,y

    doExtra=1
    for niter in range(InverMaxIter):
        # Transform (distort) current x,y guesses.
        xo,yo = xyIDCtrans(x,y,Cx,Cy)

        # Check if current distorted position agrees with specified 
        # position to within the specified error criterion.
        if(abs(xout-xo) < Err and abs(yout-yo) < Err):
            if doExtra:
                doExtra = 0
            else:
                break

        # Calculate the 4 derivatives and the determinant.
        dxx,dxy = _IDCderiv(x,y,Cx)
        dyx,dyy = _IDCderiv(x,y,Cy)
        det = dxx*dyy-dyx*dxy

        # Improve guesses using Newton-Raphson.
        x = x + ((xout-xo)*dyy-(yout-yo)*dxy)/det
        y = y + ((yout-yo)*dxx-(xout-xo)*dyx)/det
        if verb: print x,y
        niter += 1

    if niter >= InverMaxIter:
        print "xyIDCinvert: failed to converge in",niter,"iterations."
        return None,None

    if verb:
        raddist = math.sqrt(x**2 + y**2)
        if raddist > 0:
            pcerror = 100*abs(math.sqrt(xinit**2 + yinit**2) - raddist)/raddist
            print 'Linear extrapolation was off by %.4f%c\n'%(pcerror,'%')
    return x,y


def _IDCderiv(x,y,C):
    """  Evaluate derivatives in x & y of a 2-d IDC-matrix distortion.
    This is mainly for inverting the distortion using Newton-Raphson.
    Takes:   x, y, and C_i [a single IDC-Matrix, not both!]
    Returns:  dI/dx, dI/dy """

    if(C.shape[0] != C.shape[1]):
        raise ValueError, "_IDCderiv: IDC Matrix not NxN !"
    dIdx = 0.
    dIdy = 0.
    for i in range(C.shape[0]):
        for j in range(0,i+1):
            # Note: it's x^j and y^(i-j)
            xpow = j
            ypow = i-j
            if xpow > 0:
                dIdx += xpow * C[i,j] * x**(xpow-1) * y**(ypow)
            if ypow > 0:
                dIdy += ypow * C[i,j] * x**(xpow) * y**(ypow-1)

    return dIdx,dIdy

def printIDC1d(Ctuple):
    "Prints elements of each of a *tuple* of IDCTAB matrices.  "
    for i in range(Ctuple[0].shape[0]):
        for j in range(0,i+1):
            print (' C[%d,%d] ' %(i,j)),
            for C in Ctuple:
                print (" %16.7e" % C[i,j]),
            print " "

def printIDCmat(Cmat):
    "Prints an IDCTAB matrix is a nicer form than a simple 'print'. "
    for i in range(Cmat.shape[0]):
        for j in range(Cmat.shape[1]):
            print " %14.6e" % Cmat[i,j],
        print " "
    print " "
    
def _readACSIDCtab(tabname, chip=1, direction='forward', verbose=0,filter1=None,filter2=None):
    """  Function to read the IDC table and generates the two matrices
    with the geometric correction coefficients.  
            INPUT:   FITS object of open IDC table
            OUTPUT:  2 coefficient matrices (for x & y)
                     the refdict, and the order
    Function borrowed from Warren Hack's drutil module, then modified
    to work with the new pyfits.
    """
    try:
        ftab = pyfits.open(tabname)
    except:
        raise IOError,"IDC table '%s' not valid as specified!" % tabname

    # Implement default values for filters here to avoid the default
    # being overwritten by values of None passed by user.
    #pdb.set_trace()
    if filter1 == None or filter1.find('CLEAR') == 0:
        filter1 = 'CLEAR'
    if filter2 == None or filter2.find('CLEAR') == 0:
        filter2 = 'CLEAR'

    # if verbose: print "*** matutil got filters: ",filter1,filter2

    # First thing we need, is to read in the coefficients from the IDC
    # table and populate the Fx and Fy matrices.
    
    if ftab['PRIMARY'].header.has_key('DETECTOR'):
        detector = ftab['PRIMARY'].header['DETECTOR']
    else:
        if ftab['PRIMARY'].header.has_key('CAMERA'):
            detector = str(ftab['PRIMARY'].header['CAMERA'])
        else:
            detector = 1

    # Set default filters for SBC
    if detector == 'SBC':
        if filter1 == 'CLEAR':
            filter1 = 'F115LP'
            filter2 = 'N/A'
        if filter2 == 'CLEAR':
            filter2 = 'N/A'

    # Read FITS header to determine order of fit, i.e. k
    order = ftab['PRIMARY'].header['NORDER']
    
    fx = numpy.zeros((order+1,order+1),numpy.float64)
    fy = numpy.zeros((order+1,order+1),numpy.float64)
    
    #Determine row from which to get the coefficients.
    # How many rows do we have in the table...
    fshape = ftab[1].data.shape
    colnames = ftab[1].data._names
    row = -1

    # Loop over all the rows looking for the one which corresponds
    # to the value of CCDCHIP we are working on...
    for i in xrange(fshape[0]):

        try:
            # Match FILTER combo to appropriate row,
            #if there is a filter column in the IDCTAB...
            if 'FILTER1' in colnames and 'FILTER2' in colnames:

                filt1 = ftab[1].data.field('FILTER1')[i]
                if filt1.find('CLEAR') > -1: filt1 = filt1[:5]

                filt2 = ftab[1].data.field('FILTER2')[i]
                if filt2.find('CLEAR') > -1: filt2 = filt2[:5]
            else:
                if 'OPT_ELEM' in colnames:
                    filt1 = ftab[1].data.field('OPT_ELEM')
                    if filt1.find('CLEAR') > -1: filt1 = filt1[:5]
                else:
                    filt1 = filter1

                if 'FILTER' in colnames:
                    _filt = ftab[1].data.field('FILTER')[i]
                    if _filt.find('CLEAR') > -1: _filt = _filt[:5]
                    if 'OPT_ELEM' in colnames:
                        filt2 = _filt
                    else:
                        filt1 = _filt
                        filt2 = 'CLEAR'
                else:
                    filt2 = filter2
        except:
            # Otherwise assume all rows apply and compare to input filters...
            filt1 = filter1
            filt2 = filter2

        try:
            detchip = ftab[1].data[i].field('DETCHIP')
        except LookupError,NameError:
            detchip = 1

        try:			
            direct = string.lower(ftab[1].data[i].field('DIRECTION'))
        except:
            print 'matutil: _readIDCtab did not find valid DIRECTION'
            raise LookupError,'matutil: _readIDCtab did not find valid DIRECTION'

        if filt1 == filter1.strip() and filt2 == filter2.strip():
            if direct == direction.strip():
                if int(detchip) == int(chip) or int(detchip) == -999:
                    #print "change row"
                    #pdb.set_trace()
                    row = i
                    break
            
    #pdb.set_trace()
    if row < 0:
        err_str = '\nProblem finding row in IDCTAB! Could not find row matching:\n'
        err_str += '        CHIP: '+str(chip)+'\n'
        err_str += '     FILTERS: '+filter1+','+filter2+'\n'
        ftab.close()
        del ftab
        raise LookupError,err_str
    else:
        if verbose:
            print ' _readIDCtab: IDCTAB Distortion model from',tabname,'row',str(row+1),'for'
            print '      chip,direct,filt1,filt2:',detchip,direct,filter1,filter2

    refpix = {}
    refpix['XREF'] = ftab[1].data[row].field('XREF')
    refpix['YREF'] = ftab[1].data[row].field('YREF')
    refpix['XSIZE'] = ftab[1].data[row].field('XSIZE')
    refpix['YSIZE'] = ftab[1].data[row].field('YSIZE')
    refpix['PSCALE'] = ftab[1].data[row].field('SCALE')
    refpix['V2REF'] = ftab[1].data[row].field('V2REF')
    refpix['V3REF'] = ftab[1].data[row].field('V3REF')
    if ('THETA' in ftab[1].columns.names):
      refpix['THETA'] = ftab[1].data[row].field('THETA')
    else:
      refpix['THETA'] = None
    refpix['XDELTA'] = 0.0
    refpix['YDELTA'] = 0.0

    # Now that we know which row to look at, read coefficients into the
    # numeric arrays we've created.  Try to figure out which column name
    # convention the IDCTAB follows, either A,B or CX,CY.
    if 'CX10' in ftab[1].data._names:
        cxstr = 'CX'
        cystr = 'CY'
    else:
        cxstr = 'A'
        cystr = 'B'
    
    for i in xrange(1,order+1):
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

#########################

def _readIDCtab(tabname, chip=1, direction='forward', verbose=0,filter=None):
    """  Function to read the IDC table and generates the two matrices
    with the geometric correction coefficients.  
            INPUT:   FITS object of open IDC table
            OUTPUT:  2 coefficient matrices (for x & y)
                     the refdict, and the order
    Function borrowed from Warren Hack's drutil module, then modified
    to work with the new pyfits.
    """
    #pdb.set_trace()
    try:
        ftab = pyfits.open(tabname)
    except:
        raise IOError,"IDC table '%s' not valid as specified!" % tabname

    #fref = pyfits.open()
    #hd = ftab['PRIMARY'].header #WZ
    #fref.close()
    #del fref

    #flter =  hd['FILTER']

    # Implement default values for filters here to avoid the default
    # being overwritten by values of None passed by user.
    if filter == None or string.find(filter,"F105W")>-1:
        filter = 'F105W'
    #if filter2 == None or filter2.find('CLEAR') == 0:
    #    filter2 = 'CLEAR'

    # if verbose: print "*** matutil got filters: ",filter1,filter2

    # First thing we need, is to read in the coefficients from the IDC
    # table and populate the Fx and Fy matrices.

    #pdb.set_trace()
    if ftab['PRIMARY'].header.has_key('DETECTOR'):
        detector = ftab['PRIMARY'].header['DETECTOR']
        if detector == 'UVIS':    # WZ
            filter = 'F606W'
    else:
        if ftab['PRIMARY'].header.has_key('CAMERA'):
            detector = str(ftab['PRIMARY'].header['CAMERA'])
        else:
            detector = 1

    # Read FITS header to determine order of fit, i.e. k
    order = ftab['PRIMARY'].header['NORDER']
    
    fx = numpy.zeros((order+1,order+1),numpy.float64)
    fy = numpy.zeros((order+1,order+1),numpy.float64)
    
    #Determine row from which to get the coefficients.
    # How many rows do we have in the table...
    fshape = ftab[1].data.shape
    colnames = ftab[1].data._names
    row = -1

    # Loop over all the rows looking for the one which corresponds
    # to the value of CCDCHIP we are working on...
    for i in xrange(fshape[0]):

        try:
            # Match FILTER combo to appropriate row,
            #if there is a filter column in the IDCTAB...
            if 'FILTER' in colnames:
                #pdb.set_trace()
                filt = ftab[1].data.field('FILTER')[i]
                if filt.find('F200LP') > -1: filt = filt[:5]

                #filt2 = ftab[1].data.field('FILTER2')[i]
                #if filt2.find('CLEAR') > -1: filt2 = filt2[:5]
            else:
                if 'OPT_ELEM' in colnames:
                    filt = ftab[1].data.field('OPT_ELEM')
                    if filt.find('F200LP') > -1: filt = filt[:5]
                else:
                    filt = filter

                if 'FILTER' in colnames:
                    _filt = ftab[1].data.field('FILTER')[i]
                    if _filt.find('F200LP') > -1: _filt = _filt[:5]
                    if 'OPT_ELEM' in colnames:
                        pass # filt2 = _filt
                    else:
                        filt = _filt

        except:
            # Otherwise assume all rows apply and compare to input filters...
            filt = filter
        #pdb.set_trace()
        try:
            detchip = ftab[1].data[i].field('DETCHIP')
        except LookupError,NameError:
            detchip = 1

        try:			
            direct = string.lower(ftab[1].data[i].field('DIRECTION'))
        except:
            print 'matutil: _readIDCtab did not find valid DIRECTION'
            raise LookupError,'matutil: _readIDCtab did not find valid DIRECTION'

        if filt == filter: # and filt2 == filter2.strip():
            if direct == direction.strip():
                if int(detchip) == int(chip) or int(detchip) == -999:
                    row = i
                    break
    #print "Row = ", row, filt, filter
    #pdb.set_trace()
    #filter='F606W'
    #filt=filter
    #row = 29
    if row < 0:
        err_str = '\nProblem finding row in IDCTAB! Could not find row matching:\n'
        err_str += '        CHIP: '+str(chip)+'\n'
        err_str += '      FILTER: '+filter+'\n'
        #err_str += '     FILTERS: '+filter1+','+filter2+'\n'
        ftab.close()
        del ftab
        raise LookupError,err_str
    else:
        if verbose:
            print ' _readIDCtab: IDCTAB Distortion model from',tabname,'row',str(row+1),'for'
            print '      chip,direct,filt:',detchip,direct,filter

    refpix = {}
    refpix['XREF'] = ftab[1].data[row].field('XREF')
    refpix['YREF'] = ftab[1].data[row].field('YREF')
    refpix['XSIZE'] = ftab[1].data[row].field('XSIZE')
    refpix['YSIZE'] = ftab[1].data[row].field('YSIZE')
    refpix['PSCALE'] = ftab[1].data[row].field('SCALE')
    refpix['V2REF'] = ftab[1].data[row].field('V2REF')
    refpix['V3REF'] = ftab[1].data[row].field('V3REF')
    #print "THETA"
    #pdb.set_trace()
    if ('THETA' in ftab[1].columns.names):
      refpix['THETA'] = ftab[1].data[row].field('THETA')
    else:
      refpix['THETA'] = None
    refpix['XDELTA'] = 0.0
    refpix['YDELTA'] = 0.0

    # Now that we know which row to look at, read coefficients into the
    # numeric arrays we've created.  Try to figure out which column name
    # convention the IDCTAB follows, either A,B or CX,CY.
    if 'CX10' in ftab[1].data._names:
        cxstr = 'CX'
        cystr = 'CY'
    else:
        cxstr = 'A'
        cystr = 'B'
    
    for i in xrange(1,order+1):
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

def _getColIndex (extn, colname):
    """ Generic function for returning the column index value based
    on the column name.  This one also from Warren Hack's drutil."""
    found = 0
    cindex = 0

    nfields = extn.header['TFIELDS']
    for col in xrange(nfields):
        cindex = col + 1
        ttype = 'TTYPE'+str(cindex)
        if string.strip(extn.header[ttype]) == colname:
            found = 1
            return col
    if found == 0:
        raise LookupError,"Column "+colname+" not found."
    return None


#############################################################################
# Following are old Norder=3 specific functions, superceded by those above. #
#############################################################################

def _cubeval(x,y,coeff):
    """  Evaluate a cubic geometric distortion of the form given in
    the WFPC2 instrument handbook."""

    eval = coeff[0] + coeff[1]*x + coeff[2]*y + \
           coeff[3]*x*x + coeff[4]*x*y + coeff[5]*y*y + \
           coeff[6]*x*x*x + coeff[7]*x*x*y + \
           coeff[8]*x*y*y + coeff[9]*y*y*y 
    return eval

def _cubderiv(x,y,c):
    """  Evaluate derivatives in x &y of a 2-d cubic polynomial distortion.
    This is mainly for inverting the distortion using Newton-Raphson.
    Takes:   x, y, distortion coefficients
    Returns:  df/dx, df/dy """

    dx = c[1] + 2.0*c[3]*x + c[4]*y + 3.0*c[6]*x*x + 2.0*c[7]*x*y + c[8]*y*y
    dy = c[2] + c[4]*x + 2.0*c[5]*y + c[7]*x*x + 2.0*c[8]*x*y + 3.0*c[9]*y*y
    return dx,dy

def cubinvert(xout, yout, xco, yco, verb=1, Err=InversionTol):
    """    * Invert a 2-d cubic distortion using the Newton-Raphson method
    * (i.e., successive approximation).
    *  
    * Takes the following:
    *    Xout, Yout - (real) output position from the geometric distortion.
    *    Xco, Yco  - (real) arrays the cubic distortion coefficients
    *    verb    - (int) 1=> print x,y at each iteration; 0=> don't print
    *    Err     - (real) accuracy required of the inversion (in X and Y)
    * 
    * Returns the real values:
    *    x, y  -  <x,y> position [relative to XREF,YREF] which, when undistorted,
    *          ends up at the specified (xout, yout), to within <Err> tolerance.
    * 
    * This routine is similar to the one used in Richard Hook's blot.f code.
    * It's iterative and therefore slow in principle, but converges pretty fast
    * on the cube roots in practice, especially with the linear 1st-order guess.
    """
    # Make an initial guess at the (unknown) distorted position.
    # First check to see if the cubic transformation is pathological;
    # if not, use the linear order as the initial approximation.
    if(len(xco) != 10 or len(yco) != 10):
        raise ValueError, "1-d distortion matrices not of length 10."
    if(xco[1] == 0. or xco[2]*yco[1] == xco[1]*yco[2]):
        x = 2.0*xout - _cubeval(xout,yout,xco)
        y = 2.0*yout - _cubeval(xout,yout,yco)
    else:
        y0numer = yco[1]*(xout - xco[0]) - xco[1]*(yout - yco[0])
        y0denom = xco[2]*yco[1] - xco[1]*yco[2]
        y = y0numer/y0denom
        x = (xout - xco[0] - xco[2]*y)/xco[1]
    if verb:
        xinit = x; yinit=y
        print x,y

    while(1):
        # Transform (distort) current x,y guesses.
         xo = _cubeval(x,y,xco)
         yo = _cubeval(x,y,yco)

         # Check if current distorted position agrees with specified 
         # position to within the specified error criterion.
         if(abs(xout-xo) < Err and abs(yout-yo) < Err):
             break

         # Calculate the 4 derivatives and the determinant.
         dxx,dxy = _cubderiv(x,y,xco)
         dyx,dyy = _cubderiv(x,y,yco)
         det = dxx*dyy-dyx*dxy

         # Improve guesses using Newton-Raphson.
         x = x + ((xout-xo)*dyy-(yout-yo)*dxy)/det
         y = y + ((yout-yo)*dxx-(xout-xo)*dyx)/det
         if verb: print x,y

    if verb:
        raddist = math.sqrt(x**2 + y**2)
        if raddist > 0:
            pcerror = 100*abs(math.sqrt(xinit**2 + yinit**2) - raddist)/raddist
            print 'Linear extrapolation was off by %.4f%c\n'%(pcerror,'%')
    return x,y

def cubxyIDCinvert(x_c, y_c, Cx, Cy, verb=1, Err=InversionTol):
    """  Function to transform back an dedistorted <x_c, y_c> position using
    the inverse of the 'forward' IDCtab distortion solution.
    Takes:  dedistored <x,y> posn; <Cx,Cy> IDC-like matrices.
    Returns:  inverse transformed (\"redistorted\") <x,y> position, DEFINED
      RELATIVE to (XREF,YREF) pixel postion.  Use refdict get absolute (x,y).
    *************************************************************************
    ************************>   NOTE THIS WELL!   <**************************
    ***  The returned x,y "untransformed" or "redistorted" pixel positons ***
    ***  are relative to the "reference position".  THEREFORE to get the  ***
    ***  absolute (x,y) pixel positions on chip, do the following:        ***
    ***            xp,yp = xyIDCinvert(x_arcsec,y_arcsec,Cx,Cy)           ***
    ***            xp += refdict['XREF']                                  ***
    ***            yp += refdict['YREF']                                  ***
    ***  where <refdict> is returned along with Cx,Cy,Norder by the       ***
    ***  drutil.readIDCtab() method or, equivalently, getIDCinfo() above. ***
    ***  NOTE ALSO:  this converts IDC matrices to 1-d for each call.     ***
    ***  TO SAVE TIME:  Use cubinvert() directly with the 1-d arrays.     ***
    *************************************************************************"""
    if(Cx.shape != Cy.shape or Cx.shape != (4,4)):
        raise ValueError, "Can't do a cubic inversion if IDCTABs not order 3!"
    
    # Convert the IDC matrices to 1-d 10-member coefficient arrays.
    xcoeff, ycoeff = convertIDCtrans(Cx,Cy,0)
    xpix,ypix = cubinvert(x_c, y_c, xcoeff, ycoeff, verb, Err)
    return xpix,ypix

def help():
    print __doc__

def targroll(roll,dec,v2,v3,updown=0):
    """ Given the following:
         roll   -  roll angle defined up yonder at V1 Axis, i.e., header PA_V3 [deg]
         dec    -  declination of target [deg]
         v2     -  v2 position of target [arcsec]
         v3     -  v3 position of target [arcsec]
         updown -  undistorted target y-axis anti-parallel to V3? 
    Returns:
         troll  -  the roll angle evaluated at the target, i.e., the angle V3
                  makes with North at the target position [deg].
    Basic algorithm from C. Cox.
    """
    # convert all angles to radians
    m  = math
    pi = m.pi
    dr = pi / 180.
    roll *= dr
    dec  *= dr
    v2   *= (dr/3600.)
    v3   *= (dr/3600.)

    sin_rho = m.sqrt(m.sin(v2)**2 + m.sin(v3)**2 - m.sin(v2)**2 * m.sin(v3)**2)
    rho = m.asin(sin_rho)

    beta = m.asin(m.sin(v3)/sin_rho)
    gamm = m.asin(m.sin(v2)/sin_rho)

    if (v2<0):  beta = pi - beta
    if (v3<0):  gamm = pi - gamm

    A = (pi/2) + roll - beta

    abscissa = m.sin(dec)*sin_rho - m.cos(dec)*m.cos(rho)*m.cos(A)
    ordinate = m.sin(A)*m.cos(dec)
    B = m.atan2(ordinate,abscissa)
    
    troll_rad = pi - (gamm + B)

    # convert back to deg and flip upside-down if requested (as for ACS/WFC)
    if updown:
        troll = divmod((troll_rad/dr)+180,360)[1]
    else:
        troll = divmod((troll_rad/dr),360)[1]

    return troll

def makeCDdict(pa_deg,pscale):
    "Given PA(deg) and pixscale, retrun a CD matrix dictionary."
    
    cdmat = {}
    m  = math
    pa = pa_deg * m.pi / 180.
    
    cdmat['CD1_1'] = -pscale * m.cos(pa) / 3600.
    cdmat['CD1_2'] =  pscale * m.sin(pa) / 3600.
    cdmat['CD2_1'] =  pscale * m.sin(pa) / 3600.
    cdmat['CD2_2'] =  pscale * m.cos(pa) / 3600. 

    return cdmat
