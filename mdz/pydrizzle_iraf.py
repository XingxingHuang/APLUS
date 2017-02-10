import os
from pyraf import iraf

no = iraf.no
yes = iraf.yes

import xydrizzle

# Point to default parameter file for task
_parfile = 'dither$pydrizzle/pydrizzle.par'

######
# Set up Python IRAF interface here
######
def pydriz_iraf(input,output=None,kernel='square',units='cps',pixfrac=1.0,
        rotate=no,orient=None,psize=None,ra=None,dec=None,xsize=None,ysize=None,
        use_mask=yes,bits=0,idckey=None,clean=yes,save=no,build=yes,version=None):
    # Transform IRAF empty parameters to Python None when
    # PyDrizzle expects it.
    if output == '': output = None
    if kernel == '': kernel = None
    if bits == '': bits = None
    if idckey == '': idckey = None
    
    if ra == '': ra = None
    if dec == '': dec = None
    if orient == '': orient = None
    if xsize == '': xsize = None
    if ysize == '': ysize = None
    
    
    if not rotate:
        orient = None
           
    if not use_mask: 
        bits = None
    
    #print 'Running PyDrizzle Version ',version
    # Setup dictionary of parameters used by SkyField object
    _fpars = {'xsize':None,'ysize':None,'psize':None,'orient':None,'ra':None,'dec':None}
    _field=None
    # Check to see if any SkyField parameters were set by user.
    for _par in _fpars.keys():
        _eval = eval(_par)
        # If a parameter was set, then update dictionary with value
        # and set switch 
        if _eval != None:
            _fpars[_par] = float(_eval)

    # If we have user-specified parameters for the SkyField object,
    # then setup the object with those values...
    # We are counting how many _fpars values are set to None.
    # If less than full number of _fpars members, then user set something.
    if _fpars.values().count(None) < len(_fpars):
        _field = xydrizzle.SkyField(None)
        if _fpars['xsize'] != None:
            _shape = (int(_fpars['xsize']),int(_fpars['ysize']))
        else: _shape = None

        _field.set(shape=_shape,psize=_fpars['psize'],orient=_fpars['orient'],
                    ra=_fpars['ra'],dec=_fpars['dec'])
                
    drobj = xydrizzle.PyDrizzle(input,output=output,field=_field,kernel=kernel,
            units=units,pixfrac=pixfrac,bits=bits,idckey=idckey,clean=clean)
    drobj.run(save=save,build=build)


# Setup PyDrizzle as an IRAF task here 
# by setting up an absolute path to the parfile...
#_ospath = os.path
# File that gets picked up here is: iraffunctions.py
#_abspath = _ospath.split(_ospath.abspath(__file__))[0]
#parfile = os.path.join(_abspath,'pydrizzle.par')

parfile = iraf.osfn(_parfile)
pyd = iraf.IrafTaskFactory(taskname='xydrizzle', value=parfile, pkgname=PkgName, 
            pkgbinary=PkgBinary, function=pydriz_iraf)
