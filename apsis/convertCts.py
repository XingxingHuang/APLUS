#!/usr/bin/env python
import sys
import re
import numpy

import pyfits
print "Using pyfits from " + str(pyfits.__path__[0])


def convert_fits(self, file):
  """Check if the argument is a fits file having atleast one extension that has its 
  BUNIT set to counts/s or electrons/s

  Returns nothing
  """
  print 'convert cts called with ' + file  # WZ Apr 2013
  fitsfile = pyfits.open(file)

  #Check for simple fits vs mef
  if len(fitsfile) == 1:
    ext = 0
  else:
    ext = 1

  # Get the exposure time from the Primary HDU "only"
  exptime = fitsfile[0].header.get('EXPTIME')

  for ext in fitsfile:
    if ext.name == 'SCI' or ext.name == 'ERR':
      unit = ext.header.get('BUNIT')
      #TODO: Move the regular expression out of this block; also check
      #to see if we can modify the file-open permissions in transit to
      #avoid having to close it and then re-open it
      if re.search('^(counts?\/s|electrons?\/s)$', unit, re.IGNORECASE):
        fitsfile.close()
        _convert(self, file, exptime)
        break


def _convert(self, file, exptime):
  """Update all extensions and corresponding headers of the file which have their 
  BUNIT set to counts/s or electrons/s. 

  Returns nothing
  """
  print "Updating " + file + " with EXPTIME - " + str(exptime)
  fitsfile = pyfits.open(file, mode='update')
	
  for ext in fitsfile:
    if ext.name == 'SCI' or ext.name == 'ERR':
      unit = ext.header.get('BUNIT')
      if re.search('^(counts?\/s|electrons?\/s)$', unit, re.IGNORECASE):
        print "Updating extension - " + ext.name
        # Update the header
        ext.header.update('BUNIT', 'COUNTS')
        # Update the data
        ext.data = numpy.multiply(exptime, ext.data)

  if self:
    self.logfile.write('Converted file ' + file + ' from cts/s to cts')

	# Explicitly close the file to ensure write-out
  fitsfile.close()


# Main processing
if __name__ == '__main__':
  for file in sys.argv[1:]:
    convert_fits(None, file)

