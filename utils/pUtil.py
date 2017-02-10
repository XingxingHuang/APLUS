#!/usr/bin/env python

# $Id: pUtil.py,v 1.7 2003/11/11 21:50:14 anderson Exp $
# ---------------------------------------------------------------------
# pipeline utility functions for handling all kinds of par/catalog 
# files in ACS pipeline.

__version__      = '$Revision: 1.7 $ '[11:-3]
__version_date__ = '$Date: 2003/11/11 21:50:14 $ '[7:-3]
__author__       = "Ken Anderson, anderson@pha.jhu.edu"


import string,os,cPickle
import time
from   types  import *
import numpy.oldnumeric as Numeric

def deNAN(a,value=0.0):
     nans=Numeric.logical_not(Numeric.less(a,0.0)+Numeric.greater_equal(a,0.0))
     return Numeric.where(nans,value,a)

def fillHeader(infile):

    """
    This function is meant specifically to return a list of headers from ACS catalogs
    (eg Sextractor catalogs) where all gaps in column identifiers have been filled in.
    Only column number and column name are returned.  infile can be either a file 
    or a file object. Errors raised if file not found or infile is not of type FileType.
    
    """
    # Checking that infile is a file or a file object.

    if type(infile) == StringType :
        try:
            catalog = open(infile)
        except IOError,err:
            raise IOError,err

    elif type(infile) == FileType :
        catalog = infile                        # keeping names consistent catalog is the file object.
    else:
        raise TypeError, 'Object passed is not a filename, file object'

    #ok, got through that...or not.
    
    headers    = []
    newheaders = []
    
    while 1:
        newline = string.split(catalog.readline())
        if '###' in newline:
            continue
        if '##' in newline:
            continue
        elif '#' in newline:
            col,name = newline[1],newline[2]
            headers.append((int(col),name))
        else: break


    # This is a bit of a rewrite to this function which now has to deal with the 
    # possible header structure:
    #
    # 1 FLUX_APER
    # 2 FLUX_APER
    # 3 FLUX_APER
    # 4 FLUXERR_APER
    # 5 FLUXERR_APER
    # 6 FLUXERR_APER
    # 7 MAG_APER
    # 8 MAG_APER
    # 9 MAG_APER
    #10 MAGERR_APER
    #11 MAGERR_APER
    #12 MAGERR_APER
    #
    # which it probably should have been able to handle correctly anyway.
    # The usual SExtractor header would have looked like this:
    #
    # 1 FLUX_APER
    # 4 FLUXERR_APER
    # 7 MAG_APER
    #10 MAGERR_APER
    #
    # This function can now handle both properly.

    col,name = headers[0]
    same_name_counter = 0
    
    for f in range(1,len(headers)):
        nextcol,nextname = headers[f]
        params           = range(col,nextcol)
        param_count      = counter(col,nextcol)
        if len(params) == 1 and name == nextname:
            same_name_counter += 1
            newheaders.append((col,name,same_name_counter))
        elif len(params) == 1:
            if same_name_counter:
                same_name_counter += 1
                newheaders.append((col,name,same_name_counter))
            else: newheaders.append((col,name))
            same_name_counter = 0
        else:
            for f in param_count:
                newheaders.append((params[f-1],name,param_count[f-1]))
            same_name_counter = 0
        col  = nextcol
        name = nextname
    if headers[-1][1] == headers[-2][1]:
         newheaders.append((col,name,same_name_counter+1)) # append last element
    else:
         newheaders.append((col,name))
    return newheaders


def counter(n,m):

    """
    counter returns a list of values from 1 to m-n+1 for range counting.

    """

    return range(1,(m-n+1))

def makeHeaderDict(cat):
    """ Returns a dictionary of name:column pairs in the header from 
    SExtractor catalog.
    """
    headDict = {}
    for n in cat:
        newline = string.split(n)
        if '##' in newline:
            pass
        elif '#' in newline:
            fields = string.split(n)
            name = fields[2]
            col  = int(fields[1])
            headDict[name] = col
        else: break
    return headDict
     


def readInParFile(parfile):

    """ read a par file and return a dictionary of key:value pairs."""

    parFile = open(parfile)
    pars = {}
    while 1:
        newline = string.split(parFile.readline())
        if newline:
            pars[newline[0]] = newline[1]
        else:
            break
    
    return pars


def readParamFile(paramfile):

    """
    read a SExtractor .param file and return a list of key,value tuples. 
    Value is either 0 or 1. This is used instead of a dictionary
    as the order of params for output must be preserved for SExtractor 
    catalogs. 
    """

    paramFile = open(paramfile).readlines()
    params = []
    for f in range(len(paramFile)):
        line = string.split(paramFile[f])
        if '#' in line[0]:
            newline = string.split(line[0],'#')
            params.append( (newline[1],0) )
        else:
            params.append( (line[0],1) )
    return params

def rmFiles(dir,fileList):
    """Remove a list of files from directory dir."""
    curdir = os.getcwd()
    os.chdir(dir)
    nFiles = 0
    for f in fileList:
        try:
            os.remove(f)
            nFiles += 1
            print "Removed file",f
        except OSError,err:
            print "Remove file failed on",f
            print err
            pass
    os.chdir(curdir)
    return nFiles
    
class logFile:

    """ the class provides for a logfile object, which is instantiated by 
    the buildObs method of an DataSet object. buildObs passes the logFileName 
    attribute to the constructor. Write method appends a passed string
    object to the logfile. 
    """

    def __init__(self,newobspath):
        """ Initiate log file."""
        newobs = os.path.split(newobspath)[1]
        self.logFileName = newobspath + '/' + newobs + '_proc.log'
        logfile = open(self.logFileName,"w")
        logFileEntry = 'Log file started at: ' + ptime()
        logFileEntry = logFileEntry + '\n\nPipeline processing starting for observation:  ' + newobs
        logFileEntry = logFileEntry + '\n============================================\n'
        logfile.write(logFileEntry)
        logfile.close()

    def write(self,logFileEntry):
        """ write a string to the log file object.
        Method opens and closes the file for each call.
        """
        logfile = open(self.logFileName,"a")
        logfile.write('\n' + ptime() + '>  ' + logFileEntry)
        print ptime() + '>  ' + logFileEntry
        logfile.close()
        return

    
def ptime():
    """returns the current date/time in ISO 8601 standard format."""

    return time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime(time.time()))

def errorSearch(errout):
    """searches a passed list of strings for any of a set of error keywords.
    This list will usually be the result of the readlines() method
    called on a fromchild file object from popen2 having run a command.
    """
    errList = []
    errorKeys = ['Error','ERROR','WARNING','empty','Empty','keyword','no value',
          'not found','must be Y or N','unknown','too many members',
          'not enough', 'out of range','missing','at least','Not enough',
          'problem', 'Cannot','SYNTAX','NOT','overflow']


    for f in range(len(errout)):
        for SXerr in errorKeys:
            errstr = errout[f]
            if string.find(errstr,SXerr) != -1:
                newerr = string.replace(string.strip(errstr),">","")   # strip ">" from any error strings
                errList.append(newerr)

    try:
        return string_tally(errList)
    except IndexError:
        return


def string_tally(strings):
    """search a list of strings and find all unique strings and
    number of occurances for each found string.  Returns a 
    dictionary of unique strings (as keys) with the number of 
    occurrences found as the value.
    eg.,
    stringList = ['error1','error2','error3','error1','error1','error3']
    
    string_tally(stringList) returns
    {'error3': 2, 'error2': 1, 'error1': 3}
    """
    unique = {}
    unique[strings[0]] = 1
    for item in strings[1:]:
        if item in unique.keys():
            n = unique[item]
            unique[item] = n+1
        else:
            unique[item] = 1
    return unique

def sourcexy(catalogfile):
    """write out an ascii file of object number and x-y pixel coordinates
    from a SExtractor catalog.
    """
    outname = os.path.splitext(catalogfile)[0]+"_xy.list"
    ofile   = open(outname,"w")
    cat     = open(catalogfile).readlines()
    Ncol = ''
    xcol = ''
    ycol = ''
    ofile.write("N\tX_IMAGE  \tY_IMAGE\n###############################\n")
    for i in range(len(cat)):
        line = string.split(cat[i])
        # find the columns whereN,  X_IMAGE and Y_IMAGE are
        if line[0] == "#":
            if   line[2] == "NUMBER":
                Ncol = line[1]
            elif line[2] == "X_IMAGE":
                xcol = line[1]
            elif line[2] == "Y_IMAGE":
                ycol = line[1]
            else: continue
        else:
            ofile.write(line[int(Ncol)-1]+'\t'+line[int(xcol)-1]+'  \t'+line[int(ycol)-1]+'\n')
    ofile.close()
    return

def jarIt(file,obj):
    """writes a pickle file for a passed object. file can be just a file
    or a path."""
    file = open(file,'w')
    cPickle.dump(obj,file)
    print "pickled ",obj
    file.close()
    del file
    return



def _hackit(cat):
    """hack the detectionCatalog.cat file to take out a bunch of the aperture data. 
    Only include the first three apertures in the final catalog.
    This will hack the columns indicated by 

    MAG_APER
    MAGERR_APER
    FLUX_APER
    FLUXERR_APER
    """
    dir,old_file = os.path.split(cat)
    headerList = []
    headerList = fillHeader(cat)       # this returns a list of the catalog header.

    # go through the header and find the columns to keep.  We are looking for
    #'FLUX_APER', 1)
    #'FLUX_APER', 2)
    #'FLUX_APER', 3)
    #'FLUXERR_APER', 1)
    #'FLUXERR_APER', 2)
    #'FLUXERR_APER', 3)
    #'MAG_APER', 1)
    #'MAG_APER', 2)
    #'MAG_APER', 3)
    #'MAGERR_APER', 1)
    #'MAGERR_APER', 2)
    #'MAGERR_APER', 3)

    keep_apertures = [1,2,3]
    newheader = []
    for i in headerList:
        if len(i) == 2:
            newheader.append(i)
        else:
            if i[2] not in keep_apertures:
                continue
            else:
                newheader.append(i)

    #return newheader

    cols = []
    for i in newheader:
        cols.append(i[0]-1)

    new_rows = []
    for row in open(cat).readlines():
        if '#' in row: continue
        fields = row.split()
        arow = ''
        for column in cols:
            arow += '  '+fields[column]
        new_rows.append(arow)

    # OK, we have the newheader and the new data . 
    # We need to renumber the columns in the header.  the newheader
    # list has the old catalog's column identifiers and that needs to get fixed.
    
    new_newheader = []
    for i in range(len(newheader)):
        if len(newheader[i]) == 2:
            new_newheader.append((i+1,newheader[i][1]))
        else:
            new_newheader.append((i+1,newheader[i][1],newheader[i][2]))
    
    new_file=os.path.join(dir,"new_"+old_file)
    file = open(new_file,'w')
    for item in new_newheader:
        file.write('# '+str(item[0])+'\t'+str(item[1])+'\n')
    for row in new_rows:
        file.write(row+'\n')
    file.close()
    return
