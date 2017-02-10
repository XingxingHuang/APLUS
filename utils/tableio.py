#!/usr/bin/env python

# $Id: tableio.py,v 1.2 2003/11/11 21:51:02 anderson Exp $
# ---------------------------------------------------------------------
# adapted from useful.py (txitxo) to only include those ufuncs needed by
# apsis and other acs packages.
# K anderson 22-apr-2002


__version__      = '$Revision: 1.2 $ '[11:-3]
__version_date__ = '$Date: 2003/11/11 21:51:02 $ '[7:-3]

import string
import numpy

#Read/write headers

def get_header(file):
    """ Returns a string containing all the lines 
    at the top of a file which start by '#'"""
    buffer=''
    for line in open(file).readlines():
        if line[0]=='#': 
            buffer=buffer+line
        else: break
    return buffer

def put_header(file,text):
    """Adds text (starting by '#' and ending by '\n')
    to the top of a file."""
    if len(text)==0: 
        return
    if text[0]<>'#': 
        text='#'+text
    if text[-1]<>'\n':
        text=text+'\n'
    buffer=text+open(file).read()
    open(file,'w').write(buffer)

#Files containing strings

def get_str(file,cols=0,nrows='all'):
    """ 
        Reads strings from a file
        Usage: 
         x,y,z=get_str('myfile.cat',(0,1,2))
        x,y,z are returned as string lists
    """
    if type(cols)==type(0):
        cols=(cols,)
        nvar=1
    else: 
        nvar=len(cols)
    lista=[]
    for i in range(nvar): 
        lista.append([])
    buffer=open(file).readlines() 
    if nrows=='all': 
        nrows=len(buffer)
    counter=0
    for lines in buffer:
        if counter>=nrows :
            break
        pieces=string.split(lines)
        if len(pieces)==0: 
            continue
        if pieces[0][0]=='#': 
            continue 
        for j in range(nvar): 
            lista[j].append(pieces[cols[j]])
        counter=counter+1
    if nvar==1: 
        return lista[0]
    else: 
        return tuple(lista) 

def put_str(file,tupla):
    """ Writes tuple of string lists to a file
        Usage:
      put_str(file,(x,y,z))
    """    
    if type(tupla)<>type((2,)):
        raise 'Need a tuple of variables'

    f=open(file,'w')    

    for i in range(1,len(tupla)):
        if len(tupla[i])<>len(tupla[0]):
            raise 'Variable lists have different length'
    for i in range(len(tupla[0])):
        cosas=[]
        for j in range(len(tupla)):
            cosas.append(str(tupla[j][i]))
        f.write(join(cosas)+'\n')
    f.close()


#Files containing data

def get_data(file,cols=0,nrows='all'):
    """ Returns data in the columns defined by the tuple
    (or single integer) cols as a tuple of float arrays 
    (or a single float array)"""
    if type(cols)==type(0):
        cols=(cols,)
        nvar=1
    else: nvar=len(cols)

    data=get_str(file,cols,nrows)

    if nvar==1: 
        return numpy.array(map(float,data))
    else:
        data=list(data)
        for j in range(nvar): 
            data[j]=numpy.array(map(float,data[j]))
        return tuple(data) 

def put_data(file,variables,header='',format='',append='no'):
    """ Writes tuple of float variables to a file 
        Usage:
      put_data(file,(x,y,z),header,format)
    where header is any string  
        and format is a string of the type:
           '%f %f %i ' 
    """    
    if type(variables)<>type((2,)):
        raise 'Need a tuple of variables'
    if format==''   : format='%.6g  '*len(variables)
    if append=='yes': f=open(file,'a')
    else: f=open(file,'w')
    if header<>"":
        if header[0] <>'#' : header='#'+header
        if header[-1]<>'\n': header=header+'\n'
        f.write(header)
    for i in range(len(variables[0])):
        cosas=[]
        for j in range(len(variables)):
            cosas.append(variables[j][i])
        line=format % tuple(cosas)             
        f.write("\t"+line+'\n')
    f.close()

def put_data2(file,variables,variables2,header='',format='',format2='',append='no'):
    """ Writes tuple of float variables and variables2 (WCS) to a file 
        Usage:
      put_data2(file,(x,y,z),(ra,dec),header,format)
    where header is any string  
        and format is a string of the type:
           '%f %f %i ' 
    """    
    if type(variables)<>type((2,)):
        raise 'Need a tuple of variables'
    if type(variables2)<>type((2,)):
        raise 'Need a tuple of variables2'
    if format==''   : format='%.6g  '*len(variables)
    if format2==''   : format2='%s  '*len(variables2)
    if append=='yes': f=open(file,'a')
    else: f=open(file,'w')
    if header<>"":
        if header[0] <>'#' : header='#'+header
        if header[-1]<>'\n': header=header+'\n'
        f.write(header)
    for i in range(len(variables[0])):
        cosas=[]
        cosas2=[]
        for j in range(len(variables)):
            cosas.append(variables[j][i])
        #pdb.set_trace()
        line=format % tuple(cosas)
        for j in range(len(variables2)):
            cosas2.append(variables2[j][i])
        #pdb.set_trace()
        line2=format2 % tuple(cosas2)
        #Outline=line+" "+line2
        #f.write("\t"+line+'\n')
        f.write("\t"+line+" "+line2+'\n')
    f.close()
