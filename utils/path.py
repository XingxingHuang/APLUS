#!/usr/bin/env python

# $Id: path.py,v 1.5 2003/01/14 22:04:14 anderson Exp $
# ---------------------------------------------------------------------
# module provides the class Env which defines the ACS pipeline environment.

__version__      = '$Revision: 1.5 $ '[11:-3]
__version_date__ = '$Date: 2003/01/14 22:04:14 $ '[7:-3]
__author__       = "Ken Anderson, anderson@pha.jhu.edu"


import os
class Env:
    """
    This class defines an environment for the ACS pipeline.
    The methods are self-evident.
    K Anderson 22-02-01
    """
    def __init__(self):
        self.pipeline = os.environ['PIPELINE']
        self.paths = {
            'INGEST'   : os.environ['INGEST'],
            'DATASETS' : os.environ['DATASETS'],
            'CONFIGS'  : os.path.join(self.pipeline,'configs'),
            'PARS'     : os.path.join(self.pipeline,'pars'),
            'ARCHMSG'  : os.path.join(self.pipeline,'ARCHMSG')}

        try:
            self.paths['ARCH_REQ'] = os.environ['ARCHREQUESTS']
        except KeyError:
            pass

        try:
            os.environ["DUST_DIR"]
        except KeyError:
            print "DUST_DIR environment variable not defined."

        print "Testing availability of required executables..."

        test_exes = [ "match","sex","dust_getval"]

        for exe in test_exes:
            if self.which(exe): continue
            else: raise EnvironmentError,exe+": command not found."
        print "Environment locked and loaded."
        
    def getenv(self, var):
        return self.paths[var]
    
    def setenv(self, var, value):
        self.paths[var] = value
        return

    def which(self,exe):
        foundPath = None
        paths = os.environ["PATH"].split(":")
        for path in paths:
            cmdPath = os.path.join(path,exe)
            if os.path.exists(cmdPath):
                foundPath = cmdPath
            else:
                continue
        return foundPath
