# Module filter
import re

__author__ = "Wei Zheng"
__version__ = 1.0

"""
Example usage:
>>>import filter
>>>WFC3filters = filter.WFC3FilterInfo()
>>>WFC3filters.allAttrs()
['apsize', 'camera', 'position', 'description', 'cw_range', 'ftype', 'width', 'wheel']
>>>WFC3filters.searchFilters('ftype', 'broad')
['F225W', 'F275W', 'F336W', 'F390W', 'F200LP', 'F350LP','F606W','F814W','F105W', 'F110W', 'F125W', 'F140W', 'F160W','F098M']
or (see WFC3FilterInfo class for predefined WFC3 filter classifications)
>>>WFC3filters.broad
['F225W', 'F275W', 'F336W', 'F390W', 'F200LP', 'F350LP','F606W','F814W','F105W', 'F110W', 'F125W', 'F140W', 'F160W','F098M']
"""
class FilterInfo:
    def __init__(self, filterdict):
        self.filters = filterdict
    
    def allFilters(self):
        """Returns a list of all filters"""
        return self.filters.keys()

    def allAttrs(self):
        """Returns a list of filter attributes"""
        return self.filters[self.allFilters()[0]].keys()
    
    def getFilter(self, filter):
        """Returns a dictionary of a filters attrubutes"""
        return self.filters[filter]
    
    def getFilterAttr(self, filter, attr):
        """Returns value of a filters attrubute"""
        return self.filters[filter][attr]

    def searchFilters(self, attr, value):
        """Returns a list of filters that match a filters attribute value"""
        filterlist = []
        p = re.compile(value)
        for filter in self.filters.keys():
            s = self.filters[filter][attr]
            m = p.search(s)
            if m:
                filterlist.append(filter)
        return filterlist

class WFC3FilterInfo(FilterInfo):
    """An WFC3 specific filter class"""
    def __init__(self):
        WFC3Filters = {
        'F225W':{'wheel':'3', 'position':'2', 'ftype':'nuv', 'apsize':'large', 'cw_range':'2341.0', 'width':'547.3','description':'UV1', 'camera':'UVIS'},
        'F275W':{'wheel':'9', 'position':'1', 'ftype':'nuv', 'apsize':'large', 'cw_range':'2715.3', 'width':'480.8','description':'UV2', 'camera':'UVIS'},
        'F336W':{'wheel':'4', 'position':'1', 'ftype':'nuv', 'apsize':'large', 'cw_range':'3361.1', 'width':'553.8','description':'U', 'camera':'UVIS'},
        'F390W':{'wheel':'11', 'position':'1', 'ftype':'nuv', 'apsize':'large', 'cw_range':'3604.6', 'width':'953.0','description':'C', 'camera':'UVIS'},
        'F200LP':{'wheel':'11', 'position':'2', 'ftype':'broadband', 'apsize':'large', 'cw_range':'4883', 'width':'5002','description':'C', 'camera':'UVIS'}, # WZ TBD
        'F350LP':{'wheel':'11', 'position':'2', 'ftype':'broadband', 'apsize':'large', 'cw_range':'5846', 'width':'4758','description':'C', 'camera':'UVIS'}, 
        # 'F350LP':{'wheel':'11', 'position':'2', 'ftype':'broadband', 'apsize':'large', 'cw_range':'6812.0', 'width':'4000','description':'V', 'camera':'UVIS'},
        'F606W':{'wheel':'6', 'position':'2', 'ftype':'broadband', 'apsize':'large', 'cw_range':'5887.0', 'width':'2182','description':'V', 'camera':'UVIS'},
        'F814W':{'wheel':'7', 'position':'2', 'ftype':'broadband', 'apsize':'large', 'cw_range':'8024.0', 'width':'1536','description':'I', 'camera':'UVIS'},
        'F098M':{'wheel':'1', 'position':'6', 'ftype':'broadband', 'apsize':'large', 'cw_range':'9864.0', 'width':'1570.0','description':'M', 'camera':'IR'},
        'F105W':{'wheel':'1', 'position':'5', 'ftype':'broadband', 'apsize':'large', 'cw_range':'1048.9', 'width':'292.3','description':'Y', 'camera':'IR'},
        'F110W':{'wheel':'1', 'position':'4', 'ftype':'broadband', 'apsize':'small', 'cw_range':'1141.4', 'width':'503.4','description':'YJ', 'camera':'IR'},
        'F125W':{'wheel':'1', 'position':'2', 'ftype':'broadband', 'apsize':'large', 'cw_range':'1245.9', 'width':'301.5','description':'J', 'camera':'IR'},
        'F140W':{'wheel':'1', 'position':'8', 'ftype':'broadband', 'apsize':'large', 'cw_range':'1392.1', 'width':'399.0','description':'JH', 'camera':'IR'},
        'F160W':{'wheel':'1', 'position':'1', 'ftype':'broadband', 'apsize':'large', 'cw_range':'1540.52', 'width':'297.8','description':'H', 'camera':'IR'},
        'F126N':{'wheel':'1', 'position':'15', 'ftype':'narrowband', 'apsize':'large', 'cw_range':'1\
258.46', 'width':'15.2','description':'FeII', 'camera':'IR'},
        }
        FilterInfo.__init__(self, WFC3Filters)

class ACSFilterInfo(FilterInfo):
    """An ACS specific filter class"""
    def __init__(self):
        ACSFilters = {
        'CLEAR1L':{'wheel':'1', 'position':'1', 'ftype':'clear', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Clear WFC aperture', 'camera':'WFC'},
        'F555W':{'wheel':'1', 'position':'2', 'ftype':'broadband', 'apsize':'large', 'cw_range':'534.6', 'width':'119.3','description':'Johnson V', 'camera':'WFC,HRC'},
        'F775W':{'wheel':'1', 'position':'3', 'ftype':'sdss', 'apsize':'large', 'cw_range':'776.4', 'width':'152.8','description':'SDSS I', 'camera':'WFC,HRC'},
        'F625W':{'wheel':'1', 'position':'4', 'ftype':'sdss', 'apsize':'large', 'cw_range':'631.8', 'width':'144.2','description':'SDSS r', 'camera':'WFC,HRC'},
        'F550M':{'wheel':'1', 'position':'5', 'ftype':'narrowband', 'apsize':'large', 'cw_range':'558', 'width':'54.7','description':'Narrow V', 'camera':'WFC,HRC'},
        'F851LP':{'wheel':'1', 'position':'6', 'ftype':'sdss', 'apsize':'large', 'cw_range':'944.5', 'width':'122.9','description':'SDSS z', 'camera':'WFC,HRC'},
        'CLEAR1S':{'wheel':'1', 'position':'7', 'ftype':'clear', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'Clear HRC aperture', 'camera':'HRC'},
        'POL0UV':{'wheel':'1', 'position':'8', 'ftype':'polarizer', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'UV Polarizer 0', 'camera':'WFC,HRC'},
        'POL60UV':{'wheel':'1', 'position':'9', 'ftype':'polarizer', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'UV Polarizer 60', 'camera':'WFC,HRC'},
        'POL120UV':{'wheel':'1', 'position':'10', 'ftype':'polarizer', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'UV Polarizer 120', 'camera':'WFC,HRC'},
        'F892N':{'wheel':'1', 'position':'11', 'ftype':'narrowband', 'apsize':'small', 'cw_range':'891.7', 'width':'15.4','description':'Methane', 'camera':'WFC,HRC'},
        'F606W':{'wheel':'1', 'position':'12', 'ftype':'broadband', 'apsize':'large', 'cw_range':'590.7', 'width':'234.2','description':'Broad V', 'camera':'WFC,HRC'},
        'F502N':{'wheel':'1', 'position':'13', 'ftype':'narrowband', 'apsize':'large', 'cw_range':'502.2', 'width':'5.7','description':'OIII', 'camera':'WFC,HRC'},
        'G800L':{'wheel':'1', 'position':'14', 'ftype':'grism', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'GRISM', 'camera':'WFC,HRC'},
        'F658N':{'wheel':'1', 'position':'15', 'ftype':'narrowband', 'apsize':'large', 'cw_range':'658.4', 'width':'7.3','description':'H alpha', 'camera':'WFC,HRC'},
        'F475W':{'wheel':'1', 'position':'16', 'ftype':'sdss', 'apsize':'large', 'cw_range':'476.0', 'width':'145.8','description':'SDSS g', 'camera':'WFC,HRC'},
        'CLEAR2L':{'wheel':'2', 'position':'1', 'ftype':'clear', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Clear WFC aperture', 'camera':'WFC'},
        'F660N':{'wheel':'2', 'position':'2', 'ftype':'narrowband', 'apsize':'large', 'cw_range':'660.2', 'width':'3.5','description':'NII', 'camera':'WFC,HRC'},
        'F814W':{'wheel':'2', 'position':'3', 'ftype':'broadband', 'apsize':'large', 'cw_range':'833.3', 'width':'251.1','description':'Broad I', 'camera':'WFC,HRC'},
        'FR388N':{'wheel':'2', 'position':'4m', 'ftype':'ramp', 'apsize':'large', 'cw_range':'371-405', 'width':'2%','description':'OII Ramp (middle)', 'camera':'WFC,HRC'},
        'FR423N':{'wheel':'2', 'position':'4i', 'ftype':'ramp', 'apsize':'large', 'cw_range':'405-442', 'width':'2%','description':'OII Ramp (inner)', 'camera':'WFC'},
        'FR462N':{'wheel':'2', 'position':'4o', 'ftype':'ramp', 'apsize':'large', 'cw_range':'442-482', 'width':'2%','description':'OII Ramp (outer)', 'camera':'WFC'},
        'F435W':{'wheel':'2', 'position':'5', 'ftype':'broadband', 'apsize':'large', 'cw_range':'429.7', 'width':'103.8','description':'Johnson B', 'camera':'WFC,HRC'},
        'FR656N':{'wheel':'2', 'position':'6m', 'ftype':'ramp', 'apsize':'large', 'cw_range':'627-685', 'width':'2%','description':'H a  Ramp (middle)', 'camera':'WFC,HRC'},
        'FR716N':{'wheel':'2', 'position':'6i', 'ftype':'ramp', 'apsize':'large', 'cw_range':'685-747', 'width':'2%','description':'H a  Ramp (inner)', 'camera':'WFC'},
        'FR782N':{'wheel':'2', 'position':'6o', 'ftype':'ramp', 'apsize':'large', 'cw_range':'747-816', 'width':'2%','description':'H a  Ramp (outer)', 'camera':'WFC'},
        'CLEAR2S':{'wheel':'2', 'position':'7', 'ftype':'clear', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'Clear HRC aperture', 'camera':'HRC'},
        'POL0V':{'wheel':'2', 'position':'8', 'ftype':'polarizer', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'Visible Polarizer 0', 'camera':'WFC,HRC'},
        'F330W':{'wheel':'2', 'position':'9', 'ftype':'nuv', 'apsize':'small', 'cw_range':'335.4', 'width':'58.8','description':'HRC u', 'camera':'HRC'},
        'POL60V':{'wheel':'2', 'position':'10', 'ftype':'polarizer', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'Visible Polarizer 60', 'camera':'WFC,HRC'},
        'F250W':{'wheel':'2', 'position':'11', 'ftype':'nuv', 'apsize':'small', 'cw_range':'269.6', 'width':'54.9','description':'Near-UV filter', 'camera':'HRC'},
        'POL120V':{'wheel':'2', 'position':'12', 'ftype':'polarizer', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'Visible Polarizer 120', 'camera':'WFC,HRC'},
        'PR200L':{'wheel':'2', 'position':'13', 'ftype':'prism', 'apsize':'small', 'cw_range':'None', 'width':'None','description':'HRC PRISM', 'camera':'HRC'},
        'F344N':{'wheel':'2', 'position':'14', 'ftype':'narrowband', 'apsize':'small', 'cw_range':'343.4', 'width':'6','description':'NeV', 'camera':'HRC'},
        'F220W':{'wheel':'2', 'position':'15', 'ftype':'nuv', 'apsize':'small', 'cw_range':'222.8', 'width':'48.5','description':'Near-UV filter', 'camera':'HRC'},
        'FR853N':{'wheel':'2', 'position':'16i', 'ftype':'ramp', 'apsize':'large', 'cw_range':'816-891', 'width':'2%','description':'IR Ramp (inner)', 'camera':'WFC'},
        'FR914M':{'wheel':'2', 'position':'16m', 'ftype':'ramp', 'apsize':'large', 'cw_range':'757-1071', 'width':'9%','description':'Broad Ramp (middle)', 'camera':'WFC,HRC'},
        'FR931N':{'wheel':'2', 'position':'16o', 'ftype':'ramp', 'apsize':'large', 'cw_range':'891-972', 'width':'2%','description':'IR Ramp (outer)', 'camera':'WFC'},
        'FR647M':{'wheel':'2', 'position':'17i', 'ftype':'ramp', 'apsize':'large', 'cw_range':'537-757', 'width':'9%','description':'Broad Ramp (inner)', 'camera':'WFC'},
        'FR459M':{'wheel':'2', 'position':'17m', 'ftype':'ramp', 'apsize':'large', 'cw_range':'381-537', 'width':'9%','description':'Broad Ramp (middle)', 'camera':'WFC,HRC'},
        'FR1016N':{'wheel':'2', 'position':'17o', 'ftype':'ramp', 'apsize':'large', 'cw_range':'972-1061', 'width':'2%','description':'IR Ramp (outer)', 'camera':'WFC'},
        'FR505N':{'wheel':'2', 'position':'18m', 'ftype':'ramp', 'apsize':'large', 'cw_range':'482-527', 'width':'2%','description':'OIII Ramp (middle)', 'camera':'WFC,HRC'},
        'FR551N':{'wheel':'2', 'position':'18i', 'ftype':'ramp', 'apsize':'large', 'cw_range':'527-575', 'width':'2%','description':'OIII Ramp (inner)', 'camera':'WFC'},
        'FR601N':{'wheel':'2', 'position':'18o', 'ftype':'ramp', 'apsize':'large', 'cw_range':'575-627', 'width':'2%','description':'OIII Ramp (outer)', 'camera':'WFC'},
        'BLOCK1':{'wheel':'3', 'position':'1', 'ftype':'block', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Blocks path to the SBC', 'camera':'SBC'},
        'F115LP':{'wheel':'3', 'position':'2', 'ftype':'fuv', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'MgF 2', 'camera':'SBC'},
        'F125LP':{'wheel':'3', 'position':'3', 'ftype':'fuv', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'CaF 2', 'camera':'SBC'},
        'BLOCK2':{'wheel':'3', 'position':'4', 'ftype':'block', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Blocks path to the SBC', 'camera':'SBC'},
        'F140LP':{'wheel':'3', 'position':'5', 'ftype':'fuv', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'BaF 2', 'camera':'SBC'},
        'F150LP':{'wheel':'3', 'position':'6', 'ftype':'fuv', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Crystal Quartz', 'camera':'SBC'},
        'BLOCK3':{'wheel':'3', 'position':'7', 'ftype':'block', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Blocks path to the SBC', 'camera':'SBC'},
        'F165LP':{'wheel':'3', 'position':'8', 'ftype':'fuv', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Dynasil', 'camera':'SBC'},
        'F122M':{'wheel':'3', 'position':'9', 'ftype':'fuv', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Lyman alpha', 'camera':'SBC'},
        'BLOCK4':{'wheel':'3', 'position':'10', 'ftype':'block', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'Blocks path to the SBC', 'camera':'SBC'},
        'PR130L':{'wheel':'3', 'position':'11', 'ftype':'prism', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'CaF 2  Prism', 'camera':'SBC'},
        'PR110L':{'wheel':'3', 'position':'12', 'ftype':'prism', 'apsize':'large', 'cw_range':'None', 'width':'None','description':'LiF 2  Prism', 'camera':'SBC'}
        }
        FilterInfo.__init__(self, ACSFilters)


        self.all = self.allFilters()
        self.broadband = self.searchFilters('ftype', 'broadband')
        self.wfc = self.searchFilters('camera', 'WFC')
        self.uvis = self.searchFilters('camera', 'UVIS')
        self.ir   = self.searchFilters('camera', 'IR')
        self.wfc_imaging = ['F435W','F475W','F502N','F550M','F555W','F606W','F625W','F658N','F660N','F775W','F814W','F850LP']
        self.uvis_imaging = ['F225W','F275W','F336W','F390W','F200LP','F350LP','F606W','F814W']
        self.ir_imaging   = ['F098M','F105W','F110W','F125W','F140W','F160W','F126N']
        self.sdss = self.searchFilters('ftype', 'sdss')
        self.narrowband = self.searchFilters('ftype', 'narrowband')
        self.polarizer = self.searchFilters('ftype', 'polarizer')
        self.ramp = self.searchFilters('ftype', 'ramp')
        self.nuv = self.searchFilters('ftype', 'nuv')
        self.fuv = self.searchFilters('ftype', 'fuv')
        self.grism = self.searchFilters('ftype', 'grism')
        self.prism = self.searchFilters('ftype', 'prism')
        self.spectroscopic = self.grism + self.prism
        self.clear = self.searchFilters('ftype', 'clear')
        self.block = self.searchFilters('ftype', 'block')
        self.small = self.searchFilters('apsize', 'small')
        self.large = self.searchFilters('apsize', 'large')
