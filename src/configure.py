#!/usr/bin/env python
from __future__ import division,print_function
import numpy as np

"""
Collect information needed for astrometric and photometric matching of a collection of catalogs.
Put this all into various FITS tables that succeeding C++ programs can access easily.
"""


"""
FITS column names will be given in ALL CAPS.  These are also the names of "Attributes" that
are given to each extension.
In the yaml file, we'll do First Letter Cap for the 3 big categories, but expect
lowercase for all map keys otherwise.

Fields:
- name:
  coords:  
  radius:
  pmepoch:
- name:...

Files:
- glob: glob for files
  expkey: "< keyword" or _FILENAME or value for EXPOSURE to give each catalog in the file
  translation: regex to apply to the id that is derived [None]
- glob: ...

The glob for files can include numerical ranges, e.g.
glob: DECam_00(236600-236700).fits.*
will look for files with each value in this range (including 00236700).

Attributes:
- key: the column name for the attribute
  value: value or '< keyword' to look for in header
  select: regex that exposure name must match to apply this value [None]
  vtype: str, float, int - data type that the we should get for value [default=str]
  default: value to assign if the specified header keyword does not exist [None]
  translation: regex to apply to header values if found [None]
- key: 

Note that all string-valued Attributes can include "${XXX}" anywhere and this will be replaced
by the value of Attribute with name XXX.  This is *not* done for EXPOSURE, EXTENSION.

Attributes we'll require as WCSFoF needs them:
'FILENAME','EXTENSION','INSTRUMENT','DEVICE',
           'FIELD','EXPOSURE','RA','DEC','WCSIN','XKEY',
           'YKEY','IDKEY'

WCSFit requires:
'BAND'

PhotoFit requires:
'AIRMASS','MJD'

and some WCSFit and PhotoFit models will use
'EPOCH'

Special values:
EXPOSURE = _FILENAME  to pull the exposure name from the filename
FIELD = _NEAREST      to assign to nearest Field on the sky
IDKEY = _ROW          to have table row number be the object ID number
WCSIN = _ICRS         to have the input pixel coords be RA, Dec in degrees
WCSIN = _HEADER       to pull a FITS WCS out of the headers as starting WCS

translations have format <regex>=<replace>
"""

import sys
import yaml
import astropy.io.fits as pf
import astropy.coordinates as cc
import astropy.units as u
import re
import importlib
import glob
import subprocess
import gbutil

import astropy
# Major annoyance: incompatible interfaces in astropy
oldAstropy = int(astropy.__version__.split('.')[0])==0 \
  and int(astropy.__version__.split('.')[1])<3
def getDegree(a):
    if oldAstropy:
        return a.degrees
    else:
        return a.degree

def getHeader(fits, iextn, filename):
    return fits[iextn].header

# These instrument names have special meaning.  "Observations" with these instruments
# do not require a device name.
# They have special indices: **** Keep in sync with FitSubroutines.h values! ****
specialInstruments={'REFERENCE':-1,'TAG':-3,   # backwards-compatible names
                    '_REFERENCE':-1,'_PMCAT':-2,'_TAG':-3}

# These are the extension attributes that are processed in some special way before
# being written to the output table.
specialAttributes =('FILENAME','EXTENSION','EXPOSURE',\
                    'DEVICE','RA','DEC','WCSIN')

# This list enumerates all of the attributes that we will pull from extensions up to
# the Exposure table.  The first element of the tuple is the attribute name.
# Second member is the tolerance for matching of each extension to the Exposure.
#   Value of None requires an exact match. Mismatches trigger errors.
# Third member is set to True if we want to retain the attribute in the Extensions,
#   otherwise it is popped from the extension's dictionary
# Fourth member is default value if no data given for an Exposure.  If None then
#   absence of data will generate an error, i.e. required column.
exposureLevelAttributes = [ ('INSTRUMENT',None,  True,  None),\
                            ('AIRMASS',   0.001, False, 1.),\
                            ('BAND',      None,  True,  None),\
                            ('EXPTIME',   0.0002,False, 1.),\
                            ('FIELD',     None,  False, None),\
                            ('MJD',       0.0002,False, 0.),\
                            ('EPOCH',     None,  False, ""),\
                            ('WEIGHT',    0.001, False, 1.),\
                            ('MAGWEIGHT', 0.001, False, 1.),\
                            ('APCORR',    0.0001,False, 0.),\
                            ('SYSERRMMAG',0.1,   False, 0.),\
                            ('SYSERRMAS', 0.1,   False, 0.),\
                            ('SYSERRXX',  0.001, False, 0.),\
                            ('SYSERRYY',  0.001, False, 0.),\
                            ('SYSERRXY',  0.001, False, 0.),\
                            ('OBSX',      1e-7,  False, 0.),\
                            ('OBSY',      1e-7,  False, 0.),\
                            ('OBSZ',      1e-7,  False, 0.) ]
# SYSERRMMAG is photometric sys error (in mmag), SYSERR[XY][XY] are astrometric error
# covariances (in mas^2). Alternatively SYSERRMAS is a circular astrometric error (in mas).
# OBS[XYZ] are observatory position at midpoint time of the exposure, in barycentric ICRS AU
        
headerIndicator = "<"
# Character at start of a attribute value that indicates lookup in header                    
                    
def py_to_fits(val):
    """
    Return a fits type-specification string corresponding to type of first
    element in the python array.
    Use the type of the first entry in the list
    """

    vtype=type(val[0])
    
    if(vtype is int):
        return 'J'
    elif(vtype is float):
        return 'D'
    elif(vtype is str):
        size=len(max(val,key=len))
        return 'A'+str(size)
    else:
        print("ERROR: No known FITS type for python type",vtype)
        sys.exit(1)
    

class Field:
    def __init__(self, name, coords, radius, pm_epoch=0., index=-1):
        # Expecting construction from a dictionary
        self.name = name
        if oldAstropy:
            self.coords = ICRS(coords,unit=(u.hourangle,u.degree))
        else:
            self.coords = cc.SkyCoord(coords,unit=(u.hourangle,u.degree),frame='icrs')
        self.radius = float(radius)
        self.pm_epoch = float(pm_epoch)
        self.index = int(index)
        return

    def __eq__(self, other):
        return self.coords==other.coords and self.radius==other.radius \
               and self.pm_epoch==other.pm_epoch

    def distance(self,location):
        return getDegree(self.coords.separation(location))

def buildFieldTable(fields):
    """
    Make a FITS BinTableHDU containing information on all the fields,
    and assign them indices as we go
    """
    name=[]
    ra  =[]
    dec =[]
    radius = []
    pmepoch = []
    
    index = 0
    for k,v in fields.items():
        name.append(k)
        ra.append(getDegree(v.coords.ra))
        dec.append(getDegree(v.coords.dec))
        radius.append(v.radius)
        pmepoch.append(v.pm_epoch)
        v.index = index
        index += 1
        
    hdu = pf.BinTableHDU.from_columns(\
        pf.ColDefs( [pf.Column(name='NAME',format=py_to_fits(name),array=name),
                     pf.Column(name='RA',format=py_to_fits(ra),array=ra),
                     pf.Column(name='DEC',format=py_to_fits(dec),array=dec),
                     pf.Column(name='RADIUS',format=py_to_fits(radius),array=radius),
                     pf.Column(name='PM_EPOCH',format=py_to_fits(pmepoch),array=pmepoch)]),
                     name = 'Fields')
    return hdu

class Instrument:
    def __init__(self, name, index=-4, band=None):
        # Start with an empty dictionary of device names
        self.name = name
        self.band = band
        self.devices = {}
        self.index = index
        return
    def addDevice(self,devName):
        if devName not in self.devices:
            self.devices[devName] = len(self.devices)
        return

def buildInstrumentTable(inst):
    """
    Return a FITS BinTableHDU for this instrument
    """
    names = ["" for k in inst.devices]
    for k,v in inst.devices.items():
        names[v] = k
    for n in names:
        # If any of the device names is null, then the indices are screwed up
        if n=="":
            print('ERROR: Device indices screwed up for instrument',inst.name)
            sys.exit(1)
    # Make an array of zeros to build the [xy]Min/Max columns
    z = [0. for n in names]
    hdu = pf.BinTableHDU.from_columns(\
        pf.ColDefs( [pf.Column(name='NAME',format=py_to_fits(names),array=names),
                     pf.Column(name='XMIN',format=py_to_fits(z),array=z),
                     pf.Column(name='XMAX',format=py_to_fits(z),array=z),
                     pf.Column(name='YMIN',format=py_to_fits(z),array=z),
                     pf.Column(name='YMAX',format=py_to_fits(z),array=z)]),
                     name='Instrument')
    hdu.header['NAME'] = inst.name
    if inst.band is None:
        hdu.header['BAND'] = inst.name
    else:
        hdu.header['BAND'] = inst.band
    hdu.header['NUMBER'] = inst.index
    #hdu.header['EXTNAME'] = 'Instrument'
    return hdu


def extractExposureLevelAttributes(extn, expo):
    '''Extract exposure attributes from an extension, or make sure that extension
       matches previously assigned attributes, and pop attributes from the extension
       before returning.  Both arguments are dicts.
    '''
    # Process RA & DEC of extension into coordinates
    if 'RA' in extn:
        ra = extn.pop('RA')
    else:
        ra = None
    if 'DEC' in extn:
        dec = extn.pop('DEC')
    else:
        dec = None
    if ra is not None and dec is not None:
        try:
            # If ra is valid floating point, assume it's degrees
            float(ra)
            coords = cc.SkyCoord(ra=ra,dec=dec,unit=(u.degree,u.degree),\
                                 frame='icrs')
        except ValueError:
            # Otherwise it must be HH:MM:SS.sss style in hours
            coords = cc.SkyCoord(ra=ra,dec=dec,unit=(u.hourangle,u.degree),\
                                 frame='icrs')
        # Check agreement of coordinates with any previous location
        if 'coords' in expo:
            if getDegree(expo['coords'].separation(coords)) > 1.:
                print("ERROR: RA/Dec mismatch at exposure",extn['EXPOSURE'],\
                          "from file",extn['FILENAME'],
                          "extension",extn['EXTENSION'])
                sys.exit(1)
        else:
            expo['coords'] = coords

    # Extract all desired fields from extension
    for key,tolerance,keep,default in exposureLevelAttributes:
        if key not in extn:
            continue
        if key in expo:
            # Already have a value for this key in the Exposure.  Check for agreement
            if (tolerance is None and expo[key]!=extn[key]) or \
               (tolerance is not None and np.abs(expo[key]-extn[key])>tolerance):
                print("Mismatch between Exposure and Extension values of key",key,\
                         "at exposure",extn['EXPOSURE'],\
                         "vs filename",extn['FILENAME'],\
                         "extension",extn['EXTENSION'])
                print("Values: ",expo[key],"vs",extn[key])
                sys.exit(1)
        else:
            # New key for this Exposure, adopt Extension value
            expo[key] = extn[key]
        if not keep:
            extn.pop(key)
    return
    
def buildExposureTable(exposures, fields, instruments):
    """
    Return a FITS BinTableHDU with information on exposures.
    Also assign indices to each according to their row in the table.
    """

    # First do some processing on some attributes and
    # Collect union of all attributes of all Exposures
    index = 0
    allAttributes = set()
    for k,e in exposures.items():
        # Save a name and running index number
        e['NAME'] = k;
        e['index'] = index;
        index = index+1

        # Convert coords to RA and Dec
        if 'coords' not in e:
            print("Exposure",k,"does not have coordinates")
            sys.exit(1)
        coords = e.pop('coords')
        e['RA'] = getDegree(coords.ra)
        e['DEC'] = getDegree(coords.dec)
        
        # Convert field name to index number
        if 'FIELD' not in e:
            print("Exposure",k,"does not have FIELD")
            sys.exit(1)
        f = e.pop('FIELD')
        e['FIELDNUMBER'] = fields[f].index;

        # Convert instrument to an index number
        if 'INSTRUMENT' not in e:
            print("Exposure",k,"does not have INSTRUMENT")
            sys.exit(1)
        f = e.pop('INSTRUMENT')
        if f in specialInstruments:
            e['INSTRUMENTNUMBER'] = specialInstruments[f]
        else:
            e['INSTRUMENTNUMBER'] = instruments[f].index
        
        allAttributes.update(e.keys())
    
    # Remove from the list of attributes we are putting into tile
    allAttributes.discard('index')
    # Fill in defaults for all exposures, and check for missing required ones
    for key,tolerance,keep,default in exposureLevelAttributes:
        if key not in allAttributes:
            # Don't need to consider an attribute that will not go into output table
            continue
        for expName,expo in exposures.items():
            if key not in expo:
                if default is None:
                    print("Exposure",expName,"is missing required attribute",key)
                    sys.exit(1)
                else:
                    expo[key] = default

    # Make FITS columns for all attributes
    cols = []
    for key in allAttributes:
        data = [expo[key] for k,expo in exposures.items()]
        cols.append(pf.Column(name=key, format=py_to_fits(data), array=data))
    hdu = pf.BinTableHDU.from_columns(cols, name = 'Exposures')
    return hdu

class AttributeFinder:
    def __init__(self, key, value, vtype="str", default=None, translation=None, select=None):
        self.key = key
        if select != None:
            # regular expression that must match the input NAME for the Attribute to apply
            self.select = re.compile(select)
        else:
            self.select=None

        # Determine the Python builtin type that the attribute must have, defaults to str
        if sys.version_info >= (3,0):
            self.valueType = getattr(importlib.import_module('builtins'),vtype)
        else:
            self.valueType = getattr(importlib.import_module('__builtin__'),vtype)

        # Either get the value of the attribute, or information for extracting it from headers
        self.headerKey = None
        if type(value)==str and len(value.strip())>0 and value.strip()[0]==headerIndicator:
            # The given value is a string indicating a header keyword lookup
            self.headerKey = value.strip()[1:].strip()
            if default == None:
                self.default = None
            else:
                # There is a default to assign if the keyword is absent
                self.default = self.valueType(default)
            if translation==None:
                self.translation=None
            else:
                # The header value should be read as a string and run through a regex/replace
                rr = translation.split('=')
                if len(rr)!=2:
                    print("ERROR: Attribute translator is not of form <regex>=<replace>: ",\
                      initdict['Translation'])
                    sys.exit(1)
                self.translation = (re.compile(rr[0]),rr[1])
        else:
            # There is a value given, cast it to desired type
            self.value = self.valueType(value)
        return

    def __call__(self,name, extnHeader=None, primaryHeader=None):
        # Return the value of the attribute if name passes the select, otherwise return None.
        # If Attribute is to be sought in header, will do so from extnHeader first, then primaryHeader
        # if not found, then assign default if one exists.  Otherwise return None
        if self.select!=None and self.select.match(name)==None:
            # This name does not pass the select
            return None
        if self.headerKey==None:
            # Simple fixed-value attribute
            return self.value

        val = None
        # Look in both headers for a value
        if extnHeader!=None and self.headerKey in extnHeader:
            val = extnHeader[self.headerKey]
        elif primaryHeader!=None and self.headerKey in primaryHeader:
            val = primaryHeader[self.headerKey]

        if val!=None and self.translation!=None:
            # Translate the header value through the regex
            val = self.translation[0].sub(self.translation[1],str(val))
            ### Note that if regex does not match, will pass unchanged.  Confusing with empty matches possible

        if val!=None:
            # Found something; return it as desired type
            return self.valueType(val)
        elif self.default!=None:
            # return the default if nothing was found
            return self.default
        else:
            # Found nothing useful in the header(s)
            return None

    
class Attribute:
    """
    Class is a list of AttributeFinders which all encode a given keyword.
    They are added in order of increasing priority.  We will try each one (in reverse order)
    to see if it can assign a value to the attribute, and take the first one that succeeds.
    """
    def __init__(self, first):
        self.seq = [first]
        self.key = first.key
        return

    def addFinder(self, next):
        if next.key != self.key:
            print("ERROR: Attribute with wrong name",first.key,"being added to",self.key)
            sys.exit(1)
        if next.valueType != self.seq[0].valueType:
            print("ERROR: Mismatched types for Attribute",self.key)
            sys.exit(1)
        self.seq.append(next)
        return
    
    def __call__(self,name, extnHeader=None, primaryHeader=None):
        for finder in reversed(self.seq):
            val = finder(name,extnHeader=extnHeader, primaryHeader=primaryHeader)
            if val != None:
                return val
        return None   # If no finder worked..

    def nodata(self):
        # Return a value of appropriate type indicating absence of data
        if self.seq[0].valueType==str:
            return 'nodata'
        else:
            return self.seq[0].valueType(-999)

def variableSubstitution(d):
    """
    Take the dictionary d, any values in it that are string-valued and contain
    substring of form ${XXX} will be replaced with the value at key XXX.
    Throws an exception for nonexistent key or for non-string-valued keys, or
    for self-referencing variable substitution.
    Does not attempt to find proper order of evaluation!
    But will iterate until no more substitutions are needed.

    variable names will be matched case-insensitive
    """
    variable = re.compile(r"^(.*)\$\{(.*)\}(.*)")

    # translate the dictionary to lower-case keys:
    dd = {k.lower():v for k,v in d.items()}
    maxIterations=4
    
    for i in range(maxIterations):
        anyChanges=False
        for k,v in dd.items():
            if not isinstance(v,str):
                # Only operate on string-valued entries
                continue
            m = variable.match(v)
            if not m:
                continue
            anyChanges = True
            vout = str(v)
            while m:
                key = m.group(2).lower()
                if key not in dd:
                    print("ERROR: variable substitution asks for nonexistent Attribute", key, "in", v)
                    sys.exit(1)
                if key==k:
                    print("ERROR: self-reference to Attribute", key, "in", v)
                vv = dd[key]
                if not isinstance(vv,str):
                    print("ERROR: variable substitution using non-string-valued Attribute",key)
                    sys.exit(1)
                vout = m.expand(r"\g<1>"+vv+r"\g<3>")
                m = variable.match(vout)
            dd[k] = vout
        if not anyChanges:
            break   # Done
        if i==maxIterations:
            print("ERROR: Too many iterations in variableSubstitution")
            sys.exit(1)
        # restore case of original dictionary
        for k in d:
            d[k] = dd[k.lower()]
    return

if __name__=='__main__':
    # Check for input & output file names in sys.argv (multiple inputs ok?)
    if len(sys.argv)<3:
        print("Usage: configure.py <input yaml> [more yaml...] <output fits>")
        sys.exit(1)
        
    # Read the YAML configuration file(s)
    fieldInput = []
    fileInput = []
    attributeInput = []
    for yamlFile in sys.argv[1:-1]:
        stream = open(yamlFile)
        y = yaml.load(stream)
        stream.close()
        if 'Files' in y:
            fileInput += y['Files']
        if 'Fields' in y:
            fieldInput += y['Fields']
        if 'Attributes' in y:
            attributeInput += y['Attributes']
        
    outFile = sys.argv[-1]

    # Collect Fields, check that all keywords present
    fields = {}
    for d in fieldInput:
        f = Field(**d)
        if f.name in fields:
            if not f==fields[f.name]:
                # Error to have two distinct fields with same name
                print('ERROR: 2 different fields with name', name)
                sys.exit(1)
        else:
            fields[f.name] = f

    print("Number of fields: ", len(fields))

    # Make an Attribute sequence for all unique column names
    attributes = {}
    for d in attributeInput:
        af = AttributeFinder(**d)
        if af.key in attributes:
            attributes[af.key].addFinder(af)
        else:
            attributes[af.key] = Attribute(af)

    # Dictionaries for exposures, instruments
    exposures = {}
    instruments = {}
    extensions = []

    filenameSignal = 'MAGIC___'   #Phony keyword used to store the filename in primary header
    rangere = re.compile(r"^(.*)\((\d+)-(\d+)\)(.*)$")  # regex to find numerical ranges in globs
    # Loop through file choices:
    for d in fileInput:
        # Expand any numerical ranges that are in the glob
        rm=rangere.match(d['glob'])
        if rm:
            start=int(rm.group(2))
            end=int(rm.group(3))
            if end<start:
                print("ERROR: bad numerical range in file request", d)
                sys.exit(1)
            globList = [rm.group(1) + str(x) + rm.group(4) for x in range(start,end+1)]
        else:
            globList = [d['glob']]

        files = []
        for g in globList:
            files += glob.glob(g)

        if len(files)==0:
            print("WARNING: No files found matching", d['glob'])
            
        tmpdict = {'key':'EXPOSURE','value':'< EXPNUM'}
        if 'expkey' in d:
            tmpdict['value'] = d['expkey']
        if 'translation' in d:
            tmpdict['translation'] = d['translation']
        if tmpdict['value'] == '_FILENAME':
            # Special signal means to use filename, not a header keyword, as EXPOSURE value
            # I'll kludge this by adding a special keyword to the header after I read it:
            tmpdict['value'] = headerIndicator + filenameSignal
        exposureAttrib = AttributeFinder(**tmpdict)

        for fitsname in files:
            print("Reading", fitsname)
            fits = pf.open(fitsname, memmap=True)
            # Primary extension can have a useful header but no table data
            pHeader = getHeader(fits,0,fitsname)
            pHeader[filenameSignal] = fitsname
            
            eHeader = None
            for iextn in range(1,len(fits)):
                tmphead = getHeader(fits, iextn, fitsname)
                if 'EXTNAME' in tmphead and tmphead['EXTNAME'] == 'LDAC_IMHEAD':
                    # For LDAC header extension, we just read header and move on
                    # eHeader = gbutil.readLDACHeader(fitsname, iextn)
                    eHeader = gbutil.headerFromString(fits[iextn].data[0][0])
                else:
                    if eHeader == None:
                        # This should be a catalog extension.  Get its header
                        eHeader = tmphead
                    else:
                        # This should be the objects part of an LDAC pair
                        if tmphead['EXTNAME'] != 'LDAC_OBJECTS':
                            print('ERROR: Did not get expected LDAC_OBJECTS at extn',iextn,\
                                'of file',fitsname)
                            sys.exit(1)
                    # Now have our headers, start filling in attributes in a dictionary
                    extn = {'FILENAME':fitsname,
                            'EXTENSION':iextn}

                    # Determine EXPOSURE  for this extension (required)
                    expoName = exposureAttrib(fitsname, primaryHeader=pHeader, extnHeader=eHeader)
                    if expoName==None:
                        print('ERROR: no exposure ID for file',fitsname,'extension',iextn)
                        sys.exit(1)
                    extn['EXPOSURE']=expoName

                    # Read every attribute for this exposure
                    for k,a in attributes.items():
                        if k=='APCORR':
                            # For the APCORR, give priority to the primary header because
                            # we prefer the median value of all CCDs to any individual one:
                            value = a(name=expoName, primaryHeader=eHeader, extnHeader=pHeader)
                        else:
                            # Normal priority is extension header first.
                            value = a(name=expoName, primaryHeader=pHeader, extnHeader=eHeader)
                        if value is not None:
                            extn[k] = value

                    # Do any variable substitutions on the attributes:
                    variableSubstitution(extn)

                    # Fill in blank device name if needed
                    if 'INSTRUMENT' not in extn:
                        print("ERROR: Missing Instrument at file",fitsname,"extension",iextn)
                        sys.exit(1)
                    if 'DEVICE' not in extn:
                        if extn['INSTRUMENT'] in specialInstruments:
                            # Give a blank device name
                            extn['DEVICE'] = ''
                        else:
                            print("ERROR: Missing Device at file",fitsname,"exposure",iextn)
                            sys.exit(1)

                    # Add new instrument if needed
                    inst = extn['INSTRUMENT']
                    if inst not in instruments and inst not in specialInstruments:
                        instruments[inst] = Instrument(inst, band=extn['BAND'])
                    # Check that extension BAND matches the instrument's band,
                    # or set the instrument BAND name if it has been None
                    if inst not in specialInstruments and extn['BAND'] is not None:
                        if instruments[inst].band is None:
                            instruments[inst].band = extn['BAND']
                        elif instruments[inst].band != extn['BAND']:
                            print("ERROR: Band",extn['BAND'],"does not match instrument",\
                                instruments[inst].name,"band",inst.band,\
                                "in catalog file",fitsname)
                            sys.exit(1)

                    # Read the WCS if needed
                    if extn['WCSIN'].strip()=='_HEADER':
                        # Retrieve the header from FITS WCS information in headers
                        wcshdr = gbutil.extractWCS( [pHeader, eHeader])
                        extn['WCSIN'] = wcshdr.tostring(sep='\n',
                                                        padding=False,endcard=True) + '\n'
                            
                    
                    if extn['EXPOSURE'] not in exposures:
                        # Create a new exposure dictionary if needed
                        exposures[extn['EXPOSURE']] = {}

                    extractExposureLevelAttributes(extn=extn, expo=exposures[extn['EXPOSURE']])

                    # We are done with this extension.  Clear out the extension header
                    eHeader = None
                    extensions.append(extn)
                pass  ## Ends the if/else on LDAC format catalogs

            fits.close() ## End the reading of this catalog's extensions

    # Completed loop through all input files from all globs

    # Assign all exposures to a field
    for expName,expo in exposures.items():
        if 'FIELD' not in expo:
            print("Exposure",expName,"does not have a FIELD attribute")
            sys.exit(1)
        if expo['FIELD']=="_NEAREST":
            # Finding nearest field center.
            expo['FIELD']=None
            minDistance = 0.
            for k,f in fields.items():
                if expo['FIELD']==None or \
                        f.distance(expo['coords'])<minDistance:
                    expo['FIELD'] = k
                    minDistance = f.distance(expo['coords'])
    
    # Get rid of fields with only 1 exposure:
    for k in list(fields):
        count = 0
        for e,v in exposures.items():
            if v['FIELD']==k:
                count += 1
                if count>1:
                    break
        if count<=1:
            # No matches will occur in this field.  Purge.
            fields.pop(k)
            print("WARNING: Ignoring field", k, "with <=1 exposure")

    # Now get rid of exposures that are not in a usable field, and
    # extensions that are from these exposures
    killedThese = []
    for e,v in list(exposures.items()):
        if v['FIELD'] not in fields:
            exposures.pop(e)
            killedThese.append(e)
            print("WARNING: Ignoring isolated exposure",e)
    tmplist = []
    for e in extensions:
        if e['EXPOSURE'] not in killedThese:
            tmplist.append(e)
    extensions = tmplist
    del killedThese
    del tmplist
    
    # Get rid of instruments that are not used in useful exposures
    for i in list(instruments):
        used = False
        for e,v in exposures.items():
            if v['INSTRUMENT']==i:
                used = True
                break
        if not used:
            instruments.pop(i)

    # Collect device names for all Instruments from remaining extensions
    # and replace DEVICE field in Extensions with the index number
    for extn in extensions:
        inst = exposures[extn['EXPOSURE']]['INSTRUMENT']
        if inst in specialInstruments:
            extn['DEVICE'] = -1
        else:
            instruments[inst].addDevice(extn['DEVICE'])
            extn['DEVICE'] = instruments[inst].devices[extn['DEVICE']]
    
    # ------Now ready to write the output files-------
    
    # Make the output FITS object
    outfits = pf.HDUList([pf.PrimaryHDU()])
    
    # Write FIELDs, assign numerical indices
    outfits.append(buildFieldTable(fields))
    
    # Write INSTRUMENTS, assign indices
    for index,inst in enumerate(sorted(instruments.keys())):
        instruments[inst].index = index
        outfits.append(buildInstrumentTable(instruments[inst]))
        
    # Write EXPOSURES, assign indices
    outfits.append(buildExposureTable(exposures, fields, instruments))
    
    # Now build the final (large) table of extensions.
    
    # Replace EXPOSURE name with integer index
    for e in extensions:
        e['EXPOSURE'] = exposures[e['EXPOSURE']]['index']

    cols = []

    # Add the WCSIN column as a variable-length character array
    data = []
    for e in extensions:
        data.append(e['WCSIN'])
#    cols.append(pf.Column(name='WCSIN',format='PA',array=data))
# Astropy Code seems to screw up variable-length arrays, at least in 0.2.4
    cols.append(pf.Column(name='WCSIN',format=py_to_fits(data),array=data))

    # Collect union of all attributes in all extensions
    allAttributes = set()
    for e in extensions:
        allAttributes.update(e.keys())
    allAttributes.discard('WCSIN')  # Already did this one above.
    
    # Now create a column for every attribute we have
    for a in allAttributes:
        print("*** Getting",a)
        data = []
        for e in extensions:
            if a in e:
                data.append(e[a])
            else:
                data.append(attributes[a].nodata())  # Enter a nodata value for the column
        cols.append(pf.Column(name=a, format=py_to_fits(data), array=data))

    tmp = list(allAttributes)
    tmp.sort()
    print(tmp) ###
    # Add the EXTENSION table as an HDU
    thdu = pf.BinTableHDU.from_columns(pf.ColDefs(cols),
                                      name='Extensions')
    outfits.append(thdu)
    
    # Write to FITS file
    outfits.writeto(outFile, overwrite=True)
    #Done!
    
