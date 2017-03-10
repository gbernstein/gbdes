#!/usr/bin/env python
"""
Collect information needed for astrometric and photometric matching of a collection of catalogs.
Successor to Bob Armstrong's parse3.py that will take YAML input files, and also do all
collecting of information from the input files' headers, to remove some of the tedium from
the C++ programs that follow in the pipeline.
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
import gmbpy.utilities

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
specialInstruments={'REFERENCE':-1,'TAG':-2}

# These are the extension attributes that are processed in some special way before
# being written to the output table, or belong to Exposure rather than Extension
# and hence should not be written to Extension table.
specialAttributes =('FILENAME','EXTENSION','EXPOSURE','INSTRUMENT','FIELD',\
                    'DEVICE','RA','DEC','AIRMASS','EXPTIME','MJD','WCSIN','BAND',\
                    'EPOCH','APCORR')

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
        print "ERROR: No known FITS type for python type",vtype
        sys.exit(1)
    

def colForAttribute(extensions, name):
    """
    Make a FITS BinTable column out of the value under key=name in every extension
    """
    data = []
    for e in extensions:
        data.append(e[name])
    return pf.Column(name=name, format=py_to_fits(data), array=data)

class Field:
    def __init__(self, name, coords, radius, index=-1):
        # Expecting construction from a dictionary
        self.name = name
        if oldAstropy:
            self.coords = ICRS(coords,unit=(u.hourangle,u.degree))
        else:
            self.coords = cc.SkyCoord(coords,unit=(u.hourangle,u.degree),frame='icrs')
        self.radius = float(radius)
        self.index = int(index)
        return

    def __eq__(self, other):
        return self.coords==other.coords and self.radius==other.radius

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
    
    index = 0
    for k,v in fields.items():
        name.append(k)
        ra.append(getDegree(v.coords.ra))
        dec.append(getDegree(v.coords.dec))
        radius.append(v.radius)
        v.index = index
        index += 1
        
    hdu = pf.BinTableHDU.from_columns(\
        pf.ColDefs( [pf.Column(name='NAME',format=py_to_fits(name),array=name),
                     pf.Column(name='RA',format=py_to_fits(ra),array=ra),
                     pf.Column(name='DEC',format=py_to_fits(dec),array=dec),
                     pf.Column(name='RADIUS',format=py_to_fits(radius),array=radius)]),
                     name = 'Fields')
#    hdu.header['EXTNAME'] = 'Fields'
    return hdu

class Instrument:
    def __init__(self, name, index=-3, band=None):
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
            print 'ERROR: Device indices screwed up for instrument',inst.name
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

class Exposure:
    def __init__(self, name, coords, field, instrument, exptime=None, \
                 airmass=None, mjd=None, epoch=None, \
                 apcorr=None, index=-1, **kwargs):
        # Note that the constructor may contain arguments not used.
        self.name = name
        self.coords = coords
        self.field = field
        if exptime is None:
            self.exptime = 1.
        else:
            self.exptime = exptime
        if airmass is None:
            self.airmass = 1.
        else:
            self.airmass = airmass
        if mjd is None:
            self.mjd = 0.
        else:
            self.mjd = mjd
        if epoch is None:
            self.epoch = ""
        else:
            self.epoch = epoch
        if apcorr is None:
            self.apcorr = 0.
        else:
            self.apcorr = apcorr
        self.instrument = instrument
        self.index = index
        return
    
def buildExposureTable(exposures, fields, instruments):
    """
    Return a FITS BinTableHDU with information on exposures.
    Also assign indices to each according to their row in the table
    """
    name = []
    ra = []
    dec= []
    field= []
    inst = []
    airmass = []
    mjd = []
    exptime = []
    epoch = []
    apcorr = []
    index = 0
    for k,e in exposures.items():
        name.append(e.name)
        ra.append(getDegree(e.coords.ra))
        dec.append(getDegree(e.coords.dec))
        field.append(fields[e.field].index)
        if e.instrument in specialInstruments:
            inst.append(specialInstruments[e.instrument])
        else:
            inst.append(instruments[e.instrument].index)
        e.index = index
        index += 1

        airmass.append(e.airmass)
        mjd.append(e.mjd)
        exptime.append(e.exptime)
        epoch.append(e.epoch)
        apcorr.append(e.apcorr)
    hdu = pf.BinTableHDU.from_columns(\
        pf.ColDefs( [pf.Column(name='NAME',format=py_to_fits(name),array=name),
                     pf.Column(name='RA',format=py_to_fits(ra),array=ra),
                     pf.Column(name='DEC',format=py_to_fits(dec),array=dec),
                     pf.Column(name='FIELDNUMBER',format=py_to_fits(field),array=field),
                     pf.Column(name='INSTRUMENTNUMBER',format=py_to_fits(inst),\
                               array=inst),
                     pf.Column(name="MJD",format=py_to_fits(mjd),array=mjd),
                     pf.Column(name="AIRMASS",format=py_to_fits(airmass),array=airmass),
                     pf.Column(name="EXPTIME",format=py_to_fits(exptime),array=exptime),
                     pf.Column(name="EPOCH",format=py_to_fits(epoch),array=epoch) ]),
                     pf.Column(name="APCORR",format=py_to_fits(apcorr),array=apcorr) ]),
                     name = 'Exposures')
    # hdu.header['EXTNAME'] = 'Exposures'
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
        self.valueType = getattr(importlib.import_module('__builtin__'),vtype)

        # Either get the value of the attribute, or information for extracting it from headers
        self.headerKey = None
        if type(value)==str and len(value)>0 and value.strip()[0]==headerIndicator:
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
                    print "ERROR: Attribute translator is not of form <regex>=<replace>: ",\
                      initdict['Translation']
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
        if extnHeader!=None and self.headerKey in extnHeader.keys():
            val = extnHeader[self.headerKey]
        elif primaryHeader!=None and self.headerKey in primaryHeader.keys():
            val = primaryHeader[self.headerKey]

        if val!=None and self.translation!=None:
            # Translate the header value through the regex
            val = self.translation[0].sub(self.translation[1],str(val))
            ### Note that if regex does not match, will pass unchanged.  Also just one match replaced

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
            print "ERROR: Attribute with wrong name",first.key,"being added to",self.key
            sys.exit(1)
        if next.valueType != self.seq[0].valueType:
            print "ERROR: Mismatched types for Attribute",self.key
            sys.exit(1)
        self.seq.append(next)
        return
    
    def __call__(self,name, extnHeader=None, primaryHeader=None):
        for finder in reversed(self.seq):
            val = finder(name,extnHeader=extnHeader, primaryHeader=primaryHeader)
            if val != None:
                return val
        return None   # If no finder worked..

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
    dd = {k.lower():v for k,v in d.iteritems()}
    maxIterations=4
    
    for i in range(maxIterations):
        anyChanges=False
        for k,v in dd.iteritems():
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
                if key not in dd.keys():
                    print "ERROR: variable substitution asks for nonexistent Attribute", key, "in", v
                    sys.exit(1)
                if key==k:
                    print "ERROR: self-reference to Attribute", key, "in", v
                vv = dd[key]
                if not isinstance(vv,str):
                    print "ERROR: variable substitution using non-string-valued Attribute",key
                    sys.exit(1)
                vout = m.expand(r"\g<1>"+vv+r"\g<3>")
                m = variable.match(vout)
            dd[k] = vout
        if not anyChanges:
            break   # Done
        if i==maxIterations:
            print "ERROR: Too many iterations in variableSubstitution"
            sys.exit(1)
        # restore case of original dictionary
        for k in d.keys():
            d[k] = dd[k.lower()]
    return

if __name__=='__main__':
    # Check for input & output file names in sys.argv (multiple inputs ok?)
    if len(sys.argv)<3:
        print "Usage: configure.py <input yaml> [more yaml...] <output fits>"
        sys.exit(1)
        
    # Read the YAML configuration file(s)
    fieldInput = []
    fileInput = []
    attributeInput = []
    for yamlFile in sys.argv[1:-1]:
        stream = open(yamlFile)
        y = yaml.load(stream)
        stream.close()
        if 'Files' in y.keys():
            fileInput += y['Files']
        if 'Fields' in y.keys():
            fieldInput += y['Fields']
        if 'Attributes' in y.keys():
            attributeInput += y['Attributes']
        
    outFile = sys.argv[-1]

    # Collect Fields, check that all keywords present
    fields = {}
    for d in fieldInput:
        f = Field(**d)
        if f.name in fields.keys():
            if not f==fields[f.name]:
                # Error to have two distinct fields with same name
                print 'ERROR: 2 different fields with name', name
                sys.exit(1)
        else:
            fields[f.name] = f

    print "Number of fields: ", len(fields)

    # Make an Attribute sequence for all unique column names
    attributes = {}
    for d in attributeInput:
        af = AttributeFinder(**d)
        if af.key in attributes.keys():
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
                print "ERROR: bad numerical range in file request", d
                sys.exit(1)
            globList = [rm.group(1) + str(x) + rm.group(4) for x in range(start,end+1)]
        else:
            globList = [d['glob']]

        files = []
        for g in globList:
            files += glob.glob(g)

        if len(files)==0:
            print "WARNING: No files found matching", d['glob']
            
        tmpdict = {'key':'EXPOSURE','value':'< EXPNUM'}
        if 'expkey' in d.keys():
            tmpdict['value'] = d['expkey']
        if 'translation' in d.keys():
            tmpdict['translation'] = d['translation']
        if tmpdict['value'] == '_FILENAME':
            # Special signal means to use filename, not a header keyword, as EXPOSURE value
            # I'll kludge this by adding a special keyword to the header after I read it:
            tmpdict['value'] = headerIndicator + filenameSignal
        exposureAttrib = AttributeFinder(**tmpdict)

        for fitsname in files:
            print "Reading", fitsname
            fits = pf.open(fitsname, memmap=True)
            # Primary extension can have a useful header but no table data
            pHeader = getHeader(fits,0,fitsname)
            pHeader[filenameSignal] = fitsname
            
            eHeader = None
            for iextn in range(1,len(fits)):
                tmphead = getHeader(fits, iextn, fitsname)
                if 'EXTNAME' in tmphead and tmphead['EXTNAME'] == 'LDAC_IMHEAD':
                    # For LDAC header extension, we just read header and move on
                    # eHeader = gmbpy.utilities.readLDACHeader(fitsname, iextn)
                    eHeader = gmbpy.utilities.headerFromString(fits[iextn].data[0][0])
                else:
                    if eHeader == None:
                        # This should be a catalog extension.  Get its header
                        eHeader = tmphead
                    else:
                        # This should be the objects part of an LDAC pair
                        if tmphead['EXTNAME'] != 'LDAC_OBJECTS':
                            print 'ERROR: Did not get expected LDAC_OBJECTS at extn',iextn,\
                                'of file',fitsname
                            sys.exit(1)
                    # Now have our headers, start filling in attributes in a dictionary
                    extn = {'FILENAME':fitsname,
                            'EXTENSION':iextn}


                    # Determine EXPOSURE  for this extension
                    expo = exposureAttrib(fitsname, primaryHeader=pHeader, extnHeader=eHeader)
                    if expo==None:
                        print 'ERROR: no exposure ID for file',fitsname,'extension',iextn
                        sys.exit(1)
                    extn['EXPOSURE']=expo

                    # Arguments we will use for all of the Attribute calls:
                    attargs = {'name':expo,
                               'primaryHeader':pHeader,
                               'extnHeader':eHeader}

                    # Collect exposure-specific information into a dictionary
                    expoAttr = {}
                    
                    # Get INSTRUMENT
                    expoAttr['instrument'] = attributes['INSTRUMENT'](**attargs)
                    if expoAttr['instrument']==None:
                        print "ERROR: Missing Instrument at file",fitsname,"extension",iextn
                        sys.exit(1)
                    
                    # Other exposure-specific quantities:
                    ra = attributes['RA'](**attargs)
                    dec = attributes['DEC'](**attargs)
                    # ??? Check for RA already in degrees?
                    if ra!=None and dec!=None:
                        try:
                            # If ra is valid floating point, assume it's degrees
                            float(ra)
                            expoAttr['coords'] = cc.SkyCoord(ra=ra,dec=dec,unit=(u.degree,u.degree),
                                                             frame='icrs')
                        except ValueError:
                            # Otherwise it must be HH:MM:SS.sss style in hours
                            expoAttr['coords'] = cc.SkyCoord(ra=ra,dec=dec,unit=(u.hourangle,u.degree),
                                                             frame='icrs')
                    else:
                        expoAttr['coords'] = None
                    expoAttr['airmass'] = attributes['AIRMASS'](**attargs)
                    expoAttr['band'] = attributes['BAND'](**attargs)
                    expoAttr['exptime'] = attributes['EXPTIME'](**attargs)
                    expoAttr['field'] = attributes['FIELD'](**attargs)
                    expoAttr['mjd'] = attributes['MJD'](**attargs)
                    expoAttr['epoch'] = attributes['EPOCH'](**attargs)
                    expoAttr['apcorr'] = attributes['APCORR'](**attargs)

                    # Apply variable substitution to the exposure attributes
                    variableSubstitution(expoAttr)
                    # And pass some info along to the extension too
                    inst = expoAttr['instrument']
                    extn['INSTRUMENT'] = inst
                    extn['BAND'] = expoAttr['band']
                    
                    if expo in exposures:
                        # Already have an exposure with this name, make sure this is basically the same
                        e = exposures[expo]
                        if (expoAttr['airmass']!=None and abs(expoAttr['airmass']-e.airmass)>0.002) or \
                           (expoAttr['exptime']!=None and \
                            abs(expoAttr['exptime']/e.exptime-1.)>0.0002) or \
                           (expoAttr['mjd']!=None and abs(expoAttr['mjd'] - e.mjd)>0.0002) or \
                           (expoAttr['apcorr']!=None and abs(expoAttr['apcorr'] - e.apcorr)>0.0002) or \
                           (expoAttr['coords']!=None
                            and getDegree(expoAttr['coords'].separation(e.coords)) > 1.):
                            print "ERROR: info mismatch at exposure",expo, "file",fitsname, \
                                  'extension',iextn
                            sys.exit(1)
                        # Also check the field for a match, unless it is _NEAREST
                        if expoAttr['field'] != '_NEAREST' and expoAttr['field']!=e.field:
                            print "ERROR: FIELD mismatch for exposure",expo
                            print "Exposure has",e.field,"file",fitsname,"extension",iextn, \
                                  "has",expoAttr['field']
                            sys.exit(1)
                        # Check the instrument for a match
                        if inst!=e.instrument:
                            print "ERROR: INSTRUMENT mismatch for exposure",expo
                            print "Exposure has",e.instrument,"file",fitsname,"extension",iextn,"has",inst
                            sys.exit(1)
                    else:
                        # New exposure.  Need coordinates to create it
                        if expoAttr['coords']==None:
                            print "ERROR: Missing RA/DEC for new exposure",expo, \
                              "at file",fitsname,"extension",iextn
                            sys.exit(1)
                        if expoAttr['field']==None:
                            print "ERROR: Missing FIELD for new exposure",expo, \
                              "at file",fitsname,"extension",iextn
                            sys.exit(1)
                        elif expoAttr['field']=="_NEAREST":
                            # Finding nearest field center.
                            expoAttr['field']=None
                            minDistance = 0.
                            for k,f in fields.items():
                                if expoAttr['field']==None or \
                                  f.distance(expoAttr['coords'])<minDistance:
                                    expoAttr['field'] = k
                                    minDistance = f.distance(expoAttr['coords'])
                        exposures[expo] = Exposure(expo, **expoAttr)

                    # Add new instrument if needed
                    if inst not in instruments and inst not in specialInstruments:
                        instruments[inst] = Instrument(inst, band=expoAttr['band'])
                    # Check that extension BAND matches the instrument's band,
                    # or set the instrument BAND name if it has been None
                    if inst not in specialInstruments and expoAttr['band'] is not None:
                        if instruments[inst].band is None:
                            instruments[inst].band = expoAttr['band']
                        elif instruments[inst].band != expoAttr['band']:
                            print "ERROR: Band",expoAttr['band'],"does not match instrument",\
                              instruments[inst].name,"band",inst.band,\
                              "in catalog file",fitsname
                            sys.exit(1)
                        
                    # Get device name, not needed for special devices
                    dev  = attributes['DEVICE'](**attargs)
                    if dev==None:
                        if inst in specialInstruments:
                            # Give a blank device name
                            dev = ''
                        else:
                            print "ERROR: Missing Device at file",fitsname,"exposure",iextn
                            sys.exit(1)

                    # Record the device for this extension.
                    extn['DEVICE'] = dev
                        
                    # Handle the WCSIN
                    wcsin = attributes['WCSIN'](**attargs)
                    if wcsin.strip()=='_HEADER':
                        # Retrieve the header from FITS WCS information in headers
                        wcshdr = gmbpy.utilities.extractWCS( [pHeader, eHeader])
                        extn['WCSIN'] = wcshdr.tostring(sep='\n',
                                                        padding=False,endcard=True) + '\n'
                    else:
                        extn['WCSIN'] = wcsin
                        
                    # Now get a value for all other attributes and add to extension's dictionary
                    for k,a in attributes.items():
                        if k not in specialAttributes:
                            extn[k] = a(**attargs)

                    # Do any variable substitutions on the attributes:
                    variableSubstitution(extn)
                    
                    # We are done with this extension.  Clear out the extension header
                    eHeader = None
                    extensions.append(extn)
                pass  ## Ends the if/else on LDAC format catalogs

            fits.close() ## End the reading of this catalog's extensions

    # Completed loop through all input files from all globs
    # Get rid of fields with only 1 exposure:
    for k in fields.keys():
        count = 0
        for e,v in exposures.items():
            if v.field==k:
                count += 1
                if count>1:
                    break
        if count<=1:
            # No matches will occur in this field.  Purge.
            fields.pop(k)
            print "WARNING: Ignoring field", k, "with <=1 exposure"

    # Now get rid of exposures that are not in a usable field, and
    # extensions that are from these exposures
    killedThese = []
    for e,v in exposures.items():
        if v.field not in fields:
            exposures.pop(e)
            killedThese.append(e)
            print "WARNING: Ignoring isolated exposure",e
    tmplist = []
    for e in extensions:
        if e['EXPOSURE'] not in killedThese:
            tmplist.append(e)
    extensions = tmplist
    del killedThese
    del tmplist
    
    # Get rid of instruments that are not used in useful exposures
    for i in instruments.keys():
        used = False
        for e,v in exposures.items():
            if v.instrument==i:
                used = True
                break
        if not used:
            instruments.pop(i)

    # Collect device names for all Instruments from remaining extensions
    for extn in extensions:
        inst = exposures[extn['EXPOSURE']].instrument
        if inst not in specialInstruments:
            instruments[inst].addDevice(extn['DEVICE'])
    
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
    
    # Now build the final (large) table of extensions:
    cols = []
    cols.append(colForAttribute(extensions, 'FILENAME'))
    cols.append(colForAttribute(extensions, 'EXTENSION'))

    # For EXPOSURE and DEVICE: get the integer indices instead of names
    data = []
    for e in extensions:
        data.append(exposures[e['EXPOSURE']].index)
    cols.append(pf.Column(name='EXPOSURE',format=py_to_fits(data),array=data))

    # For DEVICE: get integer indices, aware that specialInstruments have none
    data = []
    for e in extensions:
        inst = exposures[e['EXPOSURE']].instrument
        if inst in specialInstruments:
            dev = -1
        else:
            dev = instruments[inst].devices[e['DEVICE']]
        data.append(dev)
    cols.append(pf.Column(name='DEVICE',format=py_to_fits(data),array=data))

    # Add the WCSIN column as a variable-length character array
    data = []
    for e in extensions:
        data.append(e['WCSIN'])
#    cols.append(pf.Column(name='WCSIN',format='PA',array=data))
# Astropy Code seems to screw up variable-length arrays, at least in 0.2.4
    cols.append(pf.Column(name='WCSIN',format=py_to_fits(data),array=data))

    # Now create a column for every other attribute we have
    for k,a in attributes.items():
        if a.key not in specialAttributes:
            print "*** Getting",a.key
            cols.append(colForAttribute(extensions, a.key))

    # Add the EXTENSION table as an HDU
    thdu = pf.BinTableHDU.from_columns(pf.ColDefs(cols),
                                      name='Extensions')
    #thdu.header['EXTNAME'] = 'Extensions'
    outfits.append(thdu)
    
    # Write to FITS file
    outfits.writeto(outFile, clobber=True)
    #Done!
    
