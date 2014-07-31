#!/usr/bin/env python
"""
Collect information needed for astrometric and photometric matching of a collection of catalogs.
Successor to Bob Armstrong's parse3.py that will take YAML input files, and also do all
collecting of information from the input files' headers, to remove some of the tedium from
the C++ programs that follow in the pipeline.
"""


"""
Fields:
- name:
  RA:
  Dec:
  Radius:
- name:...

Files:
- glob: glob for files
  expid: @keyword or _FILENAME or value for expid to give each catalog in the file
  translator: regex to apply to the expid
- glob: ...


Columns:
- name: the column name for the attribute
  type: string, float, int-valued column
  value: value or @keyword to look for in header
  select: regex that expid must match to apply this value
  default: value to assign if the specified header keyword does not exist
  translation: regex to apply to header values if found

Columns we'll require:
'FILENAME','EXTENSION','INSTRUMENT','DEVICE',
           'FIELD','EXPOSURE','RA','DEC','WCSFILE','XKEY',
           'YKEY','IDKEY','EXPID'

Special values:
EXPOSURE = _FILENAME
FIELD = _NEAREST
IDKEY = _ROW
WCS = _ICRS or _HEADER; special processing for this one.

"""

import sys
import yaml
import astropy.io.fits as pf
import astropy.coordinates as coords
import astropy.units as u
import re
import importlib

# These instrument names have special meaning.  "Observations" with these instruments
# do not require a device name.
specialInstruments=('REFERENCE','TAG')
# These are the extension attributes that are processed in some special way before
# being written to the output table.
specialAttributes =('Filename','Extension','Exposure','Instrument','Field',\
                    'Device','RA','Dec','WCSIN')

class Field:
    def __init__(self, initdict):
        # Expecting construction from a dictionary
        self.coords = coords.ICRS(ra=initdict['RA'],dec=initdict['Dec'],unit=(u.hourangle,u.degree))
        self.radius = float(initdict['Radius'])

        return

    def __eq__(self, other):
        return self.coords==other.coords and self.radius==other.radius

    def distance(self,location):
        return self.coords.separation(location).degreen

class Instrument:
    def __init__(self):
        # Start with an empty dictionary of device names
        self.devices = {}
        return
    def addDevice(self,devName):
        self.devices[devName] = len(self.devices)
        return

class Exposure:
    def __init__(self, coords, field, instrument, exptime=None, airmass=None):
        self.coords = coords
        self.field = field
        if (exptime==None):
            self.exptime = 1.
        else:
            self.exptime = exptime
        if (airmass==None):
            airmass = 1.
        else:
            self.airmass = airmass
        self.instrument = instrument
        return
    
class AttributeFinder:
    def __init__(self, initdict):
        self.filter = None
        if ('Filter') in initdict.keys():
            # regular expression that must match the input NAME for the Attribute to apply
            self.filter = re.compile(initdict['Filter'])

        # Determine the Python builtin type that the attribute must have, defaults to str
        self.valueType = str
        if 'Type' in initdict.keys():
            self.valueType = getattr(importlib.import_module('__builtin__'),initdict['Type'])

        # Either get the value of the attribute, or information for extracting it from headers
        val = initdict['Value']
        self.headerKey = None
        self.default = None
        self.translation = None
        if type(val)==str and len(val)>0 and val[0]=='@':
            # The given value is a string indicating a header keyword lookup
            self.headerKey = val.strip()[1:]
            if 'Default' in initdict.keys():
                # There is a default to assign if the keyword is absent
                self.default = self.valueType(initdict['Default'])
            if 'Translation' in initdict.keys():
                # The header value should be read as a string and run through a regex/replace
                rr = initdict['Translation'].split('=')
                if len(rr)!=2:
                    print "ERROR: Attribute translator is not of form <regex>=<replace>: ",\
                      initdict['Translation']
                    sys.exit(1)
                self.translation = (re.compile(rr[0]),rr[1])
        else:
            # There is a value given, cast it to desired type
            self.value = self.valueType(val)
        return

    def __call__(name, extnHeader=None, primaryHeader=None):
        # Return the value of the attribute if name passes the filter, otherwise return None.
        # If Attribute is to be sought in header, will do so from extnHeader first, then primaryHeader
        # if not found, then assign default if one exists.  Otherwise return None
        if self.filter!=None and self.filter.match(name)==None:
            # This name does not pass the filter
            return None
        if self.headerKey==None:
            # Simple fixed-value attribute
            return self.value

        val = None
        # Look in both headers for a value
        if extnHeader!=None and self.headerKey in extnHeader.keys():
            val = extnHeader[self.headerKey]
        elif primaryHeader!=None and self.headerKey in primaryHeader.keys():
            val = extnHeader[self.headerKey]

        if val!=None and self.translation!=None:
            # Translate the header value through the regex
            val = self.translation[0].sub(self.translation[1],str(val))
            ### Note that if regex does not match, will pass unchanged.  Also just one match replaced

        if val!=None:
            # Found something; return it as desired type
            return self.valueType(val)
        elif if self.default!=None:
            # return the default if nothing was found
            return self.default
        else:
            # Found nothing useful in the header(s)
            return None

    
class Attribute:
    """
    Class is a list of AttributeFinders which all encode a given keyword.
    They are added in order of increasing priority.  We will each one (in reverse order) to
    see if it can assign a value to the attribute, and take the first one that succeeds.
    """
    def __init__(self, first):
        self.seq = [first]
        return

    def addFinder(self, next):
        if next.valueType != self.seq[0].valueType:
            print "ERROR: Mismatched types for Attribute with value",next.value
            sys.exit(1)
        self.seq.append(next)
        return
    
    def __call__(name, extnHeader=None, primaryHeader=None):
        for finder in self.seq.reverse():
            val = finder(name,extnHeader=extnHeader, primaryHeader=primaryHeader)
            if val != None:
                return val
        return None   # If no finder worked..

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
        y = yaml.load(yamlFile)
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
        name = d['Name']
        f = Field(d)
        if name in fields.keys():
            if not f==fields[name]:
                # Error to have two distinct fields with same name
                print 'ERROR: 2 different fields with name', name
                sys.exit(1)
        else:
            fields[name] = f

    print "Number of fields: ", len(fields)

    # Make an Attribute sequence for all unique column names
    attributes = {}
    for d in attributeInput:
        name = d['Name']
        af = AttributeFinder(d)
        if name in attributes.keys():
            attributes[name].append(af)
        else:
            attributes[name] = Attribute(af)

    # Dictionaries for exposures, instruments
    exposures = {}
    instruments = {}
    extensions = []

    # Loop through file choices:
    for d in fileInput:
        # Maybe allow numerical ranges in the glob?
        files = glob.glob(d['glob'])
        tmpdict = {'Value':'@EXPNUM'}
        if 'Idkey' in d.keys():
            tmpdict['Value'] = d['Idkey']
        if 'Translation' in d.keys():
            tmpdict['Translation'] = d['Translation']
        if tmpdict['Value'] = '_FILENAME':
            # Special signal means to use filename, not a header keyword, as EXPOSURE value
            # I'll kludge this by adding a special keyword to the header after I read it:
            tmpdict['Value'] = '@MAGIC__'
        exposureAttrib = AttributeFinder(tmpdict)

        for f in files:
            fits = pf.open(f)
            # Primary extension can have a useful header but no table data
            pHeader = fits[0].header
            pHeader['@MAGIC__'] = f
            
            eHeader = None
            for iextn in range(1,len(fits)):
                if fits[iextn].header['EXTNAME'] == 'LDAC_HEADER':
                    # For LDAC header extension, we just read header and move on
                    eHeader = gmbpy.readLDACHeader(fits, iextn)
                else:
                    if eHeader == None:
                        # This should be a catalog extension.  Get its header
                        eHeader = fits[iextn].header
                    else:
                        # This should be the objects part of an LDAC pair
                        if fits[iextn].header['EXTNAME'] != 'LDAC_OBJECTS':
                            print 'ERROR: Did not get expected LDAC_OBJECTS at extn',iextn,\
                                'of file',fits
                            sys.exit(1)
                    # Now have our headers, start filling in attributes in a dictionary
                    extn = {'FILENAME':f,
                            'EXTENSION':iextn}


                    # Determine EXPOSURE  for this extension
                    expo = exposureAttrib(fits, primaryHeader=pHeader, extnHeader=eHeader)
                    if expo==None:
                        print 'ERROR: no exposure ID for file',fits', extension',iextn
                        sys.exit(1)
                    extn['Exposure']=expo

                    # Get INSTRUMENT
                    inst = attributes['Instrument'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    if inst==None:
                        print "ERROR: Missing Instrument at file",fits,"extension",iextn
                        sys.exit(1)
                    
                    # Other exposure-specific quantities:
                    ra = attributes['RA'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    dec = attributes['Dec'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    # ??? Check for RA already in degrees?
                    if ra!=None and dec!=None:
                        icrs = coords.ICRS(ra,dec,unit=(u.hourangle,u.degree))
                    else:
                        icrs = None
                    airmass = attributes['Airmass'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    exptime = attributes['Exptime'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    field = attributes['Field'](expo, primaryHeader=pHeader, extnHeader=eHeader)

                    if expo in exposures:
                        # Already have an exposure with this name, make sure this is basically the same
                        e = exposures[expo]
                        if (airmass!=None and abs(airmass-e.airmass)>0.002) or \
                           (exptime!=None and abs(exptime/e.exptime-1.)>0.0002) or \
                           (icrs!=None and icrs.separation(e.coords).degree > 1):
                            print "ERROR: info mismatch at exposure",expo, "file",fits,'extension',extn
                            sys.exit(1)
                        # Also check the field for a match, unless it is _NEAREST
                        if field != '_NEAREST' and field!=e.field:
                            print "ERROR: Field mismatch for exposure",expo
                            print "Exposure has",e.field,"file",fits,"extension",iextn,"has",field
                            sys.exit(1)
                        # Check the instrument for a match
                        if inst!=e.instrument:
                            print "ERROR: Instrument mismatch for exposure",expo
                            print "Exposure has",e.instrument,"file",fits,"extension",iextn,"has",inst
                            sys.exit(1)
                    else:
                        # New exposure.  Need coordinates to create it
                        if icrs==None:
                            print "ERROR: Missing RA/Dec for new exposure",expo \
                              "at file",fits,"extension",iextn
                            sys.exit(1)
                        if field==None:
                            print "ERROR: Missing Field for new exposure",expo \
                              "at file",fits,"extension",iextn
                            sys.exit(1)
                        elif field=="_NEAREST":
                            # Finding nearest field center.
                            field=None
                            minDistance = 0.
                            for k,f in fields.items():
                                if field==None or f.distance(icrs)<minDistance:
                                    field = k
                                    minDistance = f.distance(icrs)
                        exposures[expo] = Exposure(c, field, inst, airmass=airmass, exptime=exptime)

                    # Add new instrument if needed
                    if inst not in instruments and inst not in specialInstruments:
                        instruments[inst] = Instrument()
                        
                    # Get device name, not needed for special devices
                    dev  = attributes['Device'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    if dev==None:
                        if inst in specialInstruments:
                            # Give a blank device name
                            dev = ''
                        else:
                            print "ERROR: Missing Device at file",fits,"exposure",iextn
                            sys.exit(1)

                    # Record the device for this extension.
                    extn['Device'] = dev
                        
                    # ??? Handle the WCSIN

                    # Now get a value for all other attributes and add to extension's dictionary
                    for k,a in attributes.items():
                        if k not in specialAttributes:
                            extn[k] = a(expo, primaryHeader=pHeader, extnHeader=eHeader)

                    # We are done with this extension.  Clear out the extension header
                    extnHeader = None
                    extensions.append(extn)

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

    # Now get rid of exposures that are not in a usable field, and
    # extensions that are from these exposures
    killedThese = []
    for e,v in exposures.items():
        if v.field not in fields:
            exposures.pop(e)
            killedThese.append(e)
    tmplist = []
    for e in extensions:
        if e['Exposure'] not in killedThese:
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
        inst = exposures[extn['Exposure']].instrument
        if inst not in specialInstruments:
            instruments[inst].addDevice(extn['Device'])
    
    # Write FIELDs, assign numerical indices
    # Write INSTRUMENTS, assign indices
    # Write EXPOSURES, assign indices
    # Write Extensions

