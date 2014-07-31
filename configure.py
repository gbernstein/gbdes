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
  translator: regex to apply to header values if found

Columns we'll require:
'FILENAME','EXTENSION','INSTRUMENT','DEVICE',
           'FIELD','EXPOSURE','RA','DEC','WCSFILE','XKEY',
           'YKEY','IDKEY','EXPID'

Special values:
FIELD = _NEAREST
IDKEY = _ROW
WCS = _ICRS or _HEADER; special processing for this one.

Special processing:
FIELD - replace name with number
INSTRUMENT - number?
DEVICE - number?

"""

import sys
import yaml
import astropy.io.fits as pf
import astropy.coordinates as coords
import astropy.units as u
import re
import importlib

class Field:
    def __init__(self, **kwargs):
        # Expecting construction from a dictionary
        self.coords = coords.ICRS(ra=initdict['RA'],dec=initdict['Dec'],unit=(u.hourangle,u.degree))
        self.radius = float(initdict['Radius'])
        return

    def __eq__(self, other):
        return self.coords==other.coords and self.radius==other.radius

    def distance(self,location):
        return self.coords.separation(location).degreen

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
                    print "Attribute translator is not of form <regex>=<replace>: ",initdict['Translation']
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

    # Make a sequence of AttributeFinders for all unique column names
    attributes = {}
    for d in attributeInput:
        name = d['Name']
        af = AttributeFinder(d)
        if name in attributes.keys():
            if af.valueType != attributes[name][0].valueType:
                print "Mismatched types for Attributes", name
                sys.exit(1)
            attributes[name].append(af)
        else:
            attributes[name] = [af,]

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
            primaryHeader = fits[0].header
            primaryHeader['@MAGIC__'] = f
            
            extnHeader = None
            iextn = 1
            while iextn < len(fits):
                if fits[i].header['EXTNAME'] == 'LDAC_HEADER':
                    extnHeader = ...
                    
#  Find files with glob -> FILENAME
#  Get primary header
#  Loop through extensions
#    Check for LDAC_HEADER, if so read header and advance EXTENSION counter
#    get secondary header
#    Get EXPID and translate it
#    Loop through columns
#      Find first value entry that applies to this EXPID and yields a result.
#      Translate value if requested
#        Error if no rule can fill this column
#      Special case: INSTRUMENT.  Make new Instrument if none exist
#      Special case: DEVICE.  Add to Instrument if new.
#      Speical case: EXPOSURE.  Add new Exposure, check for consistent AIRMASS?
#      Special case: RA, Dec: turn into degrees?
#      Special case: FIELD - do _NEAREST if needed
#      Special case: WCS: might extract TPV from header(s)
# Write FIELDs that are in use with >1 exposure
# Write EXPOSURES in useful fields
# Write INSTRUMENTS
# Write Extensions

