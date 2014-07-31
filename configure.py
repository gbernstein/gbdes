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
In the yaml file, we'll do First Letter Cap for the big categories, but expect lowercase for all
map keys otherwise.

Fields:
- name:
  ra:
  dec:
  radius:
- name:...

Files:
- glob: glob for files
  expkey: @keyword or _FILENAME or value for expid to give each catalog in the file
  translation: regex to apply to the id that is derived
- glob: ...


Attributes:
- name: the column name for the attribute
  select: regex that exposure name must match to apply this value
  vtype: str, float, int - data type that the we should get for value
  value: value or '@<keyword>' to look for in header
  default: value to assign if the specified header keyword does not exist
  translation: regex to apply to header values if found

Columns we'll require:
'FILENAME','EXTENSION','INSTRUMENT','DEVICE',
           'FIELD','EXPOSURE','RA','DEC','WCSFILE','XKEY',
           'YKEY','IDKEY'

Special values:
EXPOSURE = _FILENAME
FIELD = _NEAREST
IDKEY = _ROW
WCS = _ICRS or _HEADER; special processing for this one.

translations have format <regex>=<replace>
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
# They have special indices: **** Keep in sync with FitSubroutines.h values! ****
specialInstruments={'REFERENCE':-1,'TAG':-2}
# These are the extension attributes that are processed in some special way before
# being written to the output table.
specialAttributes =('FILENAME','EXTENSION','EXPOSURE','INSTRUMENT','FIELD',\
                    'DEVICE','RA','DEC','WCSIN')

def py_to_fits(val):
    """
    Return a fits type-specification string corresponding to type of first
    element in the python array.
    Use the type of the first entry in the list
    """

    vtype=type(val[0])
    
    if(vtype is IntType):
        return 'J'
    elif(vtype is FloatType):
        return 'D'
    elif(vtype is StringType):
        size=len(max(val,key=len))
        return 'A'+str(size)

def colForAttribute(extensions, name):
    """
    Make a FITS BinTable column out of the value under key=name in every extension
    """
    data = []
    for e in extensions:
        data.append(e[name])
    return pf.Column(name=name, format=py_to_fits(data), data)

class Field:
    def __init__(self, name, ra, dec, radius, index=-1):
        # Expecting construction from a dictionary
        self.name = name
        self.coords = coords.ICRS(ra,dec,unit=(u.hourangle,u.degree))
        self.radius = float(radius)
        self.index = int(index)
        return

    def __eq__(self, other):
        return self.coords==other.coords and self.radius==other.radius

    def distance(self,location):
        return self.coords.separation(location).degree

def buildFieldTable(fields):
    """
    Make a FITS BinTableHDU containing information on all the fields,
    and assign them indices as we go
    """
    name=[]
    ra  =[]
    dec =[]
    radius = []
    
    int index = 0
    for k,v in fields.items():
        name.append(k)
        ra.append(v.coords.ra.degree)
        dec.append(v.coords.dec.degree)
        radius.append(v.radius)
        v.index = index
        index += 1
        
    return pf.new_table(pf.ColDefs( [pf.Column(name='NAME',format=py_to_fits(name),name=name),
                                     pf.Column(name='RA',format=py_to_fits(ra),name=ra),
                                     pf.Column(name='DEC',format=py_to_fits(dec),name=dec),
                                     pf.Column(name='RADIUS',format=py_to_fits(radius),name=radius)]))

class Instrument:
    def __init__(self, name, index=-3):
        # Start with an empty dictionary of device names
        self.name = name
        self.devices = {}
        self.index = index
        return
    def addDevice(self,devName):
        self.devices[devName] = len(self.devices)
        return

def buildInstrumentTable(inst):
    """
    Return a FITS BinTableHDU for this instrument
    """
    names = ["" for k in inst.devices]
    for k,v in inst.devices:
        names[v] = k
    for n in names:
        # If any of the device names is null, then the indices are screwed up
        if n="":
            print 'ERROR: Device indices screwed up for instrument',inst.name
            sys.exit(1)
    # Make an array of zeros to build the [xy]Min/Max columns
    z = np.zeros(len(names),dtype=float)
    hdu = pf.new_table(pf.ColDefs( [pf.Column(name='NAME',format=py_to_fits(names),data=names),
                                    pf.Column(name='XMIN',format=py_to_fits(z),data=z),
                                    pf.Column(name='XMAX',format=py_to_fits(z),data=z),
                                    pf.Column(name='XMIN',format=py_to_fits(z),data=z),
                                    pf.Column(name='XMAX',format=py_to_fits(z),data=z)]));
    hdu.header['NAME'] = inst.name
    hdu.header['NUMBER'] = inst.index
    return hdu

class Exposure:
    def __init__(self, name, coords, field, instrument, exptime=None, airmass=None, index=-1):
        self.name = name
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
    index = 0
    for e in exposures:
        name.append(e.name)
        ra.append(e.coords.ra.degree)
        dec.append(e.coords.dec.degree)
        field.append(fields[e.field].index)
        inst.append(instruments[e.instrument].index)
        e.index = index
        index += 1

        return pf.new_table(pf.ColDefs( [pf.Column(name='NAME',format=py_to_fits(name),name=name),
                                        pf.Column(name='RA',format=py_to_fits(ra),name=ra),
                                        pf.Column(name='DEC',format=py_to_fits(dec),name=dec),
                                        pf.Column(name='FIELDNUMBER',format=py_to_fits(field),name=field),
                                        pf.Column(name='INSTRUMENTNUMBER',format=py_to_fits(inst),name=inst)\
                                          ]))

class AttributeFinder:
    def __init__(self, key, value, vtype="str", default=None, translation=None, select=None):
        self.key = key
        if select != None
            # regular expression that must match the input NAME for the Attribute to apply
            self.select = re.compile(select)

        # Determine the Python builtin type that the attribute must have, defaults to str
        self.valueType = getattr(importlib.import_module('__builtin__'),vtype)

        # Either get the value of the attribute, or information for extracting it from headers
        self.headerKey = None
        if type(value)==str and len(value)>0 and value.strip()[0]=='@':
            # The given value is a string indicating a header keyword lookup
            self.headerKey = value.strip()[1:]
            if default == None:
                self.default = None
            else:
                # There is a default to assign if the keyword is absent
                self.default = self.valueType(default)
            if translation==None:
                self.translation==None
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

    def __call__(name, extnHeader=None, primaryHeader=None):
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
        self.name = first.name
        return

    def addFinder(self, next):
        if next.name != self.name:
            print "ERROR: Attribute with wrong name",first.name,"being added to",self.name
            sys.exit(1)
        if next.valueType != self.seq[0].valueType:
            print "ERROR: Mismatched types for Attribute",self.name
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
        if af.name in attributes.keys():
            attributes[af.name].append(af)
        else:
            attributes[af.name] = Attribute(af)

    # Dictionaries for exposures, instruments
    exposures = {}
    instruments = {}
    extensions = []

    filenameSignal = 'MAGIC___'   #Phony keyword used to store the filename in primary header
    # Loop through file choices:
    for d in fileInput:
        # Maybe allow numerical ranges in the glob?
        files = glob.glob(d['glob'])
        tmpdict = {'name':'EXPOSURE','value':'@EXPNUM'}
        if 'namekey' in d.keys():
            tmpdict['value'] = d['expkey']
        if 'Translation' in d.keys():
            tmpdict['translation'] = d['translation']
        if tmpdict['value'] = '_FILENAME':
            # Special signal means to use filename, not a header keyword, as EXPOSURE value
            # I'll kludge this by adding a special keyword to the header after I read it:
            tmpdict['value'] = '@'+filenameSignal
        exposureAttrib = AttributeFinder(tmpdict)

        for fitsname in files:
            fits = pf.open(fitsname)
            # Primary extension can have a useful header but no table data
            pHeader = fits[0].header
            pHeader[filenameSignal] = fitsname
            
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
                                'of file',fitsname
                            sys.exit(1)
                    # Now have our headers, start filling in attributes in a dictionary
                    extn = {'FILENAME':f,
                            'EXTENSION':iextn}


                    # Determine EXPOSURE  for this extension
                    expo = exposureAttrib(fits, primaryHeader=pHeader, extnHeader=eHeader)
                    if expo==None:
                        print 'ERROR: no exposure ID for file',fitsname', extension',iextn
                        sys.exit(1)
                    extn['EXPOSURE']=expo

                    # Get INSTRUMENT
                    inst = attributes['INSTRUMENT'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    if inst==None:
                        print "ERROR: Missing Instrument at file",fitsname,"extension",iextn
                        sys.exit(1)
                    
                    # Other exposure-specific quantities:
                    ra = attributes['RA'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    dec = attributes['DEC'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    # ??? Check for RA already in degrees?
                    if ra!=None and dec!=None:
                        icrs = coords.ICRS(ra,dec,unit=(u.hourangle,u.degree))
                    else:
                        icrs = None
                    airmass = attributes['AIRMASS'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    exptime = attributes['EXPTIME'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    field = attributes['FIELD'](expo, primaryHeader=pHeader, extnHeader=eHeader)

                    if expo in exposures:
                        # Already have an exposure with this name, make sure this is basically the same
                        e = exposures[expo]
                        if (airmass!=None and abs(airmass-e.airmass)>0.002) or \
                           (exptime!=None and abs(exptime/e.exptime-1.)>0.0002) or \
                           (icrs!=None and icrs.separation(e.coords).degree > 1.):
                            print "ERROR: info mismatch at exposure",expo, "file",fitsname,'extension',extn
                            sys.exit(1)
                        # Also check the field for a match, unless it is _NEAREST
                        if field != '_NEAREST' and field!=e.field:
                            print "ERROR: FIELD mismatch for exposure",expo
                            print "Exposure has",e.field,"file",fitsname,"extension",iextn,"has",field
                            sys.exit(1)
                        # Check the instrument for a match
                        if inst!=e.instrument:
                            print "ERROR: INSTRUMENT mismatch for exposure",expo
                            print "Exposure has",e.instrument,"file",fitsname,"extension",iextn,"has",inst
                            sys.exit(1)
                    else:
                        # New exposure.  Need coordinates to create it
                        if icrs==None:
                            print "ERROR: Missing RA/DEC for new exposure",expo \
                              "at file",fitsname,"extension",iextn
                            sys.exit(1)
                        if field==None:
                            print "ERROR: Missing FIELD for new exposure",expo \
                              "at file",fitsname,"extension",iextn
                            sys.exit(1)
                        elif field=="_NEAREST":
                            # Finding nearest field center.
                            field=None
                            minDistance = 0.
                            for k,f in fields.items():
                                if field==None or f.distance(icrs)<minDistance:
                                    field = k
                                    minDistance = f.distance(icrs)
                        exposures[expo] = Exposure(expo, c, field, inst, airmass=airmass, exptime=exptime)

                    # Add new instrument if needed
                    if inst not in instruments and inst not in specialInstruments:
                        instruments[inst] = Instrument(inst)
                        
                    # Get device name, not needed for special devices
                    dev  = attributes['DEVICE'](expo, primaryHeader=pHeader, extnHeader=eHeader)
                    if dev==None:
                        if inst in specialInstruments:
                            # Give a blank device name
                            dev = ''
                        else:
                            print "ERROR: Missing Device at file",fitsname,"exposure",iextn
                            sys.exit(1)

                    # Record the device for this extension.
                    extn['DEVICE'] = dev
                        
                    # ??? Handle the WCSIN

                    # Now get a value for all other attributes and add to extension's dictionary
                    for k,a in attributes.items():
                        if k not in specialAttributes:
                            extn[k] = a(expo, primaryHeader=pHeader, extnHeader=eHeader)

                    # We are done with this extension.  Clear out the extension header
                    extnHeader = None
                    extensions.append(extn)

            fits.close()

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
        inst = exposures[e['EXPOSURE'].instrument]
        if inst in specialInstruments:
            dev = -1
        else:
            dev = inst.devices[e['DEVICE']]
        data.append(dev)
    cols.append(pf.Column(name='DEVICE',format=py_to_fits(data),array=data))

    # Add the WCSIN column as a variable-length character array
    data = []
    for e in extensions:
        data.append(e['WCSIN'])
    cols.append(pf.Column(name='DEVICE','PA',array=data))
        
    # Now create a column for every other attribute we have
    for a in attributes:
        if a not in specialAttributes:
            cols.append(colForAttribute(extensions, a.name))

    # Add the EXTENSION table as an HDU
    thdu = pf.new_table(pf.ColDefs(cols))
    thdu.header['EXTNAME'] = 'EXTENSIONS'
    outfits.append(thdu)
    
    # Write to FITS file
    outfits.writeto(outFile, clobber=True)
    #Done!
    
