#!/usr/bin/env python
'''
Produce photometric prior file that breaks degeneracies arising from zeropoint shifts.
Output is YAML format suitable for input to PhotoFit.
Nominal airmass extinction coefficient is taken from table if none is given
'''
import astropy.io.fits as pf
from astropy.time import Time
import numpy as np
import sys
import yaml
import argparse
import copy

parser = argparse.ArgumentParser(description='Configure absolute and nightly relative photometric priors')
parser.add_argument("dbfile", help='Configuration database filename', type=str)
parser.add_argument('band',help='Photometric band to configure', type=str)
parser.add_argument('priorfile',help='Name of output YAML file with prior specs', type=str)

parser.add_argument('--sigma', help='Expected RMS zeropoint fluctuation (mmag)', type=float,
                        default=1.)
parser.add_argument('--extinction', help='Nominal nightly coefficient of airmass (mag/X)', type=float)
parser.add_argument('--apcorr', help='Nominal coefficient of apcorr', type=float,
                        default=1.5)
parser.add_argument('--free', help='Which nightly coefficients float? (zeropoint,airmass,color,apcorr)',
                        type=str, default="zeropoint")
parser.add_argument('--absolutes', help='Ordered, comma-separated list of exposures to '
                        'set as absolute references',
                        type=str, default="")
parser.add_argument('--airmass_span',help='Minimum span of airmass to have coefficient float',
                        type=float,default=0.03)
parser.add_argument('--apcorr_span',help='Minimum span of apcorr to have coefficient float',
                        type=float,default=0.002)
parser.add_argument('--mjd_noon',help='Fraction to add to integral MJD to get start of "night"',
                        type=float,default=0.7)
parser.add_argument('--device',help='Device holding photometric reference point',type=str)
parser.add_argument('--xpix',help='x coordinate of photometric reference point',type=float)
parser.add_argument('--ypix',help='y coordinate of photometric reference point',type=float)

args = parser.parse_args()

db = pf.open(args.dbfile)
b = args.band

# Find instruments in specified band
useInstruments = set( [h.header['NUMBER'] for h in db[1:] \
                           if 'EXTNAME' in h.header and \
                           (h.header['EXTNAME']=='Instrument' or h.header['EXTNAME']=='INSTRUMENT')\
                           and h.header['BAND']==b] )

# Get a default airmass coeff for the band
# Extinction coeffs from https://cdcvs.fnal.gov/redmine/projects/descalibration/wiki
nominalExtinction = {'u':0.45, 'g':0.20,'r':0.10,'i':0.07,'z':0.08,'Y':0.07}
if args.extinction is not None:
    extinction = args.extinction
elif b in nominalExtinction:
    extinction = nominalExtinction[b]
else:
    print "WARNING: no default extinction found for band",b
    extinction = 0.1

# Put nominal values into dictionary
baseline = {'Zeropoint':0.,
            'Airmass':-extinction,
            'Color':0.,
            'Apcorr':args.apcorr,
            'Sigma':args.sigma*0.001,
            'Free':[]}
if args.device is not None:
    baseline['Device'] = args.device
if args.xpix is not None:
    baseline['XPixel'] = args.xpix
if args.ypix is not None:
    baseline['YPixel'] = args.ypix

# Which nightly coefficients are free to vary, by default?
free = {i.strip() for i in args.free.split(',') if i.strip()}
for j in ('zeropoint','airmass','color','apcorr'):
    if j in free:
        baseline['Free'].append(j)
        free.remove(j)
if free:
    # Shouldn't be anything left in the list
    print 'Unknown arguments for --free:',free
    sys.exit(1)

absolutes = [i.strip() for i in args.absolutes.split(',') if i.strip()]

class Night:
    '''Class collecting the information we need about each night's observations
    '''
    def __init__(self, d):
        # Initialize with dictionary to be written to YAML
        self.d = copy.deepcopy(d)
        self.exposures = set()  # Exposures on this night
        self.fields = set()  # Fields observed this night
    def add(self,exposure,field,airmass,apcorr):
        # Add an exposure to this night
        if not self.exposures:
            # no data yet
            self.min_airmass = airmass
            self.max_airmass = airmass
            self.min_apcorr = apcorr
            self.max_apcorr = apcorr
        else:
            self.min_airmass = min(airmass, self.min_airmass)
            self.max_airmass = max(airmass, self.max_airmass)
            self.min_apcorr = min(apcorr, self.min_apcorr)
            self.max_apcorr = max(apcorr, self.max_apcorr)
        self.exposures.add(exposure)
        self.fields.add(field)
    def prepare(self, airmass_span, apcorr_span):
        ''' Prepare dict to fully describe prior.
        Airmass and apcorr will be held fixed if exposures do not
        span at least the specified ranges.
        '''
        if not self.exposures:
            # No data here
            return
        # Freeze airmass and apcorr coefficients if needed
        if self.max_airmass - self.min_airmass < airmass_span:
            if 'airmass' in self.d['Free']:
                self.d['Free'].remove('airmass')
        if self.max_apcorr - self.min_apcorr < apcorr_span:
            if 'apcorr' in self.d['Free']:
                self.d['Free'].remove('apcorr')
        # Add list of all exposures to dictionary
        self.d['Exposures'] = [i for i in self.exposures]
            
# Make dict of nights with exposures
nights = {}

# And the root node of the YAML file
root = {}

# Go through exposures, adding each one in the band to appropriate (night,field) list, recording
# airmass ranges per night too.
eTable = db['exposures'].data
mjds = np.floor(eTable['mjd'] - args.mjd_noon)
for i in range(len(eTable)):
    if eTable['InstrumentNumber'][i] not in useInstruments:
        continue
    mjd = int(mjds[i])  # Get integer MJD
    x = eTable['airmass'][i]
    name = ''.join(eTable['name'][i])
    if mjd not in nights:
        # Create new entry
        nights[mjd] = Night(baseline)
    nights[mjd].add(name,
                    eTable['fieldnumber'][i],
                    eTable['airmass'][i],
                    eTable['apcorr'][i])

# Produce relative constraint priors for each night's data
# Make the constraint name have a calender date instead of MJD
for mjd in nights:
    t = Time(mjd,format='mjd')
    iso = t.iso
    ymd = iso[:4] + iso[5:7] + iso[8:10]
    constraintName = 'Relative_{:s}{:s}'.format(b,ymd)
    nights[mjd].prepare(args.airmass_span,args.apcorr_span)
    # Add to YAML output
    root[constraintName] = copy.deepcopy(nights[mjd].d)
    # Add MJD to YAML file for documenation
    root[constraintName]['MJD'] = mjd

# Now start finding which exposures must be set as absolute references to
# remove any shift degeneracies from the solution.
# The assumption is that any exposure in the same night *or* in the same field
# will have its zeropoint connected to the same absolute level by either overlap
# or a relative prior

absolutes_out = []

while nights:
    # Pick one exposure from one night to fix as absolute
    mjd = None
    exposure = None
    while absolutes:
        # Start by trying the ones specified at cmd line, in order
        expo = absolutes.pop(0)
        # Which night does the exposure belong to?
        for k,v in nights.items():
            if expo in v.exposures:
                mjd = k
                exposure = expo
                break
        if mjd is not None:
            # Found it.
            break
        print 'WARNING: Did not use {:s} as absolute reference exposure'.format(expo)

    # If there are no more user suggestions, pull out an exposure from front of lists
    if mjd is None:
        mjd = nights.keys()[0]
        exposure = nights[mjd].exposures.pop()

    print 'Absolute contraint at exposure',exposure,'on',mjd
    print '...Checking',len(nights),'nights' ##
    
    # Add an absolute prior for this point
    absolutes_out.append(exposure)
    
    # Now "cross out" all nights that are tied to this exposure's night by common fields
    newFields = nights[mjd].fields.copy()

    # Iterate graph traversal until no new fields found
    while newFields:
        newerFields = set()
        for k,n in nights.items():
            if n.fields & newFields:
                # This night touches calibrated fields.
                # Record any fields it touches that are not already
                # marked as connected to our absolute calibrator
                newerFields |= n.fields - newFields
                # And get rid of this night
                nights.pop(k)
        # Now prepare for another trip through the nights
        newFields = newerFields
        
# Write the absolute constraint node
constraintName = 'Absolute_'+b
# It need only have the Exposures and sigma keywords in the node.
root[constraintName] = {'Exposures':absolutes_out,
                        'Sigma':args.sigma*0.001 }
# write the YAML file
with open(args.priorfile,'w') as f:
    yaml.dump(root,f)

sys.exit(0)



