#!/usr/bin/env python
''' Program to measure covariance ellipsoid of atmospheric turbulence (or 
generally, departures from the macro solution) for all exposures in a
WCSFit solution.  Adapting from Pedro's run_covariance.py
'''

import astropy.io.fits as pf
from astropy.table import Table
import numpy as np 
#import treecorr as tc 
from scipy.spatial import cKDTree
import sys

# Make these parameters 
inFile = 'allsf3.wcscat2.fits'
outFile = 'allsf3.turbulence.fits'
amendFile = 'tmp.fof'  # FoF file to add syserr info to.
rMin = 0.005   # Range of separation (degrees) to avg
rMax = 0.02    # to get the turbulence amplitudes
maxErr = 70.   # Maximum diagonal measurement error to retain
numThreads = 8
nbins = 10

ff = pf.open(inFile)
resids = ff['WCSOut'].data
# Map from extension to exposure numbers
expoids = ff['Extensions'].data['exposure']

# Acquire the stats of interest for all useful residuals
errorBoost = 2. / resids['chisqExpected']
use = np.logical_and(resids['Clip']!=True, errorBoost)
sig = resids['covTotalW']
use = np.logical_and(use, (sig[:,0]+sig[:,1])<2*maxErr*maxErr)

print("Using",np.count_nonzero(use),"detections")

xResid = (resids['xResW']*errorBoost)[use]
yResid = (resids['yResW']*errorBoost)[use]
xW = resids['xW'][use]
yW = resids['yW'][use]
expos = expoids[resids['extension']][use]
names = ff['Exposures'].data['name']
bands = ff['Exposures'].data['BAND']

usefulExposure = ff['Exposures'].data['INSTRUMENTNUMBER']>=0
print("Have",np.count_nonzero(usefulExposure),"potentially useful exposures")

ff.close()  # Done with direct tables

# Make a dict that will turn into our output table
out = {'name':[],'band':[],'exposure':[],'nstars':[],'npairs':[],
       'syserrxx':[],'syserryy':[],'syserrxy':[]}
doExposures,whichObjects = np.unique(expos,return_inverse=True)
print("Have data in",len(doExposures),"exposures")

# Let's do the exposures in the order of their names
order = np.argsort(names)
for expo in order:
    if not usefulExposure[expo] or expo not in doExposures:
        continue
    objects = np.where(whichObjects==np.where(doExposures==expo)[0][0])[0]

    # Extract data
    xy = np.vstack((xW[objects],yW[objects])).transpose()
    dx = xResid[objects]
    dy = yResid[objects]

    # Calculate correlations
    kdt = cKDTree(xy)
    prsmax = kdt.query_pairs(rMax,output_type='ndarray')
    prsmin = kdt.query_pairs(rMin,output_type='ndarray')
    nprs = prsmax.shape[0] - prsmin.shape[0]
    xx = np.sum(dx[prsmax[:,0]]*dx[prsmax[:,1]]) \
         -np.sum(dx[prsmin[:,0]]*dx[prsmin[:,1]])
    yy = np.sum(dy[prsmax[:,0]]*dy[prsmax[:,1]]) \
         -np.sum(dy[prsmin[:,0]]*dy[prsmin[:,1]])
    xy = np.sum(dx[prsmax[:,0]]*dy[prsmax[:,1]]) \
         -np.sum(dx[prsmin[:,0]]*dy[prsmin[:,1]])\
         +np.sum(dy[prsmin[:,0]]*dx[prsmin[:,1]])\
         -np.sum(dy[prsmin[:,0]]*dx[prsmin[:,1]])
    xx /= nprs
    yy /= nprs
    xy /= 2*nprs
    print(names[expo],bands[expo],len(objects),nprs)
    out['name'].append(names[expo])
    out['band'].append(bands[expo])
    out['exposure'].append(expo)
    out['nstars'].append(len(objects))
    out['npairs'].append(nprs)
    out['syserrxx'].append(xx)
    out['syserryy'].append(yy)
    out['syserrxy'].append(xy)

tab = Table(out)
if outFile is not None:
    tab.write(outFile,overwrite=True)

if amendFile is not None:
    # Add syserr info to Exposures extension of an existing file
    ff = pf.open(amendFile,'update')
    hdu = ff['Exposures']
    new_cols = []
    nexp = len(hdu.data)
    try:
        xxdata = hdu.data['syserrxx']
        xxdata = 0.
    except KeyError:
        xxcol = pf.Column(name='syserrxx',format='D',
                          array = np.zeros(nexp, dtype=float))
        new_cols.append(xxcol)
        xxdata = xxcol.array
    try:
        yydata = hdu.data['syserryy']
        yydata = 0.
    except KeyError:
        yycol = pf.Column(name='syserryy',format='D',
                          array = np.zeros(nexp, dtype=float))
        new_cols.append(yycol)
        yydata = yycol.array
    try:
        xydata = hdu.data['syserrxy']
        xydata = 0.
    except KeyError:
        xycol = pf.Column(name='syserrxy',format='D',
                          array = np.zeros(nexp, dtype=float))
        new_cols.append(xycol)
        xydata = xycol.array
    # Enter new info into all rows
    xxdata[out['exposure']] = out['syserrxx']
    yydata[out['exposure']] = out['syserryy']
    xydata[out['exposure']] = out['syserrxy']
    
    if len(new_cols)>0:
        # Replace hdu with a new one
        new_hdu = pf.BinTableHDU.from_columns(hdu.columns + pf.ColDefs(new_cols), 
                                              header = hdu.header,
                                              name='Exposures')
        ff['Exposures'] = new_hdu
    ff.close()

sys.exit(0)


