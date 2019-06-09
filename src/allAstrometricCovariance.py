#!/usr/bin/env python
''' Program to measure covariance ellipsoid of atmospheric turbulence (or 
generally, departures from the macro solution) for all exposures in a
WCSFit solution.  Adapting from Pedro's run_covariance.py
'''

from __future__ import division,print_function
import astropy.io.fits as pf
from astropy.table import Table
import numpy as np 
from scipy.spatial import cKDTree
import sys
import argparse


def getCovariance(inFile, rMin=0.005, rMax=0.02, maxErr=70., minPairs=500):
    '''Calculate and return a table of estimating cov matrix for
    spatially correlated astrometric errors in each exposure.
    inFile  = filename of WCSFit output file
    rMin,rMax = range of pair separations to avg for cov (degrees)
    maxErr = largest measurement error allowed for input stars
    minPairs = no. of residual pairs needed for valid covariance matrix

    Returns an astropy Table giving info on each exposure that was
    calculable
    '''

    ff = pf.open(inFile)
    resids = ff['WCSOut'].data
    # Map from extension to exposure numbers
    expoids = ff['Extensions'].data['exposure']

    # Acquire the stats of interest for all useful residuals
    errorBoost = np.sqrt(2. / resids['chisqExpected'])
    use = np.logical_and(resids['Clip']!=True, errorBoost<1.2)
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
    print("Have",np.count_nonzero(usefulExposure),
          "potentially useful exposures")

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
        if nprs < minPairs:
            # Not enough data for this one
            print("Warning: insufficient pairs for exposure",names[expo],
                  nprs)
            continue;
        xx = np.sum(dx[prsmax[:,0]]*dx[prsmax[:,1]]) \
             -np.sum(dx[prsmin[:,0]]*dx[prsmin[:,1]])
        yy = np.sum(dy[prsmax[:,0]]*dy[prsmax[:,1]]) \
             -np.sum(dy[prsmin[:,0]]*dy[prsmin[:,1]])
        xy = np.sum(dx[prsmax[:,0]]*dy[prsmax[:,1]]) \
             -np.sum(dx[prsmin[:,0]]*dy[prsmin[:,1]])\
             +np.sum(dy[prsmax[:,0]]*dx[prsmax[:,1]])\
             -np.sum(dy[prsmin[:,0]]*dx[prsmin[:,1]])
        xx /= nprs
        yy /= nprs
        xy /= 2*nprs
        out['name'].append(names[expo])
        out['band'].append(bands[expo])
        out['exposure'].append(expo)
        out['nstars'].append(len(objects))
        out['npairs'].append(nprs)
        out['syserrxx'].append(xx)
        out['syserryy'].append(yy)
        out['syserrxy'].append(xy)

    return Table(out)

def minEigenvalue(sysTable,minSigma):
    '''Change entries for syserr ellipse in table such that no
    eigenvalue is less than minSigma**2.  In particular, enforces
    pos-def ellipses.
    '''
    mm = np.array( ( (sysTable['syserrxx'],sysTable['syserrxy']),
                     (sysTable['syserrxy'],sysTable['syserryy']) ) )
    mm = np.moveaxis(mm,2,0)  # Now is Nx2x2 array
    ww, vv = np.linalg.eigh(mm)  # Get eigenvalues, eigenvecs
    # Place floor on eigenvalues
    ww = np.clip(ww, minSigma*minSigma, None)
    # Reconstruct the cov matrices
    mm = np.matmul(np.moveaxis(vv,1,2)*ww[:,np.newaxis,:],vv)
    # Put back into table
    sysTable['syserrxx'] = mm[:,0,0]
    sysTable['syserrxy'] = mm[:,0,1]
    sysTable['syserryy'] = mm[:,1,1]

def updateDB(amendFile, sysTable, defaultSigma=10., clear=False):
    '''Add/update columns to the amendFile (which is usually the
    output of a WCSFoF run) containing the syserrxx, etc., values
    stored in the sysTable.
    clear = True will zero out any syserr entries already in the
            amendFile columns, and set syswarn = True for exposures
            without valid covariance calculation.  Otherwise
            any prior values will be kept if they are not overwritten.
    defaultSigma = value of error per axis (mas) to be installed
            in exposures upon clearing (syswarn is still True for them).
    '''
    # Add syserr info to Exposures extension of an existing file
    ff = pf.open(amendFile,'update')
    hdu = ff['Exposures']
    new_cols = []
    nexp = len(hdu.data)
    defaultVar = defaultSigma**2
    try:
        xxdata = hdu.data['syserrxx']
        if clear:
            xxdata = np.ones(nexp,dtype=float) * defaultVar;
    except KeyError:
        xxcol = pf.Column(name='syserrxx',format='D',
                          array = np.ones(nexp, dtype=float)*defaultVar)
        new_cols.append(xxcol)
        xxdata = xxcol.array
    try:
        yydata = hdu.data['syserryy']
        if clear:
            yydata = np.ones(nexp,dtype=float) * defaultVar;
    except KeyError:
        yycol = pf.Column(name='syserryy',format='D',
                          array = np.ones(nexp, dtype=float) * defaultVar)
        new_cols.append(yycol)
        yydata = yycol.array
    try:
        xydata = hdu.data['syserrxy']
        if clear:
            xydata = np.zeros(nexp,dtype=float)
    except KeyError:
        xycol = pf.Column(name='syserrxy',format='D',
                          array = np.zeros(nexp, dtype=float))
        new_cols.append(xycol)
        xydata = xycol.array
    try:
        warndata = hdu.data['syserrwarn']
        if clear:
            warndata = np.ones(nexp, dtype=bool)
    except KeyError:
        warncol = pf.Column(name='syserrwarn',format='L',
                          array = np.ones(nexp, dtype=bool))
        new_cols.append(warncol)
        warndata = warncol.array
    # Enter new info into all rows
    xxdata[sysTable['exposure']] = sysTable['syserrxx']
    yydata[sysTable['exposure']] = sysTable['syserryy']
    xydata[sysTable['exposure']] = sysTable['syserrxy']
    warndata[sysTable['exposure']] = False  # Clear warning flag where we have values
    
    if len(new_cols)>0:
        # Replace hdu with a new one
        new_hdu = pf.BinTableHDU.from_columns(hdu.columns 
                                              + pf.ColDefs(new_cols), 
                                              header = hdu.header,
                                              name='Exposures')
        ff['Exposures'] = new_hdu
    ff.close()

if __name__=='__main__':

    parser = argparse.ArgumentParser(description=\
     'Calculate correlated covariance of astrometric errors for each exposure')
    parser.add_argument("input", help='Filename of WCSFit residuals', 
                        type=str)
    parser.add_argument("-o", "--output", 
                        help='Filename for output covariance table',type=str)
    parser.add_argument("-a", "--amend", 
                        help='Filename of FoF file to augment' + \
                             ' with exposure covariances', type=str)
    parser.add_argument('--rMin',help='Minimum separation for pairs (degree)',
                        type=float, default=0.001)
    parser.add_argument('--rMax',help='Maximum separation for pairs (degree)',
                        type=float, default=0.02)
    parser.add_argument('--maxErr',help='Maximum usable measurement error' + \
                        ' for point (mas)', 
                        type=float, default=100.)
    parser.add_argument('--minPairs',help='Minimum number of pairs to use',
                        type=int, default=300)
    parser.add_argument('--minSigma',
                        help='Floor to place on error ellipse axes',
                        type=float, default=0.)
    parser.add_argument('--clear',help='Erase all existing syserr entries' + \
                        ' in amend file?', action='store_true')
    parser.add_argument('--defaultSigma',
                        help='Error assigned if insufficent pairs (mas)', 
                        type=float, default=20.)
    args = parser.parse_args()

    tab = getCovariance(args.input, 
                        rMin=args.rMin, rMax=args.rMax,
                        maxErr=args.maxErr, minPairs=args.minPairs);
    
    # Clip syserr ellipse axes if requested
    if args.minSigma > 0:
        minEigenvalue(tab, args.minSigma)

    # Save the output to its own file if desired:
    if args.output is not None:
        tab.write(args.output,overwrite=True)
    
    # And insert into FoF file if desired:
    if args.amend is not None:
        updateDB(args.amend, tab, 
                 clear=args.clear, defaultSigma=args.defaultSigma)

    sys.exit(0)

