# Routines for E/B separations of astrometric residuals
# and investigation of correlation functions etc.
# Note the use of Mike Jarvis's treecorr below
from __future__ import division,print_function
import numpy as np
import astropy.io.fits as pf
import astropy.wcs as wcs
import pylab as pl
import re
import treecorr as tc
from scipy.spatial.ckdtree import cKDTree

import gbutil
from gbutil import decam
from gbutil.decam import ccdnums, plotColors
pixScale = 264.   # nominal mas per pixel

def extractExposures(residIn, residOut):
    '''
    Take a WCSFit output file as input.  Extract from it info on the unclipped residuals:
    WCS coordinates, fit residuals, which exposure they come from, measurement errors,
    if Gaia is involved in the fit, matchID of their star.
    Write this back out, along with the exposure table.
    (Austin: this is a routine I used to make your file.  You don't have
    to use it or worry about it yet!).
    '''
    ff = pf.open(residIn)
    resids = ff['WCSOut'].data
    # Skip the clipped objects and those with too much weight in their
    # fit.
    boostFactor = np.sqrt(2. / resids['chisqExpected'])
    skip = np.logical_or(boostFactor>1.15, resids['Clip'])
    use = np.logical_not(skip)
    # Get an array with exposure numbers per exposure
    matchID = resids['matchID'][use]
    u = resids['xW'][use]
    v = resids['yW'][use]
    dx = resids['xresW'][use] * boostFactor[use]
    dy = resids['yresW'][use] * boostFactor[use]
    extns = resids['extension'][use]
    exposure = ff['Extensions'].data['Exposure'][extns]
    objs = resids['object'][use]

    # Mark whether each is part of a match with Gaia
    inGaia = np.unique(ff['PMOut'].data['matchID'][np.logical_not(ff['PMOut'].data['Clip'])])
    hasGaia = np.isin(matchID, inGaia)

    # Look up the raw measurement error of every observation
    measErr = np.ones_like(u)*15.   ###
    for ee in np.unique(exposure):
        if ff['Exposures'].data['INSTRUMENTNUMBER'][ee]<0:
            continue
        # Cull the catalog to those in this exposure
        inExpo = np.where(exposure==ee)[0]
        expoExtns = extns[inExpo]
        expoObjs = objs[inExpo]

        print("Reading",ee,ff['Exposures'].data['NAME'][ee])

        
        cats = None
        # Gather the desired column from each extension
        for extn in np.unique(expoExtns):
            # Open FITS file first time
            targets = inExpo[expoExtns==extn]
            fitsfile = ff['Extensions'].data['FILENAME'][extn]
            fitsextn = ff['Extensions'].data['EXTENSION'][extn]
            if cats is None:
                cats = pf.open(fitsfile)
            objcat = cats[fitsextn].data
            measErr[targets] = objcat['ERRAWIN_IMAGE'][expoObjs[expoExtns==extn]]
        if cats is not None:
            cats.close()
    measErr *= pixScale # Convert pix to mas

    # Now create a new FITS table from this info
    outobj = pf.BinTableHDU.from_columns(\
            [pf.Column(name="exposure",array=exposure,format='I'),
             pf.Column(name="matchID",array=matchID,format='J'),
             pf.Column(name="u",array=u,format='D'),
             pf.Column(name="v",array=v,format='D'),
             pf.Column(name="dx",array=dx,format='D'),
             pf.Column(name="dy",array=dy,format='D'),
             pf.Column(name="measErr",array=measErr,format='D'),
             pf.Column(name="hasGaia",array=hasGaia,format='L')],
            name='Residuals')
    hdulist = pf.HDUList([pf.PrimaryHDU(), ff['Exposures'], outobj])
    hdulist.writeto(residOut,overwrite=True)

class Poly2d:
    # Polynomial function of 2 dimensions
    def __init__(self,order,sumOrder=True):
        # Set up order; sumOrder=True if order is max sum of powers of x and y.
        self.order = order
        mask = np.ones( (self.order+1, self.order+1), dtype=bool)
        if sumOrder:
            # Clear mask where sum of x & y orders is >order
            for i in range(1,self.order+1):
                mask[i, self.order-i+1:] = False
        self.use = mask
        self.sumOrder = sumOrder
        self.coeffs = None
        return
    def evaluate(self,x,y):
        xypow = np.ones( (x.shape[0],self.order+1), dtype=float)
        for i in range(1,self.order+1):
            xypow[:,i] = xypow[:,i-1] * x
        tmp = np.dot(xypow, self.coeffs)
        for i in range(1,self.order+1):
            xypow[:,i] = xypow[:,i-1] * y
        return np.sum(tmp * xypow, axis=1)
    def fit(self, x, y, z):
        # Least-square fit polynomial coefficients to P(x,y)=z
        # Make cofactor array
        npts = x.shape[0]
        xpow = np.ones( (npts,self.order+1), dtype=float)
        for i in range(1,self.order+1):
            xpow[:,i] = xpow[:,i-1] * x
        ypow = np.ones_like(xpow)
        for i in range(1,self.order+1):
            ypow[:,i] = ypow[:,i-1] * y
        A = xpow[:,:, np.newaxis] * ypow[:, np.newaxis, :]
        # Retain powers wanted
        A = A.reshape(npts,(self.order+1)*(self.order+1))[:,self.use.flatten()]
        b = np.linalg.lstsq(A,z,rcond=None)[0]
        self.coeffs = np.zeros( (self.order+1, self.order+1), dtype=float)
        self.coeffs[self.use] = b
        return
    def getCoeffs(self):
        # Return the current coefficients in a vector.
        return self.coeffs[self.use]
    def setCoeffs(self,c):
        # Set the coefficients to the specified vector.
        if len(c.shape)!=1 or c.shape[0]!=np.count_nonzero(self.use):
            print("Poly2d.setCoeffs did not get proper-size array", c.shape)
            sys.exit(1)
        self.coeffs[self.use] = c
        return

def getExposure(fitsobj,exposure,polyOrder=None):
    '''Extract the given exposure's information from 
    the residual file (already opened, as fitsobj).
    Will return a table.
    Optionally fit and subtract a polynomial of the given
    order in (u,v) from (dx,dy) data.  Unweighted fit.
    '''
    tab = fitsobj['Residuals'].data[fitsobj['Residuals'].data['exposure']==exposure]
    if polyOrder is not None:
        poly = Poly2d(polyOrder)
        poly.fit(tab['u'],tab['v'],tab['dx'])
        tab['dx'] -= poly.evaluate(tab['u'],tab['v'])
        poly.fit(tab['u'],tab['v'],tab['dy'])
        tab['dy'] -= poly.evaluate(tab['u'],tab['v'])
    return tab

    
def residInPixels(tab, binpix=256, scaleFudge=1.,
                  maxErr=50):
    '''
    Return a 2d vector diagram of weighted mean astrometric residual in pixelized areas.
    Residuals shown in milliarcsec, and the binned resid arrays are returned.
    binpix is the width of cells (in nominal pixel size)
    scaleFudge is multiplier of arrow length scaling (and key size) 
    to apply to default choice. By default, arrows are scaled so that 
    typical arrow is of length equal to the distance between arrows (cell size).

    maxErr is the largest error that a binned vector can have (mas) and still be plotted,
    to avoid cluttering up the plot with noisy arrows.
    '''

    noData = np.nan
    
    if len(tab) < 100:
        # Not enough residuals to even try.
        return None
    # use exposure coordinates:
    x = tab['u']  ### Problem is these are not relative to array center!!
    y = tab['v']
    residx = tab['dx']
    residy = tab['dy']
    sig = tab['measErr']
    xyBuffer = 0.05   # buffer in degrees
    cellSize = binpix * pixScale / (1000.*3600.)

    weight = np.where(sig>0, 1./(sig*sig), 0.)
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    
    xsize = int(np.ceil( (xmax-xmin) / cellSize))
    ysize = int(np.ceil( (ymax-ymin) / cellSize))

    # Figure out the image pixel for each residual point
    ix = np.array( np.floor( (x-xmin) / cellSize), dtype=int)
    ix = np.clip(ix, 0, xsize-1)
    iy = np.array( np.floor( (y-ymin) / cellSize), dtype=int)
    iy = np.clip(iy, 0, ysize-1)
    index = iy*xsize + ix
    
    sumxw = np.histogram(index, bins=xsize*ysize, range=(-0.5,xsize*ysize+0.5),
                        weights=weight*residx)[0]
    sumyw = np.histogram(index, bins=xsize*ysize, range=(-0.5,xsize*ysize+0.5),
                        weights=weight*residy)[0]
    sumw = np.histogram(index, bins=xsize*ysize, range=(-0.5,xsize*ysize+0.5),
                        weights=weight)[0]
    # Smallest weight that we'd want to plot a point for
    minWeight = maxErr**(-2)
    # Value to put in where weight is below threshold:
    sumxw = np.where( sumw > minWeight, sumxw / sumw, noData)
    sumyw = np.where( sumw > minWeight, sumyw / sumw, noData)
    sumw = np.where( sumw > minWeight, 1./sumw, noData)
    rmsx = np.std(sumxw[sumw>0.])
    rmsy = np.std(sumyw[sumw>0.])
    print('RMSx, RMSy, noise:', rmsx, rmsy, np.sqrt(np.mean(sumw[sumw>0.])))

    # Make an x and y position for each cell to be the center of its cell
    xpos = np.arange(xsize*ysize,dtype=int)
    xpos = ((xpos%xsize)+0.5)*cellSize + xmin
    ypos = np.arange(xsize*ysize,dtype=int)
    ypos = ((ypos//xsize)+0.5)*cellSize + ymin

    useful = np.logical_and(sumw!=0, sumw<(maxErr*maxErr))
    dx = sumxw[useful]
    dy = sumyw[useful]
    xpos = xpos[useful]
    ypos = ypos[useful]

    # Choose a scale which makes arrow length similar to spacing,
    # with a choice of round numbers for scale equal to spacing.
    typicalArrow = np.percentile(np.hypot(dx,dy), 90.)
    # Choose a scale at a round number, 1/2/5 x 10^N
    choices = np.array([1.,2.,5.])
    i = np.floor(np.log10(typicalArrow))
    j = np.where(typicalArrow / 10.**i >= choices)[0][-1]
    keyLength = choices[j] * 10**i
    arrowScale = keyLength / cellSize
    arrowScale /= scaleFudge   # Adjust lengths of arrows if desired
    
    # Now make the plot
    q=pl.quiver(xpos, ypos, dx, dy,
                pivot='middle',
                color='green',
                angles='xy',scale_units='xy',
                scale=arrowScale,
                #scale=0.001,
                units='x') ##width=5.,headlength=10,headwidth=8 )
            
    QK=pl.quiverkey(q, 0.5,0.5, keyLength, '{:.1f} mas'.format(keyLength),
                 coordinates='axes', color='red', labelpos='N',labelcolor='red')
    pl.xlim(xmin-xyBuffer,xmax+xyBuffer)
    pl.ylim(ymin-xyBuffer,ymax+xyBuffer)
    
    pl.gca().set_aspect('equal')
    pl.grid()

    return dx,dy,xpos,ypos, sumw

def calcEB(dx, dy, valid):
    """
    Given vector displacement (dx,dy) defined on identical 2d grids,
    and boolean 2d array valid which is True where data are valid,
    return arrays giving divergence and curl of the vector field.
    These will have NaN in pixels without useful info.
    """
    dxdx = dx[2:,1:-1] - dx[:-2,1:-1]
    dydx = dy[2:,1:-1] - dy[:-2,1:-1]
    dxdy = dx[1:-1,2:] - dx[1:-1,:-2]
    dydy = dy[1:-1,2:] - dy[1:-1,:-2]
    use = np.logical_and(valid[1:-1,:-2], valid[1:-1,2:])
    use = np.logical_and(use, valid[:-2,1:-1])
    use = np.logical_and(use, valid[2:,1:-1])
    div = np.where(use, dxdx+dydy, np.nan)
    curl= np.where(use, dydx-dxdy, np.nan)
    return div,curl

def ebPlot(dx, dy, xpos, ypos,scale=50.):
    """
    Plot div and curl of the specified vector field, assumed to be samples
    from a grid
    """
    # Rebuild the grid
    step = np.min(np.abs(np.diff(xpos)))
    ix = np.array( np.floor(xpos/step+0.1), dtype=int)
    iy = np.array( np.floor(ypos/step+0.1), dtype=int)
    ix = ix - np.min(ix)
    iy = iy - np.min(iy)
    dx2d = np.zeros( (np.max(iy)+1,np.max(ix)+1), dtype=float)
    dy2d = np.array(dx2d)
    valid = np.zeros_like(dx2d, dtype=bool)
    valid[iy,ix] = True
    dx2d[iy,ix] = dx
    dy2d[iy,ix] = dy
    div, curl = calcEB(dy2d, dx2d, valid)
    pl.subplot(121)
    pl.imshow(div, origin='lower', cmap='Spectral', vmin=-scale, vmax=scale)
    pl.colorbar()
    pl.title('Divergence')
    pl.subplot(122)
    pl.imshow(curl, origin='lower', cmap='Spectral', vmin=-scale, vmax=scale)
    pl.colorbar()
    pl.title('Curl')
    # Some stats:
    vardiv = gbutil.clippedMean(div[div==div],5.)[1]
    varcurl = gbutil.clippedMean(curl[div==div],5.)[1]
    print("RMS of div: {:.2f}; curl: {:.2f}".format(np.sqrt(vardiv),
                                                    np.sqrt(varcurl)))
    return
    
def vcorr(expotab,
          rmin=5./3600., rmax=1.5, dlogr=0.05,
          maxpts = 30000):
    """
    Produce angle-averaged 2-point correlation functions of astrometric error
    for the supplied sample of data, using brute-force pair counting.
    Output are the following functions:
    logr - mean log of radius in each bin
    xi_+ - <vr1 vr2 + vt1 vt2> = <vx1 vx2 + vy1 vy2>
    xi_- - <vr1 vr2 - vt1 vt2>
    xi_x - <vr1 vt2 + vt1 vr2>
    xi_z2 - <vx1 vx2 - vy1 vy2 + 2 i vx1 vy2>
    """

    if len(expotab) > maxpts:
        # Subsample array to get desired number of points
        rate = float(maxpts) / len(x)
        print("Subsampling rate {:5.3f}%".format(rate*100.))
        use = np.random.random(len(expotab)) <= rate
        subtab = expotab[use]
        u = subtab['u']
        v = subtab['v']
        dx = subtab['dx']
        dy = subtab['dy']
    else:
        u = expotab['u']
        v = expotab['v']
        dx = expotab['dx']
        dy = expotab['dy']

    print("Length ",len(u))
    # Get index arrays that make all unique pairs
    i1, i2 = np.triu_indices(len(u))
    # Omit self-pairs
    use = i1!=i2
    i1 = i1[use]
    i2 = i2[use]
    del use
    
    # Make complex separation vector
    dr = 1j * (v[i2]-v[i1])
    dr += u[i2]-u[i1]

    # log radius vector used to bin data
    logdr = np.log(np.absolute(dr))
    logrmin = np.log(rmin)
    bins = int(np.ceil(np.log(rmax/rmin)/dlogr))
    hrange = (logrmin, logrmin+bins*dlogr)
    counts = np.histogram(logdr, bins=bins, range=hrange)[0]
    logr = np.histogram(logdr, bins=bins, range=hrange, weights=logdr)[0] / counts

    # First accumulate un-rotated stats
    vec =  dx + 1j*dy
    vvec = dx[i1] * dx[i2] + dy[i1] * dy[i2]
    xiplus = np.histogram(logdr, bins=bins, range=hrange, weights=vvec)[0]/counts
    vvec = vec[i1] * vec[i2]
    xiz2 = np.histogram(logdr, bins=bins, range=hrange, weights=vvec)[0]/counts

    # Now rotate into radial / perp components
    print(type(vvec),type(dr)) ###
    tmp = vvec * np.conj(dr)
    vvec = tmp * np.conj(dr)
    dr = dr.real*dr.real + dr.imag*dr.imag
    vvec /= dr
    del dr
    ximinus = np.histogram(logdr, bins=bins, range=hrange, weights=vvec)[0]/counts
    xicross = np.imag(ximinus)
    ximinus = np.real(ximinus)

    return logr, xiplus, ximinus, xicross, xiz2

def xiEB(logr, xiplus, ximinus):
    """
    Return estimate of pure E- and B-mode correlation functions from xi+-
    """
    # Integral of d(log r) ximinus(r) from r to infty:
    dlogr = np.zeros_like(logr)
    dlogr[1:-1] = 0.5*(logr[2:] - logr[:-2])
    tmp = np.array(ximinus) * dlogr
    integral = np.cumsum(tmp[::-1])[::-1]
    xiB = 0.5*(xiplus-ximinus) + integral
    xiE = xiplus - xiB
    return xiE,xiB

def fastXi(x,y,dx,dy,rrange=(5./3600.,1.5), nbins=100):
    """
    Use treecorr to get xi of (dx*dx + dy*dy) (xiplus).
    Returns arrays of mean radius in each (log-spaced) bin, mean xi.
    """
    catx = tc.Catalog(x=x, y=y, k=dx)
    caty = tc.Catalog(x=x, y=y, k=dy)
    kk = tc.KKCorrelation(min_sep=rrange[0], max_sep=rrange[1], nbins=nbins)
    kk.process(catx)
    xx = np.array(kk.xi)

    kk = tc.KKCorrelation(min_sep=rrange[0], max_sep=rrange[1], nbins=nbins)
    kk.process(caty)
    yy = np.array(kk.xi)

    return np.exp(kk.meanlogr), xx+yy

def vcorr2d(expotab,
            rmax=1., bins=513,
            maxpts = 30000):
    """
    Produce 2d 2-point correlation function of total displacement power
    for the supplied sample of data, using brute-force pair counting.
    Output are 2d arrays giving the 2PCF and then the number of pairs that
    went into each bin.  The 2PCF calculated is
    xi_+ - <vr1 vr2 + vt1 vt2> = <vx1 vx2 + vy1 vy2>
    Output function will be symmetric about origin.
    """

    hrange = [ [-rmax,rmax], [-rmax,rmax] ]
    if len(expotab) > maxpts:
        # Subsample array to get desired number of points
        rate = float(maxpts) / len(x)
        print("Subsampling rate {:5.3f}%".format(rate*100.))
        use = np.random.random(len(expotab)) <= rate
        subtab = expotab[use]
        u = subtab['u']
        v = subtab['v']
        dx = subtab['dx']
        dy = subtab['dy']
    else:
        u = expotab['u']
        v = expotab['v']
        dx = expotab['dx']
        dy = expotab['dy']
    print("Length ",len(u))
    # Get index arrays that make all unique pairs
    i1, i2 = np.triu_indices(len(u))
    # Omit self-pairs
    use = i1!=i2
    i1 = i1[use]
    i2 = i2[use]
    del use
    
    # Make separation vectors and count pairs
    yshift = v[i2]-v[i1]
    xshift = u[i2]-u[i1]
    counts = np.histogram2d(xshift,yshift, bins=bins, range=hrange)[0]

    # Accumulate displacement sums
    vec =  dx + 1j*dy
    vvec = dx[i1] * dx[i2] + dy[i1] * dy[i2]
    xiplus = np.histogram2d(xshift,yshift, bins=bins, range=hrange, weights=vvec)[0]/counts

    xiplus = 0.5*(xiplus + xiplus[::-1,::-1])  # Combine pairs
    return xiplus,counts

def covAnnulus(expotab, rMin=0.005, rMax=0.02, maxErr=70.):
    '''Calculate the mean <dx dx>, <dx dy>, and <dy dy> for
    pairs within a circular annulus of separation.
    Exclude individual points with meas errors above maxErr.
    
    Returns cxx, cyy, cxy, npairs where last is number of
    object pairs used in calculation.
    '''
    # Extract data
    use = expotab['measErr'] <= maxErr
    xy = np.vstack((expotab['u'][use],expotab['v'][use])).transpose()
    dx = expotab['dx'][use]
    dy = expotab['dy'][use]

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
         +np.sum(dy[prsmax[:,0]]*dx[prsmax[:,1]])\
         -np.sum(dy[prsmin[:,0]]*dx[prsmin[:,1]])
    xx /= nprs
    yy /= nprs
    xy /= 2*nprs
    
    return xx,yy,xy,nprs
