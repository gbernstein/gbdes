# Routines for use in plotting astrometry residuals
from __future__ import division,print_function
import numpy as np
import astropy.io.fits as pf
import astropy.wcs as wcs
import pylab as pl
import re
from scipy.interpolate import UnivariateSpline
from scipy.signal import medfilt
import gbutil
from gbutil import decam
from gbutil.decam import ccdnums, plotColors
pixScale = 264.   # nominal mas per pixel

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
        b = np.linalg.lstsq(A,z)[0]
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

class Matcher:
    """
    Class is initialized with one or a sequence of strings that are treated as
    regular expressions.  The () operation returns True if the argument matches
    any of the regular expressions.
    Initializing with None or empty list means anything will match.
    """
    def __init__(self,exprs=None):
        self.relist = []
        if exprs==None: return
        if type(exprs)==str:
            # Append character to RE to force match of full string
            self.relist.append(re.compile(exprs + '$'))
        else:
            for e in exprs:
                self.relist.append(re.compile(e + '$'))
        return
    def __call__(self,arg):
        if len(self.relist)==0:
            return True
        for r in self.relist:
            if r.match(arg):
                return True
        return False
    def isNull(self):
        # Return true if the matcher will match everything
        return len(self.relist)==0
            
def exposuresUsing(fits, instrument):
    """
    Return list of names of exposures using the chosen instrument
    """
    # First find all instruments that work

    iMatch = Matcher(instrument)
    
    # Open the fits file from name given, or just assume it's already an open FITS object
    if type(fits)==str:
        ff = pf.open(fits)
    else:
        ff = fits

    expTab = ff['Exposures'].data
    useInst = []
    
    # First find instruments and devices that match criteria
    for hdu in ff[1:]:
        if hdu.header['EXTNAME']=='Instrument' and iMatch(hdu.header['NAME']):
            # Found an instrument we'll use
            useInst.append(hdu.header['NUMBER'])

    out = []
    for inst in useInst:
        for exp in expTab['Name'][expTab['InstrumentNumber']==inst]:
            if type(exp)==str:
                out.append(exp)
            else:
                # Sometimes strings come back as lists of characters...
                out.append(''.join(exp))
    return out

def select(fits, instrument=None, device=None, exposure=None,
           clip=False, reserve=None, amp=None, wtFracMax=0.35,
           fieldCoords=False):
    """
    Returns rows of the object table in the named PhotoFit output fits file that
    meet ALL of the given selection criteria.  fits can also be an opened fits object.
    A selection =None means no cut.
    
    instrument matches (one of) the string(s) given.
    device name matches string given
    exposure name matches string given
    exposures matches exposure names.
    clip selects True or False as given, same for reserve
    amp can be "A" or "B" and limits selection to appropriate side of N and S CCDs
    wtFracMax is an upper limit on wtFrac to accept
    fieldCoords says whether the "world" coords will be relative to the field center (True)
       or relative to the nominal exposure RA and Dec (False)
    """
    
    # First find all (instrument, device) pairs and instruments that work
    useInst = {}

    iMatch = Matcher(instrument)
    dMatch = Matcher(device)
    eMatch = Matcher(exposure)
    reS = re.compile('^S.*')  # Regular expression checking for South-side devices
    abSplit = 1024.   # X pixel value splitting A and B amplifiers
    
    # Open the fits file from name given, or just assume it's already an open FITS object
    if type(fits)==str:
        ff = pf.open(fits)
    else:
        ff = fits

    expTab = ff['Exposures'].data
    extnTab = ff['Extensions'].data
    objTab = ff['WCSOut'].data
    use = np.zeros(len(objTab), dtype=bool)

    maxDev = np.max(extnTab['DEVICE'])

    # First find instruments and devices that match criteria
    for hdu in ff[1:]:
        if hdu.header['EXTNAME']=='Instrument':
            if not iMatch(hdu.header['NAME']):
                continue
            # Found an instrument we'll use
            inst = hdu.header['NUMBER']
            # Make bool arrays saying whether each device is used for N or S
            useDevN = np.zeros(maxDev+1,dtype=bool)
            useDevS = np.zeros_like(useDevN)
            devnames = hdu.data['NAME']
            for i,k in enumerate(devnames):
                if dMatch(k):
                    if reS.match(k):
                        useDevS[i]=True
                    else:
                        useDevN[i]=True

            # Find exposures using selected instrument and matching spec
            useExp = np.where(expTab['InstrumentNumber']==inst)[0]
            names = expTab['Name'][useExp]
            keepExp = []
            for i,n in enumerate(names):
                if eMatch(n):
                    keepExp.append(useExp[i])
            
            useExp = np.zeros(len(expTab),dtype=bool)
            useExp[keepExp] = True

            # Now make lists of extension numbers in selected S and N exposure/device pairs
            useExtnS = np.logical_and(useExp[extnTab['EXPOSURE']],
                                      useDevS[extnTab['DEVICE']])
            useExtnN = np.logical_and(useExp[extnTab['EXPOSURE']],
                                      useDevN[extnTab['DEVICE']])

            # Collect objects from these extension, applying amp cut if any
            if amp is None:
                # Not using A/B, just have one list
                useExtnLow = np.logical_or(useExtnS,useExtnN)
                use = np.logical_or(use, useExtnLow[objTab['EXTENSION']])
            elif amp=='A':
                lowX = objTab['XPIX'] < abSplit
                use = np.logical_or(use,
                                    np.logical_and(useExtnS[objTab['EXTENSION']], lowX))
                use = np.logical_or(use,
                                    np.logical_and(useExtnN[objTab['EXTENSION']], ~lowX))
            elif amp=='B':
                lowX = objTab['XPIX'] < abSplit
                use = np.logical_or(use,
                                    np.logical_and(useExtnS[objTab['EXTENSION']], ~lowX))
                use = np.logical_or(use,
                                    np.logical_and(useExtnN[objTab['EXTENSION']], lowX))
            else:
                print('Bad amp value:', amp)
                sys.exit(1)

    if iMatch('REFERENCE') and dMatch(''):
        # include objects from reference exposures
        useExp = np.where(expTab['InstrumentNumber']==-1)[0]
        names = expTab['Name'][useExp]
        keepExp = []
        for i,n in enumerate(names):
            if eMatch(n):
                keepExp.append(useExp[i])
            
        useExp = np.zeros(len(expTab),dtype=bool)
        useExp[keepExp] = True

        # Now make lists of extension numbers 
        useExtn = useExp[extnTab['EXPOSURE']]
        use = np.logical_or(use, useExtn[objTab['EXTENSION']])

    if clip is not None:
        use = np.logical_and(use, (objTab['CLIP']>0)==clip)
    if reserve is not None:
        use = np.logical_and(use, (objTab['RESERVE']>0)==reserve)
    if wtFracMax is not None:
        use = np.logical_and(use, objTab['WTFRAC']<=wtFracMax)

    ## print('Select',np.count_nonzero(use), '/',len(use))
                
    outTable = objTab[use]

    if not fieldCoords:
        # Shift the world coordinates in exposure system, not field
        exposures = extnTab['exposure'][outTable['Extension']]
        exposuresInUse,invert = np.unique(exposures,return_inverse=True)
        dx, dy = getFieldCoords(ff, exposuresInUse)
        outTable['xw'] -= dx[invert]
        outTable['yw'] -= dy[invert]
        
    return outTable

def getFieldCoords(fits, expNums):
    """
    Update or create a new column in the exposure table giving the exposure
    axis position in gnomonic position about its field center.
    fits is the opened residual file
    expNums is an array of exposure numbers of interest

    Returns dx, dy, which are the locations (in degrees) of the exposure
    centers in the field world coordinate system.  
    """
    # check for field and instrument exposures
    try:
        exposures = fits['Exposures'].data
        fields = fits['Fields'].data
    except KeyError:
        print("File does not have necessary extensions")
        return 1
    # Add columns to the exposure table ???
    dx = []
    dy = []
    for iexpo in expNums:
        iinst = exposures['instrumentnumber'][iexpo]
        if iinst<0:
            # Do not offset reference catalogs
            dx.append(0.)
            dy.append(0.)
            continue
        #  find field RA & dec
        ifield = exposures['FIELDNUMBER'][iexpo]
        fra = fields['RA'][ifield]
        fdec = fields['Dec'][ifield]
        ra = exposures['RA'][iexpo]
        dec = exposures['Dec'][iexpo]
        # Calculate gnomonic projection (degrees here by default)
        w = wcs.WCS(naxis=2)
        w.wcs.ctype = ["RA---TAN","DEC--TAN"]
        w.wcs.crval = np.array([fra, fdec])
        dxy = w.wcs_world2pix(ra,dec,1)
        dx.append(dxy[0])
        dy.append(dxy[1])
    return np.array(dx), np.array(dy)

def residInPixels(fits, binpix=64, scaleFudge=1.,
                  maxErr=50, pixCoords=True, **kw):
    noData = np.nan
    
    # Return a 2d vector diagram of weighted mean astrometric residual in pixelized areas.
    # Binning is done in pixel coordinates if pixCoords=True (all devices will be stacked)
    # else done in exposure coordinates (projection about axis) ???.
    # Residuals shown in milliarcsec in either case.
    # binpix is the width of cells (in nominal pixel size)
    # scaleFudge is multiplier of arrow length scaling (and key size) to apply to default choice
    # By default, arrows are scaled so that typical arrow is of length equal to
    # the distance between arrows (cell size).
    # comparable to the 

    # Get the desired objects
    tab = select(fits,**kw)
    if len(tab) < 100:
        return None
    # use exposure coordinates:
    if pixCoords:
        x = tab['xPix']
        y = tab['yPix']
        residx = tab['xresPix'] * pixScale/(1-tab['wtfrac'])
        residy = tab['yresPix'] * pixScale/(1-tab['wtfrac'])
        sig = tab['sigW']  # Already in mas
        xyBuffer = 50   # 50 pixels around edges
        cellSize = binpix
    else:
        x = tab['xW']  ### Problem is these are not relative to array center!!
        y = tab['yW']
        residx = tab['xresW']/(1-tab['wtfrac'])
        residy = tab['yresW']/(1-tab['wtfrac'])
        sig = tab['sigW']
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
    #text(-0.05*ccd.SIZE_CCD_X, 1.05*ccd.SIZE_CCD_Y, 'median: %.3g pix (%.3g mas)' %(ccd.get_median_residual(), ccd.get_median_residual()*ccd.plate_scale*1000. ), size=6.5)
            
    QK=pl.quiverkey(q, 0.5,0.5, keyLength, '{:.1f} mas'.format(keyLength),
                 coordinates='axes', color='red', labelpos='N',labelcolor='red')
    pl.xlim(xmin-xyBuffer,xmax+xyBuffer)
    pl.ylim(ymin-xyBuffer,ymax+xyBuffer)
    
    pl.gca().set_aspect('equal')
    pl.grid()

    if type(fits)==str:
        # Give a title if we have the filename
        pl.title('Binned residuals for ' + fits)
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
    
def xPlot(fits, cellSize, doY=False, maxErr=2, color='g', **kwargs):
    """
    Stack residuals along one or other axis given a fits wcs output table
    cellSize is the approximate box size (in pixels) in which detections will be binned
    maxErr is the biggest error in cell mean residual that will get plotted (in pixels)
    doY = True will plot the Y displacement, default False to plot X
    kwargs is a dict giving selection criteria for the plot as per select() function
    """

    t = select(fits, **kwargs)      # Get the data to be plotted

    # Now set up grid of cells
    if doY:
        x = t['ypix']
        res = t['yrespix']
    else:
        x = t['xpix']
        res = t['xrespix']

    xmin = np.min(x)
    xmax = np.max(x)

    nx = int( np.floor( (xmax-xmin)/cellSize))
    xCell = (xmax - xmin) / nx

    # Decide which cell number each point is in
    index = np.zeros(len(t), dtype=int)
    np.clip( np.floor((x-xmin)/xCell), 0, nx-1, index)

    # Calculate weighted mean of the x and y residuals
    w = t['sigpix'] / (1-t['wtFrac'])
    w = np.where(w>0, 1/(w*w), 0.)
    sumxw = np.histogram(index, bins=nx, range=(0,nx), weights= res*w)[0]
    sumw =  np.histogram(index, bins=nx, range=(0,nx), weights= w)[0]

    dx = sumxw / sumw
    err = np.where(sumw>0. , pixScale/np.sqrt(sumw), 2*maxErr)
    print(err)

    #useful = np.logical_and(sumw!=0, sumw>1./(maxErr*maxErr))
    useful = err < maxErr

    # Make an x and y position for each cell to be the center of its cell
    xpos = np.arange(nx,dtype=int)
    xpos = (xpos+0.5)*xCell + xmin

    dx = dx[useful]
    xpos = xpos[useful]

    # Now make the plot
    #pl.plot(xpos, dx * pixScale, color+'o')
    pl.errorbar(xpos, dx*pixScale, yerr=err[useful], ecolor=color, color=color, marker='o',ls='None')
    
    pl.xlim(xmin-50,xmax+50)
    pl.plot((xmin-50,xmax+50),(0.,0.),'k-',lw=3)
    if doY:
        pl.xlabel('Y Position (pixels)')
        pl.ylabel('Mean Y displacement (mas)')
    else:
        pl.xlabel('X Position (pixels)')
        pl.ylabel('Mean X displacement (mas)')
    pl.grid()
    
    return xpos, dx, err[useful]

# 2pcf within chips?

def rPlot(fits, cellSize, xy, maxErr=0.1, color='g', **kwargs):
    """
    Stack residuals along one or other axis given a fits wcs output table
    cellSize is the approximate box size (in pixels) in which detections will be binned
    maxErr is the biggest error in cell mean residual that will get plotted (in pixels)
    xy is (x,y) of origin for radius
    kwargs is a dict giving selection criteria for the plot as per select() function
    """

    t = select(fits, **kwargs)      # Get the data to be plotted

    # Now set up grid of cells
    x = t['xpix'] - xy[0]
    y = t['ypix'] - xy[1]
    r = np.sqrt(x*x+y*y)
    
    # Get radial and tangent components of displacement
    dr = (t['xrespix']*x + t['yrespix']*y)/r
    dt = (t['xrespix']*y - t['yrespix']*x)/r

    xmin = np.min(r)
    xmax = np.max(r)

    nx = int( np.floor( (xmax-xmin)/cellSize))
    xCell = (xmax - xmin) / nx

    # Decide which cell number each point is in
    index = np.zeros(len(t), dtype=int)
    np.clip( np.floor((r-xmin)/xCell), 0, nx-1, index)

    # Calculate weighted mean of the x and y residuals
    w = t['sigpix'] / (1-t['wtFrac'])
    w = np.where(w>0, 1/(w*w), 0.)
    sumrw = np.histogram(index, bins=nx, range=(0,nx), weights= dr*w)[0]
    sumtw = np.histogram(index, bins=nx, range=(0,nx), weights= dt*w)[0]
    sumw =  np.histogram(index, bins=nx, range=(0,nx), weights= w)[0]

    dr = sumrw / sumw
    dt = sumtw / sumw

    useful = np.logical_and(sumw!=0, sumw>1./(maxErr*maxErr))

    # Make an x and y position for each cell to be the center of its cell
    xpos = np.arange(nx,dtype=int)
    xpos = (xpos+0.5)*xCell + xmin

    dr = dr[useful]
    dt = dt[useful]
    xpos = xpos[useful]
    err = 1./np.sqrt(sumw[useful])

    xmin = 0.
    xmax = 4500.

    # Now make the plot
    pl.plot((xmin,xmax),(0.,0.),'k-',lw=3)
    #    pl.plot(xpos, dr * pixScale, color+'o-', label='Radial')
    #    pl.plot(xpos, dt * pixScale, color+'o', label='Azimuthal',alpha=0.3)
    pl.errorbar(xpos, dr*pixScale, yerr=err*pixScale, ecolor=color, color=color, marker='o',ms=5, ls='None')
    pl.xlim(xmin,xmax)
    pl.xlabel('Radius from ({:.0f},{:.0f}) (pixels)'.format(xy[0],xy[1]))
    pl.ylabel('Mean displacement (mas)')
    pl.grid()
    
    return xpos, dr*pixScale

def ringplot(fits, device, photo, astro, **kwargs):
    """
    Plot the ring template, the actual residuals vs radius,
    and the difference between them.
    Also the derivative of (filtered) difference plotted vs the
    photometric template from which it came.
    photo and astro are the dictionaries holding the photometric and astrometric templates.
    
    """
    pl.subplot(311)
    x0 = astro[device]['XCenter']
    y0 = astro[device]['YCenter']
    r, dr = rPlot(fits, cellSize=8, maxErr=10./pixScale, device=device, xy=(x0,y0), clip=False, **kwargs)
    #drast = (-pixScale)*np.array(astro[device]['Values'])
    drast = (pixScale)*np.array(astro[device]['Values'])
    rast = np.arange(len(drast)) * astro[device]['ArgStep'] + astro[device]['ArgStart']
    pl.plot(rast, drast, 'm-')
    pl.xlim(0,4500.)
    pl.ylim(-40,40)
    pl.title(device)
    # Plot residual and a spline fit to it
    pl.subplot(312)
    resid = dr - np.interp(r, rast, drast)
    pl.plot(r, resid, 'm-', label='Radial residual')
    pl.grid()
    pl.ylim(-20,20)
    xmin = 0
    xmax = 4500.
    pl.xlim(xmin,xmax)
    pl.plot((xmin,xmax),(0.,0.),'k-',lw=3)
    pl.xlabel("Radius (pixels)")
    pl.ylabel("Residual (mas)")
    # Median-filter residuals to get rid of crap points
    tmp = medfilt(resid,5)
    s = UnivariateSpline(r, tmp, s=len(r)*3.)
    pl.plot(r, s(r), 'c-', label='Spline fit')
    pl.legend(loc=3)
    # Plot residual to template after fit to template and spline
    xx = np.vstack((np.interp(r,rast,drast), s(r))).transpose()
    tmp = tmp + np.interp(r,rast,drast) # put template back into data
    c = np.linalg.lstsq(xx,tmp)[0]
    print(device, c)
    resid = dr - np.dot(xx,c)
    pl.subplot(313)
    pl.plot(r, resid, 'go-', label='Resid after spline spline')
    newtemplate = drast + c[1]*s(rast)/c[0]
    pl.plot(rast, newtemplate, 'c-', label='new template')
    pl.grid()
    pl.ylim(-10,10)
    pl.xlim(xmin,xmax)
    pl.plot((xmin,xmax),(0.,0.),'k-',lw=3)
    pl.xlabel("Radius (pixels)")
    pl.ylabel("Photometric deviation")
    pl.legend(loc=3, framealpha=0.3)
    astro[device]['Values'] = (newtemplate/pixScale).tolist()
    return 

# Each-star vectors?


def allRingPlots(filename):
    """
    Plot all CCDs ring data into a big pdf
    """
    import yaml
    photo = yaml.load(open('photorings3.yaml'))
    astro = yaml.load(open('astrorings3.yaml'))
    kk = list(decam.ccdnums.keys())
    kk.sort()
    kk.pop(23) # get rid of N30
    from matplotlib.backends.backend_pdf import PdfPages
    pp= PdfPages('rrings.pdf')
    for detpos in kk:
        pl.figure(figsize=(14.,8.))
        ringplot(filename,detpos,photo,astro)
        pp.savefig()
        pl.close()
    pp.close()
    return astro

def ringFactors(wcs):
    """ Plot the multiplicative factors being used for astrometric rings in
    all filters, ccds
    The input is the YAML-loaded wcs specification file.
    """
    data = {}
    bands = ('g','r','i','z','Y')
    for b in bands:
        data[b] = np.zeros(63,dtype=float)
        num = []
        val = []
        m = re.compile('^' + b + r'/([NS]\d+)/rings$')
        for k,v in wcs.items():
            if not re.match(m,k):
                continue
            detpos = re.match(m,k).group(1)
            val.append(-v['Parameter'])
            num.append(ccdnums[detpos])
        data[b][np.array(num)] = np.array(val)

    gr = (data['g'] + data['r'])/2.
    ymin = 0.0
    ymax = 1.1
    ccd = np.arange(1,63)
    pl.plot(ccd, np.clip(gr[1:],ymin,ymax), 'ks',label='gr mean')
    for k in 'rizY':
        pl.plot(ccd, np.clip(data[k][1:]/gr[1:],ymin,ymax), plotColors[k]+'o',
                label=k+'/mean',alpha=0.7)
    pl.xlabel('CCDNUM')
    pl.ylabel('Ring template scaling factor')
    pl.xlim(0,63)
    pl.ylim(ymin-0.01,ymax+0.01)
    pl.grid()
    pl.legend(loc=4,numpoints=1,framealpha=0.3)
    #    pl.savefig('ringscales.pdf')
    
    #f = open('ringscales.dat','w')
    #for i in range(len(ccd)):
    #print(ccd[i],gri[i],file=f)
    return

def lateralColor(astro, band, scale):
    """
    Plot the color solution derived for a given filter's CCDs.  Will plot displacement
    arrows at the corners and center of each CCD using the linear-per-ccd solutions.
    astro is an already-loaded astrometry YAML structure.
    scale is the size of a large arrow
    """
    x = []  # Input position and output location vectors 
    y = []
    dx = []
    dy = []
    m = re.compile('^' + band + r'/(.*)/color')
    for k,v in astro.items():
        if not m.match(k):
            continue
        coeffs = v['Function']['Coefficients']
        detpos = m.match(k).group(1)
        x0,x1,y0,y1 = decam.ccdBounds[detpos]
        xin = np.array([x0,x0,x1,x1,0.5*(x0+x1)])
        yin = np.array([y0,y1,y0,y1,0.5*(y0+y1)])
        xout = coeffs[0] + coeffs[1]*xin + coeffs[2]*yin
        yout = coeffs[3] + coeffs[4]*xin + coeffs[5]*yin
        x.append(np.array(xin))
        y.append(np.array(yin))
        dx.append(xout-xin)
        dy.append(yout-yin)
        pl.plot( (x0,x0,x1,x1,x0), (y0,y1,y1,y0,y0), 'k:') # Plot ccd outline
    x = np.concatenate(x)
    y = np.concatenate(y)
    dx = np.concatenate(dx)
    dy = np.concatenate(dy)
    dx -= np.mean(dx)
    dy -= np.mean(dy)
    dx *= 3600. * 1000.  #change to mas
    dy *= 3600. * 1000.  #change to mas
    q = pl.quiver(x,y,dx,dy,angles='xy',scale_units='xy',scale=scale/0.3)
    pl.xlim(-1.2, 1.2)
    pl.ylim(-1.2, 1.2)
    pl.title('Chromatic distortion for band ' + band)
    pl.quiverkey(q, 0.6,-1.0,scale,"{:.0f} mas / mag".format(scale),coordinates='data',labelpos='S')
    pl.gca().set_aspect('equal')
    pl.gca().xaxis.set_visible(False)
    pl.gca().yaxis.set_visible(False)
