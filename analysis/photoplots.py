# Plots of photometric fitting results
from __future__ import division,print_function
import numpy as np
import astropy.io.fits as pf
import pylab as pl
import re
from scipy.interpolate import UnivariateSpline

bands = { 'g':'g','r':'r','i':'c','z':'m','Y':'y'}
pixScale = 264.   # nominal mas per pixel

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
            for e in exprs.split(','):
                # Append character to RE to force match of full string
                self.relist.append(re.compile(e.strip() + '$'))
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
            
    def where(self, v):
        """
        Return array of indices of objects in the list v which match
        one of the test expressions.
        """
        
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
           clip=False, reserve=None, amp=None, wtFracMax=0.35):
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
    objTab = ff['PHOTOOUT'].data
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

    if clip is not None:
        use = np.logical_and(use, (objTab['CLIP']>0)==clip)
    if reserve is not None:
        use = np.logical_and(use, (objTab['RESERVE']>0)==reserve)
    if wtFracMax is not None:
        use = np.logical_and(use, objTab['WTFRAC']<=wtFracMax)

    ## print('Select',np.count_nonzero(use), '/',len(use))
                
    return objTab[use]

def rmsVsSigma(fitsfile, color='r',label=None, **kw):
    # Plot rms photometric residual vs input noise
    # Additional keywords are given to the selection
    
    tab = select(fitsfile, **kw)
    sig = tab['sigMag']
    resid = tab['magRes'] / np.sqrt((1 - tab['wtFrac']))
    
    # take out systematic error (parameter ???)
    sysError = 0.002
    sig = np.sqrt(sig*sig-sysError*sysError)

    dsig = 0.00025;
    s = []
    rms = []
    
    for sl in np.arange(0,0.01,dsig):
        use = np.logical_and(sig>=sl, sig<sl+dsig)
        if np.count_nonzero(use)>0:
            rms.append(np.std(resid[use]))
            s.append(sl+dsig*0.5)
    s = np.array(s)
    rms = np.array(rms)

    # Report fit to model of systematic + scaled internal scatter
    A = np.array( [s*s, np.ones(len(s))] )
    w = np.linalg.lstsq(A.T, rms*rms)[0]
    print("color " , color, " Fitted error multiplier ", np.sqrt(w[0]), "RMS in first bin: ", rms[0])
    
    pl.plot(s,rms,color+'o',label=label)
    s = np.arange(0,0.01,0.0001)
    pl.plot(s, np.sqrt(w[1]+w[0]*s*s), color+':')
    pl.xlabel('SExtractor magnitude measurement error')
    pl.ylabel('Magnitude error')
    pl.title('Solution accuracy for ' + fitsfile)
    pl.xlim(xmin=0.)
    pl.ylim(ymin=0.)
    lims = pl.xlim()
    pl.plot(lims,lims,'k--')
    pl.grid()

    return

def residInPixels(fitsfile, binpix=128, plotRMS=False, scaleLimit = 0.003,
                  maxerr=0.0005, pixCoords=False, decorate=True, **kw):
    noData = np.nan
    
    # Return a 2d array (image) containing weighted mean residual in pixelized areas.
    # Done in pixel coordinates if pixCoords=True (all devices will be stacked)
    # else done in exposure coordinates (projection about axis)
    # binpix is the size of bins in nominal DECam pixels
    # or if pixCoords=True it is simply the size of superpixels.
    # plotRMS=true to plot RMS of residuals instead of mean
    # decorate = False to omit labels, titles
    tab = select(fitsfile,**kw)
    if len(tab) < 100:
        return None
    resid = tab['magRes']/(1-tab['wtfrac'])
    sig = tab['sigMag']
    weight = 1./(sig*sig)

    # use exposure coordinates:
    if pixCoords:
        x = tab['xPix']
        y = tab['yPix']
    else:
        x = tab['xExpo']
        y = tab['yExpo']
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    
    # Make 2d arrays to sum up the mean and weight per pixel
    if pixCoords:
        xyBuffer = 50   # 50 pixels around edges
        cellSize = binpix
    else:
        xyBuffer = 0.05   # buffer in degrees
        cellSize = binpix * pixScale / (1000.*3600.)

    xsize = int(np.ceil( (xmax-xmin) / cellSize))
    ysize = int(np.ceil( (ymax-ymin) / cellSize))

    # Figure out the image pixel for each residual point
    ix = np.array( np.floor( (x-xmin) / cellSize), dtype=int)
    ix = np.clip(ix, 0, xsize-1)
    iy = np.array( np.floor( (y-ymin) / cellSize), dtype=int)
    iy = np.clip(iy, 0, ysize-1)
    index = iy*xsize + ix

    if plotRMS:
        sumxw = np.histogram(index, bins=xsize*ysize, range=(-0.5,xsize*ysize+0.5),
                            weights=resid*resid)[0].reshape(ysize,xsize)
        sumw = np.histogram(index, bins=xsize*ysize,
                            range=(-0.5,xsize*ysize+0.5))[0].reshape(ysize,xsize)
        # Smallest number of points for useful RMS:
        minWeight = 20
        # Value to put in where weight is below threshold:
        sumxw = np.where( sumw > minWeight, sumxw / sumw, noData)
        sumxw = np.sqrt(sumxw);
        rms = 0.
        pl.imshow(sumxw, aspect='equal', vmin=0., vmax=scaleLimit,
                  origin='lower',cmap='Spectral',
                  interpolation='nearest', extent=(xmin,xmax,ymin,ymax))
    else:
        sumxw = np.histogram(index, bins=xsize*ysize, range=(-0.5,xsize*ysize+0.5),
                            weights=weight*resid)[0].reshape(ysize,xsize)
        sumw = np.histogram(index, bins=xsize*ysize, range=(-0.5,xsize*ysize+0.5),
                            weights=weight)[0].reshape(ysize,xsize)
        # Smallest weight that we'd want to plot a point for
        minWeight = maxerr**(-2)
        # Value to put in where weight is below threshold:
        sumxw = np.where( sumw > minWeight, sumxw / sumw, noData)
        sumw = np.where( sumw > minWeight, 1./sumw, noData)
        rms = np.std(sumxw[sumw>0.])
        print('RMS, noise:', rms, np.sqrt(np.mean(sumw[sumw>0.])))
        pl.imshow(sumxw, aspect='equal', vmin=-scaleLimit, vmax=scaleLimit, origin='lower',
                 interpolation='nearest', extent=(xmin,xmax,ymin,ymax),cmap='Spectral')

    pl.xlim(xmin-xyBuffer,xmax+xyBuffer)
    pl.ylim(ymin-xyBuffer,ymax+xyBuffer)
    if decorate:
        pl.colorbar()
        if pixCoords:
            pl.xlabel("X (pixels)")
            pl.ylabel("Y (pixels)")
        else:
            pl.xlabel("HA (degrees)")
            pl.ylabel("Dec (degrees)")

        if type(fitsfile)==str:
            pl.title('Binned residuals for ' + fitsfile)
    return rms,sumxw

def residArray(base='allsf',binpix=512,maxerr=0.01,scaleLimit=0.006):
    """
    Plot a big array of residuals to collective fit in each filter/epoch combination
    *** Needs work with formatting, filenames are hard-wired.
    """
    sfDates = ['20121120','20121223','20130221','20130829','20131115','20140118',
                '20140807','20141105','20150204','20150926','20160209',
                '20160223','20160816','20161117','20170111','20170214']
    
    f, axes = pl.subplots(len(bands), len(sfDates), sharex=True, sharey=True,
                          figsize=(16,6))
    f.subplots_adjust(hspace=0, wspace=0,
                      left=0.02, right=0.94,
                      bottom=0.02, top=0.98)
    for j,b in enumerate(('g','r','i','z','Y')):
        ff = pf.open(base+'.'+b+'.photo.fits')
        for i,date in enumerate(sfDates):
            if date=='20150204' and b=='z':
                continue  ###
            pl.sca(axes[j,i]);
            img = residInPixels(ff,binpix=binpix,maxerr=maxerr,
                                scaleLimit=scaleLimit,
                                clip=False, instrument = b + date, decorate=False)[1]
            if j==0:
                axes[j,i].text(-1.,0.8,date,fontsize=10)
            if i==0:
                axes[j,i].text(-1.,-0.8,b,fontsize=14)
            axes[j,i].tick_params(axis='both',top=False,bottom=False,
                                      left=False, right=False,
                                      labeltop=False,labelbottom=False,
                                      labelleft=False,labelright=False)
        ff.close()
    # Add colorbar
    cax = f.add_axes([0.96,0.05,0.02,0.9])
    pl.colorbar(cax=cax)
    pl.setp([a.get_xticklabels() for a in f.axes], visible=False)
    pl.setp([a.get_yticklabels() for a in f.axes], visible=False)
    f.savefig('array.pdf',dpi=600,transparent=True)
    return

def residVsR(fitsname, yaml,device, rbin=4, **kw):
    # Plot mean photometric residual vs radius from ring center, and
    # both halves of template read from yaml file
    """
    fitsname is a fits filename or object
    yaml is the decoded yaml photorings template file
    device is DETPOS string used in YAML lookup
    rbin is size of radial bins for residuals
    kw are selection criteria - device will be added.
    """

    x0 = yaml[device]['XCenter']
    y0 = yaml[device]['YCenter']
    r0 = yaml[device]['ArgStart']
    dr = yaml[device]['ArgStep']
    v = yaml[device]['Values']
    pl.plot(np.arange(len(v))*dr+r0, v, 'r-',label='Template')

    # Get the observational residuals
    kw['device']=device
    tab = select(fitsname,**kw)

    sig = tab['sigMag']
    x = tab['xPix']-x0
    y = tab['yPix']-y0
    resid = tab['magRes'] / (1 - tab['wtFrac'])
    sig = tab['sigMag']
    weight = 1./(sig*sig)

    r = np.sqrt(x*x+y*y)
    index = np.floor( (r - np.min(r))/rbin)
    nbins = np.max(index)
    sumxw = np.histogram(index, bins=nbins, range=(-0.5,nbins+0.5),
                            weights=weight*resid)[0]
    sumw = np.histogram(index, bins=nbins, range=(-0.5,nbins+0.5),
                            weights=weight)[0]
    rr = np.arange(len(sumxw))*rbin + np.min(r) + 0.5*rbin
    y = sumxw / sumw
    sig = np.sqrt(1./sumw)
    pl.errorbar(rr,y, yerr=sig, marker='o',markerfacecolor='green',linestyle='None',
                label='Fit residuals', ecolor='green',markersize=5,alpha=0.5)
    pl.xlabel('Ring radius (pixels)')
    pl.ylabel('Mean residual')
    pl.title(device + ' ring residuals')
    pl.ylim(-0.015,0.015)
    lims = pl.xlim()
    pl.axhline(0.,0.,1.,lw=2,color='k')
    pl.grid()
    pl.legend(loc=3)
    return

def resid1d(fitsfile, npix=4, plotRMS=False, scaleLimit = 3., maxerr=2., doY = False, **kw):
    
    # Make 1d plot of mean residual vs X pixel coordinate (or Y, if doY=True).
    # npix is the bin width.
    # scalelimit is the +- limit of the plot (in mmag).
    # plotRMS=True to plot std dev of photometry errors instead of mean.

    tab = select(fitsfile,**kw)
    resid = 1000.*tab['magRes']/(1-tab['wtfrac'])
    sig = 1000.*tab['sigMag']
    weight = 1./(sig*sig)

    # use exposure coordinates:
    if doY:
        x = tab['yPix']
        xmin = 0.
        xmax = 4096.
    else:
        x = tab['xPix']
        xmin = 0.
        xmax = 2048.
    
    # Figure out the bin for each residual point
    ix = np.array( np.floor( (x-xmin) / npix), dtype=int)
    nx = int(np.floor(xmax-xmin)/npix)+1
    ix = np.clip(ix, 0, nx-1)
    xbin = (np.arange(nx)+0.5)*npix

    if plotRMS:
        sumxw = np.histogram(ix, bins=nx, range=(-0.5,nx+0.5),
                            weights=resid*resid)[0]
        sumw = np.histogram(ix, bins=nx, range=(-0.5,nx+0.5))[0]
        # Smallest number of points for useful RMS:
        minWeight = 20
        # Value to put in where weight is below threshold:
        use = sumw > minWeight
        sumxw = (sumxw/sumw)[use]
        sumxw = np.sqrt(sumxw);
        xbin = xbin[use]
        pl.plot(xbin, sumxw, 'ro')
        pl.ylim(0,scaleLimit)
        pl.xlim(xmin-1, xmax+1)
        pl.grid()
    else:
        sumxw = np.histogram(ix, bins=nx, range=(-0.5,nx+0.5),
                            weights=weight*resid)[0]
        sumw = np.histogram(ix, bins=nx, range=(-0.5,nx+0.5),
                            weights=weight)[0]
        # Smallest weight that we'd want to plot a point for
        minWeight = maxerr**(-2)
        # Value to put in where weight is below threshold:
        use = sumw > minWeight
        sumxw = (sumxw/sumw)[use]
        sumw = 1./sumw[use]
        xbin = xbin[use]
        print('RMS of useful pixels: ', np.std(sumxw))
        print('Expected noise: ', np.sqrt(np.mean(sumw)))
        ##pl.errorbar(xbin, sumxw, yerr=1./np.sqrt(sumw), fmt='ro',ecolor='r')
        pl.plot(xbin, sumxw,'ro')
        pl.ylim(-scaleLimit,scaleLimit)
        pl.xlim(xmin-1, xmax+1)
        pl.grid()
        
    if doY:
        pl.xlabel("Y (pixels)")
    else:
        pl.xlabel("X (pixels)")
    if plotRMS:
        pl.ylabel("RMS (mmag)")
    else:
        pl.ylabel("Mean (mmag)")

    pl.title('Binned residuals for ' + fitsfile)
    return

def excessVariance(fitsfile, color='r',label=None, **kw):
    # Plot rms photometric residual vs input noise
    # Additional keywords are given to the selection
    
    tab = select(fitsfile, **kw)
    sig = tab['sigMag']
    resid = tab['magRes'] / np.sqrt((1 - tab['wtFrac']))
    
    # take out systematic error (parameter ???)
    sysError = 0.002
    sig = np.sqrt(sig*sig-sysError*sysError)

    dsig = 0.00025;
    s = []
    rms = []
    
    for sl in np.arange(0,0.01,dsig):
        use = np.logical_and(sig>=sl, sig<sl+dsig)
        if np.count_nonzero(use)>0:
            rms.append(np.std(resid[use]))
            s.append(sl+dsig*0.5)
    s = np.array(s)
    rms = np.array(rms)

    # Report fit to model of systematic + scaled internal scatter
    A = np.array( [s*s, np.ones(len(s))] )
    w = np.linalg.lstsq(A.T, rms*rms)[0]
    print("color " , color, " Fitted error multiplier ", np.sqrt(w[0]), "RMS in first bin: ", rms[0])
    
    pl.plot(s,rms,color+'o',label=label)
    s = np.arange(0,0.01,0.0001)
    pl.plot(s, np.sqrt(w[1]+w[0]*s*s), color+':')
    pl.xlabel('SExtractor magnitude measurement error')
    pl.ylabel('Magnitude error')
    pl.title('Solution accuracy for ' + fitsfile)
    pl.xlim(xmin=0.)
    pl.ylim(ymin=0.)
    lims = pl.xlim()
    pl.plot(lims,lims,'k--')
    pl.grid()

    return

###############################
### Old plots:
###############################

def ampVariation(table, deviceNumber, color='r'):
    # Plot mean stellar residual vs exposure number for two amps of chosen device
    good = np.logical_and(table['clip']==0, table['sigMag']<0.005)
    meanA = []
    meanB = []
    # ???
    nExposures = 90
    for i in range(nExposures):
        use = np.logical_and(good, table['Extension']==deviceNumber-1+i*nExposures)
        useA = np.logical_and(use, table['xPix']<1024)
        useB = np.logical_and(use, table['xPix']>1024)
        if np.any(useA): meanA.append(np.mean(table['magRes'][useA]))
        if np.any(useB): meanB.append(np.mean(table['magRes'][useB]))
    pl.plot(meanA, color+'o', label=str(deviceNumber)+'A')
    pl.plot(meanB, color+'^', label=str(deviceNumber)+'B')
    pl.xlabel("Exposure")
    pl.ylabel("Mean photometric residual (mags)")
    pl.legend()
    return np.std(meanA), np.std(meanB)

def colorColor(magFile, m1='g', m2='r', m3='i', m4='z', maxMagError=0.01,dc = 0.01):
    table = pf.getdata(magFile,1)
    use = table[m1+'_err'] <= maxMagError
    use = np.logical_and(use, table[m2+'_err'] <= maxMagError)
    use = np.logical_and(use, table[m3+'_err'] <= maxMagError)
    use = np.logical_and(use, table[m4+'_err'] <= maxMagError)
    use = np.logical_and(use, table[m1] > -10)
    use = np.logical_and(use, table[m2] > -10)
    use = np.logical_and(use, table[m3] > -10)
    use = np.logical_and(use, table[m4] > -10)
    #    pl.plot( (table[m1]-table[m2])[use], (table[m3]-table[m4])[use], 'r.', alpha=0.3)
    pl.xlabel( m1 + '-' + m2)
    pl.ylabel( m3 + '-' + m4)
    pl.title('Color-color diagram from ' + magFile)
    #    pl.grid()
    

    xdat = (table[m1]-table[m2])[use]
    ydat = (table[m3]-table[m4])[use]

    xmin = np.percentile(xdat, 00.5)
    xmax = np.percentile(xdat, 99.5)
    ymin = np.percentile(ydat, 00.5)
    ymax = np.percentile(ydat, 99.5)
    dx = xmax - xmin;
    xmin -= 0.05*dx
    xmax += 0.05*dx
    dy = ymax - ymin;
    ymin -= 0.05*dy
    ymax += 0.05*dy
    xyrange = [[xmin,xmax],[ymin,ymax]]
    
    bins = [int((xmax-xmin)/dc), int((ymax-ymin)/dc)] # number of bins
    thresh = 2  #density threshold

    # histogram the data
    hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)

    #select points within the histogram
    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    hhsub = hh[posx[ind] - 1, posy[ind] - 1] # values of the histogram where the points are
    xdat1 = xdat[ind][hhsub < thresh] # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = np.nan # fill the areas with low density by NaNs

    pl.imshow(np.flipud(hh.T),cmap='jet',extent=np.array(xyrange).flatten(), interpolation='nearest')
    pl.colorbar()   
    pl.plot(xdat1, ydat1, 'k.', markersize=1)

def imageChange(fits1, fits2, lowerValid=-0.5, upperValid = 0.5, colorRange=None, quadratic=False):
    # plot difference between two images, after subtraction of best-fit plane from difference
    # and also create histogram & report RMS, min, max difference.
    # Input values outside of valid range are ignored (in both images)
    # +-colorRange is range for the color scale. If None, use +-4*RMS of valid points

    # Assuming primary extensions for both images ???
    img1 = pf.getdata(fits1, 0)
    img2 = pf.getdata(fits2, 0)
    if img1.shape != img2.shape:
        print("Shape mismatch: ", img1.shape, img2.shape)
        return (0.,0.,0.)
    valid = np.logical_and(img1 >= lowerValid, img1 <= upperValid)
    valid = np.logical_and(valid, img2 >= lowerValid)
    valid = np.logical_and(valid, img2 <= upperValid)
    diff = img1 - img2
    rowx = np.arange(img1.shape[0], dtype=float).reshape(img1.shape[0],1)
    coly = np.arange(img1.shape[1], dtype=float)
    (imgx, imgy) = np.broadcast_arrays(rowx, coly)
    imgx = np.array(imgx)
    imgy = np.array(imgy)
    
    data = diff.flatten()[valid.flatten()]
    x = imgx.flatten()[valid.flatten()]
    y = imgy.flatten()[valid.flatten()]
    meanx = np.mean(x)
    meany = np.mean(y)
    x -= meanx
    y -= meany
    imgx -= meanx
    imgy -= meany
    
    # Fit linear function to the data and subtract it
    ## ?? A = np.array( [x, y, np.ones(len(data))] )
    if quadratic:
        A = np.array( [np.ones(len(data)), x, y, x*x, x*y, y*y] )
    else:
        A = np.array( [np.ones(len(data)), x, y] )
    w = np.linalg.lstsq(A.T,data)[0]
    print(A.shape, w)
    data -= np.dot(w, A)
    rms = np.std(data)
    if (colorRange==None): colorRange = 4*rms
    if quadratic:
        diff -= (w[0] + imgx * w[1] + imgy * w[2] + imgx*imgx*w[3] + imgx*imgy*w[4] + imgy*imgy*w[5])
    else:
        diff -= (w[0] + imgx * w[1] + imgy * w[2])
    diff[~valid] = np.nan
    pl.imshow(diff, interpolation='nearest', aspect='equal', origin='lower',
              vmin=-colorRange, vmax=colorRange, extent=(0,diff.shape[1],0,diff.shape[0]))

    return rms, diff

def diffgrizY(base1, base2, colorRange=None):
    # big plot with difference images and histograms for g, r, i, z, Y differences
    bands = ['g', 'r', 'i', 'z', 'Y']
    f,ax = pl.subplots(5,2, sharex='col')
    for i in range(len(bands)):
        pl.axes(ax[i,1])
        rms,data = imageChange(base1+bands[i]+".fits", base2+bands[i]+".fits",
                               colorRange = colorRange)
        rms = 1000*rms
        pl.tick_params(labelbottom=False, labelleft = False)
        pl.axes(ax[i,0])
        pl.hist(data.ravel(), bins=100, range=(-colorRange, colorRange))
        pl.xlim(-colorRange, colorRange)
        pl.text(0.1,0.9,bands[i], transform=ax[i,0].transAxes)
        pl.text(0.1,0.8,'RMS= {0:.1f} mmag'.format(rms), transform=ax[i,0].transAxes)
        pl.tick_params(labelleft=False)
        
    pl.subplots_adjust(left=0.05, right=0.82, hspace=0.02, bottom=0.05, top=0.95)
    cax = pl.axes( [0.85, 0.3, 0.05, 0.4])
    pl.colorbar(cax=cax)
    pl.axes( [0.2, 0.96, 0.47, 0.], frameon=False)
    pl.tick_params(labelbottom='off',labelleft='off',bottom='off',top='off',left='off',right='off')
    pl.title(base1 + " - " + base2)

def pixelPhase(fitsfile, device='', nphase=20, plotRMS=False, maxSigma=0.006,
               exposure=None):
    # Plot mean residual vs pixel phase.
    # nphase is the number of bins in each dimension of pixel phase
    f = pf.open(fitsfile)
    data = f['PhotoOut'].data
    use = data['Clip']==0
    use = np.logical_and(use, data['sigMag']<maxSigma)
    if (exposure != None):
        # select only points from a chosen exposure number
        exposureNumbers = f['Extensions'].data['ExposureNumber']
        use = np.logical_and(use, exposureNumbers[data['Extension']]==exposure)
    resid = data['magRes'][use]
    sig = data['sigMag'][use]
    weight = 1./(sig*sig)
    xPix = data['xPix'][use]
    yPix = data['yPix'][use]
    xPix -= np.floor(xPix)
    yPix -= np.floor(yPix)
    extn = data['Extension'][use]

    deviceNumber = -1
    if (len(device)>0):
        deviceName = device
        
        # Look for the specified device in the instrument extension
        insttab = f['Instrument'].data
        deviceNames = insttab['Name']
        for i in range(len(deviceNames)):
            if (device.lower() == ''.join(deviceNames[i]).lower()):
                deviceNumber = i;
                break
        if (deviceNumber < 0):
            print('Did not find device named ' + device)
            return

        exttab = f['Extensions'].data
        useExtensions = []
        for i in range(len(exttab)):
            if exttab['DeviceNumber'][i]==deviceNumber:
                useExtensions.append(i)

        #Now select points from one of the extensions for this device.
        extn = data['Extension'][use]
        use = np.in1d(extn, useExtensions)
        print('points in: ', len(use), ' chosen ', np.count_nonzero(use))
        xPix = xPix[use]
        yPix = yPix[use]
        resid = resid[use]
        weight = weight[use]
        extn = extn[use]
        
    else:
        deviceName = 'All'
    

    # Make 2d arrays to sum up the mean and weight per pixel
    pixelScale = 1./nphase
    sumxw = np.zeros( (nphase,nphase), dtype=float)
    sumw = np.zeros( (nphase,nphase), dtype=float)

    # Figure out the image pixel for each residual point
    ix = np.array( np.floor(xPix*nphase), dtype=int)
    iy = np.array( np.floor(yPix*nphase), dtype=int)
    
    if plotRMS:
        for i in range(len(xPix)):
            sumw[ iy[i], ix[i]] += 1.
            sumxw[ iy[i], ix[i]] += resid[i]*resid[i]
    else:
        for i in range(len(xPix)):
            sumw[ iy[i], ix[i]] += weight[i]
            sumxw[ iy[i], ix[i]] += weight[i]*resid[i]

    # Smallest weight that we'd want to plot a point for
    minWeight = 0.003**(-2)
    # Value to put in where weight is below threshold:
    noData = np.nan

    if plotRMS:
        sumxw = np.where( sumw > 10, sumxw/sumw, noData)
        sumxw = np.sqrt(sumxw);
        scaleLimit = 0.020
        pl.imshow(sumxw, aspect='equal', vmin=0., vmax=scaleLimit, origin='lower',
                  interpolation='nearest', extent=(0.,1.,0.,1.))
    else:
        sumxw = np.where( sumw > minWeight, sumxw / sumw, noData)
        sumw = np.where( sumw > minWeight, 1./sumw, noData)
        print('RMS of useful pixels: ', np.std(sumxw[sumw>0.]))
        print('Expected noise: ', np.sqrt(np.mean(sumw[sumw>0.])))
        scaleLimit = 0.002
        pl.imshow(sumxw, aspect='equal', vmin=-scaleLimit, vmax=scaleLimit, origin='lower',
                  interpolation='nearest', extent=(0.,1.,0.,1.))

    pl.colorbar()
    pl.xlabel("X pixel phase")
    pl.ylabel("Y pixel phase")

    pl.title('Binned residuals for device ' + deviceName + ' results of ' + fitsfile)
    return

def grizy(base, flux=True):
    for i in ['g', 'r', 'i', 'z', 'Y']:
        rmsVsMag(base+"."+i+".photo.fits", color = bands[i], label=i, flux=flux)
    if (flux): pl.legend(loc=1)
    else: pl.legend(loc=2)
    return

def fixRings(fits, rings, device, fudge=0.9, rbin=4, **kw):
    """
    Plot ring residual and the spline fit to it.  Then alter the template in-place
    by subtracting the residual spline (divided by the fudge factor)
    """
    
    x0 = rings[device]['XCenter']
    y0 = rings[device]['YCenter']
    r0 = rings[device]['ArgStart']
    dr = rings[device]['ArgStep']
    v = rings[device]['Values']
    rtemplate = np.arange(len(v))*dr+r0
    pl.plot(rtemplate, v, 'r-',label='Template')

    # Get the observational residuals
    kw['device']=device
    tab = select(fits,**kw)

    sig = tab['sigMag']
    x = tab['xPix']-x0
    y = tab['yPix']-y0
    resid = tab['magRes'] / (1 - tab['wtFrac'])
    sig = tab['sigMag']
    weight = 1./(sig*sig)

    r = np.sqrt(x*x+y*y)
    index = np.floor( (r - np.min(r))/rbin)
    nbins = np.max(index)
    sumxw = np.histogram(index, bins=nbins, range=(-0.5,nbins+0.5),
                            weights=weight*resid)[0]
    sumw = np.histogram(index, bins=nbins, range=(-0.5,nbins+0.5),
                            weights=weight)[0]
    rr = np.arange(len(sumxw))*rbin + np.min(r) + 0.5*rbin
    y = sumxw / sumw
    sig = np.sqrt(1./sumw)
    pl.errorbar(rr,y, yerr=sig, marker='o',markerfacecolor='green',linestyle='None',
                label='Fit residuals', ecolor='green',markersize=5,alpha=0.5)

    use = sumw>0.
    s = UnivariateSpline(rr[use],y[use],w=sumw[use], s=np.sum(sumw[use])*2)
    pl.plot(rr,s(rr),'k--',lw=2, label='Spline')

    pl.xlabel('Ring radius (pixels)')
    pl.ylabel('Mean residual')
    pl.title(device + ' ring residuals')
    pl.ylim(-0.015,0.015)
    lims = pl.xlim()
    pl.axhline(0.,0.,1.,lw=2,color='k')
    pl.grid()

    # Adjust ring template to remove residuals
    ring2 = v - s(rtemplate) / fudge
    # Truncate to five decimal places to make ASCII smaller
    ring2 = 1e-5 * np.rint(ring2*1e5);
    rings[device]['Values'] = ring2.tolist()
    v = rings[device]['Values']
    ##pl.plot(rtemplate, v, 'g-',label='New template')
    
    pl.legend(loc=3)
    return
