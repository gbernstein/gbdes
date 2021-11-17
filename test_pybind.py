
from lsst.daf.butler import Butler
from lsst.daf.persistence import Butler as ButlerGen2
import wcsfit
import wcsfit_test
import numpy as np
import sip_tpv
from astropy.io import fits

#butler = Butler('/repo/main', collections="HSC/runs/RC2/w_2021_30/DM-31182")
#butler = Butler('/repo/main', collections="HSC/runs/RC2/w_2021_34/DM-31524")
#butler = ButlerGen2('/datasets/hsc/repo/rerun/DM-23243/SFM/DEEP')
refactor_db_file = "/home/csaunder/stack_projects/gbdes_tests/refactor_tests/refactor.fof"
refactor_db = fits.open(refactor_db_file)
#import ipdb; ipdb.set_trace()

"""
print('test package:', wcsfit_test.__file__)
inputTables = '/home/csaunder/stack_projects/gbdes_tests/cosmos_pdr2_field/match_newgaia.fof'
"""
# GBDES defaults:
reserveFraction = 0.05
randomNumberSeed = 0
skipFile = ''
clipThresh = 5
maxError = 100
sysError = 2
referenceSysError = 2
freePM = True
pmEpoch = 2015.5
parallaxPrior = 10
pmPrior = 100
minMatches = 2
minFitExposures = 200
clipEntireMatch = False
chisqTolerance = 0.001
divideInPlace = False
purgeOutput = False
inputMaps = '/home/csaunder/stack_projects/gbdes_tests/polyExposure.astro'
fixMaps = ''
useInstruments = '.*'
skipExposures = ''
outCatalog = 'wcscat.fits'
outWcs = 'wcsfit.wcs'
starCatalog = 'starcat.fits'
colorExposures = ''
minColor = -10
maxColor = 10
verbose = 2

pvrange, tpvx, tpvy = sip_tpv.pvsiputils.sym_tpvexprs()


def pv_matrix_from_sip(sipx, sipy):
    # This function is a combination of sip_tpv.add_pv_keywords and ReadPV in gbdes/src/subs/TPVMap.cpp
    coeffs_size = max(pvrange) + 1

    pv1_vals = np.zeros(coeffs_size)
    pv2_vals = np.zeros(coeffs_size)
    for p in pvrange:
        # makes PV1_{p}, PV2_{p}
        pv1_vals[p] = float(sip_tpv.pvsiputils.calcpv(pvrange, 1, p, sipx, sipy, tpvx, tpvy).evalf())
        pv2_vals[p] = float(sip_tpv.pvsiputils.calcpv(pvrange, 2, p, sipx, sipy, tpvx, tpvy).evalf())
    
    coeffs_size = np.flatnonzero(pv1_vals != 0).max() + 1
    order = 1
    while ((order + 1) * (order + 2) / 2) < coeffs_size:
        order += 1

    pv1 = np.zeros((order + 1, order + 1))
    pv2 = np.zeros((order + 1, order + 1))
    i = 0
    j = 0
    for k in range(coeffs_size):
        pv1[i, j] = pv1_vals[k]
        pv2[i, j] = pv2_vals[k]
        if i == 0:
            i = j + 1
            j = 0
        else:
            i -= 1
            j += 1

    if 1 not in pvrange:
        pv1[1, 0] = 1.0
        pv2[1, 0] = 1.0

    pv1 = pv1.transpose()
    pv2 = pv2.transpose()
    return pv1, pv2

POLYSTEP = 1.0/3600.0
def convert_sip_to_tpv(butler_wcs, name=""):
    # This function is an adaptation of readTPV in gbdes/src/subs/TPVMap.cpp
    
    # Format copied from sip_tpv.sip_pv function
    try:
        cd = butler_wcs.getCdMatrix()
        fits_metadata = butler_wcs.getFitsMetadata()
    except: 
        fits_metadata = butler_wcs
        cd = np.matrix([[butler_wcs.get('CD1_1'), butler_wcs.get('CD1_2')],
                        [butler_wcs.get('CD2_1'), butler_wcs.get('CD2_2')]])
    

    
    if not ((fits_metadata.get('CTYPE1') == 'RA---TAN-SIP') and
            (fits_metadata.get('CTYPE2') == 'DEC--TAN-SIP')):
        raise ValueError(f'CTYPES {fits_metadata.get("CTYPE1")} and {fits_metadata.get("CTYPE2")}'
                         'do not match SIP convention')
    
    plin = np.zeros(6)
    crpix1 = fits_metadata.get('CRPIX1')
    crpix2 = fits_metadata.get('CRPIX2')
    plin[0] = - cd[0, 0] * crpix1 - cd[0, 1] * crpix2
    plin[1] = cd[0, 0]
    plin[2] = cd[0, 1]
    plin[3] = -cd[1, 0] * crpix1 - cd[1, 1] * crpix2
    plin[4] = cd[1, 0]
    plin[5] = cd[1, 1]
    linearMap = wcsfit.LinearMap(plin)

    pole = wcsfit.SphericalICRS(fits_metadata.get('CRVAL1') * wcsfit.DEGREE,
                                fits_metadata.get('CRVAL2') * wcsfit.DEGREE)
    orientIn = wcsfit.Orientation()
    orientIn.set(pole, 0.0)
    tp = wcsfit.Gnomonic(orientIn.getPole(), orientIn)

    a_order = int(fits_metadata.get('A_ORDER', 0))
    b_order = int(fits_metadata.get('B_ORDER', 0))
    ac = np.matrix(np.zeros((a_order+1, a_order+1), dtype=np.float64))
    bc = np.matrix(np.zeros((b_order+1, b_order+1), dtype=np.float64))
    for m in range(a_order+1):
        for n in range(0, a_order+1-m):
            ac[m, n] = fits_metadata.get('A_%d_%d' % (m, n), 0.0)
    for m in range(b_order+1):
        for n in range(0, b_order+1-m):
            bc[m, n] = fits_metadata.get('B_%d_%d' % (m, n), 0.0)
    sipx, sipy = sip_tpv.pvsiputils.real_sipexprs(cd, ac, bc)
    pv1, pv2 = pv_matrix_from_sip(sipx, sipy)
    p1 = wcsfit.Poly2d(pv1)
    p2 = wcsfit.Poly2d(pv2)

    pmlist = [linearMap]
    polyName = name + "_pv"
    if p1.nCoeffs() > 1:
        # Obtained a valid x polynomial.
        if p2.nCoeffs() > 1:
            # Also have valid y polynomial.  Make the map:
            pv = wcsfit.PolyMap(p1, p2, polyName, wcsfit.Bounds(-1.0, 1.0, -1.0, 1.0), POLYSTEP)
            
        else:
            # Did not find PV2's.  Install default:
            p2Identity = wcsfit.Poly2d(1)
            coeffs = p2Identity.getC()
            coeffs[p2Identity.vectorIndex(0, 1)] = 1.0
            p2Identity.setC(coeffs)
            pv = wcsfit.PolyMap(p1, p2Identity, polyName, wcsfit.Bounds(-1.0, 1.0, -1.0, 1.0), POLYSTEP)
        pmlist.append(pv)
    else:
        # Did not obtain any PV1's.  If there are PV2's, install identity map for x coeff:
        if p2.nCoeffs() > 1:
            p1Identity = wcsfit.Poly2d(1)
            coeffs = p1Identity.getC()
            coeffs[p1Identity.vectorIndex(1, 0)] = 1.0
            p1Identity.setC(coeffs)
            pv = wcsfit.PolyMap(p1Identity, p2, polyName, wcsfit.Bounds(-1.0, 1.0, -1.0, 1.0), POLYSTEP)
            pmlist.append(pv)

    #Create a SubMap that owns a duplicate of these one or two maps (no sharing):
    sm = wcsfit.SubMap(pmlist, name, False)
    
    # Return the Wcs, which again makes its own copy of the maps (no sharing):
    return wcsfit.Wcs(sm, tp, name, wcsfit.DEGREE, False)


def readExposures_wcsfof(refs):
    wcsList = []
    for v, visitSummaryRef in enumerate(refs):
        visitSummary = butler.get(visitSummaryRef)
        
        for row in visitSummary[:2]:
            print(f'Processing visit {row["visit"]}, detector {row["id"]}')
            calexp = butler.get('calexp', visit=row['visit'], detector=row['id'])
            butler_wcs = calexp.getWcs()
            wcs = convert_sip_to_tpv(butler_wcs)
            wcsList.append(wcs)
    return wcsList

def readExposuresExtensions(refs, input_yaml, fieldNumber=0, instrumentNumber=0, isReference=False):
    # TODO: add ExposureColorPriorities, ColorExtension / decide if necessary
    
    extensions = []
    exposures = []
    
    for v, visitSummaryRef in enumerate(refs):
        print(visitSummaryRef)        
        visitSummary = butler.get(visitSummaryRef)
        #visit = int(visitSummaryRef['NAME'][-5:])
        #print(visit)
        #visitSummary = butler.get('visitSummary', visit=int(visitSummaryRef['NAME'][-5:]))
        visInfo = visitSummary[0].getVisitInfo()
        
        # TODO: ra dec I think should be in radians based in line 482 in gbdes/src/FitSubroutines.cpp
        # -- need to double check
        ra = visInfo.getBoresightRaDec().getRa().asRadians()
        dec = visInfo.getBoresightRaDec().getDec().asRadians()
        gn = wcsfit.Gnomonic(wcsfit.Orientation(wcsfit.SphericalICRS(ra, dec)))
        # TODO: probably want a different id/name in future:
        exposure = wcsfit.Exposure(str(visInfo.getId()), gn)
        exposure.field = fieldNumber
        exposure.instrument = instrumentNumber
        airmass = visInfo.boresightAirmass
        exptime = visInfo.exposureTime
        exposure.airmass = airmass
        exposure.exptime = exptime
        exposure.mjd = visInfo.date.get(visInfo.date.DateSystem.MJD)
        if False:
            exposure.pmEpoch
        if False:
            exposure.apcorr
        # TODO: add astrometric covariance
        exposures.append(exposure)

        for row in visitSummary[:2]:
            extension = wcsfit.Extension()
            extension.exposure = v
            extension.device = row['id']
            extension.airmass = airmass
            extension.magshift = 2.5 * np.log10(exptime)
            if False:
                extension.apcorr
            # Add WCS here:
            calexpWCS = butler.get('calexp.wcs', visit=row['visit'], detector=row['id'])
            #butler_wcs = calexp.getWcs()
            wcs = convert_sip_to_tpv(calexpWCS)
            extension.startWcs = wcs
            extension.wcsName = f"{exposure.name}/{extension.device}"
            extension.mapName = f"{extension.wcsName}/base"
            print(extension.mapName)
            extensions.append(extension)
    return exposures, extensions


def readExposuresExtensions_fits(hdus, input_yaml, instruments, fieldNumber=0, instrumentNumber=0, isReference=False):
    # TODO: add ExposureColorPriorities, ColorExtension / decide if necessary
    
    extensions = []
    exposures = []
    wcss = []
    visits = []

    for v, visitSummary in enumerate(hdus[2].data):
        print(visitSummary)
        
        visit = int(visitSummary['NAME'][-5:])
        visits.append(visit)
        name = visitSummary['NAME']
        print(name, visit)
        
        # TODO: ra dec I think should be in radians based in line 482 in gbdes/src/FitSubroutines.cpp
        # -- need to double check
        ra = visitSummary['RA'] * np.pi / 180
        dec = visitSummary['DEC'] * np.pi / 180
        gn = wcsfit.Gnomonic(wcsfit.Orientation(wcsfit.SphericalICRS(ra, dec)))
        # TODO: probably want a different id/name in future:
        exposure = wcsfit.Exposure(name, gn)
        exposure.field = fieldNumber
        exposure.instrument = visitSummary['INSTRUMENTNUMBER']
    
        airmass = visitSummary['AIRMASS']
        exptime = visitSummary['EXPTIME']
        exposure.airmass = airmass
        exposure.exptime = exptime
        exposure.mjd = visitSummary['MJD']
        if False:
            exposure.pmEpoch
        if False:
            exposure.apcorr
        # TODO: add astrometric covariance
        exposures.append(exposure)

    for r, row in enumerate(hdus[4].data):
        extension = wcsfit.Extension()
        extension.exposure = row['EXPOSURE']
        visit = visits[row['EXPOSURE']]
        extension.device = row['DEVICE']
        extension.airmass = exposures[extension.exposure].airmass
        extension.magshift = 2.5 * np.log10(exposures[extension.exposure].exptime)
        if False:
            extension.apcorr
        # Add WCS here:
        #identity = wcsfit.IdentityMap()
        #icrs = wcsfit.SphericalICRS()
        #extension.startWcs = wcsfit.Wcs(identity, icrs, "ICRS_degrees", np.pi / 180.)
        
        #extension.addWcs(row["WCSIN"])
        exposure_name = exposures[row['EXPOSURE']].name
        if exposures[row['EXPOSURE']].instrument < 0:
            print("doing reference")
            """
            identity = wcsfit.IdentityMap()
            icrs = wcsfit.SphericalICRS()
            wcs = wcsfit.Wcs(identity, icrs, "ICRS_degrees", np.pi / 180.)
            """
            wcs = all_wcss[r]
            extension.startWcs = wcs
            extension.startWcs = all_wcss[r]
            extension.mapName = wcsfit.IdentityMap().getName()
        else:
            print("Extension mapname:", extension.mapName, exposure_name)
            # TODO: replace instrument name, band
            #d = {"INSTRUMENT": 'Hyper_Suprime-Cam', "DEVICE": str(extension.device), "EXPOSURE": exposure_name,
            #    "BAND": "HSC-z"}
            device_name = instruments[instrumentNumber].deviceNames.nameOf(extension.device)
            d = {"INSTRUMENT": instruments[instrumentNumber].name, 
                 "DEVICE": device_name,
                 "EXPOSURE": exposure_name,
                 "BAND": instruments[instrumentNumber].band}
            device_info = [instruments[instrumentNumber].name,
                           instruments[instrumentNumber].deviceNames.nameOf(extension.device),
                           exposure_name,
                           instruments[instrumentNumber].band]
            # TODO: UNCOMMENT BELOW (This is not getting set properly, I think startWcs when you leave this scope)
            """
            butler_wcs = butler.get('calexp', visit=int(visit), ccd=int(device_name)).getWcs()
            wcs = convert_sip_to_tpv(butler_wcs)
            """
            wcs = all_wcss[r]
            extension.startWcs = wcs  #demo_wcs
            #import ipdb; ipdb.set_trace()
            extension.startWcs.reprojectTo(exposures[extension.exposure].projection)
            extension.wcsName = f"{exposure_name}/{device_name}"
            extension.mapName = f"{extension.wcsName}/base"
            #device_info = ['Hyper_Suprime-Cam', str(extension.device), exposure_name, 'HSC-z'])
            wcsf.addMap(input_yaml, extension.mapName, device_info)
        extensions.append(extension)
        wcss.append(wcs)
    return exposures, extensions, wcss


def readInstruments(filtername):
    
    instruments = []
    instrument = wcsfit.Instrument()
    instrument.name = "Hyper_Suprime-Cam" # TODO: replace name here
    instrument.band = filtername
    # TODO: replace with real bounds (~25 pixel clipping around edges?)
    for i in list(range(9)) + list(range(10, 104)):
        instrument.addDevice(str(i), wcsfit.Bounds(0, 2047, 0, 4175))
    instruments.append(instrument)
        
    return instruments

def readInstruments_fits(hdu):
    
    instruments = []
    instrument = wcsfit.Instrument("Hyper_Suprime-Cam")
    #instrument.name = hdu.header['NAME']
    instrument.band = hdu.header['BAND']
    for row in hdu.data:
        instrument.addDevice(row['NAME'], wcsfit.Bounds(row['XMIN'], row['XMAX'],
                                                        row['YMIN'], row['YMAX']))
    instruments.append(instrument)
    return instruments


def readObjects(wcsf):
    
    sci_keys = {'xkey': 'base_SdssCentroid_x',
                'ykey': 'base_SdssCentroid_y',
                'xerrkey': 'base_SdssCentroid_xErr',
                'yerrkey': 'base_SdssCentroid_yErr'}
    # TODO: add PM keys
    ref_keys = {'xkey': 'coord_ra',
                'ykey': 'coord_dec',
                'errkey': 'coord_err'}
    
    for i, extn in enumerate(wcsf.extensions):
        print(i)
        # TODO: sort which are refs
        if False:  # extn in refs:
            keys = ref_keys
        else:
            keys = sci_keys
            visit = int(wcsf.exposures[extn.exposure].name)
            src = butler.get('src', visit=visit, detector=extn.device)
            
        ff = wcsfit.FTable(len(src))
        
        for k in keys:
            ff.addColumnDouble(src[keys[k]], keys[k])
        ff.addColumnDouble(np.zeros(len(src)), 'xyerr')
        # TODO: x and y errors need to be squared to go into covariance matrix
        xyerrkeys = ['xerrkey', 'yerrkey', 'xyerr']
        print('into readObjects_oneExt')
        wcsfit.readObjects_oneExtension(wcsf.exposures, i, ff, keys['xkey'], keys['ykey'],
                                        "", "", xyerrkeys, "", 0, "", 0, "", "", "", wcsf.extensions,
                                        wcsf.fieldProjections, True, True)
        print("one added")
            
            
# TODO: probably delete this fn
def getMatchArray(catalog, allowSelfMatches=False):
    # Probably move this back to c++
    sequence = []
    extn = []
    obj = []
    matches = 0
    # Now loop through matches in this catalog
    print("Start loop")
    
    catalogVector = catalog.toVector()
    import ipdb
    ipdb.set_trace()
    for pt in catalogVector:
        print("get point")
        point = pt.toVector()
        # Skip any Match that is below minimum match size
        if len(point) < minMatches:
            print("not enough matches")
            continue
        inExposures = []
        selfMatch = False
        if not allowSelfMatches:
            print("in aSM")
            print(len(point))
            for detection in point:
                print(detection)
                if detection.exposureNumber in inExposures:
                    selfMatch = True
                    break
                inExposures.append(detection.exposureNumber)
            if selfMatch:
                continue
        print("now read")
        #Add elements of this match to the output vectors
        seq = 0
        matches += 1
        for detection in point:
            sequence.append(seq)
            extn.append(detection.extensionNumber)
            obj.append(detection.objectNumber)
            seq += 1
        1/0

fields = [9813]
bands = ['HSC-z']
usePM = False

#skyMap = butler.get('skyMap')

fieldNames = []
fieldProjections = []
fieldObjects = []
for f, field in enumerate(fields):
    #fieldNames.append(wcsfit.spaceReplace(str(field)))
    fieldNames.append(str(field))
    # TODO: does tractinfo have a center?
    
    """
    tractInfo = skyMap.generateTract(field)    
    skyOrigin = tractInfo.getWcs().getSkyOrigin()
    ra = skyOrigin.getRa().asRadians()
    dec = skyOrigin.getDec().asRadians()
    print(skyOrigin.getRa(), skyOrigin.getDec())
    """
    ra = 2.62232
    dec = 0.0389454
    fieldProjection = wcsfit.Gnomonic(wcsfit.Orientation(wcsfit.SphericalICRS(ra, dec)))
    fieldProjections.append(fieldProjection)
    
    #expRefList = list(set(butler.registry.queryDatasets('calexp', dataId={"band": "r", "tract": field, 
    #                                                                      "patch": 25, 
    #                                                                      "skymap": "hsc_rings_v1"})))
    #visitSummaryRefs = list(set(butler.registry.queryDatasets('visitSummary', dataId={"tract": field})))
    visitSummaryRefs = []
    #exposures = readExposures(expRefList[:3], fieldNumber=f, instrumentNumber=0)
    
    ## WCSFoF setup:
    """
    fieldObject = wcsfit.Field()
    fieldObject.name = str(field)
    fieldObject.projection = fieldProjection
    # TODO: matchProjection does not need to be a property of wcsfof
    fieldObject.matchRadius = 1
    fieldObject.extent = 3.0  # TODO: set something more intelligent here
    fieldObjects.append(fieldObject)
    
    instObject = wcsfit.Instr()
    instObject.name = 'Hyper Suprime-Cam'
    for d in list(range(9)) + list(range(10, 104)):
        device = wcsfit.Device()
        device.name = str(d)
        instObject.append(device)
    instrumentObjects = [instObject]
    
    exposureObjects = []
    ## TODO: reconsider whether we want to load with calexps of visitSummary tables
    for visitSummaryRef in visitSummaryRefs[:5]:
        visitSummary = butler.get(visitSummaryRef)
        visInfo = visitSummary[0].getVisitInfo()
        expo = wcsfit.Expo()
        ra = visInfo.getBoresightRaDec().getRa().asRadians()
        dec = visInfo.getBoresightRaDec().getDec().asRadians()
        expo.name = str(visInfo.getId())
        expo.field = f
        expo.instrument = 0
        expo.pointing = wcsfit.SphericalICRS(ra, dec)
        exposureObjects.append(expo)
    
    #wcsfof = wcsfit_test.WCSFoF(fieldObjects, instrumentObjects, exposureObjects)
    wcsfof = wcsfit.FoFClass()
    wcsfof.fields = fieldObjects
    wcsfof.instruments = instrumentObjects
    wcsfof.exposures = exposureObjects
    
    
    iextn = 0
    for v, visitSummaryRef in enumerate(visitSummaryRefs[:5]):
        visitSummary = butler.get(visitSummaryRef)
        exposureNumber = v
        instrumentNumber = 0
        fieldNumber = f
        for row in visitSummary[:5]:
            thisAffinity = row['band'] # TODO: not sure if this is what is wanted for AFFINITY
            deviceNumber = row['id']
            srcCatalog = butler.get('src', dataId={"visit": row['visit'], "detector": row['id']})
            vx = srcCatalog['base_SdssCentroid_x']
            vy = srcCatalog['base_SdssCentroid_y']
            print(len(vx))
            vid = np.arange(len(vx)) # TODO: anything better here?
            isStar = np.ones(len(vx))
            calexp = butler.get('calexp', dataId={"visit": row['visit'], "detector": row['id']})
            wcs = convert_sip_to_tpv(calexp.getWcs())
            # TODO: need reprojectWCS to fields[field].projection here
            wcsfof.addCatalog(wcs, thisAffinity, exposureNumber, fieldNumber, instrumentNumber, deviceNumber,
                              iextn, isStar, vx, vy, vid)
            print('Added catalog, allPoints length: ', len(wcsfof.allPoints))
            iextn += 1  # TODO: replace with unique extn + visit number?
            
    # TODO: write out matches here, check against baseline result
    print(wcsfof.fields[0].catalogs.keys())
    catalog = wcsfof.fields[0].catalogs['STELLAR']
    print("getting matchArray")
    
    #wcsfof.writeMatches("testcat.fits")
    wcsfof.sortMatches(0)
    
    #seqs, extns, objs = getMatchArray(catalog)
    
    ## TODO: don't need next two lines anymore?
    #wcs_list = readExposures_wcsfof(visitSummaryRefs[:2])
    # TODO wcs - reprojectTo(fields[field].projection)
    """
    ## WCSFit setup:
    wcsf = wcsfit.FitClass()
    wcsf.verbose = verbose
    
    fieldNames = wcsfit.NameIndex()
    fieldNames.append(str(fields[f]))
    wcsf.fieldNames = fieldNames
    wcsf.fieldProjections = fieldProjections
    wcsf.fieldEpochs = [2015.5]
    
    #instruments = readInstruments(bands[0])
    instruments = readInstruments_fits(refactor_db[3])
    wcsf.instruments = instruments
    
    inputYAML = wcsfit.YAMLCollector("/home/csaunder/stack_projects/gbdes_tests/polyExposure.astro",
                                     "PixelMapCollection")
    inputYAML.addInput("Identity:\n  Type:  Identity\n")
    print("inputYAML made")
    
    print("reading exposures")
    ff = wcsfit.FTable(len(refactor_db[4].data))
    ff.addColumnStr(list(refactor_db[4].data["WCSIn"]), "WCSIn")
    all_wcss = wcsfit.readWCSs(ff)
    #import ipdb; ipdb.set_trace()
    #exposures, extensions = readExposuresExtensions(visitSummaryRefs, fieldNumber=f, instrumentNumber=0)
    exposures, extensions, wcss = readExposuresExtensions_fits(refactor_db, inputYAML, wcsf.instruments, fieldNumber=f, instrumentNumber=0)
    
    # TODO: Add special fn to read reference catalog
    wcsf.setExposures(exposures, sysError, referenceSysError)
    
    wcsf.setExtensions(extensions)
    wcsf.setRefWCSNames()
    print("check refset")
    
    wcsf.setupMaps(inputYAML)
    #wcsf.setupMaps()
    
    
    
    #wcsfit.readFields(refactor_db_file, "test_wcscat.fits", wcsf.fieldNames, wcsf.fieldProjections, 
    #                  wcsf.fieldEpochs, pmEpoch)
    
    skipSet = wcsfit.ExtensionObjectSet("")
    colorExtensions = [wcsfit.ColorExtension() for e in extensions]
    wcsfit.readMatches(wcsfof.sequence, wcsfof.extn, wcsfof.obj, wcsf.matches, extensions, colorExtensions,
                       skipSet, minMatches, usePM)
    import ipdb
    ipdb.set_trace()
    #print(wcsf.matches.size())
    ## TODO: fill this in
    readObjects(wcsf)
      
    wcsf.fit()
    
    # Return fit

    #fitter = wcsfit_test.WCSFit(exposures, extensions, fieldNames, fieldEpochs, fieldProjections, 
    #                            instruments, matches, fixMapsList=[1, 2])
    #print(fitter.exposures)
