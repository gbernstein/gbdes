
#import numpy as np
#import matplotlib.pyplot as plt
#from glob import glob
#from astropy.io import fits

from lsst.daf.butler import Butler
import wcsfit

#from lsst.daf.persistence import Butler as ButlerGen2

import numpy as np
import sip_tpv
from astropy.io import fits
import astropy.time
import os.path
import copy

useButler = False

if useButler:
    butler = Butler('/repo/main', collections="HSC/runs/RC2/w_2021_34/DM-31524")
    #butler = ButlerGen2('/datasets/hsc/repo/rerun/DM-23243/SFM/DEEP')
    skyMap = butler.get('skyMap', skymap='hsc_rings_v1')
    print("got skymap")
#else:
if True:
    #refactor_db_file = "/home/csaunder/stack_projects/gbdes_tests/refactor_tests/refactor.fof"
    refactor_db_file = "/home/csaunder/stack_projects/gbdes_tests/refactor_tests/gen3.fof"
    refactor_db = fits.open(refactor_db_file)

prefix = '/home/csaunder/stack_projects/gbdes_tests/refactor_tests/'

# GBDES defaults:
reserveFraction = 0.05
randomNumberSeed = 0
skipFile = ''
clipThresh = 5
maxError = 100
sysError = 2.0
referenceSysError = 2.0
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
outCatalog = 'wcscat3.fits'
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
        index = p
        if (index == 3 or index == 11 or index == 23 or index == 39):
             continue
        if (index >= 39):
            index -= 1
        if (index>=23): index -= 1
        if (index>=11): index -= 1
        if (index>=3): index -= 1
    
        pv1_vals[index] = float(sip_tpv.pvsiputils.calcpv(pvrange, 1, p, sipx, sipy, tpvx, tpvy).evalf())
        pv2_vals[index] = float(sip_tpv.pvsiputils.calcpv(pvrange, 2, p, sipx, sipy, tpvx, tpvy).evalf())
    
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

    #pv1 = pv1.transpose()
    pv2 = pv2.transpose()
    return pv1, pv2

POLYSTEP = 1.0/3600.0
def convert_sip_to_tpv(butler_wcs, name=""):
    # This function is an adaptation of readTPV in gbdes/src/subs/TPVMap.cpp
    
    # Format copied from sip_tpv.sip_pv function
    try:
        cd = np.matrix(butler_wcs.getCdMatrix())
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

    #pole = wcsfit.SphericalICRS(fits_metadata.get('CRVAL1') * wcsfit.DEGREE,
    #                            fits_metadata.get('CRVAL2') * wcsfit.DEGREE)
    #orientIn = wcsfit.Orientation()
    #orientIn.set(pole, 0.0)
    #tp = wcsfit.Gnomonic(orientIn.getPole(), orientIn)
    print("starting gnomonic")
    tp = wcsfit.Gnomonic(fits_metadata.get('CRVAL1') * wcsfit.DEGREE,
                         fits_metadata.get('CRVAL2') * wcsfit.DEGREE)
    print("gnomonic tested")


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
    print('1')
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

def getSources():
    
    blended = cal_src['base_Blendedness_abs']
    extended = cal_src['base_ClassificationExtendedness_value']
    parent = cal_src['parent']

    ind = (extended == 0) & (blended == 0) & (parent == 0)
    src_cat = src_cat[ind]

#def readExposuresExtensions(refs, input_yaml, fields, band, fieldNumber=0, instrumentNumber=0, isReference=False):
def readExposuresExtensions(fields, band, visits, visitsCcds, ccds, fieldNumber=0, instrumentNumber=0, isReference=False):
    # TODO: add ExposureColorPriorities, ColorExtension / decide if necessary
    
    #extensions = []
    #exposures = []
    expNames = []
    fieldNumbers = []
    instrumentNumbers = []
    ras = []
    decs = []
    airmasses = []
    exposureTimes = []
    mjds = []
    extExps = []
    extDevices = []
    # TODO: this has to be redone to be input to WCSFoF step!!!
    for f, field in enumerate(fields):
        visitTable = butler.get('visitTable', tract=field)
        visitTable = visitTable[visitTable['filterName'] == band]
        sortOrder = np.array([np.flatnonzero(visitTable == vis) for vis in visits]).reshape(1, -1)
        expNames.extend(visitTable['visitId'].to_numpy().astype(str)[sortOrder])
        fieldNumbers.extend(np.ones(len(visitTable), dtype=int) * f)
        instrumentNumbers.extend(np.ones(len(visitTable), dtype=int) * instrumentNumber)
        ras.extend(visitTable['ra'].to_list()[sortOrder])
        decs.extend(visitTable['decl'].to_list()[sortOrder])
        airmasses.extend(visitTable['airmass'].to_list()[sortOrder])
        exposureTimes.extend(visitTable['expTime'].to_list()[sortOrder])
        mjds.extend(astropy.time.Time(visitTable['obsStart']).mjd[sortOrder])
        
        ccdVisitTable = butler.get('CcdVisitTable', tract=field)
        for v, visit in enumerate(visitTable['visitId'].to_numpy()):
            ccds = ccdVisitTable[ccdVisitTable['visitId'] == visit]
            extExp = np.ones(len(ccds), dtype=int) * v
            extDevice = ccds.index
            extExps.extend(list(extExp))
            extDevices.extend(list(extDevice))
            
    exposuresHelper = wcsfit.ExposuresHelper((expNames),
                                             (fieldNumbers),
                                             (instrumentNumbers),
                                             (ras),
                                             (decs),
                                             (airmasses),
                                             (exposureTimes),
                                             (mjds))
    import ipdb; ipdb.set_trace()
    """
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
    return exposures, extensions"""
    return extExps, extDevices, exposuresHelper


def readExposuresExtensions_fits(hdus, input_yaml, instruments, fieldNumber=0, instrumentNumber=0, isReference=False):
    # TODO: add ExposureColorPriorities, ColorExtension / decide if necessary
    
    extensions = []
    exposures = []
    wcss = []
    visits = []
    #names, ras, decs, fieldNum, insts, airmasses, exptimes, mjds = [], [], [], [], [], [], [], []
    expTable = hdus[2].data
    exposuresHelper = wcsfit.ExposuresHelper(expTable['NAME'],
                                             np.ones(len(expTable), dtype=int) * fieldNumber,
                                             expTable['INSTRUMENTNUMBER'],
                                             expTable['RA'] * np.pi / 180,
                                             expTable['DEC'] * np.pi / 180,
                                             expTable['AIRMASS'],
                                             expTable['EXPTIME'],
                                             expTable['MJD'])
    """
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
        print( 'radec', ra, dec)
        #gn = wcsfit.Gnomonic(wcsfit.Orientation(wcsfit.SphericalICRS(ra, dec)))
        # TODO: probably want a different id/name in future:
        #exposure = wcsfit.Exposure(name, gn)
        #exposure.field = fieldNumber
        #exposure.instrument = visitSummary['INSTRUMENTNUMBER']
    
        airmass = visitSummary['AIRMASS']
        exptime = visitSummary['EXPTIME']
        #exposure.airmass = airmass
        #exposure.exptime = exptime
        #exposure.mjd = visitSummary['MJD']
        names.append(name)
        ras.append(ra)
        decs.append(dec)
        fieldNum.append(fieldNumber)
        insts.append(visitSummary['INSTRUMENTNUMBER'])
        airmasses.append(airmass)
        exptimes.append(exptime)
        mjds.append(visitSummary['MJD'])
        if False:
            exposure.pmEpoch
        if False:
            exposure.apcorr
        # TODO: add astrometric covariance
        #exposures.append(exposure)
    """
    #extExps, extDevices = [], []
    extExps = list(hdus[4].data['EXPOSURE'].astype(int))
    extDevices = list(hdus[4].data['DEVICE'].astype(int))
    print(len(extExps), len(extDevices))
    for r, row in enumerate(hdus[4].data):
        continue
        extExps.append(row['EXPOSURE'])
        extDevices.append(row['DEVICE'])
        continue
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
            print("identity: ", wcsfit.IdentityMap().getName())
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
    
    #return exposures, extensions, wcss, names, ras, decs, fieldNum, insts, airmasses, exptimes, mjds, extExps, extDevices, exposuresHelper
    return extExps, extDevices, exposuresHelper

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
    src_files = refactor_db[4].data['FILENAME']
    #for i, extn in enumerate(wcsf.extensions):
    for i, extn in enumerate(src_files):
        print(i)
        # TODO: switch to butler input
        """
        if extn.device < 0:  # extn in refs:
            keys = ref_keys
        else:
            keys = sci_keys
            # TODO: regularize visit names!
            visit = int(wcsf.exposures[extn.exposure].name[4:])
            #src = butler.get('src', visit=visit, detector=extn.device)
            src = butler.get('src', visit=visit, ccd=extn.device)
            
        ff = wcsfit.FTable(len(src))
        
        if extn.device < 0:
            keys = ref_keys
            for k in keys:
                ff.addColumnDouble(src[keys[k]], keys[k])
            xyerrkeys = ['errkey']
        else:
            for k in keys:
                if 'err' in k:
                    ff.addColumnDouble(src[keys[k]]**2, keys[k])
                else:
                    ff.addColumnDouble(src[keys[k]], keys[k])
            ff.addColumnDouble(np.zeros(len(src)), 'xyerr')
            xyerrkeys = [keys['xerrkey'], keys['yerrkey'], 'xyerr']
        """
        #wcsfit.readObjects_oneExtension(wcsf.exposures, i, ff, keys['xkey'], keys['ykey'],
        #                                "", "", xyerrkeys, "", 0, "", 0, "", "", "", wcsf.extensions,
        #                                wcsf.fieldProjections, True, True)
        #wcsf.setObjects(i, ff, keys['xkey'], keys['ykey'], "", "", xyerrkeys, "", 0, "", 0, "", "", "")
        print('into readObjects_oneExt', src_files[i])
        srcs = fits.open('/home/csaunder/stack_projects/gbdes_tests/' + src_files[i])
        objects = srcs[1].data
        ff = wcsfit.FTable(len(objects))
        
        #if extn.device < 0:
        if 'base_SdssCentroid_x' not in objects.dtype.names:
            ff.addColumnDouble(objects['coord_ra'], 'coord_ra')
            ff.addColumnDouble(objects['coord_dec'], 'coord_dec')
            ff.addColumnDouble(objects['coord_err'] * np.pi / 180, 'coord_err')
            #wcsf.setObjects(i, ff, objects['coord_ra'], 'coord_ra', 'coord_dec', '', '', ['coord_err'], '',
            #                0, '', 0, '', '', '')
            d = {'ra': objects['coord_ra'],
                 'dec': objects['coord_dec'],
                 'err': objects['coord_err'] * np.pi / 180}
            wcsf.setObjects(i, d, 'ra', 'dec', ['err'])
            
        else:
            ff.addColumnDouble(objects['base_SdssCentroid_x'], 'x')
            ff.addColumnDouble(objects['base_SdssCentroid_y'], 'y')
            ff.addColumnDouble(objects['base_SdssCentroid_Err'], 'err')
            d = {"x": objects['base_SdssCentroid_x'],
                 "y": objects['base_SdssCentroid_y'],
                 "err": objects['base_SdssCentroid_Err'],}
            #wcsf.setObjects(i, ff, objects['base_SdssCentroid_x'], 'x', 'y', '', '', ['err'], '', 0, '', 0, '', '', '')
            wcsf.setObjects(i, d, 'x', 'y', ['err'])
        
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
altbands = ['z']
bands = ['HSC-Z']
usePM = False

ras = []
decs = []
epochs = []
fieldNames = []
#fieldProjections = []
#fieldObjects = []
print("2")
for f, field in enumerate(fields):
    #fieldNames.append(wcsfit.spaceReplace(str(field)))
    fieldNames.append(str(field))
    # TODO: does tractinfo have a center?
    print("2.1")
    if useButler:
        print("2.2")
        skyMap = butler.get('skyMap', skymap='hsc_rings_v1')
        print("2.21")
        tractInfo = skyMap.generateTract(field)    
        print("2.22")
        skyOrigin = tractInfo.ctr_coord
        print("2.3")
        ra = skyOrigin.getRa().asDegrees()
        dec = skyOrigin.getDec().asDegrees()
        print(skyOrigin.getRa(), skyOrigin.getDec())
    else:
        print("2.4")
        ind = refactor_db[1].data['NAME'] == str(field)
        ra = refactor_db[1].data['RA'][ind][0]
        dec = refactor_db[1].data['DEC'][ind][0]
        #ra = 2.62232
        #dec = 0.0389454
    ras.append(ra)
    decs.append(dec)
    epochs.append(2015.5)
    
    #fieldProjection = wcsfit.Gnomonic(wcsfit.Orientation(wcsfit.SphericalICRS(ra, dec)))
    #fieldProjections.append(fieldProjection)
    
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
print("3")
## WCSFit setup:
fh = wcsfit.Fields(fieldNames, ras, decs, epochs)
print("made fields")

fieldMatchRadius = 1 * wcsfit.ARCSEC / wcsfit.DEGREE
fieldExtents = [3.0]  # TODO: set something more intelligent here

# Get instrument info:
if False:#useButler:
    # TODO: get from butler
    #instrumentName = "Hyper_Suprime-Cam"
    instrumentName = 'HSC'
    instruments = [wcsfit.Instrument(instrumentName)]
else:
    instrumentName = 'HSC'
    instruments = readInstruments_fits(refactor_db[3])
print("made instruments")

if useButler:
    visitSummaryRefs = list(butler.registry.queryDatasets('visitSummary', band=altbands[0], tract=fields[0],
                                                        instrument=instrumentName))
    visitSummaryRefs = visitSummaryRefs[:5]

inputYAML = wcsfit.YAMLCollector("/home/csaunder/stack_projects/gbdes_tests/polyExposure.astro",
                                    "PixelMapCollection")
inputYAML.addInput("Identity:\n  Type:  Identity\n")
print("inputYAML made")

# Get initial WCS info:
if useButler:
    visits = []
    visitCcds = []
    ccds = []
    all_wcss = []
    all_wcss2 = []
    #butler_wcs = butler.get('calexp', visit=95066, ccd=103).getWcs()
    #for visitSummaryRef in visitSummaryRefs:
    #    visitSummary = butler.get(visitSummaryRef)
    #    for ccdRow in visitSummary:
    #        butler_wcs = butler.get('calexp.wcs', visit=ccdRow['visit'], detector=ccdRow['id'])
    #        wcs = convert_sip_to_tpv(butler_wcs)
    #        all_wcss.append(wcs)
    for extn in refactor_db[4].data:
        if extn['DEVICE'] == -1: 
            identity = wcsfit.IdentityMap()
            icrs = wcsfit.SphericalICRS()
            #wcs = wcsfit.Wcs(identity, icrs, "ICRS_degrees", np.pi / 180.)
            wcs = wcsfit.Wcs(identity, icrs, "Identity", np.pi / 180.)
            #wcs2 = wcsfit.Wcs(identity, icrs, "ICRS_degrees", np.pi / 180.)
            wcs2 = wcsfit.Wcs(identity, icrs, "Identity", np.pi / 180.)
        else:
            visit = int(float(refactor_db[2].data[extn['EXPOSURE']]['NAME'][4:]) / 100)
            ccd = int(refactor_db[3].data[extn['DEVICE']]['NAME'])
            visitCcds.append(visit)
            ccds.append(ccd)
            print(refactor_db[2].data[extn['EXPOSURE']]['NAME'], visit, ccd)
            butler_wcs = butler.get('calexp.wcs', visit=visit, detector=ccd)
            wcs = convert_sip_to_tpv(butler_wcs)
            wcs2 = convert_sip_to_tpv(butler_wcs)
        all_wcss.append(wcs)
        all_wcss2.append(wcs2)
                
    print(len(all_wcss))
            
else:
    ff = wcsfit.FTable(len(refactor_db[4].data))
    ff.addColumnStr(list(refactor_db[4].data["WCSIn"]), "WCSIn")
    print('reading WCSs')
    all_wcss = wcsfit.readWCSs(ff)
    all_wcss2 = wcsfit.readWCSs(ff)
print("made wcss")

print("reading exposures")

if False:#useButler:
    extExps, extDevices, exposuresHelper = readExposuresExtensions(fields, bands[0], visits, visitCcds, ccds)
    for extDevice in set(extDevices):
        """deviceSum = visitSummary[visitSummary['id'] == extDevice][0]
        device_bbox = wcsfit.Bounds(deviceSum['bbox_min_x'],
                                    deviceSum['bbox_max_x'],
                                    deviceSum['bbox_min_y'],
                                    deviceSum['bbox_max_y'])"""
        # TODO: get correct bbox in
        device_bbox = wcsfit.Bounds(0, 2048, 0, 4096)
        instruments[0].addDevice(str(extDevice), device_bbox)
    print("exposures:", max(extExps), "ext devices:", len(extDevices))
else:
    #exposures, extensions, wcss = readExposuresExtensions_fits(refactor_db, inputYAML, wcsf.instruments, fieldNumber=f, instrumentNumber=0)
    extExps, extDevices, exposuresHelper = readExposuresExtensions_fits(refactor_db, inputYAML, instruments, fieldNumber=f, instrumentNumber=0)
print("made exps, exts")

#fh2 = copy.copy(fh)
#instruments2 = copy.copy(instruments)
#exposures2 = copy.copy(exposuresHelper)
#fofFit = wcsfit.FoFClass(fh2, instruments2, exposuresHelper2, fieldExtents, fieldMatchRadius)

fofFit = wcsfit.FoFClass(fh, instruments, exposuresHelper, fieldExtents, fieldMatchRadius)
print("fof built")

# Add catalogs to FoF algorithm:
for i, extn in enumerate(refactor_db[4].data):
    if extn['DEVICE'] == -1: 
        print('do ref')
        
    else:
        visit = int(float(refactor_db[2].data[extn['EXPOSURE']]['NAME'][4:]) / 100)
        ccd = int(refactor_db[3].data[extn['DEVICE']]['NAME'])
        # TODO:  add butler call here
    src_file = extn['FILENAME']
    deviceNumber = extn['DEVICE']
    exposureNumber = extn['EXPOSURE']
    fieldNumber = exposuresHelper.fieldNumbers[exposureNumber]
    instrumentNumber = exposuresHelper.instrumentNumbers[exposureNumber]
    print(src_file)
    print(fieldNumber, instrumentNumber)
    srcs = fits.open('/home/csaunder/stack_projects/gbdes_tests/' + src_file)
    if extn['DEVICE'] == -1: 
        vx = srcs[1].data['coord_ra'] - ras[0]
        vy = srcs[1].data['coord_dec'] - decs[0]
    else:
        vx = srcs[1].data['base_SdssCentroid_x']
        vy = srcs[1].data['base_SdssCentroid_y']
    vid = np.arange(len(vx))
    isStar = np.ones(len(vx))
    wcs = all_wcss2[i]
    if extn['DEVICE'] >= 0:
        fofFit.reprojectWCS(wcs, fieldNumber)
    
    # TODO: need to reproject wcs to fields->projection
    fofFit.addCatalog(wcs, extn['AFFINITY'], exposureNumber, fieldNumber, instrumentNumber, deviceNumber,
                      i, isStar, vx, vy, vid)

    # TODO: replace `i` with unique extn + visit number?

fofFit.sortMatches(0, minMatches=minMatches)

sequence = refactor_db[5].data['SequenceNumber']
extn = refactor_db[5].data['Extension']
obj = refactor_db[5].data['Object']

skipSet = wcsfit.ExtensionObjectSet("")
print(len(all_wcss), len(extExps))

print("start WCSFit")
import ipdb
ipdb.set_trace()
wcsf = wcsfit.WCSFit(fh, instruments, exposuresHelper, 
                       extExps, extDevices, inputYAML, all_wcss, list(sequence.astype(int)),
                       list(extn.astype(int)), list(obj.astype(int)),
                       sysErr=sysError, refSysErr=referenceSysError, usePM=False, verbose=2)
print("wcsfit built")


#readObjects(wcsf)
"""fieldNames = wcsfit.NameIndex()
fieldNames.append(str(fields[f]))
wcsf.fieldNames = fieldNames
wcsf.fieldProjections = fieldProjections
wcsf.fieldEpochs = [2015.5]"""
#wcsf.fields = wcsfit.Fields([fields[f]], [ra], [dec], [2015.5])

#instruments = readInstruments(bands[0])

#wcsf.instruments = instruments




#exposures, extensions = readExposuresExtensions(visitSummaryRefs, fieldNumber=f, instrumentNumber=0)
print("reading exposures")
#exposures, extensions, wcss = readExposuresExtensions_fits(refactor_db, inputYAML, wcsf.instruments, fieldNumber=f, instrumentNumber=0)

# TODO: Add special fn to read reference catalog
#wcsf.setExposures(exposures, sysError, referenceSysError)

#wcsf.setExtensions(extensions)
#wcsf.setRefWCSNames()
#print("check refset")
#wcsf.setupMaps(inputYAML)

#wcsfit.readFields(refactor_db_file, "test_wcscat.fits", wcsf.fieldNames, wcsf.fieldProjections, 
#                  wcsf.fieldEpochs, pmEpoch)

#colorExtensions = [None for e in extensions]
#wcsf.setMatches(sequence, extn, obj, skipSet)
#print(wcsf.getMatchLength())
#wcsfit.readMatches(sequence, extn, obj, wcsf.matches, wcsf.extensions, colorExtensions,
#                   skipSet, minMatches, usePM)
readObjects(wcsf)
#print(wcsf.matches.size())

### TEST WCS:
#test_wcs = all_wcss[0]
#butler_wcs = butler.get('calexp', visit=95066, ccd=103).getWcs()
#convert_sip_to_tpv(butler_wcs)


print("starting fit")
wcsf.fit()
print("fit finished")
wcsf.saveResults(os.path.join(prefix, outWcs),
                 os.path.join(prefix, outCatalog),
                 os.path.join(prefix, starCatalog))
print("save finished")
