// Code for fitting parameters of defaulted pixel maps to the staring WCS solution.
#include "Std.h"
#include "PixelMapCollection.h"
#include "Instrument.h"
#include "Match.h"
#include "WcsSubs.h"
#include "Units.h"
#include "FitSubroutines.h"
#include <random>
using namespace astrometry;

void fitDefaulted(PixelMapCollection &pmc, set<Extension *> extensions,
                  const vector<unique_ptr<Instrument>> &instruments,
                  const vector<unique_ptr<Exposure>> &exposures, bool logging) {
    // Make a new pixel map collection that will hold only the maps
    // involved in this fit.
    PixelMapCollection pmcFit;

    // Take all wcs's used by these extensions and copy them into
    // pmcFit
    for (auto extnptr : extensions) {
        auto pm = pmc.cloneMap(extnptr->mapName);
        pmcFit.learnMap(*pm);
        delete pm;
    }

    // Find all the atomic map components that are defaulted.
    // Fix the parameters of all the others
    set<string> defaultedAtoms;
    set<string> fixAtoms;
    for (auto mapname : pmcFit.allMapNames()) {
        if (!pmc.isAtomic(mapname)) continue;
        if (pmc.getDefaulted(mapname))
            defaultedAtoms.insert(mapname);
        else
            fixAtoms.insert(mapname);
    }
    // We are done if nothing is defaults
    if (defaultedAtoms.empty()) return;

    // If the exposure has any defaulted maps, initialize them
    /**/ if (logging) {
        cerr << "Initializing maps ";
        for (auto name : defaultedAtoms) cerr << name << " ";
        cerr << endl;
    }
    pmcFit.learnMap(IdentityMap());  // Make sure we know this one
    pmcFit.setFixed(fixAtoms);
    pmcFit.rebuildParameterVector();

    auto identityMap = pmcFit.issueMap(IdentityMap().getName());

    // Make Matches at a grid of points on each extension's device,
    // matching pix coords to the coordinates derived from startWCS.
    MCat matches;
    for (auto extnptr : extensions) {
        const Exposure &expo = *exposures[extnptr->exposure];
        // Get projection used for this extension, and set up
        // startWCS to reproject into this system.
        astrometry::Orientation proj = *expo.projection;
        extnptr->startWcs->reprojectTo(*expo.projection);

        // Get a realization of the extension's map
        auto map = pmcFit.issueMap(extnptr->mapName);

        // Get the boundaries of the device it uses
        Bounds<double> b = instruments[expo.instrument]->domains[extnptr->device];
        if (logging) {
            cerr << "instrument, device: " << expo.instrument << " " << extnptr->device
                 << endl;
            cerr << "bounds: " << b.getXMin() << " " << b.getXMax() << endl;
            double txw, tyw;
            extnptr->startWcs->toWorld(b.getXMin(), b.getYMin(), txw, tyw);
            cerr << "check p1 " << txw << " " << tyw << endl;
        }
        // Generate a grid of matched Detections
        const int nGridPoints = 512;  // Number of test points for map initialization

        // "errors" on world coords of test points
        const double testPointSigma = 0.01 * ARCSEC / WCS_UNIT;
        const double fitWeight = pow(testPointSigma, -2.);
        // Put smaller errors on the "reference" points.  Doesn't really matter.
        const double refWeight = 10. * fitWeight;
        // Distribute points equally in x and y, but shuffle the y coords
        // so that the points fill the rectangle
        vector<int> vx(nGridPoints);
        for (int i = 0; i < vx.size(); i++) vx[i] = i;
        vector<int> vy = vx;
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(vy.begin(), vy.end(), g);
        double xstep = (b.getXMax() - b.getXMin()) / nGridPoints;
        double ystep = (b.getYMax() - b.getYMin()) / nGridPoints;
        for (int i = 0; i < vx.size(); i++) {
            double xpix = b.getXMin() + (vx[i] + 0.5) * xstep;
            double ypix = b.getXMin() + (vy[i] + 0.5) * ystep;
            double xw, yw;
            extnptr->startWcs->toWorld(xpix, ypix, xw, yw);  // startWCS has no color!
            unique_ptr<Detection> dfit(new Detection);
            unique_ptr<Detection> dref(new Detection);
            if (logging && (i == 10)) {
                cerr << "pix check " << xpix << " " << ypix << " " << xw
                     << " " << yw << endl;
            }
            dfit->xpix = xpix;
            dfit->ypix = ypix;
            dref->xpix = xw;
            dref->ypix = yw;
            dref->xw = xw;
            dref->yw = yw;
            // Note colors are defaulted to NoData in Detection.  Set colors
            // to zero here and defaulted maps will be appropriate to
            // zero-color objects.  Hopefully zero is not far from the reference
            // color of any color terms??
            dref->color = 0.;
            dfit->color = 0.;

            map->toWorld(xpix, ypix, xw, yw, dfit->color);
            dfit->xw = xw;
            dfit->yw = yw;

            // Set up errors and maps for these matched "detections"
            Matrix22 eye(0.);
            eye(0, 0) = eye(1, 1) = 1.;
            dfit->invCov = eye;
            dfit->fitWeight = fitWeight;
            dref->invCov = eye;
            dref->fitWeight = refWeight;
            dref->map = identityMap;
            dfit->map = map;

            matches.emplace_back(new Match(std::move(dfit)));
            matches.back()->add(std::move(dref));
        }
    }

    // Build CoordAlign object and solve for defaulted parameters
    CoordAlign ca(pmcFit, matches);
    ca.setRelTolerance(0.01);  // weaker tolerance for fit convergence
    ca.fitOnce(logging);

    // Copy defaulted parameters back into the parent pmc.
    for (auto mapname : defaultedAtoms) {
        auto pm = pmcFit.cloneMap(mapname);
        pmc.copyParamsFrom(*pm);
        delete pm;
    }

    return;
}

// Define and issue WCS for each extension in use, and set projection to
// field coordinates.
void setupWCS(const vector<unique_ptr<SphericalCoords>> &fieldProjections,
              const vector<unique_ptr<Instrument>> &instruments,
              const vector<unique_ptr<Exposure>> &exposures, vector<unique_ptr<Extension>> &extensions,
              PixelMapCollection &pmc) {
    for (auto const &extnptr : extensions) {
        if (!extnptr) continue;  // Not in use
        Exposure &expo = *exposures[extnptr->exposure];
        int ifield = expo.field;
        if (expo.instrument < 0) {
            // Tag & reference exposures have no instruments and no fitting
            // being done.  Coordinates are fixed to xpix = xw.
            // Build a simple Wcs that takes its name from the field
            SphericalICRS icrs;
            if (!pmc.wcsExists(extnptr->wcsName)) {
                // If we are the first reference/tag exposure in this field:
                pmc.defineWcs(extnptr->wcsName, icrs, extnptr->mapName, WCS_UNIT);
                auto wcs = pmc.issueWcs(extnptr->wcsName);
                // And have this Wcs reproject into field coordinates, learn as map
                wcs->reprojectTo(*fieldProjections[ifield]);
                pmc.learnMap(*wcs, false, false);
            }
            extnptr->wcs = pmc.issueWcs(extnptr->wcsName);
            extnptr->map = pmc.issueMap(extnptr->wcsName);
        } else {
            // Real instrument, WCS goes into its exposure coordinates
            pmc.defineWcs(extnptr->wcsName, *expo.projection, extnptr->mapName, WCS_UNIT);
            extnptr->wcs = pmc.issueWcs(extnptr->wcsName);

            // Reproject this Wcs into the field system and get a SubMap including reprojection:
            extnptr->wcs->reprojectTo(*fieldProjections[ifield]);
            pmc.learnMap(*extnptr->wcs, false, false);
            extnptr->map = pmc.issueMap(extnptr->wcsName);
        }
    }  // end extension loop
}

list<int> pickExposuresToInitialize(const vector<unique_ptr<Instrument>> &instruments,
                                    const vector<unique_ptr<Exposure>> &exposures,
                                    const vector<unique_ptr<Extension>> &extensions,
                                    PixelMapCollection &pmc) {
    list<int> exposuresToInitialize;
    for (int iInst = 0; iInst < instruments.size(); iInst++) {
        if (!instruments[iInst]) continue;  // Not in use
        auto &instr = *instruments[iInst];
        // Classify the device maps for this instrument
        set<int> defaultedDevices;
        set<int> initializedDevices;
        set<int> unusedDevices;

        // And the exposure maps as well:
        set<int> defaultedExposures;
        set<int> initializedExposures;
        set<int> unusedExposures;
        set<int> itsExposures;  // All exposure numbers using this instrument

        for (int iDev = 0; iDev < instr.nDevices; iDev++) {
            string mapName = instr.mapNames[iDev];
            if (!pmc.mapExists(mapName)) {
                unusedDevices.insert(iDev);
                continue;
            }
            if (pmc.getDefaulted(mapName)) {
                defaultedDevices.insert(iDev);
            } else {
                initializedDevices.insert(iDev);
            }
        }

        for (int iExpo = 0; iExpo < exposures.size(); iExpo++) {
            if (!exposures[iExpo]) continue;  // Not in use
            if (exposures[iExpo]->instrument == iInst) {
                itsExposures.insert(iExpo);
                string mapName = exposures[iExpo]->name;
                if (!pmc.mapExists(mapName)) {
                    unusedExposures.insert(iExpo);
                } else if (pmc.getDefaulted(mapName)) {
                    defaultedExposures.insert(iExpo);
                } else {
                    initializedExposures.insert(iExpo);
                }
            }
        }

        /**/ cerr << "Done collecting initialized/defaulted" << endl;

        // No need for any of this if there are no defaulted devices
        if (defaultedDevices.empty()) continue;

        // Now take an inventory of all extensions to see which device
        // solutions are used in coordination with which exposure solutions
        vector<set<int>> exposuresUsingDevice(instr.nDevices);
        for (auto const &extnptr : extensions) {
            if (!extnptr) continue;  // Not in use
            int iExpo = extnptr->exposure;
            if (itsExposures.count(iExpo) > 0) {
                // Extension is from one of the instrument's exposures
                int iDev = extnptr->device;
                if (pmc.dependsOn(extnptr->mapName, exposures[iExpo]->name) &&
                    pmc.dependsOn(extnptr->mapName, instr.mapNames[iDev])) {
                    // If this extension's map uses both the exposure and device
                    // maps, then enter this dependence into our sets
                    exposuresUsingDevice[iDev].insert(iExpo);
                    if (unusedExposures.count(iExpo) > 0 || unusedDevices.count(iDev) > 0) {
                        cerr << "Logic problem: extension map " << extnptr->mapName
                             << " is using allegedly unused exposure or device map" << endl;
                        exit(1);
                    }
                }
            }
        }

        /**/ cerr << "Done building exposure/device graph" << endl;

        // Find a non-defaulted exposure using all of the defaulted devices
        int exposureForInitializing = -1;
        for (auto iExpo : initializedExposures) {
            bool hasAllDevices = true;
            for (auto iDev : defaultedDevices) {
                if (exposuresUsingDevice[iDev].count(iExpo) == 0) {
                    hasAllDevices = false;
                    break;
                }
            }
            if (hasAllDevices) {
                exposureForInitializing = iExpo;
                break;
            }
        }

        if (exposureForInitializing < 0) {
            cerr << "Could not find an exposure that can initialize defaulted devices \n"
                 << "for instrument " << instr.name << endl;
            cerr << "Write more code if you want to exploit more complex situations \n"
                 << "where some non-defaulted devices can initialize a defaulted exposure." << endl;
            exit(1);
        } else {
            exposuresToInitialize.push_back(exposureForInitializing);
            cerr << "Using exposure " << exposures[exposureForInitializing]->name
                 << " to initialize defaulted devices on instrument " << instr.name << endl;
        }
    }  // end instrument loop

    return exposuresToInitialize;
}
