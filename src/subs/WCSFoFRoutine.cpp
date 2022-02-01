// Program to match catalogs based on world coordinates.

#include "WCSFoFRoutine.h"
#include "Std.h"
#include "FoF.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "NameIndex.h"
#include "Pset.h"
#include <map>
#include <iostream>
#include "PixelMapCollection.h"
#include "TPVMap.h"
#include "TemplateMap.h"
#include "StringStuff.h"
#include "Units.h"
#include "FitSubroutines.h"

using stringstuff::stripWhite;

using astrometry::WCS_UNIT;  // we will want all WCS's to operate in same units

FoFClass::FoFClass() {
    // Teach PixelMapCollection about all types of PixelMaps it might need to deserialize
    loadPixelMapParser();

    // If WCS is coming from a serialized PixelMapCollection, I'll keep last-used one around
    // to avoid re-parsing it all the time:
    std::string currentPixelMapCollectionFileName = "";
}

FoFClass::FoFClass(Fields &fields_, std::vector<std::shared_ptr<Instrument>> instruments_,
                   ExposuresHelper exposures_, std::vector<double> fieldExtents, double matchRadius) {
    // Teach PixelMapCollection about all types of PixelMaps it might need to deserialize
    loadPixelMapParser();

    // If WCS is coming from a serialized PixelMapCollection, I'll keep last-used one around
    // to avoid re-parsing it all the time:
    std::string currentPixelMapCollectionFileName = "";

    fields.reserve(fields_.names().size());
    for (int i = 0; i < fields_.names().size(); i++) {
        std::unique_ptr<Field> f(new Field);
        f->name = fields_.names().nameOf(i);
        f->projection = std::unique_ptr<astrometry::SphericalCoords>(fields_.projections()[i]->duplicate());
        f->extent = fieldExtents[i];
        f->matchRadius = matchRadius;
        fields.emplace_back(std::move(f));
    }

    instruments.reserve(instruments_.size());
    for (int i = 0; i < instruments_.size(); i++) {
        std::unique_ptr<Instrument> inst(new Instrument(*instruments_[i]));
        instruments.emplace_back(std::move(inst));
    }

    // Set Exposures:
    exposures = std::move(exposures_.getExposuresVector());
}

void FoFClass::getExtensionInfo(long iextn, std::string &thisAffinity, int &exposureNumber,
                                int &instrumentNumber, int &fieldNumber, int &deviceNumber,
                                std::vector<bool> &isStar, std::vector<double> &vx, std::vector<double> &vy,
                                std::vector<long> &vid) {
    std::string filename;
    extensionTable.readCell(filename, "FILENAME", iextn);
    stripWhite(filename);
    int hduNumber;
    extensionTable.readCell(hduNumber, "EXTENSION", iextn);

    if (iextn % 10 == 0)
        std::cerr << "# Reading object catalog " << iextn << "/" << extensionTable.nrows() << " in "
                  << filename << " HDU #" << hduNumber << std::endl;
    img::FTable ft = FITS::FitsTable(filename, FITS::ReadOnly, hduNumber).extract();

    std::string idKey;
    extensionTable.readCell(idKey, "IDKEY", iextn);

    if (stringstuff::nocaseEqual(idKey, "_ROW")) {
        // Want the input table row number as ID column, so make such a column
        std::vector<long> vi(ft.nrows());
        for (int i = 0; i < vi.size(); i++) vi[i] = i;
        ft.addColumn(vi, "_ROW");
    }
    extensionTable.readCell(thisAffinity, "AFFINITY", iextn);
    // Get exposure number, instrument, field, and device numbers
    extensionTable.readCell(exposureNumber, "Exposure", iextn);
    instrumentNumber = exposures[exposureNumber]->instrument;
    fieldNumber = exposures[exposureNumber]->field;
    extensionTable.readCell(deviceNumber, "Device", iextn);

    // Now begin reading the input data

    // Filter the input rows and get the columns we want
    {
        std::string selectionExpression;
        extensionTable.readCell(selectionExpression, "SELECT", iextn);
        stripWhite(selectionExpression);
        if (!selectionExpression.empty()) ft.filterRows(selectionExpression);
    }

    if (ft.nrows() <= 0) {
        // No data to use in this extension.  Don't read anything.
        return;
    }

    std::string key;
    extensionTable.readCell(key, "XKEY", iextn);
    try {
        ft.readCells(vx, key);
    } catch (img::FTableError &m) {
        // Trap for using float column in source file instead of double:
        std::vector<float> vf;
        ft.readCells(vf, key);
        vx.resize(vf.size());
        for (int i = 0; i < vf.size(); i++) vx[i] = vf[i];
    }

    extensionTable.readCell(key, "YKEY", iextn);
    try {
        ft.readCells(vy, key);
    } catch (img::FTableError &m) {
        // Trap for using float column in source file instead of double:
        std::vector<float> vf;
        ft.readCells(vf, key);
        vy.resize(vf.size());
        for (int i = 0; i < vf.size(); i++) vy[i] = vf[i];
    }

    extensionTable.readCell(key, "IDKEY", iextn);
    try {
        //**/std::cerr << "Trying long read" << std::endl;
        ft.readCells(vid, key);
    } catch (img::FTableError &m) {
        // Trap for using int column in source file instead of long:
        std::vector<LONGLONG> vidint;
        //**/std::cerr << "Retrying int" << std::endl;
        ft.readCells(vidint, key);
        vid.reserve(vidint.size());
        vid.insert(vid.begin(), vidint.begin(), vidint.end());
    }

    {
        isStar = std::vector<bool>(ft.nrows(), true);
        std::string starExpression;
        extensionTable.readCell(starExpression, "STARSELECT", iextn);
        stripWhite(starExpression);
        if (!starExpression.empty()) {
            ft.evaluate(isStar, starExpression);
        }
    }
}

std::unique_ptr<astrometry::Wcs> FoFClass::getWCS(long iextn, int fieldNumber) {
    // Get the WCS for this extension
    std::unique_ptr<astrometry::Wcs> wcs = 0;

    std::string wcsin;
    extensionTable.readCell(wcsin, "WCSIN", iextn);
    std::list<std::string> wcssplit = stringstuff::split(wcsin, '@');
    std::unique_ptr<astrometry::PixelMapCollection> inputPmc = 0;
    if (wcssplit.size() == 2) {
        // @ sign signals that we are getting a map from a PMC
        std::string wcsName = wcssplit.front();
        std::string pmcFile = wcssplit.back();
        // Replace any white space in wcs name with underscore:
        wcsName = stringstuff::regexReplace("[[:space:]]+", "_", wcsName);

        if (pmcFile != currentPixelMapCollectionFileName) {
            // Need to open this PMC file
            // Try opening file as new serialized PMC
            std::ifstream wcsStream(pmcFile.c_str());
            if (!wcsStream.is_open()) {
                throw std::runtime_error("Could not open PixelMapCollection file <" + pmcFile + ">");
            }
            std::unique_ptr<astrometry::PixelMapCollection> tryPMC;
            tryPMC = std::unique_ptr<astrometry::PixelMapCollection>(new astrometry::PixelMapCollection);
            if (tryPMC->read(wcsStream)) {
                // Successfully found serialized file:
                inputPmc = std::move(tryPMC);
                currentPixelMapCollectionFileName = pmcFile;
            } else {
                throw std::runtime_error("Failed deserializing PixelMapCollection file <" + pmcFile + ">");
            }
        }  // Done reading the pmcFile

        // Read the wcs from the collection
        wcs = std::unique_ptr<astrometry::Wcs>(inputPmc->issueWcs(wcsName)->duplicate());
    } else if (stringstuff::nocaseEqual(wcsin, "_ICRS")) {
        // No mapping will be needed.  Do nothing.
        astrometry::IdentityMap identity;
        astrometry::SphericalICRS icrs;
        wcs = std::unique_ptr<astrometry::Wcs>(
                new astrometry::Wcs(&identity, icrs, "ICRS_degrees", WCS_UNIT));
    } else {
        // Assume that wcsin is a FITS WCS specification
        img::Header wcsHead;
        //	wcsin += "\nEND\n";  //*** hack to fix problem in astropy.io.fits
        std::istringstream wcsStream(wcsin);
        if (!(wcsStream >> wcsHead)) {
            throw std::runtime_error("Error reading WCS/header from WCSIN at extension " +
                                     std::to_string(iextn));
        }
        wcs = astrometry::readTPV(wcsHead);
        if (!wcs) {
            throw std::runtime_error("Failed reading TPV WCS from extension " + std::to_string(iextn));
        }
    }

    // Now replace WCSIN in table with a new serialized one,
    // and set projection to be the field coordinates
    if (wcs->getName() != "Identity") {
        std::string wcsDump;
        astrometry::PixelMapCollection pmc;
        pmc.learnWcs(*wcs);
        std::ostringstream oss;
        pmc.writeWcs(oss, wcs->getName());
        wcsDump = oss.str();
        wcs->reprojectTo(*fields[fieldNumber]->projection);
        extensionTable.writeCell(wcsDump, "WCSIN", iextn);
    }
    return wcs;
}

void FoFClass::reprojectWCS(astrometry::Wcs &wcs, int fieldNumber) {
    wcs.reprojectTo(*fields[fieldNumber]->projection);
}

void FoFClass::addCatalog(const astrometry::Wcs &wcs, std::string thisAffinity, int exposureNumber,
                          int fieldNumber, int instrumentNumber, int deviceNumber, long iextn,
                          std::vector<bool> isStar, std::vector<double> vx, std::vector<double> vy,
                          std::vector<long> vid) {
    PointCat &starCatalog = fields[fieldNumber]->catalogFor(stellarAffinity);
    PointCat &galaxyCatalog = (stringstuff::nocaseEqual(stellarAffinity, thisAffinity)
                                       ? starCatalog
                                       : fields[fieldNumber]->catalogFor(thisAffinity));

    // Now loops over objects in the catalog
    for (int iObj = 0; iObj < vx.size(); iObj++) {
        double xpix = vx[iObj];
        double ypix = vy[iObj];
        // Expand bounds of device:
        if (instrumentNumber >= 0)
            (*instruments[instrumentNumber]).domains[deviceNumber] += Position<double>(xpix, ypix);

        //  maps coords to field's tangent plane
        double xw, yw;
        if (wcs.getName() != "Identity") {
            Assert(wcs);
            wcs.toWorld(xpix, ypix, xw, yw);
        } else {
            // Already have RA, Dec, just project them
            // Assume they are in our std units (see Units.h)
            fields[fieldNumber]->projection->convertFrom(
                    astrometry::SphericalICRS(vx[iObj] * WCS_UNIT, vy[iObj] * WCS_UNIT));
            fields[fieldNumber]->projection->getLonLat(xw, yw);
            // Matching program will use our std unit (Units.h), not radians:
            xw /= WCS_UNIT;
            yw /= WCS_UNIT;
        }

        allPoints.push_back(Point(xw, yw, iextn, vid[iObj], exposureNumber));

        // Now we match!!!!
        if (isStar[iObj])
            starCatalog.add(allPoints.back());
        else
            galaxyCatalog.add(allPoints.back());
    }  // end object loop
}

void FoFClass::writeMatches(std::string outCatalogName, int minMatches, bool allowSelfMatches) {
    // Now write all the match catalogs
    long matchCount = 0;
    long pointCount = 0;
    //  Loop over fields
    for (int iField = 0; iField < fields.size(); iField++) {
        //  loop over affinities per field
        for (Field::CatMap::iterator iAffinity = fields[iField]->catalogs.begin();
             iAffinity != fields[iField]->catalogs.end(); ++iAffinity) {
            const PointCat &pcat = *iAffinity->second;
            std::cerr << "Catalog for field " << fields[iField]->name << " affinity " << iAffinity->first
                      << " with " << pcat.size() << " matches " << std::endl;
            if (pcat.empty()) continue;
            FITS::FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
            ft.setName("MatchCatalog");
            img::FTable ff = ft.use();
            ff.header()->replace("Field", fields[iField]->name, "Field name");
            ff.header()->replace("FieldNum", iField, "Field number");
            ff.header()->replace("Affinity", iAffinity->first, "Affinity name");

            std::vector<int> sequence;
            std::vector<long> extn;
            std::vector<long> obj;
            long matches = 0;
            // Now loop through matches in this catalog
            for (PointCat::const_iterator j = pcat.begin(); j != pcat.end(); ++j) {
                // Skip any Match that is below minimum match size
                if ((*j)->size() < minMatches) continue;
                bool selfMatch = false;
                if (!allowSelfMatches) {
                    std::set<long> itsExposures;
                    for (std::list<const Point *>::const_iterator k = (*j)->begin(); k != (*j)->end(); ++k)
                        if (!itsExposures.insert((*k)->exposureNumber).second) {
                            selfMatch = true;
                            break;
                        }
                }
                if (selfMatch) continue;

                // Add elements of this match to the output std::vectors
                int seq = 0;
                ++matches;
                ++matchCount;
                for (std::list<const Point *>::const_iterator k = (*j)->begin(); k != (*j)->end();
                     ++k, ++seq) {
                    sequence.push_back(seq);
                    extn.push_back((*k)->extensionNumber);
                    obj.push_back((*k)->objectNumber);
                }
            }  // end match loop
            std::cerr << "...Kept " << matches << " matches with " << sequence.size() << " points."
                      << std::endl;
            pointCount += sequence.size();
            ff.addColumn(sequence, "SequenceNumber");
            ff.addColumn(extn, "Extension");
            ff.addColumn(obj, "Object");
        }  // end Affinity loop
    }      // end Field loop

    std::cerr << "Total of " << matchCount << " matches with " << pointCount << " points." << std::endl;
}

void FoFClass::sortMatches(int fieldNumber, int minMatches, bool allowSelfMatches) {
    // Now write all the match catalogs
    long matchCount = 0;
    long pointCount = 0;

    //  loop over affinities per field
    for (Field::CatMap::iterator iAffinity = fields[fieldNumber]->catalogs.begin();
         iAffinity != fields[fieldNumber]->catalogs.end(); ++iAffinity) {
        const PointCat &pcat = *iAffinity->second;
        std::cerr << "Catalog for field " << fields[fieldNumber]->name << " affinity " << iAffinity->first
                  << " with " <<  pcat.size() << " matches " << std::endl;

        if (pcat.empty()) continue;

        long matches = 0;
        // Now loop through matches in this catalog
        for (PointCat::const_iterator j = pcat.begin(); j != pcat.end(); ++j) {
            // Skip any Match that is below minimum match size

            if ((*j)->size() < minMatches) continue;
            bool selfMatch = false;
            if (!allowSelfMatches) {
                std::set<long> itsExposures;
                for (std::list<const Point *>::const_iterator k = (*j)->begin(); k != (*j)->end(); ++k)
                    if (!itsExposures.insert((*k)->exposureNumber).second) {
                        selfMatch = true;
                        break;
                    }
            }

            if (selfMatch) continue;

            // Add elements of this match to the output std::vectors
            int seq = 0;
            ++matches;
            ++matchCount;
            for (std::list<const Point *>::const_iterator k = (*j)->begin(); k != (*j)->end(); ++k, ++seq) {
                sequence.push_back(seq);
                extn.push_back((*k)->extensionNumber);
                obj.push_back((*k)->objectNumber);
            }
        }  // end match loop

        std::cerr << "...Kept " << matches << " matches with "
                  <<  sequence.size() << " points." << std::endl;
        pointCount += sequence.size();
    }  // end Affinity loop

    std::cerr << "Total of " <<  matchCount << " matches with " <<  pointCount
              << " points." << std::endl;
}
