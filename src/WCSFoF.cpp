// Program to match catalogs based on world coordinates.

#include "Std.h"
#include "WCSFoFRoutine.h"
#include "Pset.h"
#include "FitSubroutines.h"

#include <iostream>

std::string usage =
        "Match objects using FoF algorithm on world coordinate system\n"
        "WCSFoF <configfile> [parameter file] [parameter file...]\n"
        "   [-parameter[=]value...]\n"
        "      <configfile>: FITS file with binary tables describing input catalogs\n"
        "      Program parameters specified as command-line options or read from\n"
        "          parameter file(s) specified on cmd line";

int main(int argc, char *argv[]) {
    // Read in parameters
    double matchRadius;
    std::string outCatalogName;
    int minMatches;
    bool allowSelfMatches;

    Pset parameters;
    {
        const int def = PsetMember::hasDefault;
        const int low = PsetMember::hasLowerBound;
        const int up = PsetMember::hasUpperBound;
        const int lowopen = low | PsetMember::openLowerBound;
        const int upopen = up | PsetMember::openUpperBound;

        parameters.addMember("matchRadius", &matchRadius, def | lowopen, "Matching tolerance (arcsec)", 1.,
                             0.);
        parameters.addMember("outName", &outCatalogName, def, "filename for FITS output catalog",
                             "match.fof");
        parameters.addMember("minMatch", &minMatches, def | low,
                             "Minimum number of detections for usable match", 2, 2);
        parameters.addMember("selfMatch", &allowSelfMatches, def,
                             "Retain matches that have 2 elements from same exposure?", false);
    }

    ////////////////////////////////////////////////
    // Read parameters
    ////////////////////////////////////////////////

    // Read all the command-line and parameter-file program parameters
    processParameters(parameters, usage, 1, argc, argv);
    std::string configFile = argv[1];

    // Convert matching radius to our system units for world coords
    matchRadius *= ARCSEC / astrometry::WCS_UNIT;

    // Set of Friends-of-Friends class
    FoFClass fofclass;

    try {
        ////////////////////////////////////////////////
        // Read in field table and save away orientations of tangent plane systems for each
        ////////////////////////////////////////////////
        img::FTable fieldTable;
        try {
            fieldTable = FITS::FitsTable(configFile, FITS::ReadOnly, "Fields").extract();
        } catch (std::runtime_error &e) {
            quit(e, 1);
        }

        std::vector<std::unique_ptr<Field>> fields;
        fields.reserve(fieldTable.nrows());

        for (int ifield = 0; ifield < fieldTable.nrows(); ifield++) {
            std::string name;
            double ra;
            double dec;
            double extent;
            fieldTable.readCell(name, "NAME", ifield);
            fieldTable.readCell(ra, "RA", ifield);
            fieldTable.readCell(dec, "DEC", ifield);
            fieldTable.readCell(extent, "RADIUS", ifield);

            std::unique_ptr<Field> f = std::unique_ptr<Field>(new Field);
            f->name = name;
            astrometry::Orientation orient(astrometry::SphericalICRS(ra * astrometry::WCS_UNIT, dec * astrometry::WCS_UNIT));
            f->projection = std::unique_ptr<astrometry::SphericalCoords>(new astrometry::Gnomonic(orient));
            f->extent = extent;
            f->matchRadius = matchRadius;
            fields.emplace_back(std::move(f));
        }  // Done reading fields
        Assert(!fields.empty());

        // Now read in all the instrument tables
        std::list<int> instrumentHDUs;
        {
            // Find which extensions are instrument tables
            FITS::FitsFile ff(configFile);
            for (int i = 1; i < ff.HDUCount(); i++) {
                FITS::Hdu h(configFile, FITS::HDUAny, i);
                if (stringstuff::nocaseEqual(h.getName(), "Instrument")) instrumentHDUs.push_back(i);
            }
        }

        std::vector<img::FTable> instrumentTables;
        instrumentTables.reserve(instrumentHDUs.size());
        std::vector<std::unique_ptr<Instrument>> instruments;
        instruments.reserve(instrumentHDUs.size());
        for (std::list<int>::const_iterator i = instrumentHDUs.begin(); i != instrumentHDUs.end(); ++i) {
            img::FTable ft = FITS::FitsTable(configFile, FITS::ReadOnly, *i).extract();
            std::string instrumentName;
            int instrumentNumber;
            if (!ft.header()->getValue("Name", instrumentName) ||
                !ft.header()->getValue("Number", instrumentNumber)) {
                std::cerr << "Could not read name and/or number of instrument at extension " << *i << std::endl;
            }
            // Get rid of spaces in the name
            spaceReplace(instrumentName);
            Assert(instrumentNumber < instrumentTables.size() && instrumentNumber >= 0);
            ft.header()->setValue("Name", instrumentName);

            // Now save the table and make an Instrument structure
            instrumentTables.push_back(ft);
            std::unique_ptr<Instrument> inst(new Instrument(instrumentName));
            for (int i = 0; i < ft.nrows(); i++) {
                std::string devname;
                ft.readCell(devname, "Name", i);
                double xmin, xmax, ymin, ymax;
                ft.readCell(xmin, "XMin", i);
                ft.readCell(xmax, "XMax", i);
                ft.readCell(ymin, "YMin", i);
                ft.readCell(ymax, "YMax", i);
                Bounds<double> devBounds;
                if (xmin != 0. || xmax != 0. || ymin != 0. || ymax != 0.) {
                    // initialize bounds if there are some, otherwise leave undefined
                    devBounds.setXMin(xmin);
                    devBounds.setXMax(xmax);
                    devBounds.setYMin(ymin);
                    devBounds.setYMax(ymax);
                }
                inst->addDevice(devname, devBounds);
            }
            instruments.emplace_back(std::move(inst));
        }
        // Check that all Instruments were read
        for (int i = 0; i < instruments.size(); i++)
            if (!instruments[i]) {
                std::cerr << "Failed to read instrument number " << i << std::endl;
                exit(1);
            }

        // Now get information on exposures
        img::FTable exposureTable;
        try {
            exposureTable = FITS::FitsTable(configFile, FITS::ReadOnly, "Exposures").extract();
        } catch (std::runtime_error &e) {
            quit(e, 1);
        }

        std::vector<std::unique_ptr<Exposure>> exposures;
        for (int i = 0; i < exposureTable.nrows(); i++) {
            double ra, dec;
            std::string name;
            int field, instrument;
            exposureTable.readCell(name, "Name", i);
            exposureTable.readCell(ra, "RA", i);
            exposureTable.readCell(dec, "Dec", i);
            exposureTable.readCell(field, "FieldNumber", i);
            exposureTable.readCell(instrument, "InstrumentNumber", i);
            astrometry::SphericalICRS pointing;
            pointing.setRADec(ra * astrometry::WCS_UNIT, dec * astrometry::WCS_UNIT);
            std::unique_ptr<Exposure> e(new Exposure(name, pointing));
            e->field = field;
            e->instrument = instrument;
            exposures.emplace_back(std::move(e));
        }
        Assert(exposures.size() == exposureTable.nrows());

        // Table of information about every catalog file to read:
        img::FTable extensionTable;
        try {
            extensionTable = FITS::FitsTable(configFile, FITS::ReadOnly, "Extensions").extract();
        } catch (std::runtime_error &e) {
            quit(e, 1);
        }

        fofclass.fields = std::move(fields);
        fofclass.instruments = std::move(instruments);
        fofclass.exposures = std::move(exposures);
        fofclass.extensionTable = extensionTable;

        // Loop over input catalogs
        for (long iextn = 0; iextn < extensionTable.nrows(); iextn++) {
            std::string thisAffinity;
            int exposureNumber;
            int instrumentNumber;
            int fieldNumber;
            int deviceNumber;
            std::vector<bool> isStar;
            std::vector<double> vx;
            std::vector<double> vy;
            std::vector<long> vid;
            fofclass.getExtensionInfo(iextn, thisAffinity, exposureNumber, instrumentNumber, fieldNumber,
                                      deviceNumber, isStar, vx, vy, vid);
            std::unique_ptr<astrometry::Wcs> wcs = fofclass.getWCS(iextn, fieldNumber);
            fofclass.addCatalog(*wcs, thisAffinity, exposureNumber, fieldNumber, instrumentNumber,
                                deviceNumber, iextn, isStar, vx, vy, vid);
        }

        //  Write all of our tables to output file
        {
            // Fields
            FITS::FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::OverwriteFile, "Fields");
            ft.copy(fieldTable);
        }
        {
            // Exposures
            FITS::FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Exposures");
            ft.copy(exposureTable);
        }

        for (int i = 0; i < fofclass.instruments.size(); i++) {
            // Instrument tables
            // Update the device bounds in the table first
            const Instrument &inst = *fofclass.instruments[i];
            std::vector<double> vxmin(inst.nDevices);
            std::vector<double> vxmax(inst.nDevices);
            std::vector<double> vymin(inst.nDevices);
            std::vector<double> vymax(inst.nDevices);
            for (int j = 0; j < inst.nDevices; j++) {
                vxmin[j] = floor(inst.domains[j].getXMin());
                vxmax[j] = ceil(inst.domains[j].getXMax());
                vymin[j] = floor(inst.domains[j].getYMin());
                vymax[j] = ceil(inst.domains[j].getYMax());
            }
            instrumentTables[i].writeCells(vxmin, "XMin");
            instrumentTables[i].writeCells(vxmax, "XMax");
            instrumentTables[i].writeCells(vymin, "YMin");
            instrumentTables[i].writeCells(vymax, "YMax");

            FITS::FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
            ft.setName("Instrument");
            ft.setVersion(i + 1);  // might use this later, or put into extension name???
            ft.copy(instrumentTables[i]);
        }

        {
            // Extension table
            FITS::FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Extensions");
            ft.copy(extensionTable);
        }

        fofclass.writeMatches(outCatalogName, minMatches, allowSelfMatches);

    } catch (std::runtime_error &m) {
        quit(m, 1);
    }
}
