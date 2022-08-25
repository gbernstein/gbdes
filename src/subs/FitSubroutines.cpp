// Constants and subroutines shared by the astrometric and photometric matching/fitting codes
#include "FitSubroutines.h"
#include "StringStuff.h"
#include "Bounds.h"
#include <list>
#include "Pset.h"
#include "PhotoTemplate.h"
#include "PhotoPiecewise.h"
#include "TemplateMap.h"
#include "PiecewiseMap.h"
#include "PixelMapCollection.h"
#include "PhotoMapCollection.h"
#include "FitsTable.h"
#include "Match.h"
#include "PhotoMatch.h"
#include "Random.h"
#include "Stopwatch.h"

using astrometry::AstrometryError;
using astrometry::PARALLAX_UNIT;
using astrometry::PM_UNIT;
using astrometry::PMDetection;
using astrometry::PMMatch;
using astrometry::RESIDUAL_UNIT;
using astrometry::TDB_UNIT;
using astrometry::WCS_UNIT;
using photometry::MMAG;

// A helper function that strips white space from front/back of a string and replaces
// internal white space with underscores:
void spaceReplace(string &s) {
    stripWhite(s);
    s = regexReplace("[[:space:]]+", "_", s);
}

// Another helper function to split up a string into a list of whitespace-trimmed strings.
// Get rid of any null ones.
list<string> splitArgument(string input, const char listSeperator) {
    list<string> out = split(input, listSeperator);
    for (list<string>::iterator i = out.begin(); i != out.end();) {
        stripWhite(*i);
        if (i->empty()) {
            i = out.erase(i);
        } else {
            ++i;
        }
    }
    return out;
}

// Read parameters from files on command line and from std input.
// First nRequiredArgs must be present and are not names of parameter files.
// Print the usage string (to cerr) and the default parameters (to cout) if
// there are not enough parameters.
// Then dumps final parameter values to cout.
void processParameters(Pset &parameters, const string &usage, int nRequiredArgs, int argc, char *argv[]) {
    int positionalArguments;
    try {
        // First read the command-line arguments so we know how many positional
        // arguments precede them.
        positionalArguments = parameters.setFromArguments(argc, argv);
    } catch (std::runtime_error &m) {
        // An error here might indicate someone entered "-help" or something
        cerr << usage << endl;
        cout << "#---- Parameter defaults: ----" << endl;
        parameters.dump(cerr);
        quit(m, 1);
    }
    parameters.setDefault();
    if (positionalArguments < nRequiredArgs + 1) {
        cerr << usage << endl;
        cout << "#---- Parameter defaults: ----" << endl;
        parameters.dump(cerr);
        throw std::runtime_error(
                "Number of positional arguments is less than number of required arguments + 1");
    }

    for (int i = nRequiredArgs + 1; i < positionalArguments; i++) {
        // Open & read all specified input files
        ifstream ifs(argv[i]);
        if (!ifs) {
            std::string argMessage(argv[i], 1);
            throw std::runtime_error("Can't open parameter file " + argMessage + "\n" + usage);
        }
        try {
            parameters.setStream(ifs);
        } catch (std::runtime_error &m) {
            cerr << "In file " << argv[i] << ":" << endl;
            quit(m, 1);
        }
    }
    // And now re-read the command-line arguments so they take precedence
    parameters.setFromArguments(argc, argv);

    cout << "#" << stringstuff::taggedCommandLine(argc, argv) << endl;
    parameters.dump(cout);
}

// Take a string with this format:
// <regex>=<replace> [, <regex>=<replace>] ...
// and parse it into a stringstuff::RegexReplacement object that will execute
// the regular expression replacements on an input string.
RegexReplacements parseTranslator(string specString, string errorDescription) {
    list<string> ls = split(specString, DefaultListSeperator);
    RegexReplacements translator;
    for (auto &i : ls) {
        if (i.empty()) continue;
        list<string> ls2 = split(i, '=');
        if (ls2.size() != 2) {
            std::string iMessage(i, 1);
            throw std::runtime_error(errorDescription + " has bad translation spec: " + iMessage);
        }
        string regex = ls2.front();
        string replacement = ls2.back();
        stripWhite(regex);
        stripWhite(replacement);
        translator.addRegex(regex, replacement);
    }
    return translator;
}

ExtensionObjectSet::ExtensionObjectSet(string filename) {
    if (filename.empty()) return;
    ifstream ifs(filename.c_str());
    if (!ifs) {
        throw std::runtime_error("Could not open Extension/Object pair file " + filename);
    }
    string buffer;
    while (stringstuff::getlineNoComment(ifs, buffer)) {
        int extensionNumber;
        long objectNumber;
        istringstream iss(buffer);
        if (!(iss >> extensionNumber >> objectNumber)) {
            throw std::runtime_error("Bad Extension/Object pair in " + filename + ": <" + buffer + ">");
        }
        pairs.insert(EOPair(extensionNumber, objectNumber));
    }
}

bool ExtensionObjectSet::operator()(int extensionNumber, long objectNumber) const {
    return pairs.count(EOPair(extensionNumber, objectNumber)) > 0;
}

void loadPixelMapParser() {
    // put all new kinds of PixelMap atoms here:
    astrometry::PixelMapCollection::registerMapType<astrometry::TemplateMap>();
    astrometry::PixelMapCollection::registerMapType<astrometry::PiecewiseMap>();
}

void loadPhotoMapParser() {
    // put all new kinds of PhotoMap atoms here:
    photometry::PhotoMapCollection::registerMapType<photometry::PhotoTemplate>();
    photometry::PhotoMapCollection::registerMapType<photometry::PhotoPiecewise>();
}

int elementNumber(string &key) {
    int out = -1;
    list<string> l =
            stringstuff::split(stringstuff::regexReplace("^(.*)\\[([0-9]*)\\]$", "\\1;\\2", key), ';');
    if (l.size() == 2) {
        key = l.front();
        stripWhite(key);
        out = atoi(l.back().c_str());
    }
    return out;
}

bool isDouble(img::FTable f, string key, int elementNumber) {
    bool out = true;
    try {
        if (elementNumber < 0) {
            double d;
            f.readCell(d, key, 0);
        } else {
            vector<double> d;
            f.readCell(d, key, 0);
        }
    } catch (img::FTableError &e) {
        out = false;
    }
    return out;
}

double getTableDouble(img::FTable f, string key, int elementNumber, bool isDouble, long irow) {
    double out;
    if (isDouble) {
        if (elementNumber < 0) {
            f.readCell(out, key, irow);
        } else {
            vector<double> v;
            f.readCell(v, key, irow);
            if (elementNumber >= v.size()) {
                throw std::runtime_error("Requested element " + std::to_string(elementNumber) +
                                         " of array at key " + key + " but size of array is " +
                                         std::to_string(v.size()));
            }
            out = v[elementNumber];
        }
    } else {
        if (elementNumber < 0) {
            float fl;
            f.readCell(fl, key, irow);
            out = fl;
        } else {
            vector<float> v;
            f.readCell(v, key, irow);
            if (elementNumber >= v.size()) {
                throw std::runtime_error("Requested element " + std::to_string(elementNumber) +
                                         " of array at key " + key + " but size of array is " +
                                         std::to_string(v.size()));
            }
            out = v[elementNumber];
        }
    }
    return out;
}

// This function is used to find degeneracies between exposures and device maps.
// Start with list of free & fixed devices as initial degen/ok, same for exposures.
// Will consider as "ok" any device used in an "ok" exposure and vice-versa.
// The last argument says which exposure/device pairs are used together.
void findDegeneracies(set<int> &degenerateDevices, set<int> &okDevices, set<int> &degenerateExposures,
                      set<int> &okExposures, const vector<set<int>> &exposuresUsingDevice) {
    // Device/exposures with degeneracy broken this round
    auto nowOkDevices = okDevices;
    auto nowOkExposures = okExposures;
    while (!(nowOkDevices.empty() && nowOkExposures.empty())) {
        // Mark as ok any exposures using an ok device
        for (auto iDev : nowOkDevices) {
            for (auto iExpo : exposuresUsingDevice[iDev]) {
                if (degenerateExposures.erase(iExpo) > 0) {
                    // Newly non-degenerate exposure.  Change category
                    nowOkExposures.insert(iExpo);
                    okExposures.insert(iExpo);
                }
            }
        }
        nowOkDevices.clear();

        // Now mark as ok any device that is used by a non-degenerate exposure
        for (auto iDev : degenerateDevices) {
            for (auto iExpo : exposuresUsingDevice[iDev]) {
                if (nowOkExposures.count(iExpo) > 0) {
                    // Mark device as not degenerate after all
                    nowOkDevices.insert(iDev);
                    break;
                }
            }
        }
        nowOkExposures.clear();
        for (auto iDev : nowOkDevices) {
            degenerateDevices.erase(iDev);
            okDevices.insert(iDev);
        }
    }
    return;
}

// Figure out which extensions of the FITS file inputTables
// are Instrument or MatchCatalog extensions.
void inventoryFitsTables(string inputTables, vector<int> &instrumentHDUs, vector<int> &catalogHDUs) {
    FITS::FitsFile ff(inputTables);
    for (int i = 1; i < ff.HDUCount(); i++) {
        FITS::Hdu h(inputTables, FITS::HDUAny, i);
        if (stringstuff::nocaseEqual(h.getName(), "Instrument"))
            instrumentHDUs.push_back(i);
        else if (stringstuff::nocaseEqual(h.getName(), "MatchCatalog"))
            catalogHDUs.push_back(i);
    }
}

Fields::Fields(vector<string> names, vector<double> ra, vector<double> dec, vector<double> epochs)
        : _names(), _projections(), _epochs(std::move(epochs)) {
    // TODO: check that sizes are consistent, raise if they are not.
    _projections.reserve(names.size());
    for (std::size_t i = 0; i != names.size(); ++i) {
        spaceReplace(names[i]);
        _names.append(std::move(names[i]));
        astrometry::Orientation orient(
                astrometry::SphericalICRS(std::move(ra[i]) * WCS_UNIT, std::move(dec[i]) * WCS_UNIT));
        _projections.emplace_back(new astrometry::Gnomonic(orient));
    }
}

// Read info from fields table, write to a new output file
Fields Fields::read(string inputTables, string outCatalog, double defaultEpoch) {
    const double MINIMUM_EPOCH = 1900.;

    FITS::FitsTable in(inputTables, FITS::ReadOnly, "Fields");
    FITS::FitsTable out(outCatalog, FITS::ReadWrite + FITS::OverwriteFile + FITS::Create, "Fields");
    img::FTable ft = in.extract();
    out.adopt(ft);
    vector<double> ra;
    vector<double> dec;
    vector<string> name;
    ft.readCells(name, "Name");
    ft.readCells(ra, "RA");
    ft.readCells(dec, "Dec");
    vector<double> epochs;
    if (ft.hasColumn("PM_EPOCH")) {
        ft.readCells(epochs, "PM_EPOCH");
        // Set reference epochs for each field to defaultEpoch if
        // they have a signal value like 0 or -1
        for (int i = 0; i < epochs.size(); i++) {
            if (epochs[i] < MINIMUM_EPOCH && defaultEpoch > MINIMUM_EPOCH) {
                // Update to reference epoch if it's sensible
                epochs[i] = defaultEpoch;
                ft.writeCell(defaultEpoch, "PM_EPOCH", i);
            }
        }
    } else {
        // Assign default epoch to every field
        epochs.clear();
        epochs.resize(name.size(), defaultEpoch);
        // And save back to output if this is a sensible default
        if (defaultEpoch > MINIMUM_EPOCH) {
            ft.addColumn(epochs, "PM_EPOCH");
        }
    }
    return Fields(std::move(name), std::move(ra), std::move(dec), std::move(epochs));
}

// Read in all the instrument extensions and their device info from input
// FITS file, save useful ones and write to output FITS file.
// The useInstrumentList entries are regexes, empty means use all.
// The final bool argument is set true if we have already created
// the outCatalog FITS file.
vector<unique_ptr<Instrument>> readInstruments(const vector<int> &instrumentHDUs,
                                               const list<string> &useInstrumentList, string inputTables,
                                               string outCatalog, bool &outputCatalogAlreadyOpen) {
    vector<unique_ptr<Instrument>> instruments;
    for (int iextn : instrumentHDUs) {
        FITS::FitsTable ft(inputTables, FITS::ReadOnly, iextn);
        Assert(stringstuff::nocaseEqual(ft.getName(), "Instrument"));

        string instrumentName;
        int instrumentNumber;
        if (!ft.header()->getValue("Name", instrumentName) ||
            !ft.header()->getValue("Number", instrumentNumber)) {
            cerr << "Could not read name and/or number of instrument at extension " << iextn << endl;
        }
        if (instrumentNumber >= instruments.size()) {
            int oldsize = instruments.size();
            instruments.resize(instrumentNumber + 1);
            for (; oldsize < instruments.size(); oldsize++) instruments[oldsize] = 0;
        }
        spaceReplace(instrumentName);

        // Send every instrument extension to output file
        FITS::Flags outFlags = FITS::ReadWrite + FITS::Create;
        if (!outputCatalogAlreadyOpen) {
            outFlags = outFlags + FITS::OverwriteFile;
            outputCatalogAlreadyOpen = true;
        }
        FITS::FitsTable out(outCatalog, outFlags, -1);
        out.setName("Instrument");
        img::FTable ff = ft.extract();
        out.adopt(ff);

        if (regexMatchAny(useInstrumentList, instrumentName)) {
            // This is an instrument we will use.  Make an Instrument instance
            auto instptr = unique_ptr<Instrument>(new Instrument(instrumentName));
            string band;
            if (!ff.header()->getValue("Band", band)) {
                instptr->band = instptr->name;  // Use instrument name for BAND if not found
            } else {
                spaceReplace(band);
                instptr->band = band;
            }

            vector<string> devnames;
            vector<double> vxmin;
            vector<double> vxmax;
            vector<double> vymin;
            vector<double> vymax;
            ff.readCells(devnames, "Name");
            ff.readCells(vxmin, "XMin");
            ff.readCells(vxmax, "XMax");
            ff.readCells(vymin, "YMin");
            ff.readCells(vymax, "YMax");
            for (int j = 0; j < devnames.size(); j++) {
                spaceReplace(devnames[j]);
                instptr->addDevice(devnames[j], Bounds<double>(vxmin[j], vxmax[j], vymin[j], vymax[j]));
            }
            instruments[instrumentNumber] = std::move(instptr);
        }  // Done reading an instrument.
    }

    return instruments;
}

// Read the Exposure table into an array.
vector<unique_ptr<Exposure>> readExposures(const vector<unique_ptr<Instrument>> &instruments,
                                           const vector<double> &fieldEpochs,
                                           vector<int> &exposureColorPriorities,
                                           const list<string> &useColorList, string inputTables,
                                           string outCatalog, const list<string> &skipExposureList,
                                           bool useReferenceExposures, bool &outputCatalogAlreadyOpen) {
    // Read the exposure table
    FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Exposures");
    FITS::Flags outFlags = FITS::ReadWrite + FITS::Create;
    if (!outputCatalogAlreadyOpen) {
        outFlags = outFlags + FITS::OverwriteFile;
        outputCatalogAlreadyOpen = true;
    }
    FITS::FitsTable out(outCatalog, outFlags, "Exposures");
    img::FTable ff = ft.extract();
    out.adopt(ff);
    vector<string> names;
    vector<double> ra;
    vector<double> dec;
    vector<int> fieldNumber;
    vector<int> instrumentNumber;
    vector<double> airmass;
    vector<double> exptime;
    vector<double> mjd;
    vector<double> apcorr;
    vector<string> epoch;
    vector<vector<double>> observatory;
    vector<vector<double>> astrometricCovariance;
    vector<double> photometricVariance;
    vector<double> weight;
    vector<double> magWeight;
    vector<double> syserr_mmag;  // Photometric systematic error, mmag
    vector<double> syserr_mas;   // astrometric sys err, mas (if circular)
    vector<double> syserr_xx;    // astrometric sys error ellipse (mas^2)
    vector<double> syserr_yy;
    vector<double> syserr_xy;
    vector<double> observatory_x;  // Observatory barycentric ICRS posn (AU)
    vector<double> observatory_y;
    vector<double> observatory_z;

    ff.readCells(names, "Name");
    ff.readCells(ra, "RA");
    ff.readCells(dec, "Dec");
    ff.readCells(fieldNumber, "fieldNumber");
    ff.readCells(instrumentNumber, "InstrumentNumber");
    ff.readCells(airmass, "Airmass");
    ff.readCells(exptime, "Exptime");

    // Columns that might not always be present:
    try {
        ff.readCells(mjd, "MJD");
    } catch (img::FTableError &e) {
        mjd.clear();
    }
    try {
        ff.readCells(epoch, "Epoch");
    } catch (img::FTableError &e) {
        epoch.clear();
    }
    try {
        ff.readCells(apcorr, "apcorr");
    } catch (img::FTableError &e) {
        apcorr.clear();
    }
    try {
        ff.readCells(weight, "weight");
    } catch (img::FTableError &e) {
        weight.clear();
    }
    try {
        ff.readCells(magWeight, "magWeight");
    } catch (img::FTableError &e) {
        magWeight.clear();
    }
    try {
        ff.readCells(syserr_mmag, "syserrmmag");
        photometricVariance.resize(syserr_mmag.size());
        for (int i = 0; i < syserr_mmag.size(); i++) photometricVariance[i] = pow(syserr_mmag[i] * MMAG, 2.);
        syserr_mmag.clear();
    } catch (img::FTableError &e) {
        syserr_mmag.clear();
        photometricVariance.clear();
    }
    try {
        ff.readCells(syserr_mas, "syserrmas");
        // Convert this syserr to a diagonal matrix
        vector<double> vv(3);
        astrometricCovariance.clear();
        for (int i = 0; i < syserr_mas.size(); i++) {
            double stdev = syserr_mas[i] * astrometry::RESIDUAL_UNIT / astrometry::WCS_UNIT;
            vv[0] = vv[1] = stdev * stdev;
            vv[2] = 0.;
            astrometricCovariance.push_back(vv);
        }
    } catch (img::FTableError &e) {
        astrometricCovariance.clear();
    }
    try {
        ff.readCells(syserr_xx, "syserrxx");
        ff.readCells(syserr_yy, "syserryy");
        ff.readCells(syserr_xy, "syserrxy");
        // If all three columns exist, replace
        // any existing astrometricCovariance with this one if
        // it is nonzero.
        bool fresh = astrometricCovariance.empty();
        const double factor = pow(astrometry::RESIDUAL_UNIT / astrometry::WCS_UNIT, 2.);
        vector<double> vv(3);
        if (fresh) {
            astrometricCovariance.resize(syserr_xx.size());
        }
        for (int i = 0; i < syserr_xx.size(); i++) {
            vv[0] = syserr_xx[i] * factor;
            vv[1] = syserr_yy[i] * factor;
            vv[2] = syserr_xy[i] * factor;
            bool anyInfo = vv[0] > 0. || vv[1] > 0.;
            if (anyInfo) {
                astrometricCovariance[i] = vv;
            } else if (fresh) {
                // Place zeros if there was no circular error
                vv[0] = vv[1] = vv[2] = 0.;
                astrometricCovariance[i] = vv;
            } else {
                // Do nothing, keep pre-existing circular error
            }
        }
        syserr_xx.clear();
        syserr_yy.clear();
        syserr_xy.clear();
    } catch (img::FTableError &e) {
        syserr_xx.clear();
        syserr_yy.clear();
        syserr_xy.clear();
    }
    try {
        ff.readCells(observatory_x, "obsx");
        ff.readCells(observatory_y, "obsy");
        ff.readCells(observatory_z, "obsz");
        // Have all three columns if we make it here.
        observatory.resize(observatory_x.size());
        vector<double> vv(3);
        for (int i = 0; i < observatory_x.size(); i++) {
            vv[0] = observatory_x[i];
            vv[1] = observatory_y[i];
            vv[2] = observatory_z[i];
            observatory[i] = vv;
        }
        observatory_x.clear();
        observatory_y.clear();
        observatory_z.clear();
    } catch (img::FTableError &e) {
        observatory_x.clear();
        observatory_y.clear();
        observatory_z.clear();
        observatory.clear();
    }

    // Initialize our output arrays to not-in-use values
    vector<unique_ptr<Exposure>> exposures(names.size());
    exposureColorPriorities = vector<int>(names.size(), -1);
    for (int i = 0; i < names.size(); i++) {
        spaceReplace(names[i]);

        if (regexMatchAny(skipExposureList, names[i])) continue;  // do not use exposures on the list.

        // See if this exposure name matches any of the color exposures
        int priority = 0;
        for (auto j = useColorList.begin(); j != useColorList.end(); ++j, ++priority) {
            if (stringstuff::regexMatch(*j, names[i])) {
                // Yes: give this exposure a priority according to order of the exposure
                // it first matches.
                exposureColorPriorities[i] = priority;
                break;  // First priority is given to lower numbers.
            }
        }
        // Are we going to want use data this exposure?
        // Yes, if it's a reference exposure and we are using references
        bool useThisExposure =
                (instrumentNumber[i] == REF_INSTRUMENT || instrumentNumber[i] == PM_INSTRUMENT) &&
                useReferenceExposures;
        // or if it's a normal exposure using an instrument that has been included
        useThisExposure = useThisExposure || (instrumentNumber[i] >= 0 && instruments[instrumentNumber[i]]);

        if (useThisExposure) {
            // The projection we will use for this exposure:
            astrometry::Gnomonic gn(
                    astrometry::Orientation(astrometry::SphericalICRS(ra[i] * WCS_UNIT, dec[i] * WCS_UNIT)));
            unique_ptr<Exposure> expo(new Exposure(names[i], gn));
            expo->field = fieldNumber[i];
            expo->instrument = instrumentNumber[i];
            expo->airmass = airmass[i];
            expo->exptime = exptime[i];
            if (!mjd.empty()) {
                expo->mjd = mjd[i];
                // Also calculate time in years after field's reference epoch
                astrometry::UT ut;
                ut.setMJD(expo->mjd);
                // Note that getTTyr returns years since J2000.
                expo->pmTDB = ut.getTTyr() - (fieldEpochs[expo->field] - 2000.);
            }
            if (!epoch.empty()) expo->epoch = epoch[i];
            if (!apcorr.empty()) expo->apcorr = apcorr[i];
            if (!weight.empty()) expo->weight = weight[i];
            if (!magWeight.empty()) expo->magWeight = magWeight[i];
            if (!observatory.empty())
                for (int j = 0; j < 3; j++) expo->observatory[j] = observatory[i][j];
            if (!astrometricCovariance.empty()) {
                expo->astrometricCovariance(0, 0) = astrometricCovariance[i][0];
                expo->astrometricCovariance(1, 1) = astrometricCovariance[i][1];
                expo->astrometricCovariance(0, 1) = astrometricCovariance[i][2];
                expo->astrometricCovariance(1, 0) = astrometricCovariance[i][2];
            }
            if (!photometricVariance.empty()) expo->photometricVariance = photometricVariance[i];
            exposures[i] = std::move(expo);
        }
    }  // End exposure loop
    return exposures;
}

template <class S>
void fixMapComponents(typename S::Collection &pmc, const list<string> &fixMapList,
                      const vector<unique_ptr<Instrument>> &instruments) {
    set<string> fixTheseMaps;
    for (auto iName : pmc.allMapNames()) {
        if (stringstuff::regexMatchAny(fixMapList, iName)) fixTheseMaps.insert(iName);
    }
    // Add names of all devices of instruments on the fixMapList
    for (auto const &instptr : instruments) {
        if (!instptr) continue;  // Not in use
        if (stringstuff::regexMatchAny(fixMapList, instptr->name)) {
            // Loop through all devices
            for (int i = 0; i < instptr->nDevices; i++) {
                string devMap = instptr->name + "/" + instptr->deviceNames.nameOf(i);
                // Freeze the device's map if it's in use
                if (pmc.mapExists(devMap)) fixTheseMaps.insert(devMap);
            }
        }
    }
    pmc.setFixed(fixTheseMaps);
}

vector<shared_ptr<astrometry::Wcs>> readWCSs(const img::FTable &extensionTable) {
    // cerr << inputTables << endl;
    // FITS::FitsTable ft(inputTables, FITS::ReadOnly, "Extensions");
    // cerr << "rW 1" << endl;
    // img::FTable extensionTable = ft.extract();
    vector<shared_ptr<astrometry::Wcs>> WCSs(extensionTable.nrows());

    for (int i = 0; i < extensionTable.nrows(); i++) {
        // Create the starting WCS for the exposure
        string s;
        extensionTable.readCell(s, "WCSIn", i);
        if (stringstuff::nocaseEqual(s, "_ICRS")) {
            // Create a Wcs that just takes input as RA and Dec in degrees;
            astrometry::IdentityMap identity;
            astrometry::SphericalICRS icrs;
            WCSs[i] = shared_ptr<astrometry::Wcs>(
                    new astrometry::Wcs(&identity, icrs, "ICRS_degrees", WCS_UNIT));
        } else {
            istringstream iss(s);
            astrometry::PixelMapCollection pmcTemp;
            if (!pmcTemp.read(iss)) {
                throw std::runtime_error("Could not deserialize starting WCS for extension #" +
                                         std::to_string(i));
            }
            string wcsName = pmcTemp.allWcsNames().front();
            shared_ptr<astrometry::Wcs> tmp(pmcTemp.cloneWcs(wcsName));
            WCSs[i] = shared_ptr<astrometry::Wcs>(pmcTemp.cloneWcs(wcsName));
        }
    }
    return WCSs;
}

template <class S>
vector<unique_ptr<typename S::Extension>> readExtensions(
        const img::FTable &extensionTable, const vector<unique_ptr<Instrument>> &instruments,
        const vector<unique_ptr<Exposure>> &exposures, const vector<int> &exposureColorPriorities,
        vector<unique_ptr<typename S::ColorExtension>> &colorExtensions, astrometry::YAMLCollector &inputYAML,
        bool logging) {
    vector<unique_ptr<typename S::Extension>> extensions(extensionTable.nrows());
    colorExtensions = vector<unique_ptr<typename S::ColorExtension>>(extensionTable.nrows());
    int processed = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 100)
#endif
    for (int i = 0; i < extensionTable.nrows(); i++) {
        int iExposure;
        extensionTable.readCell(iExposure, "Exposure", i);

        if (iExposure < 0 || iExposure >= exposures.size()) {
            throw std::runtime_error("Extension " + std::to_string(i) + " has invalid exposure number " +
                                     std::to_string(iExposure));
        }

        // Determine whether this extension might be used to provide colors
        int colorPriority = exposureColorPriorities[iExposure];
        if (colorPriority >= 0) {
            unique_ptr<typename S::ColorExtension> ce(new typename S::ColorExtension);
            ce->priority = colorPriority;
            colorExtensions[i] = std::move(ce);
        }

        if (!exposures[iExposure]) continue;

#ifdef _OPENMP
#pragma omp critical(processed)
#endif
        ++processed;
        unique_ptr<typename S::Extension> extn(new typename S::Extension);
        extn->exposure = iExposure;
        const Exposure &expo = *exposures[iExposure];

        bool isReference = (expo.instrument == REF_INSTRUMENT || expo.instrument == PM_INSTRUMENT);

        int iDevice;
        extensionTable.readCell(iDevice, "Device", i);
        if (logging && processed % 1000 == 0) {
            cerr << "...Extn " << i << "/" << extensions.size() << " " << expo.name << endl;
        }
        extn->device = iDevice;
        extn->airmass = expo.airmass;
        extn->apcorr = expo.apcorr;
        extn->magshift = +2.5 * log10(expo.exptime);

        // Create the starting WCS for the exposure
        string s;
        extensionTable.readCell(s, "WCSIn", i);
        if (stringstuff::nocaseEqual(s, "_ICRS")) {
            // Create a Wcs that just takes input as RA and Dec in degrees;
            astrometry::IdentityMap identity;
            astrometry::SphericalICRS icrs;
            extn->startWcs = std::unique_ptr<astrometry::Wcs>(
                    new astrometry::Wcs(&identity, icrs, "ICRS_degrees", WCS_UNIT));
        } else {
            istringstream iss(s);
            astrometry::PixelMapCollection pmcTemp;
            if (!pmcTemp.read(iss)) {
                throw std::runtime_error("Could not deserialize starting WCS for extension #" +
                                         std::to_string(i));
            }
            string wcsName = pmcTemp.allWcsNames().front();
            extn->startWcs = std::unique_ptr<astrometry::Wcs>(pmcTemp.cloneWcs(wcsName));
        }

        // destination projection for startWCS is the Exposure projection,
        // so that any exposure-level magnitude corrections use this coord system
        extn->startWcs->reprojectTo(*expo.projection);

        // Extract the map specifications for this extension from the input
        // YAML files.
        if (expo.instrument < 0) {
            // Tag & reference exposures have no instruments and no fitting
            // being done.  Coordinates are fixed to xpix = xw.
            // WCS will be the same for everything in this field
            // ???      extn->wcsName = fieldNames.nameOf(ifield);
            extn->mapName = astrometry::IdentityMap().getName();
        } else {
            // Real instrument, make a map combining its exposure with its Device map:
            astrometry::YAMLCollector::Dictionary d;
            d["INSTRUMENT"] = instruments[expo.instrument]->name;
            d["DEVICE"] = instruments[expo.instrument]->deviceNames.nameOf(extn->device);
            d["EXPOSURE"] = expo.name;
            d["BAND"] = instruments[expo.instrument]->band;
            // Add Epoch to the dictionary if it is not empty
            if (!expo.epoch.empty()) d["EPOCH"] = expo.epoch;
            // This will be the name of the extension's final WCS and photometry map:
            extn->wcsName = d["EXPOSURE"] + "/" + d["DEVICE"];
            // And this is the name of the PixelMap underlying the WCS, or
            // the name of the photometric map (for which we do not want the "base" part).
            if (S::isAstro)
                extn->mapName = extn->wcsName + "/base";
            else
                extn->mapName = extn->wcsName;
            if (!inputYAML.addMap(extn->mapName, d)) {
                throw std::runtime_error("Input YAML files do not have complete information for map " +
                                         extn->mapName);
            }
        }
        extensions[i] = std::move(extn);
    }  // End extension loop
    return extensions;
}

// Routine to find an exposure that should have exposure map set to identity
// to resolve exposure/device degeneracies in the model.  Returns <0 if
// no canonical exposure is necessary for this instrument.  Takes
// as input references to an Instrument (and its id number), the map collection
// being studied, and the exposure/extension tables.

template <class S>
int findCanonical(Instrument &instr, int iInst, vector<unique_ptr<Exposure>> &exposures,
                  vector<unique_ptr<typename S::Extension>> &extensions, typename S::Collection &pmc) {
    // Classify the device maps for this instrument
    set<int> fixedDevices;
    set<int> freeDevices;
    set<int> unusedDevices;

    // And the exposure maps as well:
    set<int> fixedExposures;
    set<int> freeExposures;
    set<int> unusedExposures;
    set<int> itsExposures;  // All exposure numbers using this instrument

    for (int iDev = 0; iDev < instr.nDevices; iDev++) {
        string mapName = instr.name + "/" + instr.deviceNames.nameOf(iDev);
        if (!pmc.mapExists(mapName)) {
            unusedDevices.insert(iDev);
            continue;
        }
        instr.mapNames[iDev] = mapName;
        if (pmc.getFixed(mapName)) {
            fixedDevices.insert(iDev);
        } else {
            freeDevices.insert(iDev);
        }
    }

    for (int iExpo = 0; iExpo < exposures.size(); iExpo++) {
        if (!exposures[iExpo]) continue;  // Skip unused
        if (exposures[iExpo]->instrument == iInst) {
            itsExposures.insert(iExpo);
            string mapName = exposures[iExpo]->name;
            if (!pmc.mapExists(mapName)) {
                unusedExposures.insert(iExpo);
            } else if (pmc.getFixed(mapName)) {
                fixedExposures.insert(iExpo);
            } else {
                freeExposures.insert(iExpo);
            }
        }
    }

    // Now take an inventory of all extensions to see which device
    // solutions are used in coordination with which exposure solutions
    vector<set<int>> exposuresUsingDevice(instr.nDevices);
    for (auto const &extnptr : extensions) {
        if (!extnptr) continue;  // Extension not in use.
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
                    throw std::runtime_error("ERROR: Logic problem: extension map " + extnptr->mapName +
                                             " is using allegedly unused exposure or device map");
                }
            }
        }
    }

    // We have a degeneracy if there is a set of exposures and devices such
    // that
    // * the device maps and exposure maps are free
    // * and the device maps are used only in exposures from this set
    // * and the exposures contain only devices from this set.
    // See here if such a subset exists, in which case we need
    // to fix one of the exposure maps as a "canonical" exposure with
    // Identity map.

    // Split in-use devices into those potentially degenerate and those
    // either fixed or tied to a fixed solution
    auto degenerateDevices = freeDevices;
    auto okDevices = fixedDevices;
    auto degenerateExposures = freeExposures;
    auto okExposures = fixedExposures;

    // propagate "ok-ness" from devices to exposures and back until no more.
    findDegeneracies(degenerateDevices, okDevices, degenerateExposures, okExposures, exposuresUsingDevice);

    // If there are no degenerate device maps, we are done!
    // Return a value indicating no canonical needed at all.
    if (degenerateDevices.empty()) return -1;

    if (degenerateExposures.empty()) {
        throw std::runtime_error("Logic problem: Instrument " + instr.name +
                                 " came up with degenerate devices but not exposures.");
    }

    int canonicalExposure = -1;
    // Find the exposure using the most degenerate devices
    {
        int maxDevices = 0;
        for (auto iExpo : degenerateExposures) {
            int nDevices = 0;
            for (auto iDev : degenerateDevices) {
                if (exposuresUsingDevice[iDev].count(iExpo) > 0) nDevices++;
            }
            if (nDevices > maxDevices) {
                // Save this exposure as best to use, no need
                // to continue if it uses all devices.
                maxDevices = nDevices;
                canonicalExposure = iExpo;
                if (maxDevices == degenerateDevices.size()) break;
            }
        }
    }

    if (canonicalExposure < 0) {
        throw std::runtime_error("Failed to locate a canonical exposure for " + instr.name);
    }

    // Check that fixing this exposure map will resolve degeneracies.
    degenerateExposures.erase(canonicalExposure);
    okExposures.insert(canonicalExposure);
    findDegeneracies(degenerateDevices, okDevices, degenerateExposures, okExposures, exposuresUsingDevice);
    if (!degenerateDevices.empty()) {
        throw std::runtime_error("But canonical did not resolve exposure/device degeneracy.");
    }
    return canonicalExposure;
}

// Add every extension's map to the YAMLCollector and then emit the
// YAML and read into the MapCollection.
// The names of all maps are already in the extension list.
// !! Returns the number of seconds spent in critical loop parts of addMap
template <class S>
void createMapCollection(const vector<unique_ptr<Instrument>> &instruments,
                         const vector<unique_ptr<Exposure>> &exposures,
                         const vector<unique_ptr<typename S::Extension>> &extensions,
                         astrometry::YAMLCollector &inputYAML, typename S::Collection &pmc) {
    double criticalTime = 0.;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 50) reduction(+ : criticalTime)
#endif
    for (int i = 0; i < extensions.size(); i++) {
        auto const &extnptr = extensions[i];
        if (!extnptr) continue;  // Not in use.
        auto &expo = *exposures[extnptr->exposure];
        // Extract the map specifications for this extension from the input
        // YAML files.
        astrometry::YAMLCollector::Dictionary d;
        if (expo.instrument >= 0) {
            // Real instrument, need the translation dictionary
            d["INSTRUMENT"] = instruments[expo.instrument]->name;
            d["DEVICE"] = instruments[expo.instrument]->deviceNames.nameOf(extnptr->device);
            d["EXPOSURE"] = expo.name;
            d["BAND"] = instruments[expo.instrument]->band;
            // Add Epoch to the dictionary if it is not empty
            if (!expo.epoch.empty()) d["EPOCH"] = expo.epoch;
        }
        if (!inputYAML.addMap(extnptr->mapName, d, &criticalTime)) {
            throw std::runtime_error("Input YAML files do not have complete information for map " +
                                     extnptr->mapName);
        }
    }  // End extension loop
    cout << "**Time in addMap critical regions:" << criticalTime << endl;
    {
        Stopwatch timer;
        timer.start();
        istringstream iss(inputYAML.dump());
        timer.stop();
        cout << "Dumping time: " << timer << endl;
        timer.reset();
        timer.start();
        if (!pmc.read(iss)) {
            throw std::runtime_error("Failure parsing the final YAML map specs");
        }
        timer.stop();
        cout << "Loading time: " << timer << endl;
    }
}

template <class S>
void whoNeedsColor(const vector<unique_ptr<typename S::Extension>> &extensions) {
    for (auto const &extnptr : extensions) {
        if (extnptr) {
            extnptr->needsColor = extnptr->map->needsColor();
        }
    }
}

// Read a MatchCatalog Extension, recording in each extension the
// objects from it that need to be read from catalog.

template <class S>
void readMatches(const vector<int> &seq, const vector<LONGLONG> &extn, const vector<LONGLONG> &obj,
                 typename S::MCat &matches, const vector<unique_ptr<typename S::Extension>> &extensions,
                 const vector<unique_ptr<typename S::ColorExtension>> &colorExtensions,
                 const ExtensionObjectSet &skipSet, int minMatches, bool usePM) {
    // Smaller collections for each match
    vector<long> matchExtns;
    vector<long> matchObjs;
    // These variables determine what the highest-priority color information available
    // in the match is so far.
    long matchColorExtension = -1;
    long matchColorObject = 0;
    int colorPriority = -1;

    for (int i = 0; i <= seq.size(); i++) {
        if (i >= seq.size() || seq[i] == 0) {
            // Processing the previous (or final) match.
            int nValid = matchExtns.size();
            if (matchColorExtension < 0) {
                // There is no color information for this match.
                // Discard any detection which requires a color to produce its map:
                nValid = 0;
                for (int j = 0; j < matchExtns.size(); j++) {
                    Assert(extensions[matchExtns[j]]);
                    if (extensions[matchExtns[j]]->needsColor) {
                        // Mark this detection as useless
                        matchExtns[j] = -1;
                    } else {
                        nValid++;
                    }
                }
            }

            if (nValid >= minMatches) {
                // Make a match from the valid entries, and note need to get data for the detections and color
                unique_ptr<typename S::Match> m;
                for (int j = 0; j < matchExtns.size(); j++) {
                    if (matchExtns[j] < 0) continue;  // Skip detections starved of their color
                    unique_ptr<typename S::Detection> d(new typename S::Detection);
                    d->catalogNumber = matchExtns[j];
                    d->objectNumber = matchObjs[j];
                    extensions[matchExtns[j]]->keepers.insert(
                            std::pair<long, typename S::Detection *>(matchObjs[j], d.get()));
                    if (m)
                        m->add(std::move(d));
                    else {
                        // Make a new Match object from this detection.
                        // It might be a PMMatch if this is an appropriate catalog.
                        m = S::makeNewMatch(std::move(d), usePM);
                    }
                }
                matches.push_back(std::move(m));

                if (matchColorExtension >= 0) {
                    // Tell the color catalog that it needs to look this guy up:
                    Assert(colorExtensions[matchColorExtension]);
                    colorExtensions[matchColorExtension]->keepers.insert(
                            std::pair<long, typename S::Match *>(matchColorObject, matches.back().get()));
                }
            }
            // Clear out previous Match:
            matchExtns.clear();
            matchObjs.clear();
            matchColorExtension = -1;
            colorPriority = -1;
        }  // Finished processing previous match

        // If we done reading entries, quit this loop
        if (i >= seq.size()) break;

        // continue loop if this is an object to ignore
        if (skipSet(extn[i], obj[i])) continue;

        // Note extn/obj number of detections with useful mag data
        if (extensions[extn[i]]) {
            // Record a Detection in a useful extension:
            matchExtns.push_back(extn[i]);
            matchObjs.push_back(obj[i]);
        }

        // Record if we've got color information here
        if (colorExtensions[extn[i]]) {
            int newPriority = colorExtensions[extn[i]]->priority;
            if (newPriority >= 0 && (colorPriority < 0 || newPriority < colorPriority)) {
                // This detection holds color info that we want
                colorPriority = newPriority;
                matchColorExtension = extn[i];
                matchColorObject = obj[i];
            }
        }
    }  // End loop of catalog entries
}

template <class S>
void readMatches(const img::FTable &table, typename S::MCat &matches,
                 const vector<unique_ptr<typename S::Extension>> &extensions,
                 const vector<unique_ptr<typename S::ColorExtension>> &colorExtensions,
                 const ExtensionObjectSet &skipSet, int minMatches, bool usePM) {
    vector<int> seq;
    vector<LONGLONG> extn;
    vector<LONGLONG> obj;
    table.readCells(seq, "sequenceNumber");
    table.readCells(extn, "extension");
    table.readCells(obj, "object");
    readMatches<S>(seq, extn, obj, matches, extensions, colorExtensions, skipSet, minMatches, usePM);
}

// Subroutine to get what we want from a catalog entry for WCS fitting
void Astro::fillDetection(Astro::Detection &d, const Exposure &e,
                          astrometry::SphericalCoords &fieldProjection, const img::FTable &table, long irow,
                          string xKey, string yKey, const vector<string> &xyErrKeys, string magKey,
                          string magErrKey, int magKeyElement, int magErrKeyElement, bool xColumnIsDouble,
                          bool yColumnIsDouble, bool errorColumnIsDouble, bool magColumnIsDouble,
                          bool magErrColumnIsDouble, double magshift, const astrometry::PixelMap *startWcs,
                          bool isTag) {
    d.xpix = getTableDouble(table, xKey, -1, xColumnIsDouble, irow);
    d.ypix = getTableDouble(table, yKey, -1, yColumnIsDouble, irow);

    // Get coordinates and transformation matrix
    startWcs->toWorld(d.xpix, d.ypix, d.xw, d.yw);  // no color in startWCS
    auto dwdp = startWcs->dWorlddPix(d.xpix, d.ypix);
    if (isTag) {
        d.invCov.setZero();
    } else {
        astrometry::Matrix22 cov(0.);
        if (xyErrKeys.size() == 1) {
            // We have a single pixel error, diagonal
            double sigma = getTableDouble(table, xyErrKeys[0], -1, errorColumnIsDouble, irow);
            cov(0, 0) = sigma * sigma;
            cov(1, 1) = sigma * sigma;
        } else if (xyErrKeys.size() == 3) {
            // We have three components of an ellipse, giving x^2, y^2, xy values:
            cov(0, 0) = getTableDouble(table, xyErrKeys[0], -1, errorColumnIsDouble, irow);
            cov(1, 1) = getTableDouble(table, xyErrKeys[1], -1, errorColumnIsDouble, irow);
            cov(0, 1) = getTableDouble(table, xyErrKeys[2], -1, errorColumnIsDouble, irow);
            cov(1, 0) = cov(0, 1);
        } else {
            for (auto s : xyErrKeys) cerr << " --- " << s;
            cerr << " ---" << endl;
            throw AstrometryError("Invalid number of xyErrKeys passed to fillDetection");
        }
        cov = dwdp * cov * dwdp.transpose();  // Note this converts to world units (degrees)
        cov += e.astrometricCovariance;
        d.invCov = cov.inverse();
        d.fitWeight = e.weight;
    }

    if (dynamic_cast<const PMMatch *>(d.itsMatch)) {
        // Build projection matrix if this Detection is being used in a PMMatch
        d.buildProjector(e.pmTDB, e.observatory, &fieldProjection);
    }
}

// This one reads a full 5d stellar solution
unique_ptr<astrometry::PMDetection> Astro::makePMDetection(astrometry::Detection const &d, const Exposure &e,
                                                           const img::FTable &table, long irow, string xKey,
                                                           string yKey, string pmRaKey, string pmDecKey,
                                                           string parallaxKey, string pmCovKey,
                                                           bool xColumnIsDouble, bool yColumnIsDouble,
                                                           bool errorColumnIsDouble,
                                                           const astrometry::PixelMap *startWcs) {
    unique_ptr<astrometry::PMDetection> out(new astrometry::PMDetection);
    out->catalogNumber = d.catalogNumber;
    out->objectNumber = d.objectNumber;
    out->map = d.map;
    out->itsMatch = d.itsMatch;

    out->xpix = getTableDouble(table, xKey, -1, xColumnIsDouble, irow);
    out->ypix = getTableDouble(table, yKey, -1, yColumnIsDouble, irow);
    // Read in the reference solution
    double pmRA, pmDec, parallax;
    table.readCell(pmRA, pmRaKey, irow);
    table.readCell(pmDec, pmDecKey, irow);
    table.readCell(parallax, parallaxKey, irow);
    // Convert from I/O units to internal units
    pmRA *= PM_UNIT / (WCS_UNIT / TDB_UNIT);
    pmDec *= PM_UNIT / (WCS_UNIT / TDB_UNIT);
    parallax *= PARALLAX_UNIT / WCS_UNIT;

    // We will want to shift the inputs to move their reference time from
    // the exposure's (i.e. catalog's) reference to that of the field,
    // which is defined as 0 here.
    double epochShift = -e.pmTDB;

    // Shift the RA and Dec according to proper motion, including
    // factors for unit difference, and cos(dec) on the RA.
    double cosdec = cos(out->ypix * WCS_UNIT);
    out->xpix += epochShift * pmRA / cosdec;
    out->ypix += epochShift * pmDec;

    // Get coordinates and transformation matrix to world coordinates
    startWcs->toWorld(out->xpix, out->ypix, out->xw, out->yw);  // no color in startWCS
    auto dwdp = startWcs->dWorlddPix(out->xpix, out->ypix);
    // Add cos(dec) factors, since all error values are assumed for ra*cos(dec).
    // and we are assuming that PM catalogs are in RA/Dec, in WCS_UNIT (degrees, see Units.h)
    dwdp.col(0) /= cosdec;

    // Fill in the PM central values
    out->pmMean[astrometry::X0] = out->xw;
    out->pmMean[astrometry::Y0] = out->yw;
    out->pmMean[astrometry::VX] = pmRA;
    out->pmMean[astrometry::VY] = pmDec;
    out->pmMean[astrometry::PAR] = parallax;

    // Read the covariance matrix
    astrometry::PMCovariance pmCov;      // Covariance from catalog
    astrometry::PMCovariance dwdp5(0.);  // 5d transformation matrix to world coords
    if (errorColumnIsDouble) {
        vector<double> tmp;
        table.readCell(tmp, pmCovKey, irow);
        int k = 0;
        for (int i = 0; i < 5; i++) {
            dwdp5(i, i) = 1.;
            for (int j = 0; j < 5; j++, k++) {
                pmCov(i, j) = tmp[k];
            }
        }
    } else {
        vector<float> tmp;
        table.readCell(tmp, pmCovKey, irow);
        int k = 0;
        for (int i = 0; i < 5; i++) {
            dwdp5(i, i) = 1.;
            for (int j = 0; j < 5; j++, k++) {
                pmCov(i, j) = tmp[k];
            }
        }
    }

    // Covariances were in I/O units, convert everything to internal WCS_UNIT
    for (int j = 0; j < 5; j++) {
        pmCov(astrometry::X0, j) *= RESIDUAL_UNIT / WCS_UNIT;
        pmCov(j, astrometry::X0) *= RESIDUAL_UNIT / WCS_UNIT;
        pmCov(astrometry::Y0, j) *= RESIDUAL_UNIT / WCS_UNIT;
        pmCov(j, astrometry::Y0) *= RESIDUAL_UNIT / WCS_UNIT;
        pmCov(astrometry::VX, j) *= PM_UNIT / (WCS_UNIT / TDB_UNIT);
        pmCov(j, astrometry::VX) *= PM_UNIT / (WCS_UNIT / TDB_UNIT);
        pmCov(astrometry::VY, j) *= PM_UNIT / (WCS_UNIT / TDB_UNIT);
        pmCov(j, astrometry::VY) *= PM_UNIT / (WCS_UNIT / TDB_UNIT);
        pmCov(astrometry::PAR, j) *= PARALLAX_UNIT / WCS_UNIT;
        pmCov(j, astrometry::PAR) *= PARALLAX_UNIT / WCS_UNIT;
    }

    if (epochShift != 0.) {
        // Transform covariance matrix for shift in reference time
        astrometry::PMCovariance shift(0.);
        for (int i = 0; i < 5; i++) shift(i, i) = 1.;
        shift(astrometry::X0, astrometry::VX) = epochShift;
        shift(astrometry::Y0, astrometry::VY) = epochShift;
        pmCov = shift * pmCov * shift.transpose();
    }

    // Now apply transformation from RA/Dec to world coordinate system.
    // Subtlety here: X0 and Y0 are in world coords.
    // The proper motion and parallax are kept in ICRS.
    // So only X0, Y0 are altered by spherical projection.
    dwdp5.subMatrix(0, 2, 0, 2) = dwdp;
    astrometry::PMCovariance wCov = (dwdp5 * pmCov * dwdp5.transpose());
    // Save inverse
#ifdef USE_EIGEN
    pmCov.setIdentity();
    out->pmInvCov = wCov.ldlt().solve(pmCov);
#else
    out->pmInvCov = wCov.inverse();
#endif
    out->fitWeight = e.weight;

    // Now fill in the Detection 2d (inverse) covariance in case we end
    // up using the PM data as a single observation at the field's reference epoch.
    astrometry::Matrix22 cov22 = wCov.subMatrix(0, 2, 0, 2);
    // Add any extra covariance associated with reference catalogs.
    cov22 += e.astrometricCovariance;
    out->invCov = cov22.inverse();

    return out;
}

void Astro::handlePMDetection(unique_ptr<astrometry::PMDetection> pmd, Astro::Detection const &d) {
    auto mm = d.itsMatch;
    if (dynamic_cast<astrometry::PMMatch *>(mm)) {
        // If this Detection is being used in a PMMatch,
        // replace plain old Detection with this PMDetection.
        mm->remove(d);
        mm->add(std::move(pmd));
    } else {
        // PMDetection is being used in non-PM Match.  Slice it down
        // to pure Detection info, replace old detection with it.
        unique_ptr<Astro::Detection> dd(new Astro::Detection(std::move(*pmd)));
        mm->remove(d);
        mm->add(std::move(dd));
    }
}

unique_ptr<astrometry::Match> Astro::makeNewMatch(unique_ptr<Astro::Detection> d, bool usePM) {
    if (usePM) {
        // Make a PMMatch
        return unique_ptr<Match>(new astrometry::PMMatch(std::move(d)));
    } else {
        // Plain old Match
        return unique_ptr<Match>(new Match(std::move(d)));
    }
}

// And a routine to get photometric information too
void Photo::fillDetection(Photo::Detection &d, const Exposure &e,
                          astrometry::SphericalCoords &fieldProjection, const img::FTable &table, long irow,
                          string xKey, string yKey, const vector<string> &xyErrKeys, string magKey,
                          string magErrKey, int magKeyElement, int magErrKeyElement, bool xColumnIsDouble,
                          bool yColumnIsDouble, bool errorColumnIsDouble, bool magColumnIsDouble,
                          bool magErrColumnIsDouble, double magshift, const astrometry::PixelMap *startWcs,
                          bool isTag) {
    d.args.xDevice = getTableDouble(table, xKey, -1, xColumnIsDouble, irow);
    d.args.yDevice = getTableDouble(table, yKey, -1, yColumnIsDouble, irow);
    startWcs->toWorld(d.args.xDevice, d.args.yDevice, d.args.xExposure, d.args.yExposure);

    // Get the mag input and its error
    d.magIn = getTableDouble(table, magKey, magKeyElement, magColumnIsDouble, irow) + magshift;
    double sigma = getTableDouble(table, magErrKey, magErrKeyElement, magErrColumnIsDouble, irow);

    if (isnan(d.magIn) || isinf(d.magIn)) {
        cerr << "** NaN input at row " << irow << endl;
        d.isClipped = true;
        d.magIn = 0.;
    }
    // Map to output and estimate output error
    d.magOut = d.map->forward(d.magIn, d.args);

    d.invVar = 1. / (e.photometricVariance + sigma * sigma);
    // Don't bother, all unity: sigma *= d.map->derivative(d.magIn, d.args);

    if (isTag) {
        d.fitWeight = 0.;
    } else {
        d.fitWeight = e.magWeight;
    }
}

// Read each Extension's objects' data from it FITS catalog
// and place into Detection structures.
template <class S>
void readObjects(const img::FTable &extensionTable, const vector<unique_ptr<Exposure>> &exposures,
                 const vector<unique_ptr<typename S::Extension>> &extensions,
                 const vector<unique_ptr<astrometry::SphericalCoords>> &fieldProjections, bool logging) {
    // Should be safe to multithread this loop as different threads write
    // only to distinct parts of memory.  Protect the FITS table read though.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 60) num_threads(CATALOG_THREADS)
#endif
    for (int iext = 0; iext < extensions.size(); iext++) {
        if (!extensions[iext]) continue;  // Skip unused

        // Relevant structures for this extension
        typename S::Extension &extn = *extensions[iext];
        Exposure &expo = *exposures[extn.exposure];
        if (extn.keepers.empty()) continue;  // or useless

        bool pmCatalog = S::isAstro && expo.instrument == PM_INSTRUMENT;
        bool isTag = expo.instrument == TAG_INSTRUMENT;

        string filename;
        extensionTable.readCell(filename, "Filename", iext);
        if (logging && iext % 50 == 0)
            cerr << "# Reading object catalog " << iext << "/" << extensions.size() << " from " << filename
                 << " seeking " << extn.keepers.size() << " objects" << endl;
        int hduNumber;
        string xKey;
        string yKey;
        vector<string> xyErrKeys;  // Holds name(s) of posn error keys
        string idKey;
        // For PM catalogs:
        string pmRaKey;
        string pmDecKey;
        string parallaxKey;
        string pmCovKey;
        // For magnitudes
        string magKey;
        string magErrKey;
        int magKeyElement;
        int magErrKeyElement;

        extensionTable.readCell(hduNumber, "extension", iext);
        extensionTable.readCell(xKey, "xKey", iext);
        extensionTable.readCell(yKey, "yKey", iext);
        extensionTable.readCell(idKey, "idKey", iext);

        // Keep track of what rows we need to read from this catalog
        vector<string> neededColumns;
        bool useRows = stringstuff::nocaseEqual(idKey, "_ROW");
        if (!useRows) neededColumns.push_back(idKey);
        neededColumns.push_back(xKey);
        neededColumns.push_back(yKey);

        if (S::isAstro) {
            if (!isTag && !pmCatalog) {
                // See if we have the keys for an error ellipse

                if (extensionTable.hasColumn("errXXKey") && extensionTable.hasColumn("errYYKey") &&
                    extensionTable.hasColumn("errXYKey")) {
                    string errXXKey, errYYKey, errXYKey;
                    extensionTable.readCell(errXXKey, "errXXKey", iext);
                    extensionTable.readCell(errYYKey, "errYYKey", iext);
                    extensionTable.readCell(errXYKey, "errXYKey", iext);
                    if (!(errXXKey.empty() || errYYKey.empty() || errXYKey.empty())) {
                        // We have all the keys for an error ellipse
                        xyErrKeys.push_back(errXXKey);
                        xyErrKeys.push_back(errYYKey);
                        xyErrKeys.push_back(errXYKey);
                    }
                }

                if (xyErrKeys.empty()) {
                    // We did not get an error ellipse.  We will need a
                    // scalar error value
                    string errKey;
                    extensionTable.readCell(errKey, "errKey", iext);
                    xyErrKeys.push_back(errKey);
                }

                for (auto s : xyErrKeys) neededColumns.push_back(s);
            } else if (pmCatalog) {
                // Need PM information
                extensionTable.readCell(pmRaKey, "pmRaKey", iext);
                extensionTable.readCell(pmDecKey, "pmDecKey", iext);
                extensionTable.readCell(pmCovKey, "pmCovKey", iext);
                extensionTable.readCell(parallaxKey, "parallaxKey", iext);
                neededColumns.push_back(pmRaKey);
                neededColumns.push_back(pmDecKey);
                neededColumns.push_back(pmCovKey);
                neededColumns.push_back(parallaxKey);
            }
        } else {
            // For photometry:
            extensionTable.readCell(magKey, "magKey", iext);
            extensionTable.readCell(magErrKey, "magErrKey", iext);
            // Be willing to get an element of array-valued bintable cell
            // for mag or magerr.  Syntax would be
            // MAGAPER[4]  to get 4th (0-indexed) element of MAGAPER column
            magKeyElement = elementNumber(magKey);
            magErrKeyElement = elementNumber(magErrKey);
            neededColumns.push_back(magKey);
            neededColumns.push_back(magErrKey);
        }

        const typename S::SubMap *sm = extn.map;

        if (!sm) cerr << "Exposure " << expo.name << " submap is null" << endl;

        astrometry::Wcs *startWcs = extn.startWcs.get();

        if (!startWcs) {
            throw std::runtime_error("Failed to find initial Wcs for exposure " + expo.name);
        }

        img::FTable ff;
#ifdef _OPENMP
#pragma omp critical(fitsio)
#endif
        {
            FITS::FitsTable ft(filename, FITS::ReadOnly, hduNumber);
            ff = ft.extract(0, -1, neededColumns);
        }
        vector<LONGLONG> id;
        if (useRows) {
            id.resize(ff.nrows());
            for (long i = 0; i < id.size(); i++) id[i] = i;
        } else {
            ff.readCells(id, idKey);
        }
        Assert(id.size() == ff.nrows());

        bool xColumnIsDouble = isDouble(ff, xKey, -1);
        bool yColumnIsDouble = isDouble(ff, yKey, -1);
        bool magColumnIsDouble;
        bool magErrColumnIsDouble;
        bool errorColumnIsDouble;
        double magshift = extn.magshift;

        // Get fieldProjection for this catalog
        unique_ptr<astrometry::SphericalCoords> fieldProjection(fieldProjections[expo.field]->duplicate());

        if (S::isAstro) {
            errorColumnIsDouble = isDouble(ff, pmCatalog ? pmCovKey : xyErrKeys[0], -1);
        } else {
            magColumnIsDouble = isDouble(ff, magKey, magKeyElement);
            magErrColumnIsDouble = isDouble(ff, magErrKey, magErrKeyElement);
        }

        for (long irow = 0; irow < ff.nrows(); irow++) {
            auto pr = extn.keepers.find(id[irow]);
            if (pr == extn.keepers.end()) continue;  // Not a desired object

            // Have a desired object now.  Fill its Detection structure
            auto d = pr->second;
            extn.keepers.erase(pr);
            d->map = sm;

            if (pmCatalog) {
                // Need to read data differently from catalog with
                // full proper motion solution.  Get PMDetection.
                auto pmd = S::makePMDetection(*d, expo, ff, irow, xKey, yKey, pmRaKey, pmDecKey, parallaxKey,
                                              pmCovKey, xColumnIsDouble, yColumnIsDouble, errorColumnIsDouble,
                                              startWcs);

                // This routine will either replace the plain old Detection in
                // a PMMatch with this PMDetection; or if it is part of a plain
                // match, it will slice the PMDetection down to a Detection.
                S::handlePMDetection(std::move(pmd), *d);
            } else {
                // Normal photometric or astrometric entry
                S::fillDetection(*d, expo, *fieldProjection, ff, irow, xKey, yKey, xyErrKeys, magKey,
                                 magErrKey, magKeyElement, magErrKeyElement, xColumnIsDouble, yColumnIsDouble,
                                 errorColumnIsDouble, magColumnIsDouble, magErrColumnIsDouble, magshift,
                                 startWcs, isTag);
            }
        }  // End loop over catalog objects

        if (!extn.keepers.empty()) {
            throw std::runtime_error("Did not find all desired objects in catalog " + filename +
                                     " extension " + std::to_string(hduNumber));
        }

    }  // end loop over catalogs to read
}

// Read each Extension's objects' data from it FITS catalog
// and place into Detection structures.
template <class S>
void readObjects_oneExtension(const vector<unique_ptr<Exposure>> &exposures, int iext, const img::FTable &ff,
                              const string &xKey, const string &yKey, const string &idKey,
                              const string &pmCovKey, const vector<string> &xyErrKeys, const string &magKey,
                              const int &magKeyElement, const string &magErrKey,
                              const int &magErrKeyElement,  // TODO: make these dictionary?
                              const string &pmRaKey, const string &pmDecKey, const string &parallaxKey,
                              const vector<unique_ptr<typename S::Extension>> &extensions,
                              const vector<unique_ptr<astrometry::SphericalCoords>> &fieldProjections,
                              bool logging, bool useRows) {
    // Relevant structures for this extension
    typename S::Extension &extn = *extensions[iext];
    Exposure &expo = *exposures[extn.exposure];
    if (extn.keepers.empty()) return;  // or useless

    bool pmCatalog = S::isAstro && expo.instrument == PM_INSTRUMENT;
    bool isTag = expo.instrument == TAG_INSTRUMENT;

    const typename S::SubMap *sm = extn.map;

    if (!sm) cerr << "Exposure " << expo.name << " submap is null" << endl;

    astrometry::Wcs *startWcs = extn.startWcs.get();

    if (!startWcs) {
        throw std::runtime_error("Failed to find initial Wcs for exposure " + expo.name);
    }

    vector<LONGLONG> id;
    if (useRows) {
        id.resize(ff.nrows());
        for (long i = 0; i < id.size(); i++) id[i] = i;
    } else {
        ff.readCells(id, idKey);
    }
    Assert(id.size() == ff.nrows());
    bool xColumnIsDouble = isDouble(ff, xKey, -1);
    bool yColumnIsDouble = isDouble(ff, yKey, -1);
    bool magColumnIsDouble;
    bool magErrColumnIsDouble;
    bool errorColumnIsDouble;
    double magshift = extn.magshift;

    // Get fieldProjection for this catalog
    astrometry::SphericalCoords *fieldProjection = fieldProjections[expo.field]->duplicate();
    if (S::isAstro) {
        if (pmCatalog) errorColumnIsDouble = isDouble(ff, pmCovKey, 0);
        else errorColumnIsDouble = isDouble(ff, xyErrKeys[0], -1);
    } else {
        magColumnIsDouble = isDouble(ff, magKey, magKeyElement);
        magErrColumnIsDouble = isDouble(ff, magErrKey, magErrKeyElement);
    }

    for (long irow = 0; irow < ff.nrows(); irow++) {
        auto pr = extn.keepers.find(id[irow]);
        if (pr == extn.keepers.end()) continue;  // Not a desired object

        // Have a desired object now.  Fill its Detection structure
        typename S::Detection *d = pr->second;
        extn.keepers.erase(pr);
        d->map = sm;
        if (pmCatalog) {
            // Need to read data differently from catalog with
            // full proper motion solution.  Get PMDetection.
            auto pmd = S::makePMDetection(*d, expo, ff, irow, xKey, yKey, pmRaKey, pmDecKey, parallaxKey,
                                          pmCovKey, xColumnIsDouble, yColumnIsDouble, errorColumnIsDouble,
                                          startWcs);

            // This routine will either replace the plain old Detection in
            // a PMMatch with this PMDetection; or if it is part of a plain
            // match, it will slice the PMDetection down to a Detection.
            S::handlePMDetection(std::move(pmd), *d);
        } else {
            // Normal photometric or astrometric entry
            S::fillDetection(*d, expo, *fieldProjection, ff, irow, xKey, yKey, xyErrKeys, magKey, magErrKey,
                             magKeyElement, magErrKeyElement, xColumnIsDouble, yColumnIsDouble,
                             errorColumnIsDouble, magColumnIsDouble, magErrColumnIsDouble, magshift, startWcs,
                             isTag);
        }
    }  // End loop over catalog objects
    if (fieldProjection) delete fieldProjection;

    if (!extn.keepers.empty()) {
        throw std::runtime_error("Did not find all desired objects extension " + iext);
    }
}

// Read color information from files marked as holding such, insert into
// relevant Matches.
template <class S>
void readColors(const img::FTable &extensionTable,
                const vector<unique_ptr<typename S::ColorExtension>> &colorExtensions, bool logging) {
    for (int iext = 0; iext < colorExtensions.size(); iext++) {
        if (!colorExtensions[iext]) continue;  // Skip unused
        auto &extn = *colorExtensions[iext];
        if (extn.keepers.empty()) continue;  // Not using any colors from this catalog
        string filename;
        extensionTable.readCell(filename, "Filename", iext);
        if (logging)
            cerr << "# Reading color catalog " << iext << "/" << colorExtensions.size() << " from "
                 << filename << endl;
        int hduNumber;
        extensionTable.readCell(hduNumber, "extension", iext);
        string idKey;
        extensionTable.readCell(idKey, "idKey", iext);
        string colorExpression;
        try {
            extensionTable.readCell(colorExpression, "colorExpression", iext);
        } catch (img::FTableError &e) {
            // If there is no colorExpression column, use a default:
            colorExpression = "COLOR";
        }
        stripWhite(colorExpression);
        if (colorExpression.empty()) {
            throw std::runtime_error("No colorExpression specified for filename " + filename + " HDU " +
                                     std::to_string(hduNumber));
        }

        // Read the entire catalog for this extension
        FITS::FitsTable ft(filename, FITS::ReadOnly, hduNumber);
        img::FTable ff = ft.use();

        bool useRows = stringstuff::nocaseEqual(idKey, "_ROW");
        vector<long> id;
        if (useRows) {
            id.resize(ff.nrows());
            for (long i = 0; i < id.size(); i++) id[i] = i;
        } else {
            ff.readCells(id, idKey);
        }
        Assert(id.size() == ff.nrows());

        vector<double> color(id.size(), 0.);
        ff.evaluate(color, colorExpression);

        for (long irow = 0; irow < ff.nrows(); irow++) {
            auto pr = extn.keepers.find(id[irow]);
            if (pr == extn.keepers.end()) continue;  // Not a desired object

            // Have a desired object. Put the color into everything it matches
            typename S::Match *m = pr->second;
            extn.keepers.erase(pr);

            for (auto const &detptr : *m) S::setColor(*detptr, color[irow]);

        }  // End loop over catalog objects

        if (!extn.keepers.empty()) {
            throw std::runtime_error("Did not find all desired objects in catalog " + filename +
                                     " extension " + std::to_string(hduNumber) + " " +
                                     std::to_string(extn.keepers.size()) + " left, first ID is " +
                                     std::to_string(extn.keepers.begin()->first));
        }
    }  // end loop over catalogs to read
}

// Find all matched Detections that exceed allowable error, then
// delete them from their Match and delete the Detection.
// Does not apply to reference/tag detections.
// Also deletes anything that is already marked as clipped
template <class S>
void purgeNoisyDetections(double maxError, typename S::MCat &matches,
                          const vector<unique_ptr<Exposure>> &exposures,
                          const vector<unique_ptr<typename S::Extension>> &extensions) {
    for (auto const &mptr : matches) {
        auto j = mptr->begin();
        int k = 0;
        while (j != mptr->end()) {
            auto const &d = *j;  // Yields pointer to each detection
            // Keep it if error is small or if it's from a tag or reference "instrument"
            if ((d->isClipped || d->getSigma() > maxError) &&
                exposures[extensions[d->catalogNumber]->exposure]->instrument >= 0) {
                j = mptr->erase(j);
            } else {
                ++j;
            }
        }
    }
}

// Get rid of Matches with too few Detections being fit: delete
// the Match and all of its Detections.
template <class S>
void purgeSparseMatches(int minMatches, typename S::MCat &matches) {
    auto im = matches.begin();
    while (im != matches.end()) {
        if ((*im)->fitSize() < minMatches) {
            // Remove entire match if it's too small, and kill its Detections too
            (*im)->clear();
            im = matches.erase(im);
        } else {
            ++im;
        }
    }
}

// Get rid of Matches with color outside of specified range.
// Color is always ok if it has NODATA value.
// Note that this even kills Detections that do not need a color for their maps.
template <class S>
void purgeBadColor(double minColor, double maxColor, typename S::MCat &matches) {
    auto im = matches.begin();
    while (im != matches.end()) {
        if ((*im)->size() > 0) {
            // See if color is in range, using color from first Detection (they should
            // all have the same color).
            // Also kills things with NaN / inf colors
            double color = S::getColor(**(*im)->begin());
            if ((color != astrometry::NODATA && (color < minColor || color > maxColor)) || !isfinite(color)) {
                // Remove entire match if it's too small, and kill its Detections too
                (*im)->clear();
                im = matches.erase(im);
                continue;
            }
        }
        ++im;
    }
}

template <class S>
void reserveMatches(typename S::MCat &matches, double reserveFraction, int randomNumberSeed) {
    ran::UniformDeviate<double> u;
    if (randomNumberSeed > 0) u.seed(randomNumberSeed);

    for (auto const &mptr : matches) mptr->setReserved(u < reserveFraction);
}

// Return a map of names of non-frozen exposure maps that
// have fewer than minFitExposure fittable Detections using
// them, also gives the number of fittable Detections they have.
template <class S>
map<string, long> findUnderpopulatedExposures(long minFitExposure, const typename S::MCat &matches,
                                              const vector<unique_ptr<Exposure>> &exposures,
                                              const vector<unique_ptr<typename S::Extension>> &extensions,
                                              const typename S::Collection &pmc) {
    // First count up useful Detections in each extension:
    vector<long> extnCounts(extensions.size(), 0);
    for (auto const &mptr : matches)
        if (!mptr->getReserved())
            for (auto const &dptr : *mptr)
                if (!dptr->isClipped) extnCounts[dptr->catalogNumber]++;

    // Now for each exposure, add up extnCounts of all extensions
    // that depend on a map with name of the exposure.
    map<string, long> bad;
    for (int iExpo = 0; iExpo < exposures.size(); iExpo++) {
        if (!exposures[iExpo]) continue;  // Not in use
        string expoName = exposures[iExpo]->name;
        if (!pmc.mapExists(expoName) || pmc.getFixed(expoName)) continue;  // No free map
        long expoCounts = 0;
        for (int iExtn = 0; iExtn < extensions.size(); iExtn++) {
            if (extnCounts[iExtn] == 0) continue;
            if (pmc.dependsOn(extensions[iExtn]->mapName, expoName)) expoCounts += extnCounts[iExtn];
        }
        if (expoCounts < minFitExposure) bad[expoName] = expoCounts;
    }
    return bad;
}

// Fix the parameters of a map, and delete all Detections making
// use of it as clipped so they will not be used in fitting
template <class S>
void freezeMap(string mapName, typename S::MCat &matches,
               vector<unique_ptr<typename S::Extension>> &extensions, typename S::Collection &pmc) {
    // Nothing to do if map is already fixed or doesn't exist
    if (!pmc.mapExists(mapName) || pmc.getFixed(mapName)) return;

    set<long> badExtensions;  // values of extensions using the map
    for (int iExtn = 0; iExtn < extensions.size(); iExtn++) {
        if (!extensions[iExtn]) continue;  // Not in use
        if (pmc.dependsOn(extensions[iExtn]->mapName, mapName)) badExtensions.insert(iExtn);
    }

    // Now freeze the relevant map
    pmc.setFixed(mapName);

    // And now delete all of the affected Detections
    if (badExtensions.empty()) return;  // nothing to do
    for (auto const &mptr : matches) {
        bool recount = false;
        for (auto dptr = mptr->begin(); dptr != mptr->end(); ++dptr)
            if (badExtensions.count((*dptr)->catalogNumber) > 0) {
                dptr = mptr->erase(dptr);
                recount = true;
            }
    }

    // And the bad extensions too while we're at it
    for (auto extn : badExtensions) {
        // Mark the extension's map as invalid for possible purging.
        if (S::isAstro)
            pmc.invalidate(extensions[extn]->wcsName);
        else
            pmc.invalidate(extensions[extn]->mapName);
        extensions[extn] = nullptr;
    }
}

template <class S>
void matchCensus(const typename S::MCat &matches, ostream &os) {
    long dcount = 0;
    int dof = 0;
    double chi = 0.;
    double maxdev = 0.;
    for (auto const &mptr : matches) {
        dcount += mptr->fitSize();
        chi += mptr->chisq(dof, maxdev);
    }
    cerr << "# Using " << matches.size() << " matches with " << dcount << " total detections." << endl;
    cerr << "#  chisq " << chi << " / " << dof << " dof maxdev " << maxdev << endl;
}

// Map and clip reserved matches
template <class S>
void clipReserved(typename S::Align &ca, double clipThresh, double minimumImprovement, bool clipEntireMatch,
                  bool reportToCerr) {
    double oldthresh = 0.;
    int nclip;
    do {
        // Report number of active Matches / Detections in each iteration:
        long int mcount = 0;
        long int dcount = 0;
        ca.count(mcount, dcount, true, 2);
        double max;
        int dof = 0;
        double chisq = ca.chisqDOF(dof, max, true);
        if (reportToCerr)
            cerr << "Clipping " << mcount << " matches with " << dcount << " detections "
                 << " chisq " << chisq << " / " << dof << " dof,  maxdev " << max << " sigma" << endl;

        double thresh = sqrt(chisq / dof) * clipThresh;  // ??? expected chisq instead of dof?
        if (reportToCerr) cerr << "  new clip threshold: " << thresh << " sigma" << endl;
        if (thresh >= max) break;
        if (oldthresh > 0. && (1 - thresh / oldthresh) < minimumImprovement) break;
        oldthresh = thresh;
        nclip = ca.sigmaClip(thresh, true, clipEntireMatch);
        if (reportToCerr) cerr << "Clipped " << nclip << " matches " << endl;
    } while (nclip > 0);
}

// Save photometric fitting results (residual) to output FITS table.
void Photo::saveResults(const MCat &matches, string outCatalog) {
    // Open the output bintable
    string tablename = "PhotoOut";
    FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, tablename);
    img::FTable outTable = ft.use();
    ;

    // Make a vector of match pointers so we can use OpenMP
    vector<Match *> vmatches;
    vmatches.reserve(matches.size());
    for (auto const &m : matches) vmatches.push_back(m.get());

    {
        // Create vectors to help type each new column
        vector<int> vint;
        vector<long> vlong;
        vector<bool> vbool;
        vector<float> vfloat;
        vector<vector<float>> vvfloat;
        outTable.addColumn(vint, "matchID");
        outTable.addColumn(vlong, "extension");
        outTable.addColumn(vlong, "object");
        outTable.addColumn(vbool, "clip");
        outTable.addColumn(vbool, "reserve");
        outTable.addColumn(vfloat, "xPix");
        outTable.addColumn(vfloat, "yPix");
        outTable.addColumn(vfloat, "xExpo");
        outTable.addColumn(vfloat, "yExpo");
        outTable.addColumn(vfloat, "color");
        outTable.addColumn(vfloat, "magOut");
        outTable.addColumn(vfloat, "magRes");
        outTable.addColumn(vfloat, "magSigma");  // Total magnitude error, incl. sys
        outTable.addColumn(vfloat, "chisq");
        outTable.addColumn(vfloat, "chisqExpected");
    }

    /** Now create a team of threads to refit all matches **/
    // Write residual vectors to table this many at a time
    const int MATCH_CHUNK = 10000;  // Matches to process at a time
    // Cumulative counter for rows written to table so far:
    long pointCount = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iChunk = 0; iChunk <= (vmatches.size() + MATCH_CHUNK - 1) / MATCH_CHUNK; iChunk++) {
        // Create vectors that will be put into output table
        vector<int> matchID;
        vector<long> catalogNumber;
        vector<long> objectNumber;
        vector<bool> clip;
        vector<bool> reserve;
        vector<float> xpix;
        vector<float> ypix;
        vector<float> xExposure;
        vector<float> yExposure;
        vector<float> color;

        vector<float> magOut;
        vector<float> magRes;
        vector<float> sigMag;
        vector<float> chisq;
        vector<float> chisqExpected;

        for (int iMatch = iChunk * MATCH_CHUNK;
             iMatch < vmatches.size() && iMatch < (iChunk + 1) * MATCH_CHUNK; iMatch++) {
            const Match *m = vmatches[iMatch];

            // Update all world coords and mean solutions
            m->remap(true);  // include remapping of clipped detections
            m->solve();

            // Don't report results for matches with no results
            if (m->getDOF() < 0) continue;

            double meanMag = m->getMean();

            for (auto const &detptr : *m) {
                // Save some basics:
                matchID.push_back(iMatch);
                catalogNumber.push_back(detptr->catalogNumber);
                objectNumber.push_back(detptr->objectNumber);
                clip.push_back(detptr->isClipped);
                reserve.push_back(m->getReserved());
                color.push_back(detptr->args.color);

                // Photometry quantities:
                xpix.push_back(detptr->args.xDevice);
                ypix.push_back(detptr->args.yDevice);
                xExposure.push_back(detptr->args.xExposure);
                yExposure.push_back(detptr->args.yExposure);
                magOut.push_back(detptr->magOut);
                magRes.push_back(detptr->residMag());
                sigMag.push_back(detptr->getSigma());
                chisq.push_back(detptr->trueChisq());
                chisqExpected.push_back(detptr->expectedTrueChisq);
            }  // End detection loop
        }      // End match loop

        // At end of each chunk, write new rows to the table.
        // Only one thread at a time writes to the table
#ifdef _OPENMP
#pragma omp critical(fits)
#endif
        {
            long nAdded = matchID.size();
            outTable.writeCells(matchID, "matchID", pointCount);
            outTable.writeCells(catalogNumber, "extension", pointCount);
            outTable.writeCells(objectNumber, "object", pointCount);
            outTable.writeCells(clip, "clip", pointCount);
            outTable.writeCells(reserve, "reserve", pointCount);
            outTable.writeCells(xpix, "xPix", pointCount);
            outTable.writeCells(ypix, "yPix", pointCount);
            outTable.writeCells(xExposure, "xExpo", pointCount);
            outTable.writeCells(yExposure, "yExpo", pointCount);
            outTable.writeCells(color, "color", pointCount);
            outTable.writeCells(magOut, "magOut", pointCount);
            outTable.writeCells(magRes, "magRes", pointCount);
            outTable.writeCells(sigMag, "magSigma", pointCount);
            outTable.writeCells(chisq, "chisq", pointCount);
            outTable.writeCells(chisqExpected, "chisqExpected", pointCount);
            pointCount += nAdded;
        }  // Done flushing the vectors to Table

    }  // End outer chunk loop and multithreaded region
}

/***** Astro version ***/
// Save fitting results (residuals) to output FITS tables.
void Astro::saveResults(const astrometry::MCat &matches, string outCatalog, string starCatalog,
                        vector<astrometry::SphericalCoords *> catalogProjections) {
    // Open the output bintable for pure position residuals
    string tableName = "WCSOut";
    string pmTableName = "PMOut";
    string starTableName = "StarCat";

    // pmTable and starTable will be created after gathering all data

    // Make a vector of match pointers so we can use OpenMP
    vector<Match *> vmatches;
    vmatches.reserve(matches.size());
    for (auto const &m : matches) vmatches.push_back(m.get());

    // Make a global of PMDetections to output
    vector<const PMDetection *> pmdets;
    // And the id's of matches they belong to
    vector<int> pmDetMatchID;
    // Make vector of units conversions to I/O units
    vector<float> units(5);
    units[astrometry::X0] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::Y0] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::PAR] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::VX] = WCS_UNIT / (RESIDUAL_UNIT / TDB_UNIT);
    units[astrometry::VY] = WCS_UNIT / (RESIDUAL_UNIT / TDB_UNIT);
    {
        // Open a scope for existence of the WCSOut table
        FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, tableName);
        img::FTable outTable = ft.use();
        {
            // Create vectors to help type each new column
            vector<int> vint;
            vector<long> vlong;
            vector<bool> vbool;
            vector<float> vfloat;
            vector<vector<float>> vvfloat;
            outTable.addColumn(vint, "matchID");
            outTable.addColumn(vlong, "extension");
            outTable.addColumn(vlong, "object");
            outTable.addColumn(vbool, "clip");
            outTable.addColumn(vbool, "reserve");
            outTable.addColumn(vbool, "hasPM");
            outTable.addColumn(vfloat, "color");
            outTable.addColumn(vfloat, "xPix");
            outTable.addColumn(vfloat, "yPix");
            outTable.addColumn(vfloat, "xresPix");
            outTable.addColumn(vfloat, "yresPix");
            outTable.addColumn(vfloat, "sigPix");
            outTable.addColumn(vfloat, "xW");
            outTable.addColumn(vfloat, "yW");
            outTable.addColumn(vfloat, "xresW");
            outTable.addColumn(vfloat, "yresW");
            outTable.addColumn(vvfloat, "covTotalW", 3);  // 3-element column
            outTable.addColumn(vfloat, "chisq");
            outTable.addColumn(vfloat, "chisqExpected");
        }
        /** Now create a team of threads to refit all matches **/
        // Write residual vectors to table this many at a time
        const int MATCH_CHUNK = 10000;  // Matches to process at a time
        // Cumulative counter for rows written to table so far:
        long pointCount = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int iChunk = 0; iChunk <= (vmatches.size() + MATCH_CHUNK - 1) / MATCH_CHUNK; iChunk++) {
            vector<int> matchID;
            vector<long> catalogNumber;
            vector<long> objectNumber;
            vector<bool> clip;
            vector<bool> reserve;
            vector<bool> hasPM;  // Was this input a full PM estimate?
            vector<float> color;

            vector<float> xresw;
            vector<float> yresw;
            vector<float> xpix;
            vector<float> ypix;
            vector<float> sigpix;  // Circularized total error in pixels
            vector<float> xrespix;
            vector<float> yrespix;
            vector<float> xworld;
            vector<float> yworld;

            // The full error ellipse
            vector<vector<float>> covTotalW;
            // The true chisq for this detection, and the expected value
            vector<float> chisq;
            vector<float> chisqExpected;

            for (int iMatch = iChunk * MATCH_CHUNK;
                 iMatch < vmatches.size() && iMatch < (iChunk + 1) * MATCH_CHUNK; iMatch++) {
                const Match *m = vmatches[iMatch];

                // Update all world coords and mean solutions
                m->remap(true);  // include remapping of clipped detections
                m->solve();

                // Don't report results for matches with no results
                if (m->getDOF() < 0) continue;

                // Create a pointer to PMMatch if this is one:
                auto pmm = dynamic_cast<const astrometry::PMMatch *>(m);

                astrometry::Vector2 xyMean;  // Match's mean position, in world coords
                if (!pmm) {
                    // We can use the same best-fit position for all detections
                    xyMean = m->predict();
                }

                double matchColor = astrometry::NODATA;
                for (auto const &detptr : *m) {
                    // Get a pointer to a PMDetection if this is one:
                    auto pmDetPtr = dynamic_cast<const PMDetection *>(detptr.get());
                    if (pmDetPtr) {
#ifdef _OPENMP
#pragma omp critical(pmdet)
#endif
                        {
                            pmdets.push_back(pmDetPtr);
                            pmDetMatchID.push_back(iMatch);
                        }
                    }

                    // Set color if this is the first detection to have one
                    if (matchColor == astrometry::NODATA && detptr->color != astrometry::NODATA)
                        matchColor = detptr->color;
                    // Save some basics:
                    matchID.push_back(iMatch);
                    catalogNumber.push_back(detptr->catalogNumber);
                    objectNumber.push_back(detptr->objectNumber);
                    clip.push_back(detptr->isClipped);
                    reserve.push_back(m->getReserved());
                    hasPM.push_back(bool(pmDetPtr));
                    color.push_back(detptr->color);
                    xpix.push_back(detptr->xpix);
                    ypix.push_back(detptr->ypix);
                    xworld.push_back(detptr->xw);
                    yworld.push_back(detptr->yw);

                    // Let's get the model prediction
                    if (pmm) xyMean = pmm->predict(detptr.get());

                    // Get world residuals (returned RESIDUAL_UNIT)
                    astrometry::Vector2 residW = detptr->residWorld();
                    xresw.push_back(residW[0]);
                    yresw.push_back(residW[1]);
                    // Get pixel residuals
                    astrometry::Vector2 residP;
                    try {
                        residP = detptr->residPix();
                    } catch (AstrometryError &e) {
                        // Do not want program to crash if an inverse fails.
                        // Just move along
#ifdef _OPENMP
#pragma omp critical(io)
#endif
                        cerr << "WARNING: Astrometry failure for catalog " << detptr->catalogNumber
                             << " object " << detptr->objectNumber << " world " << detptr->xw << ","
                             << detptr->yw << " pix " << detptr->xpix << "," << detptr->ypix << endl;
                        residP[0] = astrometry::NODATA;
                        residP[1] = astrometry::NODATA;
                    }
                    xrespix.push_back(residP[0]);
                    yrespix.push_back(residP[1]);
                    astrometry::Matrix22 cov(0.);
                    // Get error covariances, if we have any
                    double chisqThis;
                    if (detptr->invCov(0, 0) > 0.) {
                        // Save cov in RESIDUAL_UNIT, it's stored in WCS_UNIT
                        cov = detptr->invCov.inverse() * pow(WCS_UNIT / RESIDUAL_UNIT, 2.);
                        vector<float> vv(3, 0.);
                        vv[0] = cov(0, 0);
                        vv[1] = cov(1, 1);
                        vv[2] = cov(0, 1);
                        covTotalW.push_back(vv);

                        // Transform back to pixel errors to get circularized pixel error
                        double sig = 0.;
                        try {
                            cov /= pow(WCS_UNIT / RESIDUAL_UNIT, 2.);  // Back to wcs units
                            auto dpdw = detptr->map->dPixdWorld(xyMean[0], xyMean[1], detptr->color);
                            astrometry::Matrix22 covPix = dpdw * cov * dpdw.transpose();
                            double detCov = covPix(0, 0) * covPix(1, 1) - covPix(1, 0) * covPix(0, 1);
                            sig = detCov > 0. ? pow(detCov, 0.25) : 0.;
                        } catch (AstrometryError &e) {
                            // Problem here getting toPix to work in dPixdWorld
#ifdef _OPENMP
#pragma omp critical(io)
#endif
                            cerr << "WARNING: Astrometry failure for catalog " << detptr->catalogNumber
                                 << " object " << detptr->objectNumber << " world " << detptr->xw << ","
                                 << detptr->yw << " pix " << detptr->xpix << "," << detptr->ypix << endl;
                        }
                        sigpix.push_back(sig);

                        // Chisq and expected
                        chisq.push_back(detptr->trueChisq());
                        chisqExpected.push_back(detptr->expectedTrueChisq);
                    } else {
                        // Did not have a usable error matrix
                        chisqThis = 0.;
                        vector<float> vv(3, 0.);
                        covTotalW.push_back(vv);
                        sigpix.push_back(0.);
                        chisq.push_back(chisqThis);
                        chisqExpected.push_back(detptr->expectedTrueChisq);
                    }
                }  // End detection loop
            }      // End match loop

            // At end of each chunk, write new rows to the table.
            // Only one thread at a time writes to the table
#ifdef _OPENMP
#pragma omp critical(fits)
#endif
            {
                long nAdded = matchID.size();
                outTable.writeCells(matchID, "matchID", pointCount);
                outTable.writeCells(catalogNumber, "extension", pointCount);
                outTable.writeCells(objectNumber, "object", pointCount);
                outTable.writeCells(clip, "clip", pointCount);
                outTable.writeCells(reserve, "reserve", pointCount);
                outTable.writeCells(xpix, "xPix", pointCount);
                outTable.writeCells(ypix, "yPix", pointCount);
                outTable.writeCells(color, "color", pointCount);

                outTable.writeCells(xrespix, "xresPix", pointCount);
                outTable.writeCells(yrespix, "yresPix", pointCount);
                outTable.writeCells(xworld, "xW", pointCount);
                outTable.writeCells(yworld, "yW", pointCount);
                outTable.writeCells(xresw, "xresW", pointCount);
                outTable.writeCells(yresw, "yresW", pointCount);
                outTable.writeCells(sigpix, "sigPix", pointCount);
                outTable.writeCells(hasPM, "hasPM", pointCount);
                outTable.writeCells(covTotalW, "covTotalW", pointCount);
                outTable.writeCells(chisq, "chisq", pointCount);
                outTable.writeCells(chisqExpected, "chisqExpected", pointCount);

                pointCount += nAdded;
            }  // Done flushing the vectors to Table

        }  // End outer chunk loop and multithreaded region

    }  // Close scope of the WCSOut table

    /**/ cerr << "Done WCSOut table" << endl;
    {
        // Now write an extension to the output table for quality of PM detections
        // These are quantities we will want for PMDetections
        Assert(pmdets.size() == pmDetMatchID.size());
        vector<int> pmMatchID(pmdets.size());
        vector<long> pmCatalogNumber(pmdets.size());
        vector<long> pmObjectNumber(pmdets.size());
        vector<bool> pmClip(pmdets.size());
        vector<bool> pmReserve(pmdets.size());
        vector<vector<float>> pmMean(pmdets.size());
        vector<vector<float>> pmInvCov(pmdets.size());
        vector<float> pmChisq(pmdets.size());
        vector<float> pmChisqExpected(pmdets.size());

        // (Not putting residuals here since PM's will be findable
        //  by matchID in the star table)
        for (int i = 0; i < pmdets.size(); i++) {
            auto pmDetPtr = pmdets[i];
            pmMatchID[i] = pmDetMatchID[i];
            pmCatalogNumber[i] = pmDetPtr->catalogNumber;
            pmObjectNumber[i] = pmDetPtr->objectNumber;
            pmClip[i] = pmDetPtr->isClipped;
            pmReserve[i] = pmDetPtr->itsMatch->getReserved();
            vector<float> vv(5);
            for (int ii = 0; ii < 5; ii++) vv[ii] = pmDetPtr->pmMean[ii];
            // Keep position in WCS_UNITS (degrees), put others into RESIDUAL_UNITS
            vv[astrometry::VX] *= units[astrometry::VX];
            vv[astrometry::VY] *= units[astrometry::VY];
            vv[astrometry::PAR] *= units[astrometry::PAR];

            pmMean[i] = vv;
            vv.resize(25);
            int k = 0;
            for (int ii = 0; ii < 5; ii++)
                for (int j = 0; j < 5; j++, k++) vv[k] = pmDetPtr->pmInvCov(ii, j) / (units[ii] * units[j]);
            pmInvCov[i] = vv;
            pmChisq[i] = pmDetPtr->trueChisq();
            pmChisqExpected[i] = pmDetPtr->expectedTrueChisq;
        }  // End PMDetection loop

        if (!pmMatchID.empty()) {
            // Create and fill the PMDetection output table if we have any
            FITS::FitsTable ft(outCatalog, FITS::ReadWrite + FITS::Create, pmTableName);
            auto pmDetTable = ft.use();

            pmDetTable.addColumn(pmMatchID, "matchID");
            pmDetTable.addColumn(pmCatalogNumber, "extension");
            pmDetTable.addColumn(pmObjectNumber, "object");
            pmDetTable.addColumn(pmClip, "clip");
            pmDetTable.addColumn(pmReserve, "reserve");
            pmDetTable.addColumn(pmMean, "pm", 5);           // 5-element fixed-length array
            pmDetTable.addColumn(pmInvCov, "pmInvCov", 25);  // 25 elements
            pmDetTable.addColumn(pmChisq, "chisq");
            pmDetTable.addColumn(pmChisqExpected, "chisqExpected");
        }
    }  // End scope of PMDetection table info

    /**/ cerr << "Done with PMOut table" << endl;

    if (!starCatalog.empty()) {
        // Write a fresh FITS file holding the stellar catalog info for all PM matches
        // Sweep through all matches recording info for the PM ones.
        // These are quantities we will want to put into our output PM catalog
        vector<int> starMatchID;
        vector<bool> starReserve;   // Is this star reserved from fit?
        vector<float> starColor;    // Color for this star
        vector<int> starPMCount;    // How many PMDetections in it were fit?
        vector<int> starDetCount;   // How many non-PM Detections in it were fit?
        vector<int> starClipCount;  // How many detections were clipped?
        vector<int> starDOF;        // Number of DOF in PM fit
        vector<float> starChisq;    // Total chisq
        vector<float> starX;        // The solution
        vector<float> starY;
        vector<float> starPMx;
        vector<float> starPMy;
        vector<float> starParallax;
        vector<vector<float>> starInvCov;  // flattened Fisher matrix of PM

        for (int iMatch = 0; iMatch < vmatches.size(); iMatch++) {
            const Match *m = vmatches[iMatch];

            // Every match has already been remapped and solved.
            // Don't report results for matches with no results
            if (m->getDOF() < 0) continue;
            // Create a pointer to PMMatch if this is one:
            auto pmm = dynamic_cast<const astrometry::PMMatch *>(m);
            if (!pmm) continue;

            // Make an entry.  First count types of data going in
            int detCount = 0;
            int pmDetCount = 0;
            int clipCount = 0;
            double matchColor = astrometry::NODATA;

            for (auto const &detptr : *m) {
                // Get a pointer to a PMDetection if this is one:
                auto pmDetPtr = dynamic_cast<const PMDetection *>(detptr.get());

                if (m->isFit(*detptr)) {
                    if (pmDetPtr) {
                        pmDetCount++;
                    } else {
                        detCount++;
                    }
                } else {
                    clipCount++;
                }

                // Set color if this is the first detection to have one
                if (matchColor == astrometry::NODATA && detptr->color != astrometry::NODATA)
                    matchColor = detptr->color;

            }  // End detection loop

            // Add this object to the PM output catalog
            starMatchID.push_back(iMatch);
            starReserve.push_back(pmm->getReserved());
            starColor.push_back(matchColor);
            starPMCount.push_back(pmDetCount);
            starDetCount.push_back(detCount);
            starClipCount.push_back(clipCount);
            int dof = 0.;
            double dummy;
            double chisqThis = pmm->chisq(dof, dummy);
            starDOF.push_back(dof);
            starChisq.push_back(chisqThis);
            {
                // Save the solution in a table
                auto pm = pmm->getPM();
                // Convert xW, yW to RA/Dec for this table, in WCS_UNITS
                astrometry::SphericalCoords *projection;
                for (auto const &detptr : *m) {
                    projection = catalogProjections[detptr->catalogNumber];
                    if (projection) break;
                }
                if (projection) {
                    projection->setLonLat(pm[astrometry::X0] * WCS_UNIT, pm[astrometry::Y0] * WCS_UNIT);
                    double ra, dec;
                    astrometry::SphericalICRS icrs(*projection);
                    icrs.getRADec(ra, dec);
                    starX.push_back(ra / WCS_UNIT);
                    starY.push_back(dec / WCS_UNIT);
                } else {
                    cerr << "WARNING:  No projection available for matchID " << iMatch << endl;
                    starX.push_back(pm[astrometry::X0] / WCS_UNIT);
                    starY.push_back(pm[astrometry::Y0] / WCS_UNIT);
                }
                starPMx.push_back(pm[astrometry::VX] * units[astrometry::VX]);
                starPMy.push_back(pm[astrometry::VY] * units[astrometry::VY]);
                starParallax.push_back(pm[astrometry::PAR] * units[astrometry::PAR]);
            }
            {
                // And the inverse covariance
                auto fisher = pmm->getInvCovPM();
                vector<float> vv(25);
                int k = 0;
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < 5; j++, k++) vv[k] = fisher(i, j) / (units[i] * units[j]);
                starInvCov.push_back(vv);
            }
        }  // end of match loop for the star catalog.

        if (!starMatchID.empty()) {
            // Create and fill the star catalog table in the output file
            FITS::FitsTable ft(starCatalog, FITS::ReadWrite + FITS::Create + FITS::OverwriteFile,
                               starTableName);
            auto starTable = ft.use();

            starTable.addColumn(starMatchID, "matchID");
            starTable.addColumn(starReserve, "reserve");
            starTable.addColumn(starPMCount, "pmCount");
            starTable.addColumn(starDetCount, "detCount");
            starTable.addColumn(starClipCount, "clipCount");
            starTable.addColumn(starChisq, "chisq");
            starTable.addColumn(starDOF, "dof");
            starTable.addColumn(starColor, "color");
            starTable.addColumn(starX, "xW");
            starTable.addColumn(starY, "yW");
            starTable.addColumn(starPMx, "pmra");
            starTable.addColumn(starPMy, "pmdec");
            starTable.addColumn(starParallax, "parallax");
            starTable.addColumn(starInvCov, "pmInvCov", 25);
        }  // End scope of starCatalog
    }      // End starCatalog making
}

std::map<std::string, vector<float>> Astro::getOutputCatalog(const astrometry::MCat &matches) {

    std::map<std::string, vector<float>> outputDict;

    // TODO: don't need vmatches probably
    // Make a vector of match pointers so we can use OpenMP
    vector<Match *> vmatches;
    vmatches.reserve(matches.size());
    for (auto const &m : matches) vmatches.push_back(m.get());

    vector<float> matchID;
    vector<float> catalogNumber;
    vector<float> objectNumber;
    vector<float> clip;
    vector<float> reserve;
    vector<float> hasPM;  // Was this input a full PM estimate?
    vector<float> color;

    vector<float> xresw;
    vector<float> yresw;
    vector<float> xpix;
    vector<float> ypix;
    vector<float> sigpix;  // Circularized total error in pixels
    vector<float> xrespix;
    vector<float> yrespix;
    vector<float> xworld;
    vector<float> yworld;

    // The full error ellipse
    vector<float> covTotalW_00;
    vector<float> covTotalW_11;
    vector<float> covTotalW_01;
    // The true chisq for this detection, and the expected value
    vector<float> chisq;
    vector<float> chisqExpected;
    cerr << "Matches len: " << vmatches.size() << endl;
    for (int iMatch = 0; iMatch < vmatches.size(); iMatch++) {
        const Match *m = vmatches[iMatch];

        // Update all world coords and mean solutions
        m->remap(true);  // include remapping of clipped detections
        m->solve();

        // Don't report results for matches with no results
        if (m->getDOF() < 0) continue;

        // Create a pointer to PMMatch if this is one:
        auto pmm = dynamic_cast<const astrometry::PMMatch *>(m);

        astrometry::Vector2 xyMean;  // Match's mean position, in world coords
        if (!pmm) {
            // We can use the same best-fit position for all detections
            xyMean = m->predict();
        }

        double matchColor = astrometry::NODATA;
        for (auto const &detptr : *m) {
            // Get a pointer to a PMDetection if this is one:
            auto pmDetPtr = dynamic_cast<const PMDetection *>(detptr.get());

            // Set color if this is the first detection to have one
            if (matchColor == astrometry::NODATA && detptr->color != astrometry::NODATA)
                matchColor = detptr->color;
            // Save some basics:
            matchID.push_back((float)iMatch);
            catalogNumber.push_back((float)detptr->catalogNumber);
            objectNumber.push_back((float)detptr->objectNumber);
            clip.push_back((float)detptr->isClipped);
            reserve.push_back((float)m->getReserved());
            hasPM.push_back((float)(bool)pmDetPtr);
            color.push_back(detptr->color);
            xpix.push_back(detptr->xpix);
            ypix.push_back(detptr->ypix);
            xworld.push_back(detptr->xw);
            yworld.push_back(detptr->yw);

            // Let's get the model prediction
            if (pmm) xyMean = pmm->predict(detptr.get());

            // Get world residuals (returned RESIDUAL_UNIT)
            astrometry::Vector2 residW = detptr->residWorld();
            xresw.push_back(residW[0]);
            yresw.push_back(residW[1]);
            // Get pixel residuals
            astrometry::Vector2 residP;
            try {
                residP = detptr->residPix();
            } catch (AstrometryError &e) {
                // Do not want program to crash if an inverse fails.
                // Just move along
    #ifdef _OPENMP
    #pragma omp critical(io)
    #endif
                cerr << "WARNING: Astrometry failure for catalog " << detptr->catalogNumber
                        << " object " << detptr->objectNumber << " world " << detptr->xw << ","
                        << detptr->yw << " pix " << detptr->xpix << "," << detptr->ypix << endl;
                residP[0] = astrometry::NODATA;
                residP[1] = astrometry::NODATA;
            }
            xrespix.push_back(residP[0]);
            yrespix.push_back(residP[1]);
            astrometry::Matrix22 cov(0.);
            // Get error covariances, if we have any
            double chisqThis;
            if (detptr->invCov(0, 0) > 0.) {
                // Save cov in RESIDUAL_UNIT, it's stored in WCS_UNIT
                cov = detptr->invCov.inverse() * pow(WCS_UNIT / RESIDUAL_UNIT, 2.);
                covTotalW_00.push_back(cov(0,0));
                covTotalW_11.push_back(cov(1,1));
                covTotalW_01.push_back(cov(0,1));

                // Transform back to pixel errors to get circularized pixel error
                double sig = 0.;
                try {
                    cov /= pow(WCS_UNIT / RESIDUAL_UNIT, 2.);  // Back to wcs units
                    auto dpdw = detptr->map->dPixdWorld(xyMean[0], xyMean[1], detptr->color);
                    astrometry::Matrix22 covPix = dpdw * cov * dpdw.transpose();
                    double detCov = covPix(0, 0) * covPix(1, 1) - covPix(1, 0) * covPix(0, 1);
                    sig = detCov > 0. ? pow(detCov, 0.25) : 0.;
                } catch (AstrometryError &e) {
                    // Problem here getting toPix to work in dPixdWorld
    #ifdef _OPENMP
    #pragma omp critical(io)
    #endif
                    cerr << "WARNING: Astrometry failure for catalog " << detptr->catalogNumber
                            << " object " << detptr->objectNumber << " world " << detptr->xw << ","
                            << detptr->yw << " pix " << detptr->xpix << "," << detptr->ypix << endl;
                }
                sigpix.push_back(sig);

                // Chisq and expected
                chisq.push_back(detptr->trueChisq());
                chisqExpected.push_back(detptr->expectedTrueChisq);
            } else {
                // Did not have a usable error matrix
                chisqThis = 0.;
                covTotalW_00.push_back(0);
                covTotalW_11.push_back(0);
                covTotalW_01.push_back(0);
                sigpix.push_back(0.);
                chisq.push_back(chisqThis);
                chisqExpected.push_back(detptr->expectedTrueChisq);
            }
        }  // End detection loop
    }     // End match loop

    outputDict["matchID"] = matchID;
    outputDict["catalogNumber"] = catalogNumber;
    outputDict["objectNumber"] = objectNumber;
    outputDict["clip"] = clip;
    outputDict["reserve"] = reserve;
    outputDict["hasPM"] = hasPM;
    outputDict["color"] = color;
    outputDict["xresw"] = xresw;
    outputDict["yresw"] = yresw;
    outputDict["xpix"] = xpix;
    outputDict["ypix"] = ypix;
    outputDict["sigpix"] = sigpix;
    outputDict["xrespix"] = xrespix;
    outputDict["yrespix"] = yrespix;
    outputDict["xworld"] = xworld;
    outputDict["yworld"] = yworld;
    outputDict["chisq"] = chisq;
    outputDict["chisqExpected"] = chisqExpected;
    outputDict["covTotalW_00"] = covTotalW_00;
    outputDict["covTotalW_11"] = covTotalW_11;
    outputDict["covTotalW_01"] = covTotalW_01;
    return outputDict;
}

std::map<std::string, vector<float>> Astro::getPMCatalog(const astrometry::MCat &matches,
                                                         vector<vector<float>> &pmMean,
                                                         vector<vector<float>> &pmInvCov) {

    std::map<string, vector<float>> PMDict;

    // TODO: probably don't need vmatches
    vector<Match *> vmatches;
    vmatches.reserve(matches.size());

    // Make a global of PMDetections to output
    vector<const PMDetection *> pmdets;
    // And the id's of matches they belong to
    vector<int> pmDetMatchID;

    // Make vector of units conversions to I/O units
    vector<float> units(5);
    units[astrometry::X0] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::Y0] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::PAR] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::VX] = WCS_UNIT / (RESIDUAL_UNIT / TDB_UNIT);
    units[astrometry::VY] = WCS_UNIT / (RESIDUAL_UNIT / TDB_UNIT);

    for (auto const &m : matches) vmatches.push_back(m.get());

    for (int iMatch = 0; iMatch < vmatches.size(); iMatch++) {
        const Match *m = vmatches[iMatch];
        for (auto const &detptr : *m) {
            // Get a pointer to a PMDetection if this is one:
            auto pmDetPtr = dynamic_cast<const PMDetection *>(detptr.get());
            if (pmDetPtr) {
                pmdets.push_back(pmDetPtr);
                pmDetMatchID.push_back(iMatch);
            }
        }
    }

    // Now write an extension to the output table for quality of PM detections
    // These are quantities we will want for PMDetections
    Assert(pmdets.size() == pmDetMatchID.size());
    vector<float> pmMatchID(pmdets.size());
    vector<float> pmCatalogNumber(pmdets.size());
    vector<float> pmObjectNumber(pmdets.size());
    vector<float> pmClip(pmdets.size());
    vector<float> pmReserve(pmdets.size());
    //vector<vector<float>> pmMean(pmdets.size());
    //vector<vector<float>> pmInvCov(pmdets.size());
    vector<float> pmChisq(pmdets.size());
    vector<float> pmChisqExpected(pmdets.size());
    pmMean.reserve(pmdets.size());
    pmInvCov.reserve(pmdets.size());

    // (Not putting residuals here since PM's will be findable
    //  by matchID in the star table)
    for (int i = 0; i < pmdets.size(); i++) {
        auto pmDetPtr = pmdets[i];
        pmMatchID[i] = (float)pmDetMatchID[i];
        pmCatalogNumber[i] = (float)pmDetPtr->catalogNumber;
        pmObjectNumber[i] = (float)pmDetPtr->objectNumber;
        pmClip[i] = (float)pmDetPtr->isClipped;
        pmReserve[i] = (float)pmDetPtr->itsMatch->getReserved();
        vector<float> vv(5);
        for (int ii = 0; ii < 5; ii++) vv[ii] = pmDetPtr->pmMean[ii];
        // Keep position in WCS_UNITS (degrees), put others into RESIDUAL_UNITS
        vv[astrometry::VX] *= units[i];
        vv[astrometry::VY] *= units[i];
        vv[astrometry::PAR] *= units[i];

        pmMean[i] = vv;
        vv.resize(25);
        int k = 0;
        for (int ii = 0; ii < 5; ii++)
            for (int j = 0; j < 5; j++, k++) vv[k] = pmDetPtr->pmInvCov(ii, j) / (units[ii] * units[j]);
        pmInvCov[i] = vv;
        pmChisq[i] = pmDetPtr->trueChisq();
        pmChisqExpected[i] = pmDetPtr->expectedTrueChisq;
    }  // End PMDetection loop

    PMDict["pmMatchID"] = pmMatchID;
    PMDict["pmCatalogNumber"] = pmCatalogNumber;
    PMDict["pmObjectNumber"] = pmObjectNumber;
    PMDict["pmClip"] = pmClip;
    PMDict["pmReserve"] = pmReserve;
    PMDict["pmChisq"] = pmChisq;
    PMDict["pmChisqExpected"] = pmChisqExpected;
    return PMDict;
}

std::map<std::string, vector<float>> Astro::getStarCatalog(const astrometry::MCat &matches,
                                                           vector<astrometry::SphericalCoords *> catalogProjections,
                                                           vector<vector<float>> &starInvCov) {
    // Make a map holding the stellar catalog info for all PM matches
    // Sweep through all matches recording info for the PM ones.
    // These are quantities we will want to put into our output PM catalog
    vector<float> starMatchID;
    vector<float> starReserve;   // Is this star reserved from fit?
    vector<float> starColor;    // Color for this star
    vector<float> starPMCount;    // How many PMDetections in it were fit?
    vector<float> starDetCount;   // How many non-PM Detections in it were fit?
    vector<float> starClipCount;  // How many detections were clipped?
    vector<float> starDOF;        // Number of DOF in PM fit
    vector<float> starChisq;    // Total chisq
    vector<float> starX;        // The solution
    vector<float> starY;
    vector<float> starPMx;
    vector<float> starPMy;
    vector<float> starParallax;
    // We will also fill starInvCov, which is the flattened Fisher matrix of PM

    // Make vector of units conversions to I/O units
    vector<float> units(5);
    units[astrometry::X0] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::Y0] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::PAR] = WCS_UNIT / RESIDUAL_UNIT;
    units[astrometry::VX] = WCS_UNIT / (RESIDUAL_UNIT / TDB_UNIT);
    units[astrometry::VY] = WCS_UNIT / (RESIDUAL_UNIT / TDB_UNIT);

    // TODO: don't need vmatches probably
    // Make a vector of match pointers so we can use OpenMP
    vector<Match *> vmatches;
    vmatches.reserve(matches.size());
    for (auto const &m : matches) vmatches.push_back(m.get());

    for (int iMatch = 0; iMatch < vmatches.size(); iMatch++) {
        const Match *m = vmatches[iMatch];

        // Every match has already been remapped and solved.
        // Don't report results for matches with no results
        if (m->getDOF() < 0) continue;
        // Create a pointer to PMMatch if this is one:
        auto pmm = dynamic_cast<const astrometry::PMMatch *>(m);
        if (!pmm) continue;

        // Make an entry.  First count types of data going in
        int detCount = 0;
        int pmDetCount = 0;
        int clipCount = 0;
        double matchColor = astrometry::NODATA;

        for (auto const &detptr : *m) {
            // Get a pointer to a PMDetection if this is one:
            auto pmDetPtr = dynamic_cast<const PMDetection *>(detptr.get());

            if (m->isFit(*detptr)) {
                if (pmDetPtr) {
                    pmDetCount++;
                } else {
                    detCount++;
                }
            } else {
                clipCount++;
            }

            // Set color if this is the first detection to have one
            if (matchColor == astrometry::NODATA && detptr->color != astrometry::NODATA)
                matchColor = detptr->color;

        }  // End detection loop

        // Add this object to the PM output catalog
        starMatchID.push_back(iMatch);
        starReserve.push_back(pmm->getReserved());
        starColor.push_back(matchColor);
        starPMCount.push_back(pmDetCount);
        starDetCount.push_back(detCount);
        starClipCount.push_back(clipCount);
        int dof = 0.;
        double dummy;
        double chisqThis = pmm->chisq(dof, dummy);
        starDOF.push_back(dof);
        starChisq.push_back(chisqThis);
        {
            // Save the solution in a table
            auto pm = pmm->getPM();
            // Convert xW, yW to RA/Dec for this table, in WCS_UNITS
            astrometry::SphericalCoords *projection;
            for (auto const &detptr : *m) {
                projection = catalogProjections[detptr->catalogNumber];
                if (projection) break;
            }
            if (projection) {
                projection->setLonLat(pm[astrometry::X0] * WCS_UNIT, pm[astrometry::Y0] * WCS_UNIT);
                double ra, dec;
                astrometry::SphericalICRS icrs(*projection);
                icrs.getRADec(ra, dec);
                starX.push_back(ra / WCS_UNIT);
                starY.push_back(dec / WCS_UNIT);
            } else {
                cerr << "WARNING:  No projection available for matchID " << iMatch << endl;
                starX.push_back(pm[astrometry::X0] / WCS_UNIT);
                starY.push_back(pm[astrometry::Y0] / WCS_UNIT);
            }
            starPMx.push_back(pm[astrometry::VX] * units[astrometry::VX]);
            starPMy.push_back(pm[astrometry::VY] * units[astrometry::VY]);
            starParallax.push_back(pm[astrometry::PAR] * units[astrometry::PAR]);
        }
        {
            // And the inverse covariance
            auto fisher = pmm->getInvCovPM();
            vector<float> vv(25);
            int k = 0;
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++, k++) vv[k] = fisher(i, j) / (units[i] * units[j]);
            starInvCov.push_back(vv);
        }
    }  // end of match loop for the star catalog.

    std::map<std::string, vector<float>> starDict;
    starDict["starMatchID"] = starMatchID;
    starDict["starReserve"] = starReserve;
    starDict["starColor"] = starColor;
    starDict["starPMCount"] = starPMCount;
    starDict["starDetCount"] = starDetCount;
    starDict["starClipCount"] = starClipCount;
    starDict["starDOF"] = starDOF;
    starDict["starChisq"] = starChisq;
    starDict["starX"] = starX;
    starDict["starY"] = starY;
    starDict["starPMx"] = starPMx;
    starDict["starPMy"] = starPMy;
    starDict["starParallax"] = starParallax;
    return starDict;
}

void Photo::reportStatistics(const MCat &matches, const vector<unique_ptr<Exposure>> &exposures,
                             const vector<unique_ptr<Photo::Extension>> &extensions, ostream &os) {
    // Create Accum instances for fitted and reserved Detections on every
    // exposure, plus total accumulator for all reference
    // and all non-reference objects.
    typedef Accum<Photo> Acc;
    vector<Acc> vaccFit(exposures.size());
    vector<Acc> vaccReserve(exposures.size());
    Acc refAccFit;
    Acc refAccReserve;
    Acc accFit;
    Acc accReserve;

    // Accumulate stats over all detections.
    for (auto const &mptr : matches) {
        mptr->remap(true);
        mptr->solve();

        for (auto const &dptr : *mptr) {
            // Accumulate statistics for meaningful residuals
            int exposureNumber = extensions[dptr->catalogNumber]->exposure;
            Exposure *expo = exposures[exposureNumber].get();
            if (mptr->getReserved()) {
                if (expo->instrument == REF_INSTRUMENT || expo->instrument == PM_INSTRUMENT) {
                    // Treat PM like reference for photometry
                    refAccReserve.add(*dptr);
                    vaccReserve[exposureNumber].add(*dptr);
                } else if (expo->instrument == TAG_INSTRUMENT) {
                    // do nothing
                } else {
                    accReserve.add(*dptr);
                    vaccReserve[exposureNumber].add(*dptr);
                }
            } else {
                // Not a reserved object:
                if (expo->instrument == REF_INSTRUMENT || expo->instrument == PM_INSTRUMENT) {
                    // Treat PM like reference for photometry
                    refAccFit.add(*dptr);
                    vaccFit[exposureNumber].add(*dptr);
                } else if (expo->instrument == TAG_INSTRUMENT) {
                    // do nothing
                } else {
                    accFit.add(*dptr);
                    vaccFit[exposureNumber].add(*dptr);
                }
            }
        }  // end Detection loop
    }      // end match loop

    // Sort the exposures by name lexical order
    vector<int> ii;
    for (int i = 0; i < exposures.size(); i++)
        if (exposures[i]) ii.push_back(i);
    std::sort(ii.begin(), ii.end(), NameSorter<Exposure>(exposures));

    // Print summary statistics for each exposure
    os << "#   Exposure           | " << Accum<Photo>::header() << endl;
    for (int i = 0; i < ii.size(); i++) {
        int iexp = ii[i];
        if (vaccFit[iexp].ntot == 0 && vaccReserve[iexp].ntot == 0)
            // Do not report statistics if there are no detections
            continue;
        os << setw(3) << iexp << " " << setw(10) << exposures[iexp]->name << " Fit     | "
           << vaccFit[iexp].summary() << endl;
        if (vaccReserve[iexp].n > 0)
            os << setw(3) << iexp << " " << setw(10) << exposures[iexp]->name << " Reserve | "
               << vaccReserve[iexp].summary() << endl;
    }  // exposure summary loop

    // Output summary data for reference catalog and detections
    cout << "# ------------------------------------------------" << endl;
    os << "Detection fit          | " << accFit.summary() << endl;
    if (accReserve.n > 0) os << "Detection reserve      | " << accReserve.summary() << endl;
}

void Astro::reportStatistics(const MCat &matches, const vector<unique_ptr<Exposure>> &exposures,
                             const vector<unique_ptr<typename Astro::Extension>> &extensions, ostream &os) {
    // Create Accum instances for fitted and reserved Detections on every
    // exposure, plus total accumulator for all reference
    // and all non-reference objects.
    typedef Accum<Astro> Acc;
    vector<Acc> vaccFit(exposures.size());
    vector<Acc> vaccReserve(exposures.size());
    Acc refAccFit;
    Acc refAccReserve;
    Acc accFit;
    Acc accReserve;

    // Accumulate stats over all detections.
    for (auto const &mptr : matches) {
        // Make sure match is ready, world coords for all detections done
        mptr->remap(true);
        mptr->solve();

        for (auto const &dptr : *mptr) {
            // Accumulate statistics for meaningful residuals
            int exposureNumber = extensions[dptr->catalogNumber]->exposure;
            Exposure *expo = exposures[exposureNumber].get();
            if (mptr->getReserved()) {
                if (expo->instrument == REF_INSTRUMENT || expo->instrument == PM_INSTRUMENT) {
                    refAccReserve.add(*dptr);
                    vaccReserve[exposureNumber].add(*dptr);
                } else if (expo->instrument == TAG_INSTRUMENT) {
                    // do nothing
                } else {
                    accReserve.add(*dptr);
                    vaccReserve[exposureNumber].add(*dptr);
                }
            } else {
                // Not a reserved object:
                if (expo->instrument == REF_INSTRUMENT || expo->instrument == PM_INSTRUMENT) {
                    refAccFit.add(*dptr);
                    vaccFit[exposureNumber].add(*dptr);
                } else if (expo->instrument == TAG_INSTRUMENT) {
                    // do nothing
                } else {
                    accFit.add(*dptr);
                    vaccFit[exposureNumber].add(*dptr);
                }
            }
        }  // end Detection loop
    }      // end match loop

    // Sort the exposures by name lexical order
    vector<int> ii;
    for (int i = 0; i < exposures.size(); i++)
        if (exposures[i]) ii.push_back(i);
    std::sort(ii.begin(), ii.end(), NameSorter<Exposure>(exposures));

    // Print summary statistics for each exposure
    os << "#   Exposure           | " << Accum<Astro>::header() << endl;
    for (int i = 0; i < ii.size(); i++) {
        int iexp = ii[i];
        if (vaccFit[iexp].ntot == 0 && vaccReserve[iexp].ntot == 0)
            // Do not report statistics if there are no detections
            continue;
        os << setw(3) << iexp << " " << setw(10) << exposures[iexp]->name << " Fit     | "
           << vaccFit[iexp].summary() << endl;
        if (vaccReserve[iexp].n > 0)
            os << setw(3) << iexp << " " << setw(10) << exposures[iexp]->name << " Reserve | "
               << vaccReserve[iexp].summary() << endl;
    }  // exposure summary loop

    // Output summary data for reference catalog and detections

    cout << "# ------------------------------------------------" << endl;

    os << "Reference fit          | " << refAccFit.summary() << endl;
    if (refAccReserve.n > 0) os << "Reference reserve      | " << refAccReserve.summary() << endl;

    os << "Detection fit          | " << accFit.summary() << endl;
    if (accReserve.n > 0) os << "Detection reserve      | " << accReserve.summary() << endl;
}

//////////////////////////////////////////////////////////////////
// For those routines differing for Astro/Photo, here is where
// we instantiate the two cases.
//////////////////////////////////////////////////////////////////

#define INSTANTIATE(AP)                                                                                    \
    template vector<unique_ptr<AP::Extension>> readExtensions<AP>(                                         \
            const img::FTable &extensionTable, const vector<unique_ptr<Instrument>> &instruments,          \
            const vector<unique_ptr<Exposure>> &exposures, const vector<int> &exposureColorPriorities,     \
            vector<unique_ptr<AP::ColorExtension>> &colorExtensions, astrometry::YAMLCollector &inputYAML, \
            bool logging);                                                                                 \
    template int findCanonical<AP>(Instrument & instr, int iInst, vector<unique_ptr<Exposure>> &exposures, \
                                   vector<unique_ptr<AP::Extension>> &extensions, AP::Collection &pmc);    \
    template void fixMapComponents<AP>(AP::Collection &, const list<string> &,                             \
                                       const vector<unique_ptr<Instrument>> &);                            \
    template void createMapCollection<AP>(const vector<unique_ptr<Instrument>> &instruments,               \
                                          const vector<unique_ptr<Exposure>> &exposures,                   \
                                          const vector<unique_ptr<AP::Extension>> &extensions,             \
                                          astrometry::YAMLCollector &inputYAML, AP::Collection &pmc);      \
    template void whoNeedsColor<AP>(const vector<unique_ptr<AP::Extension>> &extensions);                  \
    template void readObjects<AP>(                                                                         \
            const img::FTable &extensionTable, const vector<unique_ptr<Exposure>> &exposures,              \
            const vector<unique_ptr<AP::Extension>> &extensions,                                           \
            const vector<unique_ptr<astrometry::SphericalCoords>> &fieldProjections, bool logging);        \
    template void readObjects_oneExtension<AP>(                                                            \
            const vector<unique_ptr<Exposure>> &exposures, int iext, const img::FTable &ff,                \
            const string &xKey, const string &yKey, const string &idKey, const string &pmCovKey,           \
            const vector<string> &xyErrKeys, const string &magKey, const int &magKeyElement,               \
            const string &magErrKey, const int &magErrKeyElement, const string &pmRaKey,                   \
            const string &pmDecKey, const string &parallaxKey,                                             \
            const vector<unique_ptr<typename AP::Extension>> &extensions,                                  \
            const vector<unique_ptr<astrometry::SphericalCoords>> &fieldProjections, bool logging,         \
            bool useRows);                                                                                 \
    template void readMatches<AP>(const img::FTable &table, typename AP::MCat &matches,                    \
                                  const vector<unique_ptr<AP::Extension>> &extensions,                     \
                                  const vector<unique_ptr<AP::ColorExtension>> &colorExtensions,           \
                                  const ExtensionObjectSet &skipSet, int minMatches, bool usePM);          \
    template void readColors<AP>(const img::FTable &extensionTable,                                        \
                                 const vector<unique_ptr<AP::ColorExtension>> &colorExtensions,            \
                                 bool logging);                                                            \
    template void purgeNoisyDetections<AP>(double maxError, AP::MCat &matches,                             \
                                           const vector<unique_ptr<Exposure>> &exposures,                  \
                                           const vector<unique_ptr<AP::Extension>> &extensions);           \
    template void purgeSparseMatches<AP>(int minMatches, AP::MCat &matches);                               \
    template void purgeBadColor<AP>(double minColor, double maxColor, AP::MCat &matches);                  \
                                                                                                           \
    template void reserveMatches<AP>(AP::MCat & matches, double reserveFraction, int randomNumberSeed);    \
                                                                                                           \
    template map<string, long> findUnderpopulatedExposures<AP>(                                            \
            long minFitExposure, const AP::MCat &matches, const vector<unique_ptr<Exposure>> &exposures,   \
            const vector<unique_ptr<AP::Extension>> &extensions, const AP::Collection &pmc);               \
    template void freezeMap<AP>(string mapName, AP::MCat & matches,                                        \
                                vector<unique_ptr<AP::Extension>> & extensions, AP::Collection & pmc);     \
    template void matchCensus<AP>(const AP::MCat &matches, ostream &os);                                   \
    template void clipReserved<AP>(AP::Align & ca, double clipThresh, double minimumImprovement,           \
                                   bool clipEntireMatch, bool reportToCerr);

INSTANTIATE(Astro)
INSTANTIATE(Photo)
