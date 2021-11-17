// Program to match catalogs based on world coordinates.

#include "Std.h"
#include "WCSFoF_match.h"
#include "Pset.h"
#include "FitSubroutines.h"

#include <iostream>

using namespace std;
using namespace FITS;

string usage=
  "Match objects using FoF algorithm on world coordinate system\n"
  "WCSFoF <configfile> [parameter file] [parameter file...]\n"
  "   [-parameter[=]value...]\n"
  "      <configfile>: FITS file with binary tables describing input catalogs\n"
  "      Program parameters specified as command-line options or read from\n"
  "          parameter file(s) specified on cmd line";
  

int
main(int argc,
     char *argv[])
{
  // Read in parameters
  double matchRadius;
  bool useAffinities;
  string outCatalogName;
  int minMatches;
  bool allowSelfMatches;

  Pset parameters;
  {
    const int def=PsetMember::hasDefault;
    const int low=PsetMember::hasLowerBound;
    const int up=PsetMember::hasUpperBound;
    const int lowopen = low | PsetMember::openLowerBound;
    const int upopen = up | PsetMember::openUpperBound;

    parameters.addMember("matchRadius",&matchRadius, def | lowopen,
			 "Matching tolerance (arcsec)", 1., 0.);
    parameters.addMember("useAffinities",&useAffinities, def,
			 "Disallow star-galaxy matches & inter-filter galaxy matches?", true);
    parameters.addMember("outName",&outCatalogName, def,
			 "filename for FITS output catalog", "match.fof");
    parameters.addMember("minMatch",&minMatches, def | low,
			 "Minimum number of detections for usable match", 2, 2);
    parameters.addMember("selfMatch",&allowSelfMatches, def,
			 "Retain matches that have 2 elements from same exposure?", false);
  }

  
  ////////////////////////////////////////////////
  // Read parameters
  ////////////////////////////////////////////////

  // Read all the command-line and parameter-file program parameters
  processParameters(parameters, usage, 1, argc, argv);
  string configFile = argv[1];

  // Set of Friends-of-Friends class
  FoFClass fofclass;
  fofclass.matchRadius = matchRadius;
  fofclass.useAffinities = useAffinities;
  fofclass.minMatches = minMatches;
  fofclass.allowSelfMatches = allowSelfMatches;

  try {
    ////////////////////////////////////////////////
    // Read in field table and save away orientations of tangent plane systems for each
    ////////////////////////////////////////////////
    FTable fieldTable;
    try {
    fieldTable = FitsTable(configFile, FITS::ReadOnly, "Fields").extract();
    } catch (std::runtime_error& e) {
    quit(e,1);
    }

    vector<Field*> fields;

    for (int ifield = 0; ifield < fieldTable.nrows(); ifield++) {
      string name;
      double ra;
      double dec;
      double extent;
      fieldTable.readCell(name, "NAME", ifield);
      fieldTable.readCell(ra, "RA", ifield);
      fieldTable.readCell(dec, "DEC", ifield);
      fieldTable.readCell(extent, "RADIUS", ifield);

      Field* f = new Field;
      f->name = name;
      astrometry::Orientation orient(astrometry::SphericalICRS(ra*WCS_UNIT, dec*WCS_UNIT));
      f->projection = new astrometry::Gnomonic(orient);
      f->extent = extent;
      f->matchRadius = matchRadius;
      fields.push_back(f);
    } // Done reading fields
    Assert(!fields.empty());

    // Now read in all the instrument tables
    list<int> instrumentHDUs;
    {
      // Find which extensions are instrument tables
      FITS::FitsFile ff(configFile);
      for (int i=1; i<ff.HDUCount(); i++) {
        FITS::Hdu h(configFile, FITS::HDUAny, i);
        if (stringstuff::nocaseEqual(h.getName(), "Instrument"))
	    instrumentHDUs.push_back(i);
      }
    }
    vector<FTable> instrumentTables(instrumentHDUs.size());
    vector<Instr*> instruments(instrumentHDUs.size(),0);

    for (list<int>::const_iterator i=instrumentHDUs.begin();
        i != instrumentHDUs.end();
	    ++i) {
      FTable ft = FitsTable(configFile, FITS::ReadOnly, *i).extract();
      string instrumentName;
      int instrumentNumber;
      if (!ft.header()->getValue("Name", instrumentName)
          || !ft.header()->getValue("Number", instrumentNumber)) {
        cerr << "Could not read name and/or number of instrument at extension "
	         << *i << endl;
      }
      // Get rid of spaces in the name
      spaceReplace(instrumentName);
      Assert(instrumentNumber < instrumentTables.size() && instrumentNumber>=0);
      ft.header()->setValue("Name",instrumentName);

      // Now save the table and make an Instrument structure
      instrumentTables[instrumentNumber] = ft;
      instruments[instrumentNumber] = new Instr(ft);
    }
    // Check that all Instruments were read
    for (int i=0; i<instruments.size(); i++)
      if (!instruments[i]) {
        cerr << "Failed to read instrument number " << i << endl;
        exit(1);
    }

    // Now get information on exposures
    FTable exposureTable;
    try {
      exposureTable = FitsTable(configFile, FITS::ReadOnly, "Exposures").extract();
    } catch (std::runtime_error& e) {
      quit(e,1);
    }
      
    vector<Expo*> exposures;
    for (int i=0; i<exposureTable.nrows(); i++) {
      Expo* e = new Expo;
      e->read(exposureTable, i);
      exposures.push_back(e);
    }
    Assert(exposures.size() == exposureTable.nrows());

    // Table of information about every catalog file to read:
    FTable extensionTable;
    try {
      extensionTable = FitsTable(configFile, FITS::ReadOnly, "Extensions").extract();
    } catch (std::runtime_error& e) {
      quit(e,1);
    }

    fofclass.fields = fields;
    fofclass.instruments = instruments;
    fofclass.exposures = exposures;
    fofclass.extensionTable = extensionTable;
    
    // Loop over input catalogs
    for (long iextn = 0; iextn < extensionTable.nrows(); iextn++) {
      astrometry::Wcs* wcs = 0;
      string thisAffinity;
      int exposureNumber;
      int instrumentNumber;
      int fieldNumber;
      int deviceNumber;
      vector<bool> isStar;
      vector<double> vx;
      vector<double> vy;
      vector<long> vid;
      fofclass.getExtensionInfo(iextn, thisAffinity, exposureNumber, instrumentNumber, fieldNumber,
                                deviceNumber, isStar, vx, vy, vid);
      fofclass.getWCS(iextn, fieldNumber, wcs);
      // TODO: take out fieldNumber, use fields[fieldNumber]
      fofclass.addCatalog(wcs, thisAffinity, exposureNumber, fieldNumber, instrumentNumber, deviceNumber, iextn, isStar, vx, vy, vid);
    }
    fofclass.writeMatches(outCatalogName);

    cerr << "Passed fof class setup" << endl;
    friendOfFriend(matchRadius, useAffinities, outCatalogName, minMatches, allowSelfMatches, fields,
                   instruments, exposures, extensionTable);

    //  Write all of our tables to output file
    {
      // Fields
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::OverwriteFile, "Fields");
      ft.copy(fieldTable);
    }
    {
      // Exposures
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Exposures");
      ft.copy(exposureTable);
    }

    for (int i=0; i<instruments.size(); i++) {
      // Instrument tables
      // Update the device bounds in the table first
      const Instr& inst = *instruments[i];
      int nDevices = inst.size();
      vector<double> vxmin(nDevices);
      vector<double> vxmax(nDevices);
      vector<double> vymin(nDevices);
      vector<double> vymax(nDevices);
      for (int j=0; j<nDevices; j++) {
        vxmin[j] = floor(inst[j].getXMin());
        vxmax[j] = ceil(inst[j].getXMax());
        vymin[j] = floor(inst[j].getYMin());
        vymax[j] = ceil(inst[j].getYMax());
      }
      instrumentTables[i].writeCells(vxmin, "XMin");
      instrumentTables[i].writeCells(vxmax, "XMax");
      instrumentTables[i].writeCells(vymin, "YMin");
      instrumentTables[i].writeCells(vymax, "YMax");

      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
      ft.setName("Instrument");
      ft.setVersion(i+1); // might use this later, or put into extension name???
      ft.copy(instrumentTables[i]);
    }
    
    {
      // Extension table
      FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, "Extensions");
      ft.copy(extensionTable);
    }

    } catch (std::runtime_error &m) {
    quit(m,1);
  }
}