// Program to match catalogs based on world coordinates.

#include "WCSFoF_match.h"
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

using namespace std;
using namespace img;
using namespace FITS;
using stringstuff::stripWhite;

using astrometry::WCS_UNIT;  // we will want all WCS's to operate in same units


FoFClass::FoFClass(){
  // Teach PixelMapCollection about all types of PixelMaps it might need to deserialize
  loadPixelMapParser();

  // And a list holding all Points being matched.  Note that catalogs do not make copies
  // of the Points, they just hold pointers to them.  So this is a warehouse.
  // If we wanted to save the memory of the list links we could just let all the Points
  // be zombies that are destroyed when the program ends.
  list<Point> allPoints;

  // If WCS is coming from a serialized PixelMapCollection, I'll keep last-used one around
  // to avoid re-parsing it all the time:
  string currentPixelMapCollectionFileName = "";
  astrometry::PixelMapCollection* inputPmc=0;

}

void FoFClass::getExtensionInfo(long iextn, string thisAffinity, int exposureNumber, int instrumentNumber,
                                int fieldNumber, int deviceNumber, vector<bool> isStar, vector<double> vx, vector<double> vy, vector<long> vid){
  string filename;
  extensionTable.readCell(filename, "FILENAME", iextn);
  stripWhite(filename);
  int hduNumber;
  extensionTable.readCell(hduNumber, "EXTENSION", iextn);

  if (iextn%10==0) cerr << "# Reading object catalog " << iextn
      << "/" << extensionTable.nrows()
      << " in " << filename 
      << " HDU #" << hduNumber
      << endl;
  FTable ft = FitsTable(filename, FITS::ReadOnly, hduNumber).extract();

  string idKey;
  extensionTable.readCell(idKey, "IDKEY",iextn);
  
  if (stringstuff::nocaseEqual(idKey, "_ROW")) {
    // Want the input table row number as ID column, so make such a column
    vector<long> vi(ft.nrows());
    for (int i=0; i<vi.size(); i++) vi[i] = i;
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
    string selectionExpression;
    extensionTable.readCell(selectionExpression, "SELECT", iextn);
    stripWhite(selectionExpression);
    if (!selectionExpression.empty())
      ft.filterRows(selectionExpression);
  }

  if (ft.nrows()<=0) {
    // No data to use in this extension.  Don't read anything.
    return;
  }
  
  //vector<double> vx;
  string key;
  extensionTable.readCell(key, "XKEY", iextn);
  try {
    ft.readCells(vx, key);
  } catch (FTableError& m) {
    // Trap for using float column in source file instead of double:
    vector<float> vf;
    ft.readCells(vf, key);
    vx.resize(vf.size());
    for (int i=0; i<vf.size(); i++) vx[i]=vf[i];
  }
  //vector<double> vy;
  extensionTable.readCell(key, "YKEY", iextn);
  try {
    ft.readCells(vy, key);
  } catch (FTableError& m) {
    // Trap for using float column in source file instead of double:
    vector<float> vf;
    ft.readCells(vf, key);
    vy.resize(vf.size());
    for (int i=0; i<vf.size(); i++) vy[i]=vf[i];
  }
  //vector<long> vid;
  extensionTable.readCell(key, "IDKEY", iextn);
  try {
    //**/cerr << "Trying long read" << endl;
    ft.readCells(vid, key);
        } catch (FTableError& m) {
    // Trap for using int column in source file instead of long:
    vector<LONGLONG> vidint;
    //**/cerr << "Retrying int" << endl;
    ft.readCells(vidint, key);
    vid.reserve(vidint.size());
    vid.insert(vid.begin(), vidint.begin(), vidint.end());
  }

  {
    string starExpression;
    extensionTable.readCell(starExpression, "STARSELECT", iextn);
    stripWhite(starExpression);
    if (!starExpression.empty())
      ft.evaluate(isStar, starExpression);
  }
}

void FoFClass::getWCS(long iextn, int fieldNumber, astrometry::Wcs* wcs) {
  // Get the WCS for this extension
  string wcsin;
  extensionTable.readCell(wcsin, "WCSIN", iextn);

  list<string> wcssplit = stringstuff::split(wcsin,'@');
  if (wcssplit.size()==2) {
    // @ sign signals that we are getting a map from a PMC
    string wcsName = wcssplit.front();
    string pmcFile = wcssplit.back();
    // Replace any white space in wcs name with underscore:
    wcsName = stringstuff::regexReplace("[[:space:]]+","_",wcsName);

    if (pmcFile!=currentPixelMapCollectionFileName) {
      // Need to open this PMC file
      // Try opening file as new serialized PMC
      ifstream wcsStream(pmcFile.c_str());
      if (!wcsStream.is_open()) {
        cerr << "Could not open PixelMapCollection file <" << pmcFile 
        << ">" << endl;
        exit(1);
      }
      astrometry::PixelMapCollection* tryPMC = new astrometry::PixelMapCollection;
      if (tryPMC->read(wcsStream)) {
        // Successfully found serialized file:
        if (inputPmc) delete inputPmc;
        inputPmc = tryPMC;
        currentPixelMapCollectionFileName = pmcFile;
      } else {
        cerr << "Failed deserializing PixelMapCollection file <" << pmcFile 
        << ">" << endl;
        exit(1);
      }
    } // Done reading the pmcFile

    // Read the wcs from the collection
    wcs = inputPmc->issueWcs(wcsName)->duplicate();
  }
  else if (stringstuff::nocaseEqual(wcsin, "_ICRS")) {
    // No mapping will be needed.  Do nothing.
  } 
  else {
    // Assume that wcsin is a FITS WCS specification
    Header wcsHead;
    //	wcsin += "\nEND\n";  //*** hack to fix problem in astropy.io.fits
    istringstream wcsStream(wcsin);
    if (!(wcsStream >> wcsHead)) {
      cerr << "Error reading WCS/header from WCSIN at extension " << iextn << endl;
      exit(1);
    }
    wcs = astrometry::readTPV(wcsHead);
    if (!wcs) {
      cerr << "Failed reading TPV WCS from extension " << iextn << endl;
      exit(1);
    }
  }

  // Now replace WCSIN in table with a new serialized one, 
  // and set projection to be the field coordinates
  if (wcs) {
    string wcsDump;
    astrometry::PixelMapCollection pmc;
    pmc.learnWcs(*wcs);
    ostringstream oss;
    pmc.writeWcs(oss,wcs->getName());
    wcsDump = oss.str();
    wcs->reprojectTo(*fields[fieldNumber]->projection);
    extensionTable.writeCell(wcsDump, "WCSIN", iextn);
  }

}

void FoFClass::reprojectWCS(astrometry::Wcs* wcs, int fieldNumber){
  wcs->reprojectTo(*fields[fieldNumber]->projection);
}

void FoFClass::addCatalog(astrometry::Wcs* wcs, string thisAffinity, int exposureNumber, int fieldNumber,
                          int instrumentNumber, int deviceNumber, long iextn, vector<bool> isStar, vector<double> vx, vector<double> vy, vector<long> vid) 
{
  //cerr << matchRadius << endl;
  //cerr << "In addCatalog: " << matchRadius << endl;
  //cerr << "test2 " << endl;
  PointCat* starCatalog = fields[fieldNumber]->catalogFor(stellarAffinity);
  PointCat* galaxyCatalog = (stringstuff::nocaseEqual(stellarAffinity,
                                                      thisAffinity) ?
                             starCatalog : 
                             fields[fieldNumber]->catalogFor(thisAffinity));

  // Now loops over objects in the catalog
  for (int iObj = 0; iObj < vx.size(); iObj++) {
    double xpix = vx[iObj];
    double ypix = vy[iObj];
    // Expand bounds of device:
    if (instrumentNumber >=0) 
      (*instruments[instrumentNumber])[deviceNumber]+=Position<double>(xpix,ypix);

    //  maps coords to field's tangent plane
    double xw, yw;
    if (wcs) {
      Assert(wcs);
      wcs->toWorld(xpix, ypix, xw, yw);
    } else {
      // Already have RA, Dec, just project them
      // Assume they are in our std units (see Units.h)
      fields[fieldNumber]->projection->convertFrom(astrometry::SphericalICRS(vx[iObj]*WCS_UNIT, 
                                                                              vy[iObj]*WCS_UNIT ));
      fields[fieldNumber]->projection->getLonLat(xw, yw);
      // Matching program will use our std unit (Units.h), not radians:
      xw /= WCS_UNIT;
      yw /= WCS_UNIT;
    }

    allPoints.push_back(Point(xw, yw, iextn, vid[iObj], exposureNumber));

    // Now we match!!!!
    if (isStar[iObj])
      starCatalog->add(allPoints.back());
    else
      galaxyCatalog->add(allPoints.back());
  } // end object loop
  cerr << "AllPoints size: " << allPoints.size() << endl;
}

void FoFClass::addTest() {
  c = a + b;
}

void FoFClass::writeMatches(string outCatalogName) {

    // Now write all the match catalogs
    long matchCount = 0;
    long pointCount = 0;
    //  Loop over fields
    for (int iField=0; iField<fields.size(); iField++) {
      //  loop over affinities per field
      for (Field::CatMap::iterator iAffinity=fields[iField]->catalogs.begin();
            iAffinity != fields[iField]->catalogs.end();
            ++iAffinity) {
        const PointCat* pcat= iAffinity->second;
        cerr << "Catalog for field " << fields[iField]->name
              << " affinity " << iAffinity->first
              << " with " << pcat->size() << " matches "
              << endl;
        if (pcat->empty()) continue;
        FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
        ft.setName("MatchCatalog");
        FTable ff=ft.use();
        ff.header()->replace("Field", fields[iField]->name, "Field name");
        ff.header()->replace("FieldNum", iField, "Field number");
        ff.header()->replace("Affinity", iAffinity->first, "Affinity name");

        vector<int> sequence;
        vector<long> extn;
        vector<long> obj;
        long matches=0;
        // Now loop through matches in this catalog
        for (PointCat::const_iterator j=pcat->begin();
              j != pcat->end();
              ++j) {
          // Skip any Match that is below minimum match size
          if ((*j)->size() < minMatches) continue;
          bool selfMatch = false;
          if (!allowSelfMatches) {
            set<long> itsExposures;
            for (list<const Point*>::const_iterator k=(*j)->begin();
            k != (*j)->end();
            ++k) 
              if ( !itsExposures.insert((*k)->exposureNumber).second) {
          selfMatch = true;
          break;
              }
          }
          if (selfMatch) continue;

          // Add elements of this match to the output vectors
          int seq=0;
          ++matches;
          ++matchCount;
          for (list<const Point*>::const_iterator k=(*j)->begin();
                k != (*j)->end();
                ++k, ++seq) {
            sequence.push_back(seq);
            extn.push_back((*k)->extensionNumber);
            obj.push_back((*k)->objectNumber);
	  }
	} // end match loop
	cerr << "...Kept " << matches << " matches with " << sequence.size()
	     << " points." << endl;
	pointCount += sequence.size();
	ff.addColumn(sequence, "SequenceNumber");
	ff.addColumn(extn, "Extension");
	ff.addColumn(obj, "Object");
      } // end Affinity loop
    } // end Field loop
    
    cerr << "Total of " << matchCount
	 << " matches with " << pointCount
	 << " points." << endl;
}

void FoFClass::sortMatches(int fieldNumber){//, vector<int> sequence, vector<long> extn,
                             //vector<long> obj) {

    // Now write all the match catalogs
    long matchCount = 0;
    long pointCount = 0;
    //vector<int> sequence;
    //vector<long> extn;
    //vector<long> obj;
    
    //  loop over affinities per field
    for (Field::CatMap::iterator iAffinity=fields[fieldNumber]->catalogs.begin();
          iAffinity != fields[fieldNumber]->catalogs.end();
          ++iAffinity) {
      
      const PointCat* pcat= iAffinity->second;
      cerr << "Catalog for field " << fields[fieldNumber]->name
            << " affinity " << iAffinity->first
            << " with " << pcat->size() << " matches "
            << endl;
      
      if (pcat->empty()) continue;

      //vector<int> sequence;
      //vector<long> extn;
      //vector<long> obj;
      long matches=0;
      // Now loop through matches in this catalog
      for (PointCat::const_iterator j=pcat->begin();
            j != pcat->end();
            ++j) {
        // Skip any Match that is below minimum match size
        
        if ((*j)->size() < minMatches) continue;
        bool selfMatch = false;
        if (!allowSelfMatches) {
          set<long> itsExposures;
          for (list<const Point*>::const_iterator k=(*j)->begin();
          k != (*j)->end();
          ++k) 
            if ( !itsExposures.insert((*k)->exposureNumber).second) {
        selfMatch = true;
        break;
            }
        }
        
        if (selfMatch) continue;

        // Add elements of this match to the output vectors
        int seq=0;
        ++matches;
        ++matchCount;
        for (list<const Point*>::const_iterator k=(*j)->begin();
              k != (*j)->end();
              ++k, ++seq) {
          sequence.push_back(seq);
          extn.push_back((*k)->extensionNumber);
          obj.push_back((*k)->objectNumber);
        }
      } // end match loop
      
      cerr << "...Kept " << matches << " matches with " << sequence.size()
            << " points." << endl;
      pointCount += sequence.size();
    } // end Affinity loop
    
    cerr << "Total of " << matchCount
      << " matches with " << pointCount
      << " points." << endl;
}


void friendOfFriend(double matchRadius, bool useAffinities, string outCatalogName, int minMatches,
                    bool allowSelfMatches, vector<Field*> fields, vector<Instr*> instruments,
                    vector<Expo*> exposures, FTable extensionTable)
{
  
  // Convert matching radius to our system units for world coords
  matchRadius *= ARCSEC/WCS_UNIT;

  
    // Teach PixelMapCollection about all types of PixelMaps it might need to deserialize
    loadPixelMapParser();


    // And a list holding all Points being matched.  Note that catalogs do not make copies
    // of the Points, they just hold pointers to them.  So this is a warehouse.
    // If we wanted to save the memory of the list links we could just let all the Points
    // be zombies that are destroyed when the program ends.
    list<Point> allPoints;

    // If WCS is coming from a serialized PixelMapCollection, I'll keep last-used one around
    // to avoid re-parsing it all the time:
    string currentPixelMapCollectionFileName = "";
    astrometry::PixelMapCollection* inputPmc=0;

    // Loop over input catalogs
    for (long iextn = 0; iextn < extensionTable.nrows(); iextn++) {
      string filename;
      extensionTable.readCell(filename, "FILENAME", iextn);
      stripWhite(filename);
      int hduNumber;
      extensionTable.readCell(hduNumber, "EXTENSION", iextn);

      if (iextn%10==0) cerr << "# Reading object catalog " << iextn
			    << "/" << extensionTable.nrows()
			    << " in " << filename 
			    << " HDU #" << hduNumber
			    << endl;
      FTable ft = FitsTable(filename, FITS::ReadOnly, hduNumber).extract();

      string idKey;
      extensionTable.readCell(idKey, "IDKEY",iextn);
      
      if (stringstuff::nocaseEqual(idKey, "_ROW")) {
        // Want the input table row number as ID column, so make such a column
        vector<long> vi(ft.nrows());
        for (int i=0; i<vi.size(); i++) vi[i] = i;
	      ft.addColumn(vi, "_ROW");
      }

      // Get exposure number, instrument, field, and device numbers
      int exposureNumber;
      extensionTable.readCell(exposureNumber, "Exposure", iextn);
      int instrumentNumber = exposures[exposureNumber]->instrument;
      int fieldNumber = exposures[exposureNumber]->field;
      int deviceNumber;
      extensionTable.readCell(deviceNumber, "Device", iextn);

      // Now get the WCS for this extension
      astrometry::Wcs* wcs = 0;
      string wcsin;
      extensionTable.readCell(wcsin, "WCSIN", iextn);
      list<string> wcssplit = stringstuff::split(wcsin,'@');
      if (wcssplit.size()==2) {
	// @ sign signals that we are getting a map from a PMC
	string wcsName = wcssplit.front();
	string pmcFile = wcssplit.back();
	// Replace any white space in wcs name with underscore:
	wcsName = stringstuff::regexReplace("[[:space:]]+","_",wcsName);

	if (pmcFile!=currentPixelMapCollectionFileName) {
	  // Need to open this PMC file
	  // Try opening file as new serialized PMC
	  ifstream wcsStream(pmcFile.c_str());
	  if (!wcsStream.is_open()) {
	    cerr << "Could not open PixelMapCollection file <" << pmcFile 
		 << ">" << endl;
	    exit(1);
	  }
	  astrometry::PixelMapCollection* tryPMC = new astrometry::PixelMapCollection;
	  if (tryPMC->read(wcsStream)) {
	    // Successfully found serialized file:
	    if (inputPmc) delete inputPmc;
	    inputPmc = tryPMC;
	    currentPixelMapCollectionFileName = pmcFile;
	  } else {
	    cerr << "Failed deserializing PixelMapCollection file <" << pmcFile 
		 << ">" << endl;
	    exit(1);
	  }
	} // Done reading the pmcFile

	// Read the wcs from the collection
	wcs = inputPmc->issueWcs(wcsName)->duplicate();
      } else if (stringstuff::nocaseEqual(wcsin, "_ICRS")) {
	// No mapping will be needed.  Do nothing.
      } else {
	// Assume that wcsin is a FITS WCS specification
	Header wcsHead;
	//	wcsin += "\nEND\n";  //*** hack to fix problem in astropy.io.fits
	istringstream wcsStream(wcsin);
	if (!(wcsStream >> wcsHead)) {
	  cerr << "Error reading WCS/header from WCSIN at extension " << iextn << endl;
	  exit(1);
	}
	wcs = astrometry::readTPV(wcsHead);
	if (!wcs) {
	  cerr << "Failed reading TPV WCS from extension " << iextn << endl;
	  exit(1);
	}
      }

      // Now replace WCSIN in table with a new serialized one, 
      // and set projection to be the field coordinates
      if (wcs) {
        string wcsDump;
        astrometry::PixelMapCollection pmc;
        pmc.learnWcs(*wcs);
        ostringstream oss;
        pmc.writeWcs(oss,wcs->getName());
        wcsDump = oss.str();
        wcs->reprojectTo(*fields[fieldNumber]->projection);
        extensionTable.writeCell(wcsDump, "WCSIN", iextn);
      }

      // Now begin reading the input data

      // Filter the input rows and get the columns we want
      {
        string selectionExpression;
        extensionTable.readCell(selectionExpression, "SELECT", iextn);
        stripWhite(selectionExpression);
        if (!selectionExpression.empty())
          ft.filterRows(selectionExpression);
      }

      if (ft.nrows()<=0) {
        // No data to use in this extension.  Don't read anything.
        continue;
      }
      
      vector<double> vx;
      string key;
      extensionTable.readCell(key, "XKEY", iextn);
      try {
        ft.readCells(vx, key);
      } catch (FTableError& m) {
        // Trap for using float column in source file instead of double:
        vector<float> vf;
        ft.readCells(vf, key);
        vx.resize(vf.size());
        for (int i=0; i<vf.size(); i++) vx[i]=vf[i];
      }
      vector<double> vy;
      extensionTable.readCell(key, "YKEY", iextn);
      try {
        ft.readCells(vy, key);
      } catch (FTableError& m) {
        // Trap for using float column in source file instead of double:
        vector<float> vf;
        ft.readCells(vf, key);
        vy.resize(vf.size());
        for (int i=0; i<vf.size(); i++) vy[i]=vf[i];
      }
      vector<long> vid;
      extensionTable.readCell(key, "IDKEY", iextn);
      try {
        //**/cerr << "Trying long read" << endl;
        ft.readCells(vid, key);
            } catch (FTableError& m) {
        // Trap for using int column in source file instead of long:
        vector<LONGLONG> vidint;
        //**/cerr << "Retrying int" << endl;
        ft.readCells(vidint, key);
        vid.reserve(vidint.size());
        vid.insert(vid.begin(), vidint.begin(), vidint.end());
      }
      //**/cerr << "Done" << endl;
      vector<bool> isStar(ft.nrows(), true);
      {
        string starExpression;
        extensionTable.readCell(starExpression, "STARSELECT", iextn);
        stripWhite(starExpression);
        if (!starExpression.empty())
          ft.evaluate(isStar, starExpression);
      }

      // Select the Catalog(s) to which points from this extension will be added
      string thisAffinity;
      extensionTable.readCell(thisAffinity, "AFFINITY", iextn);

      PointCat* starCatalog = fields[fieldNumber]->catalogFor(stellarAffinity);
      PointCat* galaxyCatalog = (stringstuff::nocaseEqual(stellarAffinity,
                                                          thisAffinity) ?
                                 starCatalog : 
                                 fields[fieldNumber]->catalogFor(thisAffinity));

      // Now loops over objects in the catalog
      for (int iObj = 0; iObj < vx.size(); iObj++) {
        double xpix = vx[iObj];
        double ypix = vy[iObj];
        // Expand bounds of device:
        if (instrumentNumber >=0) 
          (*instruments[instrumentNumber])[deviceNumber]+=Position<double>(xpix,ypix);

        //  maps coords to field's tangent plane
        double xw, yw;
        if (wcs) {
          Assert(wcs);
          wcs->toWorld(xpix, ypix, xw, yw);
        } else {
          // Already have RA, Dec, just project them
          // Assume they are in our std units (see Units.h)
          fields[fieldNumber]->projection->convertFrom(astrometry::SphericalICRS(vx[iObj]*WCS_UNIT, 
                                                                                 vy[iObj]*WCS_UNIT ));
          fields[fieldNumber]->projection->getLonLat(xw, yw);
          // Matching program will use our std unit (Units.h), not radians:
          xw /= WCS_UNIT;
          yw /= WCS_UNIT;
        }

        allPoints.push_back(Point(xw, yw, iextn, vid[iObj], exposureNumber));

        // Now we match!!!!
        if (isStar[iObj])
          starCatalog->add(allPoints.back());
        else
          galaxyCatalog->add(allPoints.back());
      } // end object loop

      if (wcs) delete wcs;
    } // end input extension loop

    if (inputPmc) delete inputPmc;

    cerr << "*** Read " << allPoints.size() << " objects" << endl;



    // Now write all the match catalogs
    long matchCount = 0;
    long pointCount = 0;
    //  Loop over fields
    for (int iField=0; iField<fields.size(); iField++) {
      //  loop over affinities per field
      for (Field::CatMap::iterator iAffinity=fields[iField]->catalogs.begin();
            iAffinity != fields[iField]->catalogs.end();
            ++iAffinity) {
        const PointCat* pcat= iAffinity->second;
        cerr << "Catalog for field " << fields[iField]->name
              << " affinity " << iAffinity->first
              << " with " << pcat->size() << " matches "
              << endl;
        if (pcat->empty()) continue;
        FitsTable ft(outCatalogName, FITS::ReadWrite + FITS::Create, -1);
        ft.setName("MatchCatalog");
        FTable ff=ft.use();
        ff.header()->replace("Field", fields[iField]->name, "Field name");
        ff.header()->replace("FieldNum", iField, "Field number");
        ff.header()->replace("Affinity", iAffinity->first, "Affinity name");

        vector<int> sequence;
        vector<long> extn;
        vector<long> obj;
        long matches=0;
        // Now loop through matches in this catalog
        for (PointCat::const_iterator j=pcat->begin();
              j != pcat->end();
              ++j) {
          // Skip any Match that is below minimum match size
          if ((*j)->size() < minMatches) continue;
          bool selfMatch = false;
          if (!allowSelfMatches) {
            set<long> itsExposures;
            for (list<const Point*>::const_iterator k=(*j)->begin();
              k != (*j)->end();
              ++k) 
            if ( !itsExposures.insert((*k)->exposureNumber).second) {
              selfMatch = true;
              break;
              }
          }
          if (selfMatch) continue;

          // Add elements of this match to the output vectors
          int seq=0;
          ++matches;
          ++matchCount;
          for (list<const Point*>::const_iterator k=(*j)->begin();
                k != (*j)->end();
                ++k, ++seq) {
            sequence.push_back(seq);
            extn.push_back((*k)->extensionNumber);
            obj.push_back((*k)->objectNumber);
          }
        } // end match loop
        cerr << "...Kept " << matches << " matches with " << sequence.size()
              << " points." << endl;
        pointCount += sequence.size();
        ff.addColumn(sequence, "SequenceNumber");
        ff.addColumn(extn, "Extension");
        ff.addColumn(obj, "Object");
      } // end Affinity loop
    } // end Field loop
    
    cerr << "Total of " << matchCount
        << " matches with " << pointCount
        << " points." << endl;

    // Clean up
    cerr << "Cleaning fields: " << endl;
    for (int i=0; i<fields.size(); i++)
      delete fields[i];
    cerr << "Cleaning instruments: " << endl;
    for (int i=0; i<instruments.size(); i++)
      delete instruments[i];
    cerr << "Cleaning exposures: " << endl;
    for (int i=0; i<exposures.size(); i++)
      delete exposures[i];
    cerr << "Done!: " << endl;

  
}
