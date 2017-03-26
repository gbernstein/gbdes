// Saving some code here so can later have a standalone program that will output
// a TPV-format FITS WCS that is fit to an arbitrary device's map.

// Choose exposure name, WCS YAML file, and output file(s) on command line
// Get the center coordinates of the exposure
// Get WCS solutions for each device from this exposure
// For each device:
//    * Fit TPV to it
//    * Make into header
//    * Calculate per-CCD deviations from TPV; convert back into pixel shifts?
//    * Save these shifts, plus pixel area corrections, to FITS extensions (truncate accuracy)
//    * Put the TPVs into the header extensions, save in FITS.
//    * Report RMS and max of the corrections.


    // Save the new header files
    
    // Keep track of those we've already written to:
    set<string> alreadyOpened;

    // accuracy of SCAMP-format fit to real solution:
    const double SCAMPTolerance=0.0001*ARCSEC/DEGREE;
  
    for (int iext=0; iext<extensions.size(); iext++) {
      Extension& extn = *extensions[iext];
      Exposure& expo = *exposures[extn.exposure];
      if (expo.instrument < 0) continue;
      Instrument& inst = *instruments[expo.instrument];
      // Open file for new ASCII headers  ??? Save serialized per exposure too??
      string newHeadName = extn.tpvOutFile;
      bool overwrite = false;
      if ( alreadyOpened.find(newHeadName) == alreadyOpened.end()) {
	overwrite = true;
	alreadyOpened.insert(newHeadName);
      }
      ofstream newHeads(newHeadName.c_str(), overwrite ? 
			ios::out : (ios::out | ios::in | ios::ate));
      if (!newHeads) {
	cerr << "WARNING: could not open new header file " << newHeadName << endl;
	cerr << "...Continuing with program" << endl;
	continue;
      }
      // Fit the TPV system to a projection about the field's pointing center:
      expo.projection->setLonLat(0.,0.);
      Wcs* tpv = fitTPV(inst.domains[extn.device],
			*extn.wcs,
			*expo.projection,
			"NoName",
			SCAMPTolerance);
      newHeads << writeTPV(*tpv);
      delete tpv;
    }
