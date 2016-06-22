// Code for fitting parameters of defaulted pixel maps to the staring WCS solution.

    


	// Initialize the new PixelMap so that it is best fit to the canonical exposures'
	// map from pixel coordinates to world coords=projection in exposure frame
	extn.startWcs->reprojectTo(*expo.projection);

	{
	  const int nGridPoints=512;	// Number of test points for map initialization
	  cerr << "Initializing pixel map " << inst->mapNames[idev] << endl;

	  Bounds<double> b=inst->domains[idev];
	  // Distribute points equally in x and y, but shuffle the y coords
	  // so that the points fill the rectangle
	  vector<int> vx(nGridPoints);
	  for (int i=0; i<vx.size(); i++) vx[i]=i;
	  vector<int> vy = vx;
	  std::random_shuffle(vy.begin(), vy.end());
	  double xstep = (b.getXMax()-b.getXMin())/nGridPoints;
	  double ystep = (b.getYMax()-b.getYMin())/nGridPoints;
	  list<Detection*> testPoints;

	  for (int i=0; i<vx.size(); i++) {
	    double xpix = b.getXMin() + (vx[i]+0.5)*xstep;
	    double ypix = b.getXMin() + (vy[i]+0.5)*ystep;
	    double xw, yw;
	    extn.startWcs->toWorld(xpix, ypix, xw, yw);
	    Detection* d = new Detection;
	    d->xpix = xpix;
	    d->ypix = ypix;
	    d->xw = xw;
	    d->yw = yw;
	    testPoints.push_back(d);
	  }

	  const double testPointSigma = 0.01*ARCSEC/DEGREE; // "errors" on world coords of test points
	  mapFit(testPoints, pm, testPointSigma);

	  // Clear out the testPoints:
	  for (list<Detection*>::iterator i = testPoints.begin();
	       i != testPoints.end();
	       ++i)
	    delete *i;
	} // Done initializing parameters of new PixelMap

	// Save this map into the collection
	mapCollection.learnMap(*pm);
	delete pm;
      } // end loop creating PixelMaps for all Devices.

      *****/
      // Now create an exposure map for all exposures with this instrument,
      // and set initial parameters by fitting to input WCS map.

      for (list<long>::const_iterator i= exposuresWithInstrument.begin();
	   i != exposuresWithInstrument.end();
	   ++i) {
	int iexp = *i;
	Exposure& expo = *exposures[iexp];
	PixelMap* warp=0;
	// First we check to see if we are going to use an existing exposure map:
	string loadFile = expo.name;
	if (true) { //****existingMapFinder(loadFile)) {
	  // Adopt the existing map for this exposure:
	  PixelMapCollection pmcExpo;
	  ifstream ifs(loadFile.c_str());
	  if (!ifs) {
	    cerr << "Could not open serialized map file " 
		 << loadFile
		 << endl;
	    exit(1);
	  }
	  pmcExpo.read(ifs);
	  warp = pmcExpo.issueMap(expo.name);
	  mapCollection.learnMap(*warp);
	  expo.name = warp->getName();
	  // If this imported map is to be held fixed, add it to the list:
	  if (regexMatchAny(fixMapList, expo.name))
	    fixedMapNames.insert(expo.name);

	  // We're done with this exposure if we have pre-set map for it.
	  continue;
	}

	if (iexp == canonicalExposure && noDevicesFixed) {
	  /**/cerr << "Giving identity map to canonical exposure " << expo.name <<endl;
	  // Give this canonical exposure identity map to avoid degeneracy with Instrument
	  expo.name = IdentityMap().getName();
	  continue;
	}
	
	// We will create a new exposure map and then initialize it
	// ??? Potentially more sophisticated here
	/***/int exposureOrder=1;
	if (exposureOrder==1) 
	  warp = new LinearMap(expo.name);
	else
	  warp = new PolyMap(exposureOrder, expo.name, worldTolerance);
	expo.name = expo.name;
	// Set numerical-derivative step to 1 arcsec, polynomial is working in degrees:
	warp->setPixelStep(ARCSEC/DEGREE);

	// Build test points
	list<Detection*> testPoints;

	const int nGridPoints=800;	// Number of test points for map initialization
	int pointsPerDevice = nGridPoints / inst->nDevices + 1;

	// Get points from each extension that is part of this exposure
	for (long iextn=0; iextn<extensions.size(); iextn++) {
	  if (extensions[iextn]->exposure != iexp) continue;
	  Extension& extn = *extensions[iextn];
	  int idev = extn.device;
	  Bounds<double> b=inst->domains[idev];
	  double step = sqrt(b.area()/pointsPerDevice);
	  int nx = static_cast<int> (ceil((b.getXMax()-b.getXMin())/step));
	  int ny = static_cast<int> (ceil((b.getYMax()-b.getYMin())/step));
	  double xstep = (b.getXMax()-b.getXMin())/nx;
	  double ystep = (b.getYMax()-b.getYMin())/ny;
	  // starting pixel map for the first exposure with this instrument
	  if (!extn.startWcs) {
	    cerr << "Failed to get starting Wcs for device "
		 << inst->deviceNames.nameOf(idev)
		 << " in exposure " << expo.name
		 << endl;
	    exit(1);
	  }
	  extn.startWcs->reprojectTo(*expo.projection);
	  PixelMap* devpm = mapCollection.issueMap(inst->mapNames[idev]);
	
	  for (int ix=0; ix<=nx; ix++) {
	    for (int iy=0; iy<=ny; iy++) {
	      double xpix = b.getXMin() + ix*xstep;
	      double ypix = b.getYMin() + iy*ystep;
	      double xw, yw;
	      // Map coordinates into the exposure projection = world:
	      extn.startWcs->toWorld(xpix, ypix, xw, yw);
	      // And process the pixel coordinates through Instrument map:
	      devpm->toWorld(xpix,ypix,xpix,ypix);
	      Detection* d = new Detection;
	      d->xpix = xpix;
	      d->ypix = ypix;
	      d->xw = xw;
	      d->yw = yw;
	      testPoints.push_back(d);
	    }
	  } // finish collecting test points for this extension
	} // finish extension loop for this exposure
	
	double testPointSigma = 0.01*ARCSEC/DEGREE; // 10 mas scale for "errors" on the test point positions
	mapFit(testPoints, warp, testPointSigma);
	
	// Insert the Exposure map into the collection
	mapCollection.learnMap(*warp);

	// Check for fixed maps
	if (regexMatchAny(fixMapList, warp->getName())) {
	  cerr << "WARNING: Requested fixed parameters for Exposure map "
	       << warp->getName()
	       << " that does not have initial values.  Ignoring request to fix params."
	       << endl;
	}

	delete warp;

	// Clear out the testPoints:
	for (list<Detection*>::iterator i = testPoints.begin();
	     i != testPoints.end();
	     ++i)
	  delete *i;
      } // end solution for this exposure map

      /**/cerr << "Done with map initialization for instrument " << inst->name << endl;
    } // End instrument loop

    // Freeze the parameters of all the fixed maps
    mapCollection.setFixed(fixedMapNames, true);

    // Now construct Wcs and a reprojected SubMap for every extension

    // and save SubMap for it
    for (int iext=0; iext<extensions.size(); iext++) {
    } // Extension loop
