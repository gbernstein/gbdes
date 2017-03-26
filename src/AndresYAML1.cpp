// Quick program to reformat Andres's tree rings files into YAML template format
// For case of a single ring pattern per CCD.
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include <cmath>

#include "yaml-cpp/yaml.h"

using namespace std;

int
main(int argc, char *argv[])
{
  ifstream centers(argv[1]);

  YAML::Emitter out;
  out << YAML::BeginMap;

  string buffer1;
  while (getline(centers, buffer1)) {
    // Loop for CCDs.
    // Skip comments
    {
      string junk;
      istringstream iss(buffer1);
      iss >> junk;
      if (junk.empty() || junk[0]=='#') continue;
    }
    // Read CCD name, number, and ring center
    string device;
    int ccdno;
    double xcen;
    double ycen;
    istringstream iss(buffer1);
    iss >> device >> ccdno >> xcen >> ycen;

    for (int ii = 0; ii < 2; ii++) {
      // Do left, then right *** not this time.
      if (ii>0) continue;
      //      string side = (ii==0) ? "left" : "right";
      string inFile = "PHOTO1/rings_pattern_flat_" + device + "_g.dat";
      // _" + side + ".dat";
      //      string tableName = device + (ii==0 ? "_low" : "_high");
      string tableName = device;
      ifstream ifs(inFile.c_str());

      double nodeStart;
      double nodeStep;
      double val0;
      double val1;

      string buffer;
    
      vector<double> r;
      vector<double> v;

      // ** Skip two lines
      getline(ifs,buffer);
      getline(ifs,buffer);

      while (getline(ifs, buffer)) {
	// Skip comments
	{
	  string junk;
	  istringstream iss(buffer);
	  iss >> junk;
	  if (junk.empty() || junk[0]=='#') continue;
	}
	istringstream iss(buffer);
	double rr, vv;
	iss >> rr >> vv;
	r.push_back(rr);
	v.push_back(vv);
      }
      nodeStart = r.front();
      nodeStep = (r.back() - r.front())/(r.size()-1);

      for (int i = 0; i<r.size(); i++) {
	// Check that nodes are equally spaced
	if (abs( r[i] - nodeStart - i*nodeStep) > 0.01) {
	  cerr << "not equally spaced: file " << inFile << endl;
	  cerr << "r = " << r[i]
	       << " n = " << i
	       << " nodes at " << nodeStart << " " << nodeStep
	       << " expect " << nodeStart + i*nodeStep
	       << " error " << r[i] - (nodeStart + i*nodeStep)
	       << endl;
	  exit(1);
	}
      }
      out << YAML::Key << tableName 
	  << YAML::Value << YAML::BeginMap
	  << YAML::Key << "Axis" << YAML::Value << "R"
	  << YAML::Key << "XCenter" << YAML::Value << xcen
	  << YAML::Key << "YCenter" << YAML::Value << ycen
	  << YAML::Key << "ArgStart" << YAML::Value << nodeStart
	  << YAML::Key << "ArgStep" << YAML::Value << nodeStep
	  << YAML::Key << "Values" << YAML::Flow << YAML::Value << v
	  << YAML::EndMap;
    } // end left/right loop
  } // end ccd loop
  out << YAML::EndMap;
  cout << out.c_str() << endl;
}

