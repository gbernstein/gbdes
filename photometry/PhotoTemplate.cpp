// Radial and linear lookup-table photometric adjustments
#include <fstream>
#include <sstream>

#include "PhotoTemplate.h"
#include "StringStuff.h"

using namespace photometry;
using namespace stringstuff;


void
LookupTable::read(istream& is) {
  values.clear();
  string buffer;
  if (!getlineNoComment(is, buffer))
    throw PhotometryError("LookupTable::read finds no node start/step specs: ");
  {
    istringstream iss(buffer);
    if (!(iss >> nodeStart >> nodeStep)) 
      throw PhotometryError("LookupTable::read has bad node start/step specs: " + buffer);
  }
  while (getlineNoComment(is, buffer)) {
    double value;
    istringstream iss(buffer);
    if (!(iss >> value))
      throw PhotometryError("Bad data line reading LookupTable: " + buffer);
    values.push_back(value);
  }
}

double
LookupTable::lookup(double arg) const {
  double nodeCount = (arg-nodeStart) / nodeStep;
  int index = static_cast<int> (std::floor(nodeCount));
  if (index < 0) {
    return values.front();
  } else if (index >= values.size()-1) {
    return values.back();
  } else {
    double f = nodeCount - index;
    return values[index] * (1.-f) + values[index+1] * f;
  }
}

PhotoTemplate1d::PhotoTemplate1d(string sourceFile, 
				 string name): PhotoMap(name),
  filename(sourceFile) {
  read();
}

PhotoMap*
PhotoTemplate1d::create(std::istream& is, string name) {
  // We just read a filename and read the file to deserialize
  string sourceFile;
  double factor;
  if (!(is >> sourceFile >> factor))
    throw PhotometryError("Missing filename/scaling in PhotoTemplate1d::create()");
  PhotoMap* retval = new PhotoTemplate1d(sourceFile);
  retval->setParams(DVector(1,factor));
  return retval;
}

void
PhotoTemplate1d::write(std::ostream& os, int precision) const {
  // Serialization of this PhotoMap is just the name of the source file and the scaling factor
  os << filename << " " << scaling << endl;
}

void
PhotoTemplate1d::read() {
  ifstream ifs(filename.c_str());
  if (!ifs)
    throw PhotometryError("Could not open PhotoTemplate1d source file <" + filename + ">");
  string buffer;
  if (!getlineNoComment(ifs, buffer))
    throw PhotometryError("Missing data in PhotoTemplate1d source file <" + filename + ">");
  string axis;
  istringstream iss(buffer);
  if (!(iss >> axis))
      throw PhotometryError("Bad axis line in PhotoTemplate1d source file <" + filename + ">"
			    ": " + buffer);
  if (nocaseEqual(axis, "X"))
    applyToX = true;
  else if (nocaseEqual(axis, "Y"))
    applyToX = false;
  else 
    throw PhotometryError("Bad PhotoTemplate1d axis specification '" + axis + "' in source file <"
			  + filename + ">");

  // Read rest of table through the LookupTable:
  table.read(ifs);
}

PhotoRings::PhotoRings(string sourceLeft, 
		       string sourceRight, 
		       string name): PhotoMap(name),
				     xSplit(1024.),
				     fileLeft(sourceLeft),
				     fileRight(sourceRight) {
  read();
}

PhotoMap*
PhotoRings::create(std::istream& is, string name) {
  // reconstruct from the 2 file names and the scaling factor
  string sourceLeft;
  string sourceRight;
  double factor;
  if (!(is >> sourceLeft >> sourceRight >> factor))
    throw PhotometryError("Missing filename(s) or scaling in PhotoRings::create()");
  PhotoMap* retval = new PhotoRings(sourceLeft,sourceRight);
  retval->setParams(DVector(1,factor));
  return retval;
}

void
PhotoRings::write(std::ostream& os, int precision) const {
  // Simply write the names of the two template files and the scaling factor
  os << fileLeft << " " << fileRight << " " << scaling << endl;
}

void
PhotoRings::read() {
  for (int i=0; i<2; i++) {
    // Do left, then right
    bool isLeft = (i==0);
    string filename = isLeft ? fileLeft : fileRight;
    ifstream ifs(filename.c_str());
    if (!ifs)
      throw PhotometryError("Could not open PhotoRings source file <" + filename + ">");
    string buffer;
    if (!getlineNoComment(ifs, buffer))
      throw PhotometryError("Missing data in PhotoRings source file <" + filename + ">");
    string axis;
    double x0, y0;
    istringstream iss(buffer);
    if (!(iss >> axis >> x0 >> y0))
      throw PhotometryError("Bad axis line in PhotoRings source file <" + filename + ">"
			    ": " + buffer);
    if (!nocaseEqual(axis, "R"))
      throw PhotometryError("PhotoRings file must have axis specification 'R', has '" + axis +
			    "' in file <" + filename + ">");
    if (isLeft) {
      xCenterLeft = x0;
      yCenterLeft = y0;
      tableLeft.read(ifs);
    } else {
      xCenterRight = x0;
      yCenterRight = y0;
      tableRight.read(ifs);
    }
  }
}

double
PhotoRings::lookup(double xPix, double yPix) const {
  if (xPix < xSplit) {
    double radius = sqrt( (xPix - xCenterLeft)*(xPix - xCenterLeft) +
			  (yPix - yCenterLeft)*(yPix - yCenterLeft) );
    /**/cerr << "left, radius " << radius << " x " << xPix << endl;
    return tableLeft.lookup(radius);
  } else {
    double radius = sqrt( (xPix - xCenterRight)*(xPix - xCenterRight) +
			  (yPix - yCenterRight)*(yPix - yCenterRight) );
    /**/cerr << "right, radius " << radius << " x " << xPix << endl;
    return tableRight.lookup(radius);
  }
}

