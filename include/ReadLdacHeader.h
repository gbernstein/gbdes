// Function to read the header stuffed into character array cell of an LDAC file
#include "FitsTable.h"
#include <sstream>

inline img::Header ReadLdacHeader(const string& filename, int hduNumber) {
  FITS::FitsTable readit(filename, FITS::ReadOnly, hduNumber);
  img::FTable ft = readit.extract();
  std::vector<std::string> vs;
  ft.readCell(vs, "Field Header Card",0);
  string buffer;
  for (int i=0; i<vs.size(); i++) buffer+= vs[i] + "\n";
  buffer += "END\n";
  istringstream iss(buffer);
  img::Header result;
  iss >> result;
  return result;
}
