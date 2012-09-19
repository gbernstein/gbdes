// Read/Write FTables from FITS binary tables (ASCII works too?)
// Right now limited interface with no headers
#ifndef FITSTABLE_H
#define FITSTABLE_H

#include "FTable.h"
#include "FITS.h"

namespace img {
  class FITSTable {
  private:
    FITS::fitsfile* fptr;
    string fname;
    int hduNumber;
    long nrows;
    int ncols;
    // Find columns matching template and add them to vector (no duplicates)
    void findColumns(FITSTable& ft,
		     string template,
		     vector<string>& names,
		     vector<int>& numbers) const;
    void createColumns(FITSTable& ft, 
		       const vector<string>& names,
		       const vector<int>& numbers) const;
    void readData(FITSTable& ft, 
		  const vector<string>& names,
		  const vector<int>& numbers,
		  long rowStart,
		  long rowEnd) const;
    template <class T>
    void 
    getFitsColumnData(FTable ft, string colName, int icol, 
		      long rowStart, long rowEnd) const;
    // Specialize for string
    template <>
    void 
    getFitsColumnData<string>(FTable ft, string colName, int icol, 
			      long rowStart, long rowEnd) const;
    // Specialize for logical
    template <>
    void 
    getFitsColumnData<bool>(FTable ft, string colName, int icol, 
			    long rowStart, long rowEnd) const;
    // Make this the current HDU:
    void moveTo() const;
  public:
    // Open table at given extn, filename (read only right now***)
    FITSTable(string filename, int hduNumber_=2);
    // Close on destruction
    ~FITSTable();
    // Extract all columns and range of rows (rowEnd=1-past-end, zero-indexed, -1=go to end)
    FTable extract(long rowStart=0; long rowEnd=-1) const;
    // Extract columns matching any of the input strings (FITSIO wildcards OK)
    FTable extract(const vector<string>& templates, Start=0; long rowEnd=-1) const;
  };

} // end namespace img
#endif
