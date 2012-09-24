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
    void findColumns(FTable& ft,
		     string matchMe,
		     vector<string>& names,
		     vector<int>& numbers) const;
    void createColumns(FTable& ft, 
		       const vector<string>& names,
		       const vector<int>& numbers) const;
    void readData(FTable& ft, 
		  const vector<string>& names,
		  const vector<int>& numbers,
		  long rowStart,
		  long rowEnd) const;
    template <class T>
    void getFitsColumnData(FTable ft, string colName, int icol, 
			 long rowStart, long rowEnd) const;

    // Make this the current HDU, return status:
    int moveTo() const;
  public:
    // Open table at given extn, filename (read only right now***)
    FITSTable(string filename, int hduNumber_=2);
    // Close on destruction
    ~FITSTable();
    string getFilename() const {return fname;}
    // Extract all columns and range of rows (rowEnd=1-past-end, zero-indexed, -1=go to end)
    FTable extract(long rowStart=0, long rowEnd=-1) const;
    // Extract columns matching any of the input strings (FITSIO wildcards OK)
    FTable extract(const vector<string>& templates, long rowStart=0, long rowEnd=-1) const;
  };

    // Specialize FITS-reading for string and for logical
    template <>
    void FITSTable::getFitsColumnData<string>(FTable ft, string colName, int icol, 
				   long rowStart, long rowEnd) const;
    template <>
    void FITSTable::getFitsColumnData<bool>(FTable ft, string colName, int icol, 
				 long rowStart, long rowEnd) const;
} // end namespace img
#endif
