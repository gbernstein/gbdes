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
    const FITS::Flags flags;

    // Make this the current HDU, return status:
    int moveTo() const;

    // Find columns matching template and add them to vector (no duplicates)
    void findColumns(FTable& ft,
		     string matchMe,
		     vector<string>& names,
		     vector<int>& numbers) const;
    // Create new column(s) in FTable to hold data in FITSTable's column with given names and no's
    void createColumns(FTable& ft, 
		       const vector<string>& names,
		       const vector<int>& numbers) const;
    // Read data from FITS table file into the FTable, colunns with specified names/numbers.
    void readData(FTable& ft, 
		  const vector<string>& names,
		  const vector<int>& numbers,
		  long rowStart,
		  long rowEnd) const;
    // Template function used in above to read data for a single column
    template <class T>
    void getFitsColumnData(FTable ft, string colName, int icol, 
			 long rowStart, long rowEnd) const;

    // For writing FTable to FITS
    // Add all of ft's columns to end of table, and accumulate name/number pairs
    void addFitsColumns(FTable ft, vector<string>& colNames, vector<int>& colNums);
    // Write data from FTable data to FITS file, for requested row range
    template <class DT>
    void writeFitsColumn(FTable ft, string colName, int colNums,
			 long rowStart=0, long rowEnd=-1);

  public:
    // Open table at given extn, filename (read only right now***)
    FITSTable(string filename,
	      FITS::Flags f=FITS::ReadOnly,
	      int hduNumber_=1);
    // Close on destruction
    ~FITSTable();
    string getFilename() const {return fname;}
    // Extract all columns and range of rows (rowEnd=1-past-end, zero-indexed, -1=go to end)
    FTable extract(long rowStart=0, long rowEnd=-1) const;
    // Extract columns matching any of the input strings (FITSIO wildcards OK)
    FTable extract(const vector<string>& templates, long rowStart=0, long rowEnd=-1) const;

    // Empty the FITS table:
    void clear();

    // Change an FTable to be able to store as FITS bintable.  Namely, 
    // change any string cells with variable-length arrays into fixed-length at 
    // the maximum size:
    static void makeFitsWritable(FTable ft);

    // Write contents of FTable into this FITSTable:
    void replaceWith(FTable ft);
  };


    // Specialize FITS-reading for string and for logical
    template <>
    void FITSTable::getFitsColumnData<string>(FTable ft, string colName, int icol, 
				   long rowStart, long rowEnd) const;
    template <>
    void FITSTable::getFitsColumnData<bool>(FTable ft, string colName, int icol, 
				 long rowStart, long rowEnd) const;
    template<>
    void FITSTable::writeFitsColumn<string>(FTable ft, string colName, int colNum,
					    long rowStart, long rowEnd);
    template<>
    void FITSTable::writeFitsColumn<bool>(FTable ft, string colName, int colNum,
					  long rowStart, long rowEnd);

} // end namespace img
#endif
