// Read/Write FTables from FITS binary tables (ASCII works too?)
// Right now limited interface with no headers
#ifndef FITSTABLE_H
#define FITSTABLE_H

#include "FTable.h"
#include "FITS.h"

namespace fits {
  class FITSTable {
  private:
    FITS::FitsFile parent;
    int hduNumber;
    long nrows;
    int ncols;
    const bool writeable;

    // Make this the current HDU, return status:
    int moveTo() const;

    // Find columns matching template and add them to vector (no duplicates)
    void findColumns(img::FTable& ft,
		     string matchMe,
		     vector<string>& names,
		     vector<int>& numbers) const;
    // Create new column(s) in img::FTable to hold data in FITSTable's column with given names and no's
    void createColumns(img::FTable& ft, 
		       const vector<string>& names,
		       const vector<int>& numbers) const;
    // Read data from FITS table file into the FTable, colunns with specified names/numbers.
    void readData(img::FTable& ft, 
		  const vector<string>& names,
		  const vector<int>& numbers,
		  long rowStart,
		  long rowEnd) const;
    // Template function used in above to read data for a single column
    template <class T>
    void getFitsColumnData(img::FTable ft, string colName, int icol, 
			 long rowStart, long rowEnd) const;

    // For writing FTable to FITS
    // Add all of ft's columns to end of table, and accumulate name/number pairs
    void addFitsColumns(img::FTable ft, vector<string>& colNames, vector<int>& colNums);
    // Write data from FTable data to FITS file, for requested row range
    template <class DT>
    void writeFitsColumn(img::FTable ft, string colName, int colNums,
			 long rowStart=0, long rowEnd=-1);

  public:
    // Open table at given extn, filename
    FITSTable(string filename,
	      FITS::Flags f=FITS::ReadOnly,
	      int hduNumber_=1);
    // Close on destruction
    ~FITSTable();
    string getFilename() const {return parent.getFilename();}
    // Extract all columns and range of rows (rowEnd=1-past-end, zero-indexed, -1=go to end)
    bool isWriteable() const {return writeable;}

    img::FTable extract(long rowStart=0, long rowEnd=-1) const;
    // Extract columns matching any of the input strings (FITSIO wildcards OK)
    img::FTable extract(const vector<string>& templates, long rowStart=0, long rowEnd=-1) const;

    // Empty the FITS table:
    void clear();

    // Change an FTable to be able to store as FITS bintable.  Namely, 
    // change any string cells with variable-length arrays into fixed-length at 
    // the maximum size:
    static void makeFitsCompatible(img::FTable ft);

    // Write contents of FTable into this FITSTable:
    void replaceWith(img::FTable ft);
  };


    // Specialize FITS-reading for string and for logical
    template <>
    void FITSTable::getFitsColumnData<string>(img::FTable ft, string colName, int icol, 
				   long rowStart, long rowEnd) const;
    template <>
    void FITSTable::getFitsColumnData<bool>(img::FTable ft, string colName, int icol, 
				 long rowStart, long rowEnd) const;
    template<>
    void FITSTable::writeFitsColumn<string>(img::FTable ft, string colName, int colNum,
					    long rowStart, long rowEnd);
    template<>
    void FITSTable::writeFitsColumn<bool>(img::FTable ft, string colName, int colNum,
					  long rowStart, long rowEnd);

} // end namespace img
#endif
