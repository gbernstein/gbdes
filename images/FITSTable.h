// Read/Write FTables from FITS binary tables (ASCII works too?)
// Right now limited interface with no headers
#ifndef FITSTABLE_H
#define FITSTABLE_H

#include "FTable.h"
#include "FITS.h"

namespace fits {
  class FitsTable: public Hdu {
  public:
    // Open table at given extn, filename
    FitsTable(string filename,
	      FITS::Flags f=FITS::ReadOnly,
	      int hduNumber_=1);
    FitsTable(string filename,
	      FITS::Flags f=FITS::ReadOnly,
	      string hduName);
    // Close on destruction
    virtual ~FitsTable();

    img::FTable extract(long rowStart=0, long rowEnd=-1) const;
    // Extract columns matching any of the input strings (FITSIO wildcards OK)
    img::FTable extract(const vector<string>& templates, long rowStart=0, long rowEnd=-1) const;

    // Issue FTable that is mirrored to FITS file data
    img::FTable use();

    // Write any changes in memory version (mirror) back to FITS file
    // [happens automatically on object destruction]
    virtual void flush() const;

    // Empty the FITS table and header:
    virtual void clear();

    // Replace current header & data with that in external file.
    void copy(FTable ft);

    // Change an FTable to be able to store as FITS bintable.  Namely, 
    // change any string cells with variable-length arrays into fixed-length at 
    // the maximum size:
    static void makeFitsCompatible(img::FTable ft);

  private:
    // Have a copy of the whole table here if we have anyone using it.
    mutable img::FTable* mirror;

    // Find columns in FITS table matching template and add them to vector (no duplicates)
    void findFitsColumns(string matchMe,
			 vector<string>& names,
			 vector<int>& numbers) const;
    // Create new column(s) in img::FTable to hold data in 
    // FitsTable's column with given names and no's
    void createColumns(img::FTable& ft, 
		       const vector<string>& names,
		       const vector<int>& numbers) const;
    // Read data from FITS table file into the FTable, columns with specified names/numbers.
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

    void writeFitsTableData(img::FTable ft);
    void clearFitsTableData(img::FTable ft);
  };


  // Specialize FITS-reading for string and for logical
  template <>
  void FitsTable::getFitsColumnData<string>(img::FTable ft, string colName, int icol, 
					    long rowStart, long rowEnd) const;
  template <>
  void FitsTable::getFitsColumnData<bool>(img::FTable ft, string colName, int icol, 
					  long rowStart, long rowEnd) const;
  template<>
  void FitsTable::writeFitsColumn<string>(img::FTable ft, string colName, int colNum,
					  long rowStart, long rowEnd);
  template<>
  void FitsTable::writeFitsColumn<bool>(img::FTable ft, string colName, int colNum,
					long rowStart, long rowEnd);

} // end namespace FITS
#endif
