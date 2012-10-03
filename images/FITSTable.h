// Read/Write FTables from FITS binary tables (ASCII works too?)
// Right now limited interface with no headers
#ifndef FITSTABLE_H
#define FITSTABLE_H

#include "FTable.h"
#include "Hdu.h"

namespace FITS {
  class FitsTable: public Hdu {
  private:
    // Adds a TableData to the HDU
    mutable img::TableData* dptr;
    // And its link counter
    mutable int* dcount;

  public:
    // Open table at given extn, filename
    FitsTable(string filename,
	      FITS::Flags f=FITS::ReadOnly,
	      int hduNumber_=1);
    FitsTable(string filename,
	      FITS::Flags f,
	      string hduName);
    // Close on destruction
    virtual ~FitsTable();

    // Write any changes in memory version (mirror) back to FITS file
    // [happens automatically on object destruction]
    virtual void flush() const {}
    virtual void flush() {if (isWriteable()) {flushData(); Hdu::flush();}}
    virtual bool isChanged() const {return Hdu::isChanged() || (dptr && dptr->isChanged());}
    // Empty the FITS table and header:
    virtual void clear() {
      checkWriteable("FitsTable::clear()");
      Hdu::clear();	// clear header
      if (dptr) {
	// Clear out mirror if we've got one
	dptr->clear();
      } else {
	// Make ourselves an empty mirror (which will be flushed later)
	dptr = new img::TableData();
	dcount = new int(1);
      }
      dptr->touch();
    }

    img::FTable extract(long rowStart=0, long rowEnd=-1) const;
    // Extract columns matching any of the input strings (FITSIO wildcards OK)
    img::FTable extract(const vector<string>& templates, long rowStart=0, long rowEnd=-1) const;

    // Issue FTable that is mirrored to FITS file data.  
    // If FITS HDU is read-only, FTable will be locked.
    img::FTable use();

    // Copy values of external FTable into current header & data
    void copy(img::FTable ft);

    // Link the FitsTable to an external FTable.  Note that this break connection
    // to any FTable that was obtained via use() command, and also
    // there will be an exception if you destroy this FitsTable before destroying
    // all other FTables that are sharing these data.
    void adopt(img::FTable ft);

    // Change an FTable to be able to store as FITS bintable.  Namely, 
    // change any string cells with variable-length arrays into fixed-length at 
    // the maximum size:
    static void makeFitsCompatible(img::FTable ft);

  private:
    void flushData();	// Push Table data back to FITS file
    img::TableData* loadData() const; // Turn the FITS file into a TableData structure
    // Find columns in FITS table matching template and add them to vector (no duplicates)
    void findFitsColumns(string matchMe,
			 vector<string>& names,
			 vector<int>& numbers) const;
    // Create new column(s) in TableData to hold data in 
    // FitsTable's column with given names and no's
    void createColumns(img::TableData* tptr, 
		       const vector<string>& names,
		       const vector<int>& numbers) const;
    // Read data from FITS table file into the FTable, columns with specified names/numbers.
    void readFitsData(img::TableData* tptr, 
		      const vector<string>& names,
		      const vector<int>& numbers,
		      long rowStart,
		      long rowEnd) const;
    // Template function used in above to read data for a single column
    template <class T>
    void getFitsColumnData(img::TableData* tptr, string colName, int icol, 
			   long rowStart, long rowEnd) const;

    // For writing FTable to FITS
    // Add all of TableData's columns to end of FITS table, and accumulate name/number pairs
    void addFitsColumns(const img::TableData* tptr, 
			vector<string>& colNames, vector<int>& colNums);
    // Write data from FTable data to FITS file, for requested row range
    template <class DT>
    void writeFitsColumn(const img::TableData* tptr, string colName, int colNums,
			 long rowStart=0, long rowEnd=-1);

    void writeFitsTableData(const img::TableData* tptr);
    void clearFitsTableData();
  };


  // Specialize FITS-reading for string and for logical
  template <>
  void FitsTable::getFitsColumnData<string>(img::TableData* tptr, string colName, int icol, 
					    long rowStart, long rowEnd) const;
  template <>
  void FitsTable::getFitsColumnData<bool>(img::TableData* tptr, string colName, int icol, 
					  long rowStart, long rowEnd) const;
  template<>
  void FitsTable::writeFitsColumn<string>(const img::TableData* tptr, string colName, int colNum,
					  long rowStart, long rowEnd);
  template<>
  void FitsTable::writeFitsColumn<bool>(const img::TableData* tptr, string colName, int colNum,
					long rowStart, long rowEnd);

} // end namespace FITS
#endif
