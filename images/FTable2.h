// $Id: FTable2.h,v 1.4 2012/09/08 20:10:46 garyb Exp $
// In-memory version of the FITS table structure.
// This file contains classes / declarations that should not be needed
// by the user.

// Define RANGE_CHECK to have row arguments always checked:
//#define RANGE_CHECK

#ifndef FTABLE2_H
#define FTABLE2_H

#include <vector>
#include <map>
#include <typeinfo>

#include "StringStuff.h"

// sorts & filters

namespace img {

  class ColumnBase {
  public:
    ColumnBase(string name): colName(name) {}
    virtual ~ColumnBase() {}
    virtual ColumnBase* duplicate() const=0;
    const string& name() const {return colName;}
    // Number of rows of data
    virtual long size() const =0;
    virtual void reserve(long rowsReserved) =0;
    virtual void resize(long newSize) =0;

    // Return 1 (or 0) for scalar cells, >1 for fixed-length array cell, <0 for var-length
    virtual long repeat() const=0;
    // Return max allowed length of string-valued elements of column, or -1 to be free
    virtual long stringLength() const=0;
    virtual FITS::DataType elementType() const =0;  // FITS DataType for scalar or array element
    virtual char columnCode() const =0;  // FITSIO code for appropriate column data storage

    virtual void eraseRows(long rowStart, long rowEnd) =0;
    // insertBeforeRow = -1 will add row(s) to the end.
    // Return value is number of rows after insertion.
    virtual long insertRows(long insertBeforeRow, long insertNumber) =0;
  private:
    const string colName;
    // Ancillary information ???
  };

  template <class DT>
  class ScalarColumn: public ColumnBase {
  protected:
    vector<DT> v;  // Holds the column's data
    // Protected function to check length of input strings
    void checkLength(const string& input) {
      if (input.size() > length) 
	throw FTableError("Input string <" + input + "> too long at column " + name());
    }
  private:
    const FITS::DataType dType;
    const char colCode;
    long length;   // max length of string, if this is holding strings
  public:
    ScalarColumn(string name, long length_=-1): 
      ColumnBase(name), dType(FITS::FITSTypeOf<DT>()), 
      colCode(FITS::ColumnCode<DT>()), length(length_) {}
    ScalarColumn(const vector<DT>& in, string name, long length_=-1):
      ColumnBase(name), dType(FITS::FITSTypeOf<DT>()), 
      colCode(FITS::ColumnCode<DT>()), v(in), length(length_) {}
    virtual ~ScalarColumn() {}
    virtual ColumnBase* duplicate() const {return new ScalarColumn(*this);}
    virtual long size() const {return v.size();}
    virtual void reserve(long rowsReserved) {v.reserve(rowsReserved);}
    virtual void resize(long newSize) {v.resize(newSize);}
    virtual long repeat() const {return 1;}  // scalar default
    virtual long stringLength() const {return length;} 
    virtual FITS::DataType elementType() const {return dType;}
    virtual char columnCode() const {return colCode;}

    // Range checking will be done at TableData level, so none here
    virtual void readCell(DT& value, long row) const {value = v[row];}
    virtual void readCells(vector<DT>& values, long rowStart=0, long rowEnd=-1) const {
      int nCopy = (rowEnd >=0 ? rowEnd : v.size()) - rowStart;
      if (nCopy < 0) nCopy = 0;
      values.resize(nCopy);
      for (int i=0; i<nCopy; i++)
	values[i] = v[i+rowStart];
    }

    // Writes return the length of column after writing.
    virtual long writeCell(const DT& value, long row) {
      // extend column to hold data
      if (row >= v.size()) v.resize(row+1);
      v[row] = value;
      return v.size();
    }
    // Or a range: length of input vector determines rows altered; grow table if needed.
    virtual long writeCells(const vector<DT>& values, long rowStart=0) {
      long requiredSize = rowStart + values.size();
      if (requiredSize > v.size()) v.resize(requiredSize);
      for (int i=0; i<values.size(); i++)
	v[i+rowStart] = values[i];
      return v.size();
    }

    virtual void eraseRows(long rowStart, long rowEnd) {
      Assert(rowStart > 0);
      Assert(rowEnd > rowStart);
      Assert(rowEnd <= v.size());
      v.erase(v.begin()+rowStart, v.begin()+rowEnd);
    }

    // insertBeforeRow = -1 will add row(s) to the end.
    virtual long insertRows(long insertBeforeRow, long insertNumber) {
      if (insertBeforeRow < 0 || insertBeforeRow >= v.size())
	v.insert(v.end(), insertNumber, DT());
      else
	v.insert(v.begin()+insertBeforeRow, insertNumber, DT());
      return v.size();
    }
  };

  // Specialization for string-valued scalar arrays, check length as needed
  template<>
  inline ScalarColumn<string>::ScalarColumn(const vector<string>& in, string name, 
				     long length_): 
    ColumnBase(name), dType(FITS::FITSTypeOf<string>()), 
    colCode(FITS::ColumnCode<string>()), length(length_) {
    if (stringLength()>=0) 
      for (int i=0; i<in.size(); i++)
	checkLength(in[i]);
    v = in;
  }

  template<>
  long
  inline ScalarColumn<string>::writeCell(const string& value, long row) {
    if (stringLength()>=0) checkLength(value);
    // extend column to hold data
    if (row >= v.size()) v.resize(row+1);
    v[row] = value;
    return v.size();
  }

  template<>
  long
  inline ScalarColumn<string>::writeCells(const vector<string>& values, long rowStart) {
    if (stringLength()>=0) 
      for (int i=0; i<values.size(); i++)
	checkLength(values[i]);
    long requiredSize = rowStart + values.size();
    if (requiredSize > v.size()) v.resize(requiredSize);
    for (int i=0; i<values.size(); i++)
      v[i+rowStart] = values[i];
    return v.size();
  }
  
  template <class DT>
  class ArrayColumn: public ScalarColumn<vector<DT> > {
  private:
    typedef ScalarColumn<vector<DT> > Base;
    const FITS::DataType dType2;
    const char colCode2;
  protected:
    using ScalarColumn<vector<DT> >::v;
  public:
    ArrayColumn(string name, long length_=-1): 
      Base(name, length_), dType2(FITS::FITSTypeOf<DT>()),
      colCode2(FITS::ColumnCode<DT>()) {}
    ArrayColumn(const vector<vector<DT> >& in, string name, long length_=-1): 
      Base(in, name, length_), dType2(FITS::FITSTypeOf<DT>()),
      colCode2(FITS::ColumnCode<DT>()) {}
    virtual ColumnBase* duplicate() const {return new ArrayColumn(*this);}
    virtual ~ArrayColumn() {}

    // For variable-length array, the only thing we need to override is 
    // the repeat() count and the elementType should be DT, not the vector<DT>:
    virtual long repeat() const {return -1;}
    virtual FITS::DataType elementType() const {return dType2;}
    virtual char columnCode() const {return colCode2;}

    // Need to declare these so they can be specialized for string:
    long writeCell(const vector<DT>& value, long row) {
      return Base::writeCell(value,row);}
    long writeCells(const vector<vector<DT> >& values, long rowStart=0) {
      return Base::writeCells(values, rowStart);}

  };

  // Override the constructor & writing routines for Array column to check length
  // of strings in arrays of strings
  template<>
  inline ArrayColumn<string>::ArrayColumn(const vector<vector<string> >& in, string name,
				   long length_):
    Base(name, length_), dType2(FITS::FITSTypeOf<string>()),
    colCode2(FITS::ColumnCode<string>()) {
    if (stringLength()>=0)
      for (int i=0; i<in.size(); i++)
	for (int j=0; j<in[i].size(); j++)
	  checkLength(in[i][j]);
    v = in;
  }

  template<>
  long
  inline ArrayColumn<string>::writeCell(const vector<string>& value, long row) {
    if (stringLength()>=0) 
      for (int i=0; i<value.size(); i++)
	checkLength(value[i]);
    return Base::writeCell(value, row);
  }

  template <>
  long
  inline ArrayColumn<string>::writeCells(const vector<vector<string> >& values, long rowStart) {
    if (stringLength()>=0) {
      for (int i=0; i<values.size(); i++) 
	for (int j=0; j<values[i].size(); j++) 
	  checkLength(values[i][j]);
    }
    return Base::writeCells(values, rowStart);
  }
  
  template <class DT>
  class FixedArrayColumn: public ArrayColumn<DT> {
  private:
    typedef ArrayColumn<DT> Base;
    using ArrayColumn<DT>::v;
    const int width;
    virtual ColumnBase* duplicate() const {return new FixedArrayColumn(*this);}
    // Function to pad input array to required width; exception if input too long
    void checkWidth(const vector<DT>& in) {
      if (in.size() > width)
	FormatAndThrow<FTableError>() 
	  << "Input of fixed-length array cell is too long for column "
	  << ColumnBase::name()
	  << " width " << width
	  << " input size " << in.size();
    }
  public:
    FixedArrayColumn(string name, long repeat, long length_=-1): 
      Base(name, length_), width(repeat) {}
    // Initialize with data assumes that all elements have size()<=repeat
    FixedArrayColumn(const vector<vector<DT> >& in, string name, long repeat, long length_=-1): 
      Base(name, length_), width(repeat) {
      v.resize(in.size());
      for (long i=0; i<in.size(); i++)
	writeCell(in[i], i);
      }

    virtual long repeat() const {return width;}
  
    virtual long writeCell(const vector<DT>& value, long row) {
      checkWidth(value);
      v[row].reserve(width);
      // Use normal variable-length array writer
      long out = Base::writeCell(value, row);
      // Fill out cell with null values as needed to meet fixed array length
      v[row].resize(width);
      return out;
    }

    // Override variable-array version with this one to make sure each row
    // is not too long, and is padded out
    virtual long writeCells(const vector<vector<DT> >& values, long rowStart=0) {
      // extend column to hold data
      if (rowStart + values.size() > v.size()) 
	v.resize(rowStart+values.size(), vector<DT>(width));
      for (int i=0; i<values.size(); i++)
	writeCell(values[i], i+rowStart);
      return v.size();
    }

    // override insert to put in proper-length vectors to cells
    // insertBeforeRow = -1 will add row(s) to the end.
    virtual long insertRows(long insertBeforeRow, long insertNumber) {
      if (insertBeforeRow < 0 || insertBeforeRow >= v.size())
	v.insert(v.end(), insertNumber, vector<DT>(width));
      else
	v.insert(v.begin()+insertBeforeRow, insertNumber, vector<DT>(width));
      return v.size();
    }
  };

  class TableData {
  private:
    typedef std::map<string,ColumnBase*,stringstuff::LessNoCase> Index;
    Index columns;	// This is where the columns are stored.
    long rowReserve;	// New columns reserve at least this much space
    long rowCount;	// Keep this number and all column lengths in synch 
    bool isAltered;
    bool lock;
  public:
    TableData(long rowReserve_=0): rowReserve(rowReserve_), rowCount(0),
				   isAltered(false), lock(false)  {}
    ~TableData() {
      for (iterator i=begin(); i!=end(); ++i) delete *i;
    }
    TableData* duplicate() const {
      TableData* dup = new TableData(rowReserve);
      for (const_iterator i=begin(); i!=end(); ++i) dup->add((*i)->duplicate());
      return dup;
    }

    void copyFrom(const TableData& rhs) {
      checkLock("TableData::copyFrom()");
      clear();
      rowReserve = rhs.rowReserve;
      rowCount = rhs.rowCount;
      isAltered = true;
      for (const_iterator i=rhs.begin(); i!=rhs.end(); ++i) add((*i)->duplicate());
    }      
      
    long nrows() const {return rowCount;}
    int ncols() const {return columns.size();}

    bool isChanged() const {return isAltered;}
    void clearChanged() {isAltered = false;}
    void touch() {isAltered = true;}
    bool isLocked() const {return lock;}
    void setLock() {lock = true;}

    // **** Make iterator classes to run over columns ***
    class iterator {
      friend class TableData;
    private:
      explicit iterator(Index::iterator i_): i(i_) {}
      Index::iterator i;
    public:
      iterator& operator++() {++i; return *this;}
      iterator& operator--() {--i; return *this;}
      bool operator==(const iterator& rhs) const {return i==rhs.i;}
      bool operator!=(const iterator& rhs) const {return i!=rhs.i;}
      ColumnBase* operator*() const {return i->second;}
    };
    class const_iterator {
      friend class TableData;
    private:
      explicit const_iterator(Index::const_iterator i_): i(i_) {}
      Index::const_iterator i;
    public:
      const_iterator(const iterator& rhs): i(rhs.i) {}
      const_iterator& operator++() {++i; return *this;}
      const_iterator& operator--() {--i; return *this;}
      bool operator==(const const_iterator& rhs) const {return i==rhs.i;}
      bool operator!=(const const_iterator& rhs) const {return i!=rhs.i;}
      const ColumnBase* operator*() const {return i->second;}
    };

    const_iterator begin() const {return const_iterator(columns.begin());}
    const_iterator end() const {return const_iterator(columns.end());}

    const_iterator const_begin() {return const_iterator(columns.begin());}
    const_iterator const_end()   {return const_iterator(columns.end());}

    iterator begin() {checkLock(); touch(); return iterator(columns.begin());}
    iterator end() {checkLock(); touch(); return iterator(columns.end());}

    const ColumnBase* constColumn(string colname) const {
      Index::const_iterator i=columns.find(colname);
      if (i==columns.end()) throw FTableNonExistentColumn(colname);
      return i->second;
    }

    const ColumnBase* operator[](string colname) const {
      return constColumn(colname);
    }

    ColumnBase* operator[](string colname) {
      checkLock(); touch();
      Index::iterator i=columns.find(colname);
      if (i==columns.end()) throw FTableNonExistentColumn(colname);
      return i->second;
    }

    // ***** Add / remove columns:

    // Add new column given the data.
    template<class T>
    void addColumn(const vector<T>& values, string columnName, 
		   long repeat=-1, long stringLength=-1 ) {
      checkLock();   touch();
      add(new ScalarColumn<T>(values, columnName, stringLength));
    }
    // Specialization to recognize vector<vector> as array column
    template<class T>
    void addColumn(const vector<vector<T> >& values, string columnName, 
		   long repeat=-1, long stringLength=-1) {
      checkLock();   touch();
      if (repeat<0) add(new ArrayColumn<T>(values, columnName, stringLength));
      else add(new FixedArrayColumn<T>(values, columnName, repeat, stringLength));
    }
        
    void erase(string columnName) {
      checkLock();   touch();
      Index::iterator i=columns.find(columnName);
      if (i==columns.end())
	throw FTableNonExistentColumn(columnName);
      // Destroy the column:
      delete i->second;
      columns.erase(i);
    }
    void erase(iterator i) {
      checkLock();   touch();
      delete *i;
      columns.erase(i.i);
    }

    void clear() {
      checkLock();   touch();
      for (Index::iterator i=columns.begin();
	   i != columns.end();
	   ++i)
	delete i->second;
      columns.clear();
      rowCount = 0;
      rowReserve=0; // Choose to eliminate reservations too. ??
    }

    // ***** Row erase / insert:
    void eraseRows(long rowStart, long rowEnd=-1) {
      checkLock();
      touch();
      Assert(rowStart>=0);
      if (rowStart>=nrows()) return; // don't erase past end
      if (columns.empty()) return;
      if (rowEnd<0) rowEnd = nrows();
      if (rowStart>=rowEnd) return;
      for (iterator i=begin(); i!=end(); ++i) (*i)->eraseRows(rowStart, rowEnd);
      // Update row count from first column:
      rowCount = (*begin())->size();
    }
    void insertRows(long insertBeforeRow, long insertNumber) {
      checkLock();
      touch();
      if (insertBeforeRow > nrows()) throw FTableError("insertRows beyond end of table");
      for (iterator i=begin(); i!=end(); ++i) (*i)->insertRows(insertBeforeRow, insertNumber);
      rowCount += insertNumber;
      Assert(rowCount==(*begin())->size());
    }

    // ******* Data access:
    // Access single element (copy created)
    template <class T>
    void readCell(T& value, string columnName, long row) const {
      rangeCheck(row);
      const ScalarColumn<T>* col = dynamic_cast<const ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch reading column " + columnName);
      col->readCell(value, row);
    }
    // Or a range: rowEnd is one-past-last row; -1 means go to end.
    template <class T>
    void readCells(vector<T>& values, string columnName, 
		   long rowStart=0, long rowEnd=-1) const {
      rangeCheck(rowStart);
      if (rowEnd>0) rangeCheck(rowEnd-1);
      const ScalarColumn<T>* col = dynamic_cast<const ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch reading column " + columnName);
      col->readCells(values, rowStart, rowEnd);
    }

    // Writing single or multiple cells:
    template <class T>
    void writeCell(const T& value, string columnName, long row) {
      if (row<0) rangeCheck(row);
      checkLock();
      touch();
      ScalarColumn<T>* col = dynamic_cast<ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCell(value, row);
      growRows(doneRows);
    }

    template <class T>
    void writeCells(const vector<T>& values, string columnName, long rowStart=0) {
      if (rowStart<0) rangeCheck(rowStart);
      checkLock();
      touch();
      ScalarColumn<T>* col = dynamic_cast<ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCells(values, rowStart);
      growRows(doneRows);
    }

    // Specialization to recognize array-valued Cells
    template <class T>
    void writeCell(const vector<T>& value, string columnName, long row) {
      if (row<0) rangeCheck(row);
      checkLock();
      touch();
      ArrayColumn<T>* col = dynamic_cast<ArrayColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCell(value, row);
      growRows(doneRows);
    }

    template <class T>
    void writeCells(const vector<vector<T> >& values, string columnName, long rowStart=0) {
      if (rowStart<0) rangeCheck(rowStart);
      checkLock();
      touch();
      ArrayColumn<T>* col = dynamic_cast<ArrayColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      touch();
      int doneRows = col->writeCells(values, rowStart);
      growRows(doneRows);
    }

  private:
    // No copy or assignment:
    TableData(const TableData& rhs);
    void operator=(const TableData& rhs);

    // call for anything that can alter data
    // Note not a const routine since anything that is forbidden by lock should
    // not be coded in a const routine.
    void checkLock(const string& s="") {
      if (isLocked()) throw FTableLockError(s);
    }
    // call for anything that can alter data
    void rangeCheck(long row, const string& s="") const {
#ifdef RANGE_CHECK
      if (row<0 || row >= nrows())
	FormatAndThrow<FTableRangeError>() << " accessing " << row << " of " << nrows()
					   << " " << s;
#endif
    }

    // Expand all columns to some new size:
    void growRows(long targetRows) {
      if (targetRows <= nrows()) return;
      rowCount = targetRows;
      for (iterator i = begin(); i!=end(); ++i)
	(*i)->resize(targetRows);
    }

    void add(ColumnBase* newColumn) {
      string addname = newColumn->name();
      if (columns.find(addname) != columns.end())
	throw FTableError("Adding column with duplicate name <" + addname + ">");
      // Pad column to current table length, or vice-versa, to keep all columns same length
      if (newColumn->size() < nrows()) newColumn->resize(nrows());
      else if (newColumn->size() > nrows()) growRows(newColumn->size());
      // Reserve requested space
      newColumn->reserve(rowReserve);
      columns.insert(std::pair<string, ColumnBase*>(addname, newColumn));
    }
  };

} // namespace img
#endif // FTABLE2_H
