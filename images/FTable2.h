// $Id: FTable2.h,v 1.4 2012/09/08 20:10:46 garyb Exp $
// In-memory version of the FITS table structure.
// This file contains classes / declarations that should not be needed
// by the user.

#ifndef FTABLE2_H
#define FTABLE2_H

#include <vector>
#include <map>
#include <typeinfo>

#include "Image.h"


// ??? TODO:
// Case-insensitive column names
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
    long length;   // max length of string, if this is holding strings
  public:
    ScalarColumn(string name, long length_=-1): 
      ColumnBase(name), dType(FITS::FITSTypeOf<DT>()), length(length_) {}
    ScalarColumn(const vector<DT>& in, string name, long length_=-1):
      ColumnBase(name), dType(FITS::FITSTypeOf<DT>()), v(in), length(length_) {}
    virtual ~ScalarColumn() {}
    virtual ColumnBase* duplicate() const {return new ScalarColumn(*this);}
    virtual long size() const {return v.size();}
    virtual void reserve(long rowsReserved) {v.reserve(rowsReserved);}
    virtual void resize(long newSize) {v.resize(newSize);}
    virtual long repeat() const {return 1;}  // scalar default
    virtual long stringLength() const {return length;} 
    virtual FITS::DataType elementType() const {return dType;}

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
  ScalarColumn<string>::StringColumn(const vector<string>& in, string name, 
				     long length_=-1): 
    Base(name, length_), dType(FITS::FITSTypeOf<string>()) {
    if (stringLength()>=0) 
      for (int i=0; i<in.size(); i++)
	checkLength(in[i]);
    v = in;
  }

  template<>
  virtual long
  ScalarColumn<string>::writeCell(const string& value, long row) {
    if (stringLength()>=0) checkLength(value);
    // extend column to hold data
    if (row >= v.size()) v.resize(row+1);
    v[row] = value;
    return v.size();
  }

  template<>
  virtual long
  ScalarColumn<string>::writeCells(const vector<string>& values, long rowStart=0) {
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
    const FITS::DataType dType;
  public:
    ArrayColumn(string name, long length_=-1): Base(name, length_), dType(FITS::FITSTypeOf<DT>()) {}
    ArrayColumn(const vector<vector<DT> >& in, string name, long length=-1): 
      Base(in, name, length_), dType(FITS::FITSTypeOf<DT>()) {}
    virtual ColumnBase* duplicate() const {return new ArrayColumn(*this);}
    virtual ~ArrayColumn() {}
    // For variable-length array, the only thing we need to override is 
    // the repeat() count and the elementType should be DT, not the vector<DT>:
    virtual long repeat() const {return -1;}
    virtual FITS::DataType elementType() const {return dType;}
  };

  // Override the constructor & writing routines for Array column to check length
  // of strings in arrays of strings
  template<>
  ArrayColumn<string>::ArrayColumn(const vector<vector<string> >& in, string name,
				   long length_=-1):
    Base(name, length_), dType(FITS::FITSTypeOf<string>()) {
    if (stringLength()>=0)
      for (int i=0; i<in.size(); i++)
	for (int j=0; j<in[i].size(); j++)
	  checkLength(in[i][j]);
    v = in;
  }

  template<>
  ArrayColumn<string>::WriteCell(const vector<string>& value, long row) {
    if (stringLength()>=0) 
      for (int i=0; i<value.size(); i++)
	checkLength(value[i]);
    // extend column to hold data
    if (row >= v.size()) v.resize(row+1);
    v[row] = value;
  }

  template <>
  ArrayColumn<string>::writeCells(const vector<vector<string> >& values, long rowStart=0) {
    if (stringLength()>=0) {
      for (int i=0; i<values.size(); i++) 
	for (int j=0; j<values[i].size(); j++) 
	  checkLength(values[i][j]);
    long requiredSize = rowStart + values.size();
    if (requiredSize > v.size()) v.resize(requiredSize);
    for (int i=0; i<values.size(); i++)
      v[i+rowStart] = values[i];
    return v.size();
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
	WriteCell(in[i], i);
    }
    virtual long repeat() const {return width;}

    virtual writeCell(const vector<string>& value, long row) {
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
    virtual writeCells(const vector<vector<string> >& values, long rowStart=0) {
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
    typedef std::map<string,ColumnBase*> Index;
    Index columns;	// This is where the columns are stored.
    long rowReserve;	// New columns reserve at least this much space
    long rowCount;	// Keep this number and all column lengths in synch 
  public:
    TableData(long rowReserve_=0): rowReserve(rowReserve_), rowCount(0) {}
    ~TableData() {
      for (iterator i=begin(); i!=end(); ++i) delete *i;
    }
    TableData* duplicate() const {
      TableData* dup = new TableData(rowReserve);
      for (const_iterator i=begin(); i!=end(); ++i) dup->add((*i)->duplicate());
      return dup;
    }
    long nrows() const {return rowCount;}
    int ncols() const {return columns.size();}

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
    iterator begin() {return iterator(columns.begin());}
    iterator end() {return iterator(columns.end());}

    const ColumnBase* operator[](string colname) const {
      Index::const_iterator i=columns.find(colname);
      if (i==columns.end()) throw FTableNonExistentColumn(colname);
      return i->second;
    }

    ColumnBase* operator[](string colname) {
      Index::iterator i=columns.find(colname);
      if (i==columns.end()) throw FTableNonExistentColumn(colname);
      return i->second;
    }


    // ***** Add / remove columns:

    // Add new column given the data.
    template<class T>
    void addColumn(const vector<T>& values, string columnName, long repeat=-1) {
      add(new ScalarColumn<T>(values, columnName, repeat));
    }
    // Specialization to recognize vector<vector> as array column
    template<class T>
    void addColumn(const vector<vector<T> >& values, string columnName, long repeat=-1) {
      if (repeat<0) add(new ArrayColumn<T>(values, columnName));
      else add(new FixedArrayColumn<T>(values, columnName, repeat));
    }
    // Another specialization to recognize string input
    //    template<>
    void addColumn(const vector<string>& values, string columnName, long repeat=-1) {
      add(new StringColumn(values, columnName, repeat));
    }
        
    void erase(string columnName) {
      Index::iterator i=columns.find(columnName);
      if (i==columns.end())
	throw FTableNonExistentColumn(columnName);
      // Destroy the column:
      delete i->second;
      columns.erase(i);
    }
    void erase(iterator i) {
      delete *i;
      columns.erase(i.i);
    }

    void clear() {
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
      if (insertBeforeRow > nrows()) throw FTableError("insertRows beyond end of table");
      for (iterator i=begin(); i!=end(); ++i) (*i)->insertRows(insertBeforeRow, insertNumber);
      rowCount += insertNumber;
      Assert(rowCount==(*begin())->size());
    }

    // ******* Data access:
    // Access single element (copy created)
    template <class T>
    void readCell(T& value, string columnName, long row) const {
#ifdef RANGE_CHECK
      if (row<0 || row >= rowCount)
	  FormatAndThrow<FTableRangeError>() << row << " for rowCount= " << rowCount;
#endif
      const ScalarColumn<T>* col = dynamic_cast<const ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch reading column " + columnName);
      col->readCell(value, row);
    }
    // Or a range: rowEnd is one-past-last row; -1 means go to end.
    template <class T>
    void readCells(vector<T>& values, string columnName, long rowStart=0, long rowEnd=-1) const {
#ifdef RANGE_CHECK
      if (rowStart<0 || rowStart>=rowCount)
	  FormatAndThrow<FTableRangeError>() << rowStart << " for rowCount= " << rowCount;
      if (rowEnd>rowCount)
	  FormatAndThrow<FTableRangeError>() << rowEnd-1 << " for rowCount= " << rowCount;
#endif
      const ScalarColumn<T>* col = dynamic_cast<const ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch reading column " + columnName);
      col->readCells(values, rowStart, rowEnd);
    }

    // Writing single or multiple cells:
    template <class T>
    void writeCell(const T& value, string columnName, long row) {
#ifdef RANGE_CHECK
      if (row<0)
	  FormatAndThrow<FTableRangeError>() << row << " for rowCount= " << rowCount;
#endif
      ScalarColumn<T>* col = dynamic_cast<ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCell(value, row);
      growRows(doneRows);
    }

    template <class T>
    void writeCells(const vector<T>& values, string columnName, long rowStart=0) {
#ifdef RANGE_CHECK
      if (rowStart<0)
	  FormatAndThrow<FTableRangeError>() << row << " for rowCount= " << rowCount;
#endif
      ScalarColumn<T>* col = dynamic_cast<ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCells(values, rowStart);
      growRows(doneRows);
    }

    // Specialization to recognize array-valued Cells
    template <class T>
    void writeCell(const vector<T>& value, string columnName, long row) {
#ifdef RANGE_CHECK
      if (row<0)
	  FormatAndThrow<FTableRangeError>() << row << " for rowCount= " << rowCount;
#endif
      ArrayColumn<T>* col = dynamic_cast<ArrayColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCell(value, row);
      growRows(doneRows);
    }

    template <class T>
    void writeCells(const vector<vector<T> >& values, string columnName, long rowStart=0) {
#ifdef RANGE_CHECK
      if (rowStart<0)
	  FormatAndThrow<FTableRangeError>() << row << " for rowCount= " << rowCount;
#endif
      ArrayColumn<T>* col = dynamic_cast<ArrayColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCells(values, rowStart);
      growRows(doneRows);
    }

    // More specialization to recognize string-valued Cells
    void writeCell(const string& value, string columnName, long row) {
#ifdef RANGE_CHECK
      if (row<0)
	  FormatAndThrow<FTableRangeError>() << row << " for rowCount= " << rowCount;
#endif
      StringColumn* col = dynamic_cast<StringColumn*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCell(value, row);
      growRows(doneRows);
    }

    void writeCells(const vector<string>& values, string columnName, long rowStart=0) {
#ifdef RANGE_CHECK
      if (rowStart<0)
	  FormatAndThrow<FTableRangeError>() << row << " for rowCount= " << rowCount;
#endif
      StringColumn* col = dynamic_cast<StringColumn*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCells(values, rowStart);
      growRows(doneRows);
    }

  private:
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
