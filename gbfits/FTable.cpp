//  More code for the FTable

#include "FTable.h"

using namespace img;

void 
FTable::filter(long rowStart, long rowEnd,
	       const vector<string>& regexps) {
  D->filter(D, rowStart, rowEnd, regexps);
}

void 
FTable::filterRows(vector<bool>& vb) {
  D->filterRows(D, vb);
}

void 
FTable::filterRows(const string& expression) {
  vector<bool> vb;
  D->evaluate(vb, expression, H);
  D->filterRows(D, vb);
}

FTable
FTable::extract(long rowStart, long rowEnd,
		const vector<string>& regexps) const {
  FTable result;
  D->filter(result.D, rowStart, rowEnd, regexps);
  result.D->clearChanged();
  result.H->copyFrom(*H);
  return result;
}

FTable
FTable::extractRows(vector<bool>& vb) const {
  FTable result;
  D->filterRows(result.D, vb);
  result.D->clearChanged();
  result.H->copyFrom(*H);
  return result;
}

FTable
FTable::extractRows(const string& expression) const {
  FTable result;
  vector<bool> vb;
  D->evaluate(vb, expression, H);
  D->filterRows(result.D, vb);
  result.D->clearChanged();
  result.H->copyFrom(*H);
  return result;
}


/////////////////////////////////////////////////
/// TableData implementations
/////////////////////////////////////////////////

TableData::~TableData() {
  // Unlock for destruction:
  lock = false;
  for (iterator i=begin(); i!=end(); ++i) delete *i;
}

TableData* 
TableData::duplicate() const {
  TableData* dup = new TableData(rowReserve);
  for (const_iterator i=begin(); i!=end(); ++i) dup->add((*i)->duplicate());
  return dup;
}

void 
TableData::copyFrom(const TableData& rhs) {
  checkLock("TableData::copyFrom()");
  clear();
  rowReserve = rhs.rowReserve;
  rowCount = rhs.rowCount;
  isAltered = true;
  for (const_iterator i=rhs.begin(); i!=rhs.end(); ++i) {
      add((*i)->duplicate());
  }
}      

vector<string> 
TableData::listColumns() const {
  vector<string> out;
  for (const_iterator i = begin();
       i != end();
       ++i) 
    out.push_back((*i)->name());
  return out;
}

bool
TableData::hasColumn(const string& colname) const {
  for (const_iterator i = begin();
       i != end();
       ++i) {
    if (stringstuff::nocaseEqual((*i)->name(),colname))
      return true;
  }
  return false;
}

const ColumnBase* 
TableData::constColumn(string colname) const {
  Index::const_iterator i=columns.find(colname);
  if (i==columns.end()) throw FTableNonExistentColumn(colname);
  return i->second;
}

ColumnBase* 
TableData::operator[](string colname) {
  checkLock("operator[]");
  Index::iterator i=columns.find(colname);
  if (i==columns.end()) throw FTableNonExistentColumn(colname);
  return i->second;
}

void 
TableData::erase(string columnName) {
  checkLock("erase()");
  Index::iterator i=columns.find(columnName);
  if (i==columns.end())
    throw FTableNonExistentColumn(columnName);
  // Destroy the column:
  delete i->second;
  columns.erase(i);
}

void 
TableData::erase(iterator i) {
  checkLock("erase()");
  delete *i;
  columns.erase(i.i);
}

void 
TableData::clear() {
  checkLock("clear()");
  for (Index::iterator i=columns.begin();
       i != columns.end();
       ++i)
    delete i->second;
  columns.clear();
  rowCount = 0;
  rowReserve=0; // Choose to eliminate reservations too. ??
}

void 
TableData::eraseRows(long rowStart, long rowEnd) {
  checkLock("eraseRows()");
  Assert(rowStart>=0);
  if (rowStart>=nrows()) return; // don't erase past end
  if (columns.empty()) return;
  if (rowEnd<0) rowEnd = nrows();
  if (rowStart>=rowEnd) return;
  for (iterator i=begin(); i!=end(); ++i) (*i)->eraseRows(rowStart, rowEnd);
  // Update row count from first column:
  rowCount = (*begin())->size();
}

void 
TableData::insertRows(long insertBeforeRow, long insertNumber) {
  checkLock("insertRows()");
  if (insertBeforeRow > nrows()) throw FTableError("insertRows beyond end of table");
  for (iterator i=begin(); i!=end(); ++i) (*i)->insertRows(insertBeforeRow, insertNumber);
  rowCount += insertNumber;
  Assert(rowCount==(*begin())->size());
}

// Expand all columns to some new size:
void 
TableData::growRows(long targetRows) {
  checkLock("growRows()");
  if (targetRows <= nrows()) return;
  rowCount = targetRows;
  for (iterator i = begin(); i!=end(); ++i)
    (*i)->resize(targetRows);
}

void 
TableData::add(ColumnBase* newColumn) {
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

TableData*
TableData::extract(long rowStart, long rowEnd,
		   const vector<string>& regexps) const {
  if (rowEnd < 0) rowEnd = nrows();
  if (rowEnd < rowStart) 
    throw FTableError("TableData::extract() with rowEnd < rowStart");
  TableData* td = new TableData(rowEnd - rowStart);
  filter(td, rowStart, rowEnd, regexps);
  return td;
}

void
TableData::filter(TableData* td, long rowStart, long rowEnd,
		  const vector<string>& regexps) const {
  if (rowEnd < 0) rowEnd = nrows();
  // rowStart and rowEnd checked before entry here.
  td->checkLock("filter()");
  // Get set of all columns for output:

  bool writeSelf = (this==td);
  if (!writeSelf) td->clear();
  td->touch();

  // Tell destination table the size that all columns will end up with:
  td->rowCount = rowEnd - rowStart;

  vector<string> allColumns = listColumns();
  set<string> getThese = stringstuff::findMatches(regexps, allColumns);

  for (vector<string>::iterator i = allColumns.begin();
       i != allColumns.end();
       ++i) {
    bool keeper = getThese.find(*i) != getThese.end();
    if ( keeper) {
      ColumnBase* shorter = (*this)[*i]->copyRows(rowStart, rowEnd);
      if (writeSelf) td->erase(*i);
      td->add(shorter);
    } else {
      // Unwanted column
      if (writeSelf) td->erase(*i);
    }
  }
}

void
TableData::filterRows(TableData* td, vector<bool>& vb) const {
  td->checkLock("filter()");
  // Get set of all columns for output:
  vector<string> allColumns = listColumns();

  bool writeSelf = (this==td);
  if (!writeSelf) td->clear();
  td->touch();

  // Count up the number of output rows:
  int outRows = 0;
  for (int i=0; i<vb.size(); i++)
    if (vb[i]) outRows++;

  // Tell destination table the size that all columns will end up with:
  td->rowCount = outRows;

  for (vector<string>::iterator i = allColumns.begin();
       i != allColumns.end();
       ++i) {
    ColumnBase* shorter = (*this)[*i]->copyRows(vb, outRows);
    if (writeSelf) td->erase(*i);
    td->add(shorter);
  }
}
