// Code to evaluate expressions from FTable columns

#include "FTable.h"
#include "Expressions.h"

using namespace img;

class TableToken: public expression::Token;

// Converter function.  Will specialize those that are possible.  Return false if conversion
// would not be wise.
template <class T>
bool ValueConverter(vector<T>& vec, const expressions::Value* val) {return false;}

// Expression evaluation, convert to type T.  Use Header hh to evaluate scalars.
template <class T>
void 
TableData::evaluate(vector<T>& result,
		    const string& expression,
		    const Header* hh) const {

  TableToken tokenReader(this, hh);

  // tokenize
  std::list<expressions::Token*> 
    tokenList = tokenize(expression, tokenReader);

  // parse
  expression::Evaluable* eval = parse(tokenList,
				      tokenList.begin(),
				      tokenList.end(),
				      expression);

  expression::Value* val = eval.evaluate();

  if (!ValueConverter(result, val)) 
    throw FTableError("Non-convertible datatype result from expression");

}

class ColumnToken: public expression::Token {
public:
  ColumnToken(const TableData* tptr_, const string& columnName):
    tptr(tptr_), colName(columnName) {}
  virtual expression::Evaluable* createEvaluable() const {
    // read data from the column, whichever type is possible
  }
};

  

class TableToken: public expression::Token {
public:
  TableToken(const TableData* tptr_, const Header* hptr_):
    tptr(tptr_), hptr(hptr_) {}
  const set<string>& columnNames() const {return columns;}
  virtual expression::Token* createFromString(const std::string& input,
					      size_t& begin, size_t& end,
					      bool lastTokenWasOperator) const {
    Assert(end>begin);
    if (input[begin]=='@') {
      // Read a header keyword as a constant
      string keyword = input.substr(begin+1, end);
      {
	bool val;
	if (hptr->getValue(keyword, val)) 
	  return new ConstantEvaluable(val);
      }
      {
	string val;
	if (hptr->getValue(keyword, val)) 
	  return new ConstantEvaluable(val);
      }
      {
	int val;
	if (hptr->getValue(keyword, val)) 
	  return new ConstantEvaluable((long) val);
      }
      {
	double val;
	if (hptr->getValue(keyword, val)) 
	  return new ConstantEvaluable(val);
      }
      throw FTableError("FTable expression evaluation cannot find Header keyword "
			+ keyword);
    } else {
      string colName = input.substr(begin,end);
      columns.insert(colName);
      return new ColumnToken(tptr, colName);
    }
  }
private:
  const TableData* tptr;
  const Header* hptr;
  std::set<std::string> columns;
};
