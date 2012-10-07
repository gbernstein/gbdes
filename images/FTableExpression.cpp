// Code to evaluate expressions from FTable columns

#include "FTable.h"
#include "Expressions.h"

using namespace img;
using namespace expressions;

class TableToken;

// Converter function.  Will specialize those that are possible.  Return false if conversion
// would not be wise.
// I will allow conversion of any type to bool following usual rule that empty/zero = false.
// Also allow up-conversion of long to double, but no down-conversions.
template <class T>
bool ValueConverter(vector<T>& vec, const expressions::Value* val) {return false;}

class ColumnEvaluable: public expressions::Evaluable {
public:
  ColumnEvaluable(const TableData* tptr_, const string& columnName):
    tptr(tptr_), colName(columnName) {}
  virtual expressions::Value* evaluate() const {
    const ColumnBase* col = (*tptr)[colName];
    vector<bool> vb;
    vector<string> vs;
    vector<long> vl;
    vector<double> vd;
    if (col->makeVector(vb)) 
      return new expressions::VectorValue<bool>(vb);
    else if (col->makeVector(vs)) 
      return new expressions::VectorValue<string>(vs);
    else if (col->makeVector(vl)) 
      return new expressions::VectorValue<long>(vl);
    else if (col->makeVector(vd))
      return new expressions::VectorValue<double>(vd);
    else
      throw FTableError("Could not convert FTable column " + col->name()
			+ " to VectorValue for expression evaluation");
  }
  virtual expressions::Value* returnEmptyEvaluation() const {
    const ColumnBase* col = (*tptr)[colName];
    vector<bool> vb;
    vector<string> vs;
    vector<long> vl;
    vector<double> vd;
    if (col->makeVector(vb)) 
      return new expressions::VectorValue<bool>();
    else if (col->makeVector(vs)) 
      return new expressions::VectorValue<string>();
    else if (col->makeVector(vl)) 
      return new expressions::VectorValue<long>();
    else if (col->makeVector(vd)) 
      return new expressions::VectorValue<double>();
    else
      throw FTableError("Could not convert FTable column " + col->name()
			+ " to VectorValue for expression evaluation");
  }
private:
  const TableData* tptr;
  const string colName;
};

class TableToken: public expressions::Token {
public:
  TableToken(const TableData* tptr_, const Header* hptr_):
    tptr(tptr_), hptr(hptr_) {}
  const set<string>& columnNames() const {return columns;}
  virtual expressions::Token* createFromString(const std::string& input,
					      size_t& begin, size_t& end,
					      bool lastTokenWasOperator) const {
    Assert(end>begin);
    if (input[begin]=='@') {
      // Read a header keyword as a constant
      string keyword = input.substr(begin+1, end-(begin+1));
      Evaluable* eval=0;
      {
	bool val;
	if (hptr->getValue(keyword, val)) 
	  eval =  new ConstantEvaluable<ScalarValue<bool> >(val);
      }
      {
	string val;
	if (hptr->getValue(keyword, val)) 
	  eval =  new ConstantEvaluable<ScalarValue<string> >(val);
      }
      {
	int val;
	if (hptr->getValue(keyword, val)) 
	  eval =  new ConstantEvaluable<ScalarValue<long> >((long) val);
      }
      {
	double val;
	if (hptr->getValue(keyword, val)) 
	  eval =  new ConstantEvaluable<ScalarValue<double> >(val);
      }
      if (!eval) 
	throw FTableError("FTable expression evaluation cannot find Header keyword "
			  + keyword);
      return new EvaluableToken(eval);
    } else {
      string colName = input.substr(begin,end-begin);
      columns.insert(colName);
      return new EvaluableToken( new ColumnEvaluable(tptr, colName));
    }
  }
private:
  const TableData* tptr;
  const Header* hptr;
  mutable std::set<std::string> columns;
};

// Expression evaluation, convert to type T.  Use Header hh to evaluate scalars.
template <class T>
void 
TableData::evaluate(vector<T>& result,
		    const string& expression,
		    const Header* hh) const {

  TableToken tokenReader(this, hh);

  // tokenize
  std::list<expressions::Token*> 
    tokenList = expressions::tokenize(expression, tokenReader);

  // parse
  std::list<expressions::Token*>::iterator b = tokenList.begin();
  expressions::Evaluable* eval = expressions::parse(tokenList,
						    b,
						    tokenList.end(),
						    expression);

  expressions::Value* val = eval->evaluate();

  if (!ValueConverter(result, val)) 
    throw FTableError("Non-convertible datatype result from expression");

}

// Special the conversion from Expression Value to a vector
template <>
bool ValueConverter(vector<bool>& vec, const expressions::Value* val) {
  if ( const VectorValue<bool>* vptr 
       = dynamic_cast<const VectorValue<bool>*> (val) )  {
    vec = vptr->values;
    return true;
  } else if ( const VectorValue<string>* vptr
	      = dynamic_cast<const VectorValue<string>*> (val) )  {
    vec.resize(vptr->values.size());
    for (int i=0; i<vec.size(); i++) vec[i] = !vptr->values[i].empty();
    return true;
  } else if ( const VectorValue<long>* vptr 
	      = dynamic_cast<const VectorValue<long>*> (val) )  {
    vec.resize(vptr->values.size());
    for (int i=0; i<vec.size(); i++) vec[i] = ( vptr->values[i] != 0);
    return true;
  } else if ( const VectorValue<double>* vptr 
	      = dynamic_cast<const VectorValue<double>*> (val) )  {
    vec.resize(vptr->values.size());
    for (int i=0; i<vec.size(); i++) vec[i] = ( vptr->values[i] != 0.);
    return true;
  } else {
    return false;
  }
}

template <>
bool ValueConverter(vector<string>& vec, const expressions::Value* val) {
  if ( const VectorValue<string>* vptr 
       = dynamic_cast<const VectorValue<string>*> (val) )  {
    vec = vptr->values;
    return true;
  } else {
    return false;
  }
}

template <>
bool ValueConverter(vector<long>& vec, const expressions::Value* val) {
  if ( const VectorValue<long>* vptr 
       = dynamic_cast<const VectorValue<long>*> (val) )  {
    vec = vptr->values;
    return true;
  } else {
    return false;
  }
}

template <>
bool ValueConverter(vector<double>& vec, const expressions::Value* val) {
  if ( const VectorValue<double>* vptr 
       = dynamic_cast<const VectorValue<double>*> (val) )  {
    vec = vptr->values;
    return true;
  } else if ( const VectorValue<long>* vptr 
	      = dynamic_cast<const VectorValue<long>*> (val) )  {
    vec.resize(vptr->values.size());
    for (int i=0; i<vec.size(); i++) vec[i] = vptr->values[i];
    return true;
  } else {
    return false;
  }
}

template 
void TableData::evaluate(vector<bool>&, const string&, const Header* hh) const;
