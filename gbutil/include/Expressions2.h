// More classes to run the expression evaluator
#ifndef EXPRESSIONS2_H
#define EXPRESSIONS2_H

#include <functional>

namespace expressions {

// useful inline for seeing if our token is present at front
inline bool matchesThis(const std::string& key, 
			const std::string& input,
			size_t& begin, size_t& end) {
    
  size_t keyIndex=0;
  size_t inputIndex=begin;
  while ( keyIndex < key.size()) {
    if (inputIndex >= input.size() || inputIndex >= end) 
      return false; // Not enough input characters to match
    if (input[inputIndex] != key[keyIndex]) return false;
    keyIndex++;
    inputIndex++;
  }
  // Matched all key characters if here:
  begin = inputIndex;
  return true;
}

class OpenParenthesis: public Token {
public:    
  OpenParenthesis(size_t begin=0): Token(begin) {}
  virtual Token* createFromString(const std::string& input, size_t& begin, size_t& end,
				  bool lastTokenWasOperator) const {
    if (matchesThis("(", input, begin, end)) {
      return new OpenParenthesis(begin-1);
    } else {
      return 0;
    }
  }
};

class CloseParenthesis: public Token {
public:    
  CloseParenthesis(size_t begin=0): Token(begin) {}
  virtual Token* createFromString(const std::string& input, size_t& begin, size_t& end,
				  bool lastTokenWasOperator) const {
    if (matchesThis(")", input, begin, end)) {
      return new CloseParenthesis(begin-1);
    } else {
      return 0;
    }
  }
};

class StringConstantToken: public Token {
private:
  std::string value;
public:
  StringConstantToken(const std::string& vv="", size_t begin=0): Token(begin), value(vv) {}
  virtual Token* createFromString(const std::string& input, size_t& begin, size_t& end,
				  bool lastTokenWasOperator) const {
    char delim;
    if (input[begin]=='\"')
      delim = '\"';
    else
      if (input[begin]=='\'')
	delim = '\'';
      else
	return 0;
    // Found a string, continue to matching delimiter
    std::string vv;
    size_t firstDelim = begin;
    while (true) {
      begin++;
      if (begin == end || begin >= input.size()) {
	throw SyntaxError("Unmatched string delimiter", firstDelim);
      } else if (input[begin]==delim) {
	// done.
	begin++;
	return new StringConstantToken(vv);
      } else {
	vv.push_back(input[begin]);
      }
    }
  }
  virtual Evaluable* createEvaluable() const {
    return new ConstantEvaluable<ScalarValue<std::string> >(ScalarValue<std::string>(value));
  }
};

/////////////////////////////////////////////////
// Token class that reads numbers
/////////////////////////////////////////////////
class NumberToken: public Token {
private:
  long longValue;
  double doubleValue;
  bool isDouble;
public:
  NumberToken(size_t begin=0): Token(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    if (!std::isdigit(input[begin])) return 0;
    // read the number
    bool foundDouble = false;
    long tmpLong;
    double tmpDouble;
    size_t inChar = begin;
    std::string buffer;
    while (begin != end) {
      char c=input[begin];
      if (std::isdigit(c)) {
	buffer.push_back(c);
	begin++;
      } else if (c=='.') {
	buffer.push_back(c);
	foundDouble = true;
	begin++;
      } else if (std::toupper(c)=='E' || std::toupper(c)=='D') {
	foundDouble = true;
	// Next char should be digit, +, or -
	size_t next = begin+1;
	if (next == end) 
	  throw SyntaxError("Malformed number", begin);
	char nextc = input[next];
	if ( std::isdigit(nextc) || nextc=='+' || nextc=='-') {
	  buffer.push_back(c);
	  buffer.push_back(nextc);
	  begin += 2;
	}
      } else {
	// Non-numerical character here, end of input (or force read until whitespace???)
	break;
      }
    } // character loop
    std::istringstream iss(buffer);
    if (foundDouble)
      iss >> tmpDouble;
    else 
      iss >> tmpLong;
    if (iss.fail() || !iss.eof()) 
      throw SyntaxError("Malformed number", inChar);

    NumberToken* result = new NumberToken(inChar);
    result->isDouble = foundDouble;
    if (foundDouble)
      result->doubleValue = tmpDouble;
    else
      result->longValue = tmpLong;

    return result;
  }
  virtual Evaluable* createEvaluable() const {
    if (isDouble) 
      return new ConstantEvaluable<ScalarValue<double> >(ScalarValue<double>(doubleValue));
    else
      return new ConstantEvaluable<ScalarValue<long> >(ScalarValue<long>(longValue));
  }
};

/////////////////////////////////////////////////
// Unary Tokens and Evaluables
/////////////////////////////////////////////////

class NoOpEvaluable: public UnaryOpEvaluable {
public:
  NoOpEvaluable(Evaluable* right): UnaryOpEvaluable(right) {}
  ~NoOpEvaluable() {}
  virtual Value* returnEmptyEvaluation() const {return right->returnEmptyEvaluation();}
  virtual Value* evaluate() const {return right->evaluate();}
};

// Op is derived from std::binary_function so it has the typedef of result_type
template <class Op, class Arg>
class GenericUnaryEvaluableS: public UnaryOpEvaluable {
public:
  GenericUnaryEvaluableS(Evaluable* right): UnaryOpEvaluable(right) {}
  ~GenericUnaryEvaluableS() {}
  virtual Value* returnEmptyEvaluation() const {
    return new ScalarValue<typename Op::result_type>();
  }
  virtual Value* evaluate() const {
    ScalarValue<Arg>* rval = dynamic_cast<ScalarValue<Arg>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericUnaryEvaluable::right type");
    Op f;
    Value* retval = new ScalarValue<typename Op::result_type>( f(rval->value));
    delete rval;
    return retval;
  }
};

template <class Op, class Arg>
class GenericUnaryEvaluableV: public UnaryOpEvaluable {
public:
  GenericUnaryEvaluableV(Evaluable* right): UnaryOpEvaluable(right) {}
  ~GenericUnaryEvaluableV() {}
  virtual Value* returnEmptyEvaluation() const {
    return new VectorValue<typename Op::result_type>();
  }
  virtual Value* evaluate() const {
    VectorValue<Arg>* rval = dynamic_cast<VectorValue<Arg>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericUnaryEvaluable::right type");
    std::vector<typename Op::result_type> vv(rval->values.size());
    Op f;
    for (size_t i=0; i<vv.size(); i++) {
      vv[i] = f( rval->values[i] );
    }
    delete rval;
    return new VectorValue<typename Op::result_type>(vv);
  }
};


#define USTEST(OP,Type)							\
  if ( dynamic_cast<ScalarValue<Type>*> (rVal))                         \
     return new GenericUnaryEvaluableS<OP<Type>, Type>(right); 
#define UVTEST(OP,Type)							\
  if ( dynamic_cast<VectorValue<Type>*> (rVal))                         \
      return new GenericUnaryEvaluableV<OP<Type>, Type>(right); 

class UnaryPlusToken: public UnaryOpToken {
public:
  UnaryPlusToken(size_t begin=0): UnaryOpToken(begin) {}
  ~UnaryPlusToken() {}
  virtual Evaluable* createEvaluableUnary(Evaluable* right) const {
    Value* rVal = right->returnEmptyEvaluation();
    if (dynamic_cast<ScalarValue<long>*> (rVal)
	|| dynamic_cast<ScalarValue<double>*> (rVal)
	|| dynamic_cast<VectorValue<long>*> (rVal)
	|| dynamic_cast<VectorValue<double>*> (rVal)) return new NoOpEvaluable(right);
    throwSyntaxError("Type mismatch");
    return 0; //Can't get here
  }
};

class UnaryMinusToken: public UnaryOpToken {
public:
  UnaryMinusToken(size_t begin=0): UnaryOpToken(begin) {}
  ~UnaryMinusToken() {}
  virtual Evaluable* createEvaluableUnary(Evaluable* right) const {
    Value* rVal = right->returnEmptyEvaluation();
    USTEST(std::negate, double);
    USTEST(std::negate, long);
    UVTEST(std::negate, double);
    UVTEST(std::negate, long);
    throwSyntaxError("Type mismatch");
    return 0; //Can't get here
  }
};

class NotToken: public UnaryOpToken {
public:
  NotToken(size_t begin=0): UnaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    return matchesThis("!",input, begin, end) ? new NotToken(begin-1) : 0;
  }
  virtual Evaluable* createEvaluableUnary(Evaluable* right) const {
    Value* rVal = right->returnEmptyEvaluation();
    USTEST(std::logical_not, bool);
    USTEST(std::logical_not, long);
    USTEST(std::logical_not, double);
    UVTEST(std::logical_not, bool);
    UVTEST(std::logical_not, long);
    UVTEST(std::logical_not, double);
    throwSyntaxError("Type mismatch");
    return 0; //Can't get here
  }
};

  // Redefine the macros to no longer put template type onto operator:
#undef USTEST
#undef UVTEST

#define USTEST(OP,Type)					          \
  if ( dynamic_cast<ScalarValue<Type>*> (rVal))                   \
    return new GenericUnaryEvaluableS<OP, Type>(right); 
#define UVTEST(OP,Type)						  \
  if ( dynamic_cast<VectorValue<Type>*> (rVal))                   \
    return new GenericUnaryEvaluableV<OP, Type>(right); 

// Example of a unary math function:
class SinFunction: public std::unary_function<double,double> {
public:
  double operator()(double x) const {return std::sin(x);}
};

class SinToken: public UnaryOpToken { 
public:
  SinToken(size_t begin=0): UnaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    size_t inChar = begin;
    const std::string name="sin";
    if (name == input.substr(begin,end)) {
      begin = end;
      return new SinToken(inChar);
    } else {
      return 0;
    }
  }
  virtual Evaluable* createEvaluableUnary(Evaluable* right) const {
    Value* rVal = right->returnEmptyEvaluation();
    USTEST(SinFunction, long);
    USTEST(SinFunction, double);
    UVTEST(SinFunction, long);
    UVTEST(SinFunction, double);
    throwSyntaxError("Type mismatch");
    return 0; //Can't get here
  }
};

#undef USTEST
#undef UVTEST


/////////////////////////////////////////////////
// Binary Tokens and Evaluables
/////////////////////////////////////////////////

// Op is derived from std::binary_function so it has the typedef of result_type
// Scalar Op Scalar
template <class Op, class Arg1, class Arg2>
class GenericBinaryEvaluableSS: public BinaryOpEvaluable {
public:
  GenericBinaryEvaluableSS(Evaluable* left, Evaluable* right): 
    BinaryOpEvaluable(left, right) {}
  ~GenericBinaryEvaluableSS() {}
  virtual Value* returnEmptyEvaluation() const {
    return new ScalarValue<typename Op::result_type>();
  }
  virtual Value* evaluate() const {
    ScalarValue<Arg1>* lval = dynamic_cast<ScalarValue<Arg1>*> (left->evaluate());
    if (!lval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    ScalarValue<Arg2>* rval = dynamic_cast<ScalarValue<Arg2>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    Op f;
    ScalarValue<typename Op::result_type>* retval = 
      new ScalarValue<typename Op::result_type>( f(lval->value, rval->value));
    delete lval;
    delete rval;
    return retval;
  }
};

// Op is derived from std::binary_function so it has the typedef of result_type
// Scalar Op Vector
template <class Op, class Arg1, class Arg2>
class GenericBinaryEvaluableSV: public BinaryOpEvaluable {
public:
  GenericBinaryEvaluableSV(Evaluable* left, Evaluable* right): 
    BinaryOpEvaluable(left, right) {}
  ~GenericBinaryEvaluableSV() {}
  virtual Value* returnEmptyEvaluation() const {
    return new VectorValue<typename Op::result_type>();
  }
  virtual Value* evaluate() const {
    ScalarValue<Arg1>* lval = dynamic_cast<ScalarValue<Arg1>*> (left->evaluate());
    if (!lval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    VectorValue<Arg2>* rval = dynamic_cast<VectorValue<Arg2>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    std::vector<typename Op::result_type> vv(rval->values.size());
    Op f;
    for (size_t i=0; i<vv.size(); i++) {
      vv[i] = f(lval->value, rval->values[i] );
    }
    delete lval;
    delete rval;
    return new VectorValue<typename Op::result_type>(vv);
  }
};

// Op is derived from std::binary_function so it has the typedef of result_type
// Vector Op Scalar
template <class Op, class Arg1, class Arg2>
class GenericBinaryEvaluableVS: public BinaryOpEvaluable {
public:
  GenericBinaryEvaluableVS(Evaluable* left, Evaluable* right): 
    BinaryOpEvaluable(left, right) {}
  ~GenericBinaryEvaluableVS() {}
  virtual Value* returnEmptyEvaluation() const {
    return new VectorValue<typename Op::result_type>();
  }
  virtual Value* evaluate() const {
    VectorValue<Arg1>* lval = dynamic_cast<VectorValue<Arg1>*> (left->evaluate());
    if (!lval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    ScalarValue<Arg2>* rval = dynamic_cast<ScalarValue<Arg2>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    std::vector<typename Op::result_type> vv(lval->values.size());
    Op f;
    for (size_t i=0; i<vv.size(); i++) {
      vv[i] = f(lval->values[i], rval->value);
    }
    delete lval;
    delete rval;
    return new VectorValue<typename Op::result_type>(vv);
  }
};

// Vector Op Vector
template <class Op, class Arg1, class Arg2>
class GenericBinaryEvaluableVV: public BinaryOpEvaluable {
public:
  GenericBinaryEvaluableVV(Evaluable* left, Evaluable* right): 
    BinaryOpEvaluable(left, right) {}
  ~GenericBinaryEvaluableVV() {}
  virtual Value* returnEmptyEvaluation() const {
    return new VectorValue<typename Op::result_type>();
  }
  virtual Value* evaluate() const {
    VectorValue<Arg1>* lval = dynamic_cast<VectorValue<Arg1>*> (left->evaluate());
    if (!lval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    VectorValue<Arg2>* rval = dynamic_cast<VectorValue<Arg2>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    std::vector<typename Op::result_type> vv(lval->values.size());
    Assert(lval->values.size() == rval->values.size());
    Op f;
    for (size_t i=0; i<vv.size(); i++) {
      vv[i] = f(lval->values[i], rval->values[i]);
    }
    delete lval;
    delete rval;
    return new VectorValue<typename Op::result_type>(vv);
  }
};

#define BSSTEST(OP,Type,Type1,Type2)					\
  if ( dynamic_cast<ScalarValue<Type1>*> (lVal)                         \
      && dynamic_cast<ScalarValue<Type2>*> (rVal))                      \
    return new GenericBinaryEvaluableSS<OP<Type>, Type1,Type2>(left,right);
#define BSVTEST(OP,Type,Type1,Type2)					\
  if ( dynamic_cast<ScalarValue<Type1>*> (lVal)                         \
      && dynamic_cast<VectorValue<Type2>*> (rVal))                      \
    return new GenericBinaryEvaluableSV<OP<Type>, Type1,Type2>(left,right);
#define BVSTEST(OP,Type,Type1,Type2)					\
  if ( dynamic_cast<VectorValue<Type1>*> (lVal)                         \
      && dynamic_cast<ScalarValue<Type2>*> (rVal))                      \
    return new GenericBinaryEvaluableVS<OP<Type>, Type1,Type2>(left,right);
#define BVVTEST(OP,Type,Type1,Type2)					\
  if ( dynamic_cast<VectorValue<Type1>*> (lVal)                         \
      && dynamic_cast<VectorValue<Type2>*> (rVal))                      \
    return new GenericBinaryEvaluableVV<OP<Type>, Type1,Type2>(left,right);


class BinaryPlusToken: public BinaryOpToken {
public:
  BinaryPlusToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(std::plus, std::string, std::string, std::string);
    BSSTEST(std::plus, long, long, long);
    BSSTEST(std::plus, double, long, double);
    BSSTEST(std::plus, double, double, long);
    BSSTEST(std::plus, double, double, double);
    BSVTEST(std::plus, long, long, long);
    BSVTEST(std::plus, double, long, double);
    BSVTEST(std::plus, double, double, long);
    BSVTEST(std::plus, double, double, double);
    BVSTEST(std::plus, long, long, long);
    BVSTEST(std::plus, double, long, double);
    BVSTEST(std::plus, double, double, long);
    BVSTEST(std::plus, double, double, double);
    BVVTEST(std::plus, long, long, long);
    BVVTEST(std::plus, double, long, double);
    BVVTEST(std::plus, double, double, long);
    BVVTEST(std::plus, double, double, double);
    throwSyntaxError("Type mismatch");
    return 0; //Cannot get here;
  }
};

class BinaryMinusToken: public BinaryOpToken {
public:
  BinaryMinusToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(std::minus, long, long, long);
    BSSTEST(std::minus, double, long, double);
    BSSTEST(std::minus, double, double, long);
    BSSTEST(std::minus, double, double, double);
    BSVTEST(std::minus, long, long, long);
    BSVTEST(std::minus, double, long, double);
    BSVTEST(std::minus, double, double, long);
    BSVTEST(std::minus, double, double, double);
    BVSTEST(std::minus, long, long, long);
    BVSTEST(std::minus, double, long, double);
    BVSTEST(std::minus, double, double, long);
    BVSTEST(std::minus, double, double, double);
    BVVTEST(std::minus, long, long, long);
    BVVTEST(std::minus, double, long, double);
    BVVTEST(std::minus, double, double, long);
    BVVTEST(std::minus, double, double, double);
    throwSyntaxError("Type mismatch");
    return 0; // Can't get here
  }
};

/*
class ModulusToken: public BinaryOpToken {
public:
  ModulusToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    return matchesThis("%",input, begin, end) ? new ModulusToken(begin-1) : 0;
  }
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(std::modulus, long, long, long);
    BSVTEST(std::modulus, long, long, long);
    BVSTEST(std::modulus, long, long, long);
    BVVTEST(std::modulus, long, long, long);
    throwSyntaxError("Type mismatch");
    return 0; // Can't get here
  }
};
*/
class MultipliesToken: public BinaryOpToken {
public:
  MultipliesToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    return matchesThis("*",input, begin, end) ? new MultipliesToken(begin-1) : 0;
  }
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(std::multiplies, long, long, long);
    BSSTEST(std::multiplies, double, long, double);
    BSSTEST(std::multiplies, double, double, long);
    BSSTEST(std::multiplies, double, double, double);
    BSVTEST(std::multiplies, long, long, long);
    BSVTEST(std::multiplies, double, long, double);
    BSVTEST(std::multiplies, double, double, long);
    BSVTEST(std::multiplies, double, double, double);
    BVSTEST(std::multiplies, long, long, long);
    BVSTEST(std::multiplies, double, long, double);
    BVSTEST(std::multiplies, double, double, long);
    BVSTEST(std::multiplies, double, double, double);
    BVVTEST(std::multiplies, long, long, long);
    BVVTEST(std::multiplies, double, long, double);
    BVVTEST(std::multiplies, double, double, long);
    BVVTEST(std::multiplies, double, double, double);
    throwSyntaxError("Type mismatch");
    return 0; // Can't get here
  }
};

class DividesToken: public BinaryOpToken {
public:
  DividesToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    return matchesThis("/",input, begin, end) ? new DividesToken(begin-1) : 0;
  }
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(std::divides, long, long, long);
    BSSTEST(std::divides, double, long, double);
    BSSTEST(std::divides, double, double, long);
    BSSTEST(std::divides, double, double, double);
    BSVTEST(std::divides, long, long, long);
    BSVTEST(std::divides, double, long, double);
    BSVTEST(std::divides, double, double, long);
    BSVTEST(std::divides, double, double, double);
    BVSTEST(std::divides, long, long, long);
    BVSTEST(std::divides, double, long, double);
    BVSTEST(std::divides, double, double, long);
    BVSTEST(std::divides, double, double, double);
    BVVTEST(std::divides, long, long, long);
    BVVTEST(std::divides, double, long, double);
    BVVTEST(std::divides, double, double, long);
    BVVTEST(std::divides, double, double, double);
    throwSyntaxError("Type mismatch");
    return 0; // Can't get here
  }
};

#define BINARY_LOGICAL(NAME,FUNC,CODE)                                         \
class NAME: public BinaryOpToken {					       \
public:									       \
  NAME(size_t begin=0): BinaryOpToken(begin) {}				       \
  virtual Token* createFromString(const std::string& input,		       \
				  size_t& begin,			       \
				  size_t& end,				       \
				  bool lastTokenWasOperator) const {	       \
    long inChar = begin;						       \
    return matchesThis(CODE,input, begin, end) ? new NAME(inChar) : 0;	       \
  }									       \
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const { \
    Value* lVal = left->returnEmptyEvaluation();			       \
    Value* rVal = right->returnEmptyEvaluation();			       \
    BSSTEST(FUNC, bool, bool, bool)					       \
    BSSTEST(FUNC, bool, bool, long)					       \
    BSSTEST(FUNC, bool, bool, double)					       \
    BSSTEST(FUNC, bool, long, bool)					       \
    BSSTEST(FUNC, bool, long, long)					       \
    BSSTEST(FUNC, bool, long, double)					       \
    BSSTEST(FUNC, bool, double, bool)					       \
    BSSTEST(FUNC, bool, double, long)					       \
    BSSTEST(FUNC, bool, double, double)					       \
									       \
    BSVTEST(FUNC, bool, bool, bool)					       \
    BSVTEST(FUNC, bool, bool, long)					       \
    BSVTEST(FUNC, bool, bool, double)					       \
    BSVTEST(FUNC, bool, long, bool)					       \
    BSVTEST(FUNC, bool, long, long)					       \
    BSVTEST(FUNC, bool, long, double)					       \
    BSVTEST(FUNC, bool, double, bool)					       \
    BSVTEST(FUNC, bool, double, long)					       \
    BSVTEST(FUNC, bool, double, double)					       \
									       \
    BVSTEST(FUNC, bool, bool, bool)					       \
    BVSTEST(FUNC, bool, bool, long)					       \
    BVSTEST(FUNC, bool, bool, double)					       \
    BVSTEST(FUNC, bool, long, bool)					       \
    BVSTEST(FUNC, bool, long, long)					       \
    BVSTEST(FUNC, bool, long, double)					       \
    BVSTEST(FUNC, bool, double, bool)					       \
    BVSTEST(FUNC, bool, double, long)					       \
    BVSTEST(FUNC, bool, double, double)					       \
									       \
    BVVTEST(FUNC, bool, bool, bool)					       \
    BVVTEST(FUNC, bool, bool, long)					       \
    BVVTEST(FUNC, bool, bool, double)					       \
    BVVTEST(FUNC, bool, long, bool)					       \
    BVVTEST(FUNC, bool, long, long)					       \
    BVVTEST(FUNC, bool, long, double)					       \
    BVVTEST(FUNC, bool, double, bool)					       \
    BVVTEST(FUNC, bool, double, long)					       \
    BVVTEST(FUNC, bool, double, double)					       \
    throwSyntaxError("Type mismatch");					       \
    return 0;								       \
  }									       \
}


  BINARY_LOGICAL(AndToken, std::logical_and, "&&");
  BINARY_LOGICAL(OrToken, std::logical_or, "||");

#undef BINARY_LOGICAL

#define BINARY_INTEGER(NAME,FUNC,CODE)                                         \
class NAME: public BinaryOpToken {					       \
public:									       \
  NAME(size_t begin=0): BinaryOpToken(begin) {}				       \
  virtual Token* createFromString(const std::string& input,		       \
				  size_t& begin,			       \
				  size_t& end,				       \
				  bool lastTokenWasOperator) const {	       \
    long inChar = begin;						       \
    return matchesThis(CODE,input, begin, end) ? new NAME(inChar) : 0;         \
  }									       \
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const { \
    Value* lVal = left->returnEmptyEvaluation();			       \
    Value* rVal = right->returnEmptyEvaluation();			       \
    BSSTEST(FUNC, long, long, long)					\
    BVSTEST(FUNC, long, long, long)					\
    BSVTEST(FUNC, long, long, long)					\
    BVVTEST(FUNC, long, long, long)					\
    throwSyntaxError("Type mismatch");					       \
    return 0;								       \
  }									       \
}

  using std::bit_and; 
  using std::bit_or;

  // Define binary operators which only work with integer arguments
  BINARY_INTEGER(ModulusToken, std::modulus, "%");
  BINARY_INTEGER(BitAndToken, bit_and, "&");
  BINARY_INTEGER(BitOrToken, bit_or, "|");

#undef BINARY_INTEGER
  
#define COMPARISON(NAME,FUNC,CODE)                                          \
class NAME: public BinaryOpToken {					    \
public:									    \
  NAME(size_t begin=0): BinaryOpToken(begin) {}				    \
  virtual Token* createFromString(const std::string& input,		    \
				  size_t& begin,			    \
				  size_t& end,				    \
				  bool lastTokenWasOperator) const {	    \
    long inChar = begin;						    \
    return ( matchesThis(CODE,input, begin, end)) ? new NAME(inChar) : 0;   \
  }     								    \
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const { \
    Value* lVal = left->returnEmptyEvaluation();			    \
    Value* rVal = right->returnEmptyEvaluation();			    \
    BSSTEST(FUNC, std::string, std::string, std::string);                   \
    BSSTEST(FUNC, long, long, long);					    \
    BSSTEST(FUNC, double, long, double);				    \
    BSSTEST(FUNC, double, double, long);				    \
    BSSTEST(FUNC, double, double, double);				    \
									    \
    BSVTEST(FUNC, std::string, std::string, std::string);                   \
    BSVTEST(FUNC, long, long, long);					    \
    BSVTEST(FUNC, double, long, double);				    \
    BSVTEST(FUNC, double, double, long);				    \
    BSVTEST(FUNC, double, double, double);				    \
									    \
    BVSTEST(FUNC, std::string, std::string, std::string);                   \
    BVSTEST(FUNC, long, long, long);					    \
    BVSTEST(FUNC, double, long, double);				    \
    BVSTEST(FUNC, double, double, long);				    \
    BVSTEST(FUNC, double, double, double);				    \
									    \
    BVVTEST(FUNC, std::string, std::string, std::string);                   \
    BVVTEST(FUNC, long, long, long);					    \
    BVVTEST(FUNC, double, long, double);				    \
    BVVTEST(FUNC, double, double, long);				    \
    BVVTEST(FUNC, double, double, double);				    \
									    \
    throwSyntaxError("Type mismatch");   				    \
    return 0;                 				                    \
  }									    \
}

COMPARISON(GreaterToken,std::greater,">");
COMPARISON(LessToken,std::less,"<");
COMPARISON(EqualToken,std::equal_to,"==");
COMPARISON(GreaterEqualToken,std::greater_equal,">=");
COMPARISON(LessEqualToken,std::less_equal,"<=");
COMPARISON(NotEqualToken,std::not_equal_to,"!=");

#undef COMPARISON

template <class T>
class PowerOf {
public:
  typedef T result_type;
  T operator()(const T& v1, const T& v2) const {return std::pow(v1,v2);}
};

class PowerToken: public BinaryOpToken {
public:
  PowerToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    size_t inChar = begin;
    return (matchesThis("**",input, begin, end)
	    || matchesThis("^",input, begin, end)) ? new PowerToken(inChar) : 0;
  }
  virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(PowerOf, double, long, long);
    BSSTEST(PowerOf, double, long, double);
    BSSTEST(PowerOf, double, double, long);
    BSSTEST(PowerOf, double, double, double);
    BSVTEST(PowerOf, double, long, long);
    BSVTEST(PowerOf, double, long, double);
    BSVTEST(PowerOf, double, double, long);
    BSVTEST(PowerOf, double, double, double);
    BVSTEST(PowerOf, double, long, long);
    BVSTEST(PowerOf, double, long, double);
    BVSTEST(PowerOf, double, double, long);
    BVSTEST(PowerOf, double, double, double);
    BVVTEST(PowerOf, double, long, long);
    BVVTEST(PowerOf, double, long, double);
    BVVTEST(PowerOf, double, double, long);
    BVVTEST(PowerOf, double, double, double);
    throwSyntaxError("Type mismatch");
    return 0;  // Can't get here
  }
};

#undef BSSTEST
#undef BSVTEST
#undef BVSTEST
#undef BVVTEST

////////////////////////////////
// Tokens that really just decide whether unary or binary +/- are needed
////////////////////////////////


class PlusToken: public Token {
private:
  bool isBinary;
public:
  PlusToken()  {}
  virtual ~PlusToken() {}
  virtual bool isOperator() const {return true;}
  virtual Token* createFromString(const std::string& input, size_t& begin, size_t& end,
				  bool lastTokenWasOperator) const {
    if (matchesThis("+", input, begin, end)) {
      size_t opIndex = begin-1;
      if (lastTokenWasOperator)
	return new UnaryPlusToken(opIndex);
      else
	return new BinaryPlusToken(opIndex);
    } else {
      return 0;
    }
  }
};

class MinusToken: public Token {
private:
  bool isBinary;
public:
  MinusToken()  {}
  virtual ~MinusToken() {}
  virtual bool isOperator() const {return true;}
  virtual Token* createFromString(const std::string& input, size_t& begin, size_t& end,
				  bool lastTokenWasOperator) const {
    if (matchesThis("-", input, begin, end)) {
      size_t opIndex = begin-1;
      if (lastTokenWasOperator)
	return new UnaryMinusToken(opIndex);
      else
	return new BinaryMinusToken(opIndex);
    } else {
      return 0;
    }
  }
};

} // end namespace expressions
#endif // EXPRESSIONS2_H
