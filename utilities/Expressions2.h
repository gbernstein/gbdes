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
      return new OpenParenthesis(begin);
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
      return new CloseParenthesis(begin);
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
      std::istringstream iss(buffer);
      if (foundDouble)
	iss >> tmpDouble;
      else 
	iss >> tmpLong;
      if (iss.fail() || !iss.eof()) 
	throw SyntaxError("Malformed number", inChar);
    } // character loop

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

// Op is derived from std::binary_function so it has the typedef of Result
template <class Op, class Arg>
class GenericUnaryEvaluableS: public UnaryOpEvaluable {
public:
  GenericUnaryEvaluableS(Evaluable* right): UnaryOpEvaluable(right) {}
  ~GenericUnaryEvaluableS() {}
  virtual Value* returnEmptyEvaluation() const {
    return ScalarValue<typename Op::Result>();
  }
  virtual Value* evaluate() const {
    ScalarValue<Arg>* rval = dynamic_cast<ScalarValue<Arg>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericUnaryEvaluable::right type");
    Op f;
    Value* retval = new ScalarValue<typename Op::Result>( f(rval->value));
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
    return new VectorValue<typename Op::Result>();
  }
  virtual Value* evaluate() const {
    VectorValue<Arg>* rval = dynamic_cast<VectorValue<Arg>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericUnaryEvaluable::right type");
    std::vector<typename Op::Result> vv(rval->values.size());
    Op f;
    for (size_t i=0; i<vv.size(); i++) {
      vv[i] = f( rval->values[i] );
    }
    delete rval;
    return new VectorValue<typename Op::Result>(vv);
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
  virtual Evaluable* createEvaluable(Evaluable* right) const {
    Value* rVal = right->returnEmptyEvaluation();
    if (dynamic_cast<ScalarValue<long>*> (rVal)
	|| dynamic_cast<ScalarValue<double>*> (rVal)
	|| dynamic_cast<VectorValue<long>*> (rVal)
	|| dynamic_cast<VectorValue<double>*> (rVal)) return new NoOpEvaluable(right);
    throwSyntaxError("Type mismatch");
  }
};

class UnaryMinusToken: public UnaryOpToken {
public:
  UnaryMinusToken(size_t begin=0): UnaryOpToken(begin) {}
  ~UnaryMinusToken() {}
  virtual Evaluable* createEvaluable(Evaluable* right) const {
    Value* rVal = right->returnEmptyEvaluation();
    USTEST(std::negate, double);
    USTEST(std::negate, long);
    UVTEST(std::negate, double);
    UVTEST(std::negate, long);
    throwSyntaxError("Type mismatch");
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
  virtual Evaluable* createEvaluable(Evaluable* right) const {
    Value* rVal = right->returnEmptyEvaluation();
    USTEST(std::logical_not, bool);
    USTEST(std::logical_not, long);
    USTEST(std::logical_not, double);
    UVTEST(std::logical_not, bool);
    UVTEST(std::logical_not, long);
    UVTEST(std::logical_not, double);
    throwSyntaxError("Type mismatch");
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
class SinFunction: std::unary_function<double,double> {
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
  virtual Evaluable* createEvaluable(Evaluable* right) const {
    Value* rVal = right->returnEmptyEvaluation();
    USTEST(SinFunction, long);
    USTEST(SinFunction, double);
    UVTEST(SinFunction, long);
    UVTEST(SinFunction, double);
    throwSyntaxError("Type mismatch");
  }
};

#undef USTEST
#undef UVTEST


/////////////////////////////////////////////////
// Binary Tokens and Evaluables
/////////////////////////////////////////////////

// Op is derived from std::binary_function so it has the typedef of Result
// Scalar Op Scalar
template <class Op, class Arg1, class Arg2>
class GenericBinaryEvaluableSS: public BinaryOpEvaluable {
public:
  GenericBinaryEvaluableSS(Evaluable* left, Evaluable* right): 
    BinaryOpEvaluable(left, right) {}
  ~GenericBinaryEvaluableSS() {}
  virtual Value* returnEmptyEvaluation() const {
    return new ScalarValue<typename Op::Result>();
  }
  virtual Value* evaluate() const {
    ScalarValue<Arg1>* lval = dynamic_cast<ScalarValue<Arg1>*> (left->evaluate());
    if (!lval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    ScalarValue<Arg2>* rval = dynamic_cast<ScalarValue<Arg2>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    Op f;
    ScalarValue<typename Op::Result>* retval = new ScalarValue<typename Op::Result>( f(lval->value, rval->value));
    delete lval;
    delete rval;
    return retval;
  }
};

// Op is derived from std::binary_function so it has the typedef of Result
// Scalar Op Vector
template <class Op, class Arg1, class Arg2>
class GenericBinaryEvaluableSV: public BinaryOpEvaluable {
public:
  GenericBinaryEvaluableSV(Evaluable* left, Evaluable* right): 
    BinaryOpEvaluable(left, right) {}
  ~GenericBinaryEvaluableSV() {}
  virtual Value* returnEmptyEvaluation() const {
    return new VectorValue<typename Op::Result>();
  }
  virtual Value* evaluate() const {
    ScalarValue<Arg1>* lval = dynamic_cast<ScalarValue<Arg1>*> (left->evaluate());
    if (!lval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    VectorValue<Arg2>* rval = dynamic_cast<VectorValue<Arg2>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    std::vector<typename Op::Result> vv(rval->values.size());
    Op f;
    for (size_t i=0; i<vv.size(); i++) {
      vv[i] = f(lval->value, rval->values[i] );
    }
    delete lval;
    delete rval;
    return new VectorValue<typename Op::Result>(vv);
  }
};

// Op is derived from std::binary_function so it has the typedef of Result
// Vector Op Scalar
template <class Op, class Arg1, class Arg2>
class GenericBinaryEvaluableVS: public BinaryOpEvaluable {
public:
  GenericBinaryEvaluableVS(Evaluable* left, Evaluable* right): 
    BinaryOpEvaluable(left, right) {}
  ~GenericBinaryEvaluableVS() {}
  virtual Value* returnEmptyEvaluation() const {
    return new VectorValue<typename Op::Result>();
  }
  virtual Value* evaluate() const {
    VectorValue<Arg1>* lval = dynamic_cast<VectorValue<Arg1>*> (left->evaluate());
    if (!lval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    ScalarValue<Arg2>* rval = dynamic_cast<ScalarValue<Arg2>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    std::vector<typename Op::Result> vv(lval->values.size());
    Op f;
    for (size_t i=0; i<vv.size(); i++) {
      vv[i] = f(lval->values[i], rval->value);
    }
    delete lval;
    delete rval;
    return new VectorValue<typename Op::Result>(vv);
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
    return new VectorValue<typename Op::Result>();
  }
  virtual Value* evaluate() const {
    VectorValue<Arg1>* lval = dynamic_cast<VectorValue<Arg1>*> (left->evaluate());
    if (!lval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    VectorValue<Arg2>* rval = dynamic_cast<VectorValue<Arg2>*> (right->evaluate());
    if (!rval) throw ExpressionError("Bad GenericBinaryEvaluable::left type");
    std::vector<typename Op::Result> vv(lval->values.size());
    Assert(lval->values.size() == rval->values.size());
    Op f;
    for (size_t i=0; i<vv.size(); i++) {
      vv[i] = f(lval->values[i], rval->values[i]);
    }
    delete lval;
    delete rval;
    return new VectorValue<typename Op::Result>(vv);
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
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {
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
  }
};

class BinaryMinusToken: public BinaryOpToken {
public:
  BinaryMinusToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {
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
  }
};

class ModulusToken: public BinaryOpToken {
public:
  ModulusToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    return matchesThis("%",input, begin, end) ? new ModulusToken(begin-1) : 0;
  }
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(std::modulus, long, long, long);
    BSVTEST(std::modulus, long, long, long);
    BVSTEST(std::modulus, long, long, long);
    BVVTEST(std::modulus, long, long, long);
    throwSyntaxError("Type mismatch");
  }
};

class MultipliesToken: public BinaryOpToken {
public:
  MultipliesToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    return matchesThis("*",input, begin, end) ? new MultipliesToken(begin-1) : 0;
  }
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {
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
  }
};

class DividesToken: public BinaryOpToken {
public:
  DividesToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    return matchesThis("*",input, begin, end) ? new DividesToken(begin-1) : 0;
  }
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {
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
  }
};

class AndToken: public BinaryOpToken {
public:
  AndToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    long inChar = begin;
    return ( matchesThis("&&",input, begin, end)
	     || matchesThis("&",input, begin, end)) ? new AndToken(inChar) : 0;
  }
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(std::logical_and, bool, bool, bool)
    BSSTEST(std::logical_and, bool, bool, long)
    BSSTEST(std::logical_and, bool, bool, double)
    BSSTEST(std::logical_and, bool, long, bool)
    BSSTEST(std::logical_and, bool, long, long)
    BSSTEST(std::logical_and, bool, long, double)
    BSSTEST(std::logical_and, bool, double, bool)
    BSSTEST(std::logical_and, bool, double, long)
    BSSTEST(std::logical_and, bool, double, double)

    BSVTEST(std::logical_and, bool, bool, bool)
    BSVTEST(std::logical_and, bool, bool, long)
    BSVTEST(std::logical_and, bool, bool, double)
    BSVTEST(std::logical_and, bool, long, bool)
    BSVTEST(std::logical_and, bool, long, long)
    BSVTEST(std::logical_and, bool, long, double)
    BSVTEST(std::logical_and, bool, double, bool)
    BSVTEST(std::logical_and, bool, double, long)
    BSVTEST(std::logical_and, bool, double, double)

    BVSTEST(std::logical_and, bool, bool, bool)
    BVSTEST(std::logical_and, bool, bool, long)
    BVSTEST(std::logical_and, bool, bool, double)
    BVSTEST(std::logical_and, bool, long, bool)
    BVSTEST(std::logical_and, bool, long, long)
    BVSTEST(std::logical_and, bool, long, double)
    BVSTEST(std::logical_and, bool, double, bool)
    BVSTEST(std::logical_and, bool, double, long)
    BVSTEST(std::logical_and, bool, double, double)

    BVVTEST(std::logical_and, bool, bool, bool)
    BVVTEST(std::logical_and, bool, bool, long)
    BVVTEST(std::logical_and, bool, bool, double)
    BVVTEST(std::logical_and, bool, long, bool)
    BVVTEST(std::logical_and, bool, long, long)
    BVVTEST(std::logical_and, bool, long, double)
    BVVTEST(std::logical_and, bool, double, bool)
    BVVTEST(std::logical_and, bool, double, long)
    BVVTEST(std::logical_and, bool, double, double)
    throwSyntaxError("Type mismatch");
  }
};

class OrToken: public BinaryOpToken {
public:
  OrToken(size_t begin=0): BinaryOpToken(begin) {}
  virtual Token* createFromString(const std::string& input,
				  size_t& begin,
				  size_t& end,
				  bool lastTokenWasOperator) const {
    long inChar = begin;
    return ( matchesThis("||",input, begin, end)
	     || matchesThis("|",input, begin, end)) ? new OrToken(inChar) : 0;
  }
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {
    Value* lVal = left->returnEmptyEvaluation();
    Value* rVal = right->returnEmptyEvaluation();
    BSSTEST(std::logical_or, bool, bool, bool)
    BSSTEST(std::logical_or, bool, bool, long)
    BSSTEST(std::logical_or, bool, bool, double)
    BSSTEST(std::logical_or, bool, long, bool)
    BSSTEST(std::logical_or, bool, long, long)
    BSSTEST(std::logical_or, bool, long, double)
    BSSTEST(std::logical_or, bool, double, bool)
    BSSTEST(std::logical_or, bool, double, long)
    BSSTEST(std::logical_or, bool, double, double)

    BSVTEST(std::logical_or, bool, bool, bool)
    BSVTEST(std::logical_or, bool, bool, long)
    BSVTEST(std::logical_or, bool, bool, double)
    BSVTEST(std::logical_or, bool, long, bool)
    BSVTEST(std::logical_or, bool, long, long)
    BSVTEST(std::logical_or, bool, long, double)
    BSVTEST(std::logical_or, bool, double, bool)
    BSVTEST(std::logical_or, bool, double, long)
    BSVTEST(std::logical_or, bool, double, double)

    BVSTEST(std::logical_or, bool, bool, bool)
    BVSTEST(std::logical_or, bool, bool, long)
    BVSTEST(std::logical_or, bool, bool, double)
    BVSTEST(std::logical_or, bool, long, bool)
    BVSTEST(std::logical_or, bool, long, long)
    BVSTEST(std::logical_or, bool, long, double)
    BVSTEST(std::logical_or, bool, double, bool)
    BVSTEST(std::logical_or, bool, double, long)
    BVSTEST(std::logical_or, bool, double, double)

    BVVTEST(std::logical_or, bool, bool, bool)
    BVVTEST(std::logical_or, bool, bool, long)
    BVVTEST(std::logical_or, bool, bool, double)
    BVVTEST(std::logical_or, bool, long, bool)
    BVVTEST(std::logical_or, bool, long, long)
    BVVTEST(std::logical_or, bool, long, double)
    BVVTEST(std::logical_or, bool, double, bool)
    BVVTEST(std::logical_or, bool, double, long)
    BVVTEST(std::logical_or, bool, double, double)
    throwSyntaxError("Type mismatch");
  }
};

#define COMPARISON(NAME,FUNC,CODE)                                          \
class NAME: public BinaryOpToken {					    \
public:									    \
  NAME(size_t begin=0): BinaryOpToken(begin) {}				    \
  virtual Token* createFromString(const std::string& input,		    \
				  size_t& begin,			    \
				  size_t& end,				    \
				  bool lastTokenWasOperator) const {	    \
    long inChar = begin;						    \
    return ( matchesThis(CODE,input, begin, end)) ? new NAME(inChar) : 0;	    \
  }									    \
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {   \
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
    throwSyntaxError("Type mismatch");				    \
  }									    \
}


  ///  Some kludges needed because the comparison operators do not appear to 
  /// meet the binary_function standard of having Result typecast
  template <class T>
  struct _NE: public std::not_equal_to<T> {
    typedef bool Result;
  };
  template <class T>
  struct _GT: public std::greater<T> {
    typedef bool Result;
  };
  template <class T>
  struct _LT: public std::less<T> {
    typedef bool Result;
  };
  template <class T>
  struct _EQ: public std::equal_to<T> {
    typedef bool Result;
  };
  template <class T>
  struct _GE: public std::greater_equal<T> {
    typedef bool Result;
  };
  template <class T>
  struct _LE: public std::less_equal<T> {
    typedef bool Result;
  };

  /**
COMPARISON(GreaterToken,std::greater,">");
COMPARISON(LessToken,std::less,"<");
COMPARISON(EqualToken,std::equal_to,"==");
COMPARISON(GreaterEqualToken,std::greater_equal,">=");
COMPARISON(LessEqualToken,std::less_equal,"<=");
COMPARISON(NotEqualToken,std::not_equal_to,"!=");
  **/
COMPARISON(GreaterToken,_GT,">");
COMPARISON(LessToken,_LT,"<");
COMPARISON(EqualToken,_EQ,"==");
COMPARISON(GreaterEqualToken,_GE,">=");
COMPARISON(LessEqualToken,_LE,"<=");
COMPARISON(NotEqualToken,_NE,"!=");


#undef COMPARISON

template <class T>
class PowerOf {
public:
  typedef T Result;
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
	    || matchesThis("^",input, begin, end)) ? new DividesToken(inChar) : 0;
  }
  virtual Evaluable* createEvaluable(Evaluable* left, Evaluable* right) {
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
      size_t opIndex = begin++;
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
      size_t opIndex = begin++;
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
