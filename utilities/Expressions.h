// Classes to parse and evaluate string expressions.
// Can operate on either scalar or vector data types.

#ifndef EXPRESSIONS_H
#define EXPRESSIONS_H

#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <list>
#include <ios>
#include "Std.h"
#include <typeinfo>

namespace expressions {

  // Exception classes
  class ExpressionError: public std::runtime_error {
  public:
    ExpressionError(const std::string& m): std::runtime_error("Expression error: " +m) {}
    virtual ~ExpressionError() throw() {}
  };

  // Exception class for syntax error, includes index of error location
  class SyntaxError: public ExpressionError {
  public:
    SyntaxError(const std::string& m, size_t charBegin_): ExpressionError(m),
							  charBegin(charBegin_),
							  why(m) {}
    
    virtual ~SyntaxError() throw() {}
    size_t charBegin;
    std::string why;
    std::string longWhat(const std::string& fullInput) {
      std::ostringstream oss;
      oss << why << " at character " << charBegin;
      oss << "\n" << fullInput << "\n";
      for (size_t i=0; i<charBegin; i++) oss << " ";
      oss << "^";
      return oss.str();
    }
  };
      
  //////  First are classes that be the result of evaluation of expressions:
  class Value {
  public:
    virtual ~Value() {};
    virtual Value* duplicate() const {return new Value;}
  };

  // Derived class for null return:
  class VoidValue: public Value {
  public:
    virtual ~VoidValue() {};
    virtual Value* duplicate() const {return new VoidValue;}
  };

  template <class T>
  class ScalarValue: public Value {
  public:
    ScalarValue(const T& v=T()): value(v) {}
    virtual ~ScalarValue() {};
    virtual Value* duplicate() const {return new ScalarValue(value);}
    T value;
    typedef T ValueType;
  };

  template <class T>
  class VectorValue: public Value {
  public:
    VectorValue(const std::vector<T>& v = std::vector<T>()): values(v) {}
    virtual ~VectorValue() {};
    virtual Value* duplicate() const {return new VectorValue(values);}
    std::vector<T> values;
    typedef T ValueType;
  };

  ////// Now class that can be evaluated to give a result:
  class Evaluable {
  public:
    virtual ~Evaluable() {};
    virtual Value* evaluate() const {return new VoidValue;}
    // Following function used to get return type of evaluation
    virtual Value* returnEmptyEvaluation() const {return new VoidValue;}
  };

  // V shoud be a Value class
  template <class V>
  class ConstantEvaluable: public Evaluable {
  private:
    V value;
  public:
    ConstantEvaluable(const V& v_): value(v_) {}
    virtual ~ConstantEvaluable() {}
    virtual Value* evaluate() const {return value.duplicate();}
    virtual Value* returnEmptyEvaluation() const {return new V();}
  };

  class UnaryOpEvaluable: public Evaluable {
  protected:
    Evaluable* right;
  public:
    UnaryOpEvaluable(Evaluable* right_): right(right_) {}
    virtual ~UnaryOpEvaluable() {delete right;}
  };

  class BinaryOpEvaluable: public Evaluable {
  protected:
    Evaluable* left;
    Evaluable* right;
  public:
    BinaryOpEvaluable(Evaluable* left_, Evaluable* right_): 
      left(left_), right(right_) {}
    virtual ~BinaryOpEvaluable() {delete left; delete right;}
  };
    

  ///// Now base class for a token that can be found in expression string
  class Token {
  protected:
    size_t beginChar;
  public:
    Token(size_t begin=0): beginChar(begin) {}
    virtual ~Token() {};
    // Function to create a token of this type from the front of specified string.
    // Updates begin index as characters are used.
    // Return null pointer if string does not match token.
    virtual Token* createFromString(const std::string& input,
				    size_t& begin, size_t& end,
				    bool lastTokenWasOperator=false) const {return 0;}
    // Return true if this token will parse into an operator
    virtual bool isOperator() const {return false;}
    virtual bool isUnaryOperator() const {return false;}
    virtual bool isBinaryOperator() const {return false;}
    
    virtual Evaluable* createEvaluable() const {
      throw ExpressionError("called unimplemented Token::createEvaluable() ");
    }
     
    void throwSyntaxError(const std::string& m) const {
      throw SyntaxError(m, beginChar);
    }
  };

  // First derived class is something that is directly Evaluable.
  // Should always be a value, not an operator.
  class EvaluableToken: public Token {
  private:
    Evaluable* contents;
  public:
    EvaluableToken(Evaluable* contents_): Token(0), contents(contents_) {}
    virtual ~EvaluableToken() {} // Don't delete contents, it will be passed on
    virtual Evaluable* createEvaluable() const {return contents;}
  };

  // Derived Token that represents a unary operator
  class UnaryOpToken: public  Token {
  public:
    UnaryOpToken(size_t begin=0): Token(begin) {}
    virtual bool isOperator() const {return true;}
    virtual bool isUnaryOperator() const {return true;}
    virtual Evaluable* createEvaluableUnary(Evaluable* right) const {
      throw ExpressionError("called unimplemented UnaryOpToken::createEvaluable;");
    }
  };

  // ...and a binary operator
  class BinaryOpToken: public  Token {
  public:
    BinaryOpToken(size_t begin=0): Token(begin) {}
    virtual bool isOperator() const {return true;}
    virtual bool isBinaryOperator() const {return true;}
    virtual Evaluable* createEvaluableBinary(Evaluable* left, Evaluable* right) const {
      throw ExpressionError("called unimplemented BinaryOpToken::createEvaluable;");
    }
  };

  // Tokenizer function: takes as argument a Token which will be used to
  // tokenize strings that don't match any standard ones
  std::list<Token*> tokenize(const std::string& input,
			     const Token& lastToken = Token());

  // Parser function: takes iterators in a container of tokens and turns them into 
  // an Evaluable tree.
  // Returns ptr to the Evaluable at root of tree.  Delete it to delete full tree.

  extern Evaluable* parse(std::list<Token*>& tokenList,
			  std::list<Token*>::iterator& begin,
			  std::list<Token*>::iterator end,
			  const std::string& input);
    
  ///// Above is the framework.  Below here are the implementations of the
  ///// standard tokens and Evaluables

} // end namespace expressions
#endif // EXPRESSIONS_H
