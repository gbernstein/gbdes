// tokenizer and parser for expressions

#include "Expressions.h"
#include "Expressions2.h"

#include <cctype>
#include <typeinfo>

using namespace expressions;
using std::list;

std::list<Token*> 
expressions::tokenize(const std::string& input,
		      const Token& lastToken) {
  // Static list of all special-character tokens to attempt:
  static list<Token*> trialTokens;
  if (trialTokens.empty()) {
    // Initialize list of all token types we want to look for:
    trialTokens.push_back( new StringConstantToken());
    trialTokens.push_back( new OpenParenthesis());
    trialTokens.push_back( new CloseParenthesis());
    trialTokens.push_back( new PlusToken());  // ** note binary included here ahead of */
    trialTokens.push_back( new MinusToken());
    trialTokens.push_back( new PowerToken()); 
    trialTokens.push_back( new MultipliesToken()); 
    trialTokens.push_back( new DividesToken());
    trialTokens.push_back( new ModulusToken());
    trialTokens.push_back( new GreaterEqualToken()); // >= before >
    trialTokens.push_back( new LessEqualToken());    // <= before <
    trialTokens.push_back( new NotEqualToken());     // != before !
    trialTokens.push_back( new NotToken()); 
    trialTokens.push_back( new EqualToken());        // == or =
    trialTokens.push_back( new GreaterToken());
    trialTokens.push_back( new LessToken()); 
    trialTokens.push_back( new GreaterEqualToken()); // >= before >
    trialTokens.push_back( new LessEqualToken());    // <= before <
    trialTokens.push_back( new NotEqualToken());     // != before !
    trialTokens.push_back( new GreaterToken());
    trialTokens.push_back( new LessToken()); 
    trialTokens.push_back( new EqualToken());        // == or =
    trialTokens.push_back( new AndToken());	     // Check for && before &
    trialTokens.push_back( new OrToken());           // And || before |
    trialTokens.push_back( new BitAndToken()); 
    trialTokens.push_back( new BitOrToken()); 
  }
  // Static list of function names (all start with letters):
  static list<Token*> functionTokens;
  if (functionTokens.empty()) {
    functionTokens.push_back( new SinToken());
    // ??? more...
  }

  size_t begin = 0;
  size_t end = input.size();
  std::list<Token*> output;
  bool lastWasOp = false;
  try {
    while (begin < end) {
      // strip leading white spaces
      if (std::isspace(input[begin])) {
	++begin;
	continue;
      }
      Token* next = 0;
      for (std::list<Token*>::const_iterator i = trialTokens.begin();
	   !next && (i != trialTokens.end());
	   ++i) {
	next = (*i)->createFromString(input, begin, end, lastWasOp);
      }
      // Look for a numerical constant:
      if (!next) next = NumberToken().createFromString(input, begin, end, lastWasOp);

      if (next) {
	output.push_back(next);
	lastWasOp = next->isOperator();
	continue;
      }

      // This token starts with a non-digit and matches none of the operators.
      // Collect characters until whitespace or match an operator
      size_t tokenEnd = begin;
      Token* following = 0;
      while (tokenEnd < end && !following) {
	tokenEnd++;
	// End token at white space
	if (std::isspace(input[tokenEnd])) break;
	// Or if any other token matches:
	size_t nextBegin = tokenEnd;
	for (std::list<Token*>::const_iterator i = trialTokens.begin();
	     !following && (i != trialTokens.end());
	     ++i) {
	  following = (*i)->createFromString(input, nextBegin, end, false);
	}
      }
      if (following) delete following;

      // At this point we have terminated the token.  Try function names first:
      for (std::list<Token*>::const_iterator i = functionTokens.begin();
	   !next && (i != functionTokens.end());
	   ++i) {
	next = (*i)->createFromString(input, begin, tokenEnd, lastWasOp);
      }
      if (!next) {
	// Last resort is to ask the lastToken to parse:
	next = lastToken.createFromString(input, begin, tokenEnd, lastWasOp);
      }
      if (!next) {
	// We have an undecipherable token:
	throw SyntaxError("Unknown token", begin);
      }
      // Save the token found
      output.push_back(next);
      lastWasOp = next->isOperator();
      // and continue tokenizing where it ended
      begin = tokenEnd;
    } // End character loop

  } catch (SyntaxError& e) {
    // Make a longer explanation text from full expression string:
    std::string err = e.longWhat(input);
    throw SyntaxError(err, e.charBegin);
  }
  return output;
}


Evaluable*
expressions::parse(std::list<Token*>& tokenList,
		   std::list<Token*>::iterator& begin,
		   std::list<Token*>::iterator end,
		   const std::string& input) {
  typedef list<Token*>::iterator iter;

  // This table gives precedence groups for all binary operators.
  static list<list<BinaryOpToken*> > precedenceTable;
  if (precedenceTable.empty()) {
    // Exponentiation is above multiplies, below unaries.
    precedenceTable.push_back( list<BinaryOpToken*>(1, new PowerToken()));
    {
      list<BinaryOpToken*> multiplicative;
      multiplicative.push_back(new MultipliesToken());
      multiplicative.push_back(new DividesToken());
      multiplicative.push_back(new ModulusToken());
      precedenceTable.push_back(multiplicative);
    }
    {
      list<BinaryOpToken*> additive;
      additive.push_back(new BinaryPlusToken());
      additive.push_back(new BinaryMinusToken());
      precedenceTable.push_back(additive);
    }
    {
      list<BinaryOpToken*> relational;
      relational.push_back(new LessToken());
      relational.push_back(new GreaterToken());
      relational.push_back(new LessEqualToken());
      relational.push_back(new GreaterEqualToken());
      precedenceTable.push_back(relational);
    }
    {
      list<BinaryOpToken*> equality;
      equality.push_back(new EqualToken());
      equality.push_back(new NotEqualToken());
      precedenceTable.push_back(equality);
    }
    precedenceTable.push_back(std::list<BinaryOpToken*>(1,new BitAndToken()));
    precedenceTable.push_back(std::list<BinaryOpToken*>(1,new BitOrToken()));
    precedenceTable.push_back(std::list<BinaryOpToken*>(1,new AndToken()));
    precedenceTable.push_back(std::list<BinaryOpToken*>(1,new OrToken()));
  } // End preparing precedence table

  try {
    if (begin==end) {
      if (begin==tokenList.end())
	throw SyntaxError("Empty expression", tokenList.size());
      else
	(*begin)->throwSyntaxError("Empty expression");
    }
    for (iter i=begin; i!=end; ++i) {
      // First look for any open-parentheses
      if ( dynamic_cast<OpenParenthesis*> (*i)) {
	// Look for close-parenthesis to match
	iter j = i;
	int pLevel = 1;
	for (++j; j!=end; ++j) {
	  if ( dynamic_cast<OpenParenthesis*> (*j))  ++pLevel;
	  if ( dynamic_cast<CloseParenthesis*> (*j)) --pLevel;
	  if (pLevel==0) break;
	}
	if (j == end) 
	  (*i)->throwSyntaxError("Mismatched open parenthesis");

	// Replace stretch of Tokens between parentheses (inclusive)
	// with a placeholder having its Evaluable
	iter afterOpen = i;
	++afterOpen;
	Evaluable* contents = parse(tokenList, afterOpen, j, input);
	delete *j;
	tokenList.erase(j);
	delete *i;
	*i = new EvaluableToken(contents);
      }
    }

    // Should not encounter a close-parentheses first.
    if ( dynamic_cast<CloseParenthesis*> (*begin)) {
      (*begin)->throwSyntaxError("Mismatched closed parenthesis");
    }

    // search for all unary operators, back to front
    iter i=end;
    do {
      --i;
      if (!(*i)->isUnaryOperator()) continue;
      iter rightIter = i;
      ++rightIter;
      if (rightIter==end || (*rightIter)->isOperator() ) 
	(*i)->throwSyntaxError("Missing argument for unary operator");
      Evaluable* right = (*rightIter)->createEvaluable();
      delete *rightIter;
      tokenList.erase(rightIter);
      UnaryOpToken* uot = dynamic_cast<UnaryOpToken*> (*i);
      Assert(uot);
      Evaluable* contents = uot->createEvaluableUnary(right);
      delete *i;
      *i = new EvaluableToken(contents);
    } while (i != begin);

    // Now descend binary-operator precedence list
    for ( list<list<BinaryOpToken*> >::iterator i1= precedenceTable.begin();
	  i1 != precedenceTable.end();
	  ++i1) {
      // Look for tokens of this precedence level, from left to right
      for ( iter j = begin;
	    j != end;
	    ++j) {
	bool foundOne = false;
	for (list<BinaryOpToken*>::iterator i2 = i1->begin();
	     !foundOne && i2 != i1->end();
	     ++i2) { 
	  foundOne = (typeid(**i2) == typeid(**j));
	}
	if (foundOne) {
	  // Build this binary op from (everything to left) Op (immediately to right)

	  // Get left operand and kill its Token
	  if (j == begin)
	    (*j)->throwSyntaxError("Missing left operand");
	  iter leftIter = j;
	  if ( leftIter==begin || (*(--leftIter))->isOperator())
	    (*j)->throwSyntaxError("Missing left operand");
	  Evaluable* left = (*leftIter)->createEvaluable();
	  delete *leftIter;
	  if (leftIter == begin) begin = j;
	  tokenList.erase(leftIter);

	  // Get right operand and kill its Token
	  iter rightIter = j;
	  ++rightIter;
	  if (rightIter==end || (*rightIter)->isOperator())
	    (*j)->throwSyntaxError("Missing right operand");
	  Evaluable* right = (*rightIter)->createEvaluable();
	  delete *rightIter;
	  tokenList.erase(rightIter);

	  BinaryOpToken* bot = dynamic_cast<BinaryOpToken*> (*j);
	  Assert(bot);
	  Evaluable* contents = bot->createEvaluableBinary(left, right);
	  // delete the BinOpToken and replace by EvaluableToken
	  delete (*j);
	  *j = new EvaluableToken(contents);
	}
      } // End iteration through tokenList
    } // End descent of precedence table

    // At this point we should be down to a single non-operator element
    if ( begin == end) 
      throw SyntaxError("parsed down to nothing!!",0);  // where to point??
    if ( (*begin)->isOperator()) 
      (*begin)->throwSyntaxError("parsed down to an operator!!");
    Evaluable* result = (*begin)->createEvaluable();
    delete *begin;
    begin = tokenList.erase(begin);
    if (begin != end) 
      throw SyntaxError("Missing operator", 0); // Not clear where to point

    return result;
  } catch (SyntaxError& e) {
    std::string err = e.longWhat(input);
    throw SyntaxError(err, e.charBegin);
  }
}
