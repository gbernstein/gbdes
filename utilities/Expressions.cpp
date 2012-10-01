// tokenizer and parser for expressions

#include "Expressions.h"
#include <cctype>

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
    trialTokens.push_back( new PlusToken());
    trialTokens.push_back( new MinusToken());
    trialTokens.push_back( new PowerToken()); // ** before *
    trialTokens.push_back( new MultipliesToken()); 
    trialTokens.push_back( new DividesToken());
    trialTokens.push_back( new ModulusToken());
    trialTokens.push_back( new AndToken());
    trialTokens.push_back( new OrToken());
    trialTokens.push_back( new GreaterEqualToken()); // >= before >
    trialTokens.push_back( new LessEqualToken());    // <= before <
    trialTokens.push_back( new NotEqualToken());     // != before !
    trialTokens.push_back( new EqualToken());        // == or =
    trialTokens.push_back( new GreaterToken());
    trialTokens.push_back( new LessToken()); 
    trialTokens.push_back( new NotToken()); 
    trialTokens.push_back( new AndToken()); 
    trialTokens.push_back( new OrToken()); 
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
	next = i->createFromString(input, begin, end);
      }
      // Look for a numerical constant:
      if (!next) next = NumericalToken().createFromString(input, begin, end);

      if (next) {
	output.push_back(next);
	continue;
      }

      // This token starts with a letter or non-special character.  Record all such:
      size_t tokenEnd = begin+1;
      Token* following = 0;
      size_t nextBegin;
      while (tokenEnd < end && !following) {
	// End token at white space
	if (std::isspace(input[tokenEnd])) break;
	// Or if any other token matches:
	nextBegin = tokenEnd;
	for (std::list<Token*>::const_iterator i = trialTokens.begin();
	     !following && (i != trialTokens.end());
	     ++i) {
	  following = i->createFromString(input, nextBegin, end);
	}
      }
      if (following) delete following;

      // At this point we have terminated the token.  Try function names first:
      for (std::list<Token*>::const_iterator i = functionTokens.begin();
	   !next && (i != functionTokens.end());
	   ++i) {
	next = i->createFromString(input, begin, tokenEnd);
      }
      if (!next) {
	// Last resort is to ask the lastToken to parse:
	next = lastToken.createFromString(input, begin, tokenEnd);
      }
      if (!next) {
	// We have an undecipherable token:
	throw ExpressionSyntaxError("Unknown token", begin);
      }
      // Save the token found
      output.push_back(next);
      // and continue tokenizing where it ended
      begin = tokenEnd;
    } // End character loop

  } catch (ExpressionSyntaxError& e) {
    // Make a longer explanation text from full expression string:
    string err = e.longWhat(input);
    throw ExpressionSyntaxError(err, e.charBegin);
  }
  return output;
}

// For the parser, we'll need something that can sit in list<Token*> but
// is already parsed.
// Should always be a value, not an operator.
class ParsedToken: public Token {
private:
  Evaluable* contents;
  ParsedToken(Evaluable* contents_): Token(0), contents(contents_) {}
  virtual ~ParsedToken() {} // Don't delete contents, it will be passed on
  virtual Evaluable* createEvaluable() const {return contents;}
};

// Parse will delete all the Tokens in the range it parses
// and remove them from the tokenList.
Evaluable*
expressions::parse(std::list<Token>& tokenList,
		   std::list<Token*>::iterator& begin,
		   std::list<Token*>::iterator end,
		   const std::string& input) {
  typedef list<Token*>::iterator iter;

  // This table gives precedence groups for all binary operators.
  static list<list<BinaryOpToken*> > precedenceTable;
  if (precedenceTable.empty()) {
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
    precedenceTable.push_back(List<BinaryOpToken*>(1,new AndToken()));
    precedenceTable.push_back(List<BinaryOpToken*>(1,new OrToken()));
  } // End preparing precedence table

  try {
    if (begin==end) {
      if (begin==tokenList.end())
	throw SyntaxError("Empty expression", tokenList.size());
      else
	begin->throwSyntaxError("Empty expression");
    }
    for (iter i=begin; i!=end; ++i) {
      // First look for any open-parentheses
      if ( dynamic_cast<OpenParenthesis*> (*i)) {
	// Look for close-parenthesis to match
	iter j = end;
	for (--j; j!= i; --j)
	  if ( dynamic_cast<CloseParenthesis*> (*i)) break;

	if (j == i) 
	  (*i)->throwSyntaxError("Mismatched open parenthesis");

	// Replace stretch of Tokens between parentheses (inclusive)
	// with a placeholder having its Evaluable
	iter afterOpen = i;
	++afterOpen;
	iter afterClose = j;
	afterClose++;
	Evaluable* contents = parse(tokenList, afterOpen, j, input);
	delete *j;
	tokenList.erase(j);
	delete *i;
	*i = new ParsedToken(contents);
      }
    }
    // No open-parens.  Should therefore be no close-parens left:
    for (iter i=begin; i!=end; ++i) {
      // First look for any open-parentheses
      if ( dynamic_cast<CloseParenthesis*> (*i)) {
	(*i)->throwSyntaxError("Mismatched closed parenthesis");
      }
    }

    // search for all unary operators, back to front
    iter i=end;
    do {
      --i;
      if (!(*i)->isUnaryOperator()) continue;
      iter rightIter = i;
      ++rightIter;
      if (rightIter==end || (*rightIter)->isOperator() ) 
	(*i)->throwSyntaxError("Missing argument for unanry operator");
      Evaluable* right = (*rightIter)->createEvaluable();
      delete *rightIter;
      tokenList.erase(rightIter);
      UnaryOpToken* uot = dynamic_cast<UnaryOpToken*> (*i);
      Assert(uot);
      Evaluable* contents = uot->createEvaluable(right);
      delete *i;
      *i = new ParsedToken(contents);
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
	  foundOne = std::typeid(*i2) == std::typeid(*j);
	}
	if (foundOne) {
	  // Build this binary op from (everything to left) Op (immediately to right)

	  // Get left operand and kill its Token
	  if (j == begin)
	    (*j)->throwSyntaxError("Missing left operand");
	  iter leftIter = j;
	  if ( leftIter==begin || (*(--leftIter))->isOperand())
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
	  Evaluable* right = (*righIter)->createEvaluable();
	  delete *rightIter;
	  tokenList.erase(rightIter);

	  BinaryOpToken* bot = dynamic_cast<BinaryOpToken*> (*j);
	  Assert(bot);
	  Evaluable* contents = bot->createEvaluable(left, right);
	  // delete the BinOpToken and replace by ParsedToken
	  delete (*j);
	  *j = new ParsedToken(contents);
	}
      } // End iteration through tokenList
    } // End descent of precedence table

    // At this point we should be down to a single non-operator element
    if ( begin == end) 
      throw SyntaxError("parsed down to nothing!!",0);  // where to point??
    if ( (*begin)->isOperator()) 
      (*begin)->throwsyntaxError("parsed down to an operator!!");
    Evaluable result = (*begin)->createEvaluable();
    delete *begin;
    begin = tokenList.erase(begin);
    if (begin != end) 
      throw SyntaxError("Missing operator", 0); // Not clear where to point

    return result;
  } catch (SyntaxError& e) {
    string err = e.longWhat(input);
    throw SyntaxError(err, e.charBegin);
  }
}

    
	
      
    
