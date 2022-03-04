// Routines for parameter sets
#include <string>
#include <cctype>
#include "Pset.h"


void Pset::addMemberNoValue (const char *k,
			     const int _f,
			     const char *c) {
  l.push_back( new PsetMember(k, _f, c) ); }


//??? discern null strings from default string request?

// Take an input line, pull out strings representing keyword & value.
// Returns true if this line has information.  value string may be empty
// if no value is specified, but keyword should always have something.
// Throw exception for an unterminated quoted string.
bool
Pset::read_keyvalue(const string &in, string &keyword, string &value) {

  const char quote='"';	//will bound a quoted string keyword or value
  const string comments="#;";	//starts comment (unless quoted)
  const char eq='=';	//can be present after keyword

  keyword.resize(0);
  value.resize(0);

  string::size_type i=0;
  string::size_type l=in.length();

  // skip opening white space:
  while (i<l && isspace(in[i])) ++i;
  if (i==l) return false;


  if (in[i]==quote) {
    //If keyword is quoted, get everything until next quote
    ++i;	//skip the quote
    while (i<l && in[i]!=quote) keyword+=in[i++];
    if (i==l) throw PsetUnboundedQuote(in);
    ++i;	//skip the quote
  } else {
    //Keyword ends at next white, =, comment, or end:
    while (i<l 
	   && in[i]!=quote
	   && in[i]!=eq
	   && !isspace(in[i])
	   && comments.find(in[i])==string::npos ) keyword+=in[i++];
  }
  if (keyword.length()==0) return false;	//no useful info.

  // Skip whitespace or equals; done for end or comment
  while (i<l && (in[i]==eq || isspace(in[i])) ) ++i;

  if (i==l) return true;	//end - keyword with no value.
  if (comments.find(in[i])!=string::npos) return true;	//comment

  // Get the value string:
  if (in[i]==quote) {
    //If value is quoted, get everything until next quote
    ++i;	//skip the quote
    while (i<l && in[i]!=quote) value+=in[i++];
    if (i==l) throw PsetUnboundedQuote(in);
    ++i;	//skip the closing quote
  } else {
    //Value ends at comment, or end:
    while (i<l 
	   && comments.find(in[i])==string::npos ) value+=in[i++];
    // Strip trailing whitespace if was not a quoted response:
    string::size_type pos = value.size();
    while (pos > 0 && isspace(value[pos - 1])) pos--;
    value.erase(pos);
  }

  return true;	//could warn here if there are extra characters...
}

int 
Pset::setFromArguments(int argc, char **argv) {
  int firstKeyword = 0;
  while (firstKeyword < argc && argv[firstKeyword][0]!='-') 
    firstKeyword++;

  for (int iarg = firstKeyword; iarg < argc; iarg++) {
    if (argv[iarg][0] != '-') {
      // Should have a keyword specified here
      throw PsetError("Expected cmd-line argument beginning with '-', got: " 
		      + string(argv[iarg]));
    }

    string arg(argv[iarg]);
    arg.erase(0,1);	// Pop the '-' off the front of the argument
    string keyword;
    string value;

    size_t equalsPosn = arg.find('=');
    bool gotEquals = true;
    if (equalsPosn==string::npos) {
      gotEquals = false;
      stringstuff::stripWhite(arg);
      keyword = arg;
    } else {
      gotEquals = true;
      keyword = arg.substr(0, equalsPosn);
      stringstuff::stripWhite(keyword);
      if (keyword.empty())
	throw PsetError("Missing keyword at cmd-line argument: " 
			+ string(argv[iarg]));
      if (equalsPosn < arg.size()-1) {
	value = arg.substr(equalsPosn+1);
	stringstuff::stripWhite(value);
	if (!value.empty()) {
	  // Both keyword and value were in this argument:
	  setKeyValue(keyword,value);
	  continue;
	}
      }
    }

    // At this point we are expecting the value in the next keyword.
    iarg++;
    if (!gotEquals && iarg < argc) {
      // Although it might be a free-floating '=' sign we want to ignore
      string next = argv[iarg];
      stringstuff::stripWhite(next);
      if (next=="=") {
	gotEquals = true;
	iarg++;
      }
    }

    if (iarg >= argc)
      throw PsetError("Ran out of command-line arguments awaiting a value");

    value = argv[iarg];

    // Drop a leading '=' if we did not get one already
    if (!gotEquals && !value.empty() && value[0]=='=') {
      gotEquals = true;
      value.erase(0,1);
    }

    setKeyValue(keyword,value);
  }
  return firstKeyword;
}

