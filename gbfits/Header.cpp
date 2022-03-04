// Class that mirrors the FITS header data structure.

#include "Header.h"
#include "StringStuff.h"
#include <cstring>
#include <cstdio>
#include <cctype>
#include <ios>
#include <limits>

using namespace img;
using namespace std;

// Case-raising & blank-stripping function for keywords.
string img::KeyFormat(const string input) {
  string output;
  //strip leading blanks
  int i=0;
  while (i<input.size() && isspace(input[i])) i++;
  //Convert to upper case, stop at whitespace
  while (i<input.size() && !isspace(input[i]))
    output += toupper(input[i++]);
  return output;
}

string
HdrRecordBase::writeCard() const {
  string vv=getValueString();
  string card= (keyword.size() <=8) ? keyword : "HIERARCH " + keyword;
  for (int i=card.size(); i<8; i++) card += " ";
  if (!vv.empty()) {
    card += "= ";
    for (int i=vv.size(); i<20; i++) card += " ";
    card += vv;
    if (!comment.empty() || !units.empty())
      card += " /";
  }
  card += " ";
  if (!units.empty()) card+= "[" + units + "] ";
  card += comment;
  return card;
}

namespace img {
  //specializations for bool
  template<>
  string
  HdrRecord<bool>::getValueString() const {
    if (val) return "T";
    else return "F";
  }

  template<>
  bool 
  HdrRecord<bool>::setValueString(const string _v) {
      istringstream iss(_v.c_str());
      string s;
      string leftover;
      // Should be only one word, T or F:
      if (!(iss >> s) || (iss >> leftover)) return true;
      if (s=="T" || s=="t" || s=="true") val=true;
      else if (s=="F" || s=="f" || s=="false") val=false;
      else return true;
      return false;
  };

  //and for string - enclose in quotes
  template<>
  string 
  HdrRecord<string>::getValueString() const {
    return "'" + val + "'";
  }

  //and for double - use capital E, many digits...
  template<>
  string 
  HdrRecord<double>::getValueString() const {
    ostringstream oss;
    oss << right << setw(20) << uppercase 
	<< showpoint << setprecision(12) << val;
    // Above puts lots of trailing zeroes... get rid of them:
    string work=oss.str();
    size_t eSpot = work.find('E');
    int nPad=0;  // Number of zeros removed
    bool hasExp = (eSpot != string::npos);
    // where is last non-zero part of mantissa?
    int last = (hasExp ? eSpot : work.size() ) -1;
    while (last>0 && work[last]=='0') {
      last--;
      nPad++;
    }
    string out;
    for (int i=0; i<nPad; i++) out += " ";
    out += work.substr(0,last+1);
    if (hasExp)
      out += work.substr(eSpot);
    return out;
  }

  // ??? Force decimal point on float, double?

  ///////////////////////////////////////////////
  // Headers to/from ASCII streams
  ///////////////////////////////////////////////
  HdrRecordBase*
  ReadASCIIHeader(string in) {

    const char quote='\'';	//will bound a quoted string keyword or value
    const string comments="/";	//starts comment (unless quoted)
    const char eq='=';	//can be present after keyword

    string keyword;
    string vstring;
    string comment;

    string::size_type i=0;
    string::size_type l=in.length();

    // skip opening white space:
    while (i<l && std::isspace(in[i])) ++i;
    if (i==l) return 0;

    //Keyword ends at next white, =, comment, or end:
    while (i<l 
	   && in[i]!=quote
	   && in[i]!=eq
	   && !std::isspace(in[i])
	   && comments.find(in[i])==string::npos ) keyword+=in[i++];

    // HIERARCH means an extension whereby keyword keeps going until = sign:
    if (keyword=="HIERARCH") {
      keyword.clear();
      while (i<l && std::isspace(in[i])) ++i;
      if (i==l) return 0;
      while (i<l && in[i]!=eq) keyword+=in[i++];
    }
    if (keyword.length()==0) return 0;	//no useful info.

    // Skip whitespace or equals; done for end or comment
    while (i<l && (in[i]==eq || std::isspace(in[i])) ) ++i;

    // Null record if we have nothing or comment left:
    if (i==l) return new HdrRecordNull(keyword);
    if (comments.find(in[i])!=string::npos) return new HdrRecordNull(keyword, in.substr(i+1));

    // A keyword named "HISTORY" or "COMMENT" is all comment, really
    if (keyword=="COMMENT" || keyword=="HISTORY")
      return new HdrRecordNull(keyword, in.substr(i));

    // A quoted value is string:
    if (in[i]==quote) {
      //If value is quoted, get everything until next quote
      ++i;	//skip the quote
      while (i<l && in[i]!=quote) vstring+=in[i++];
      if (i==l) return 0;	// unbounded quote, failure!!!
      ++i;	//skip the closing quote
      while (i<l && std::isspace(in[i])) i++;	// skip whitespace
      if (i==l) return new HdrRecord<string>(keyword, vstring);
      else if (comments.find(in[i])!=string::npos) // Comment left?
	return new HdrRecord<string>(keyword, vstring, in.substr(i+1));
      else return 0; // ??? failure - something other than comment after string
    }

    if (in[i]=='T' || in[i]=='F') {
      // Boolean valued:
      bool value= (in[i]=='T');
      i++;
      while (i<l && isspace(in[i])) i++;	// skip whitespace
      if (i==l) return new HdrRecord<bool>(keyword, value);
      else if (comments.find(in[i])!=string::npos) // Comment left?
	return new HdrRecord<bool>(keyword, value, in.substr(i+1));
      else return 0; // ??? failure - something other than comment after T/F
    }

    // Otherwise we are getting either an integer or a float (ignore complex)
    while (i<l && comments.find(in[i])==string::npos ) vstring+=in[i++];
    // Strip trailing whitespace if was not a quoted response:
    string::size_type pos = vstring.size();
    while (pos > 0 && isspace(vstring[pos - 1])) pos--;
    vstring.erase(pos);

    // Collect comment
    if (comments.find(in[i])!=string::npos) // Comment left?
      comment = in.substr(i+1);
    HdrRecord<int>* hi = new HdrRecord<int>(keyword, 0, comment);
    if (!hi->setValueString(vstring)) return hi;
    // If that failed, try a double
    delete hi;
    HdrRecord<double>* hd = new HdrRecord<double>(keyword, 0., comment);
    if (!hd->setValueString(vstring)) return hd;
    delete hd;
    // Last try: some DECam headers coming back with NAN.0 listed. 
    if (stringstuff::nocaseEqual(vstring.substr(0,3),"NAN")) {
      return new HdrRecord<double>(keyword,
				   std::numeric_limits<double>::quiet_NaN(),
				   comment);
    }
    /**/cerr << "Got ASCII header with bad value <" << vstring << "> keyword " << keyword << endl;
    return 0;	// Formatting error
  }

  istream& 
  operator>>(istream& is, Header& h) {
    string buffer;
    while (stringstuff::getlineNoComment(is, buffer)) {
      HdrRecordBase* hrb = ReadASCIIHeader(buffer);
      if (!hrb) {
	is.setstate(ios::failbit);  // ??? do we want to throw here?
	continue;
      }
      if (hrb->getKeyword()=="END") {
	delete hrb;
	return is;
      } else if (hrb->getKeyword()=="COMMENT") {
	h.addComment(hrb->getComment());
	delete hrb;
      } else if (hrb->getKeyword()=="HISTORY") {
	h.addHistory(hrb->getComment());
	delete hrb;
      } else {
	h.append(hrb);
      }
    }
    // Get here is we exhaust the istream before END.  eof should already be set.
    return is;
  }

  ostream&
  operator<<(ostream& os,
	     const img::Header& h) {
    for (h.rewind(); !h.atEnd(); h.incr())
      os << h.current()->writeCard() << endl;
    os << "END     " << endl;
    return os;
  }

}

