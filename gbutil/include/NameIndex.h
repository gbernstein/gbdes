#ifndef NAMEINDEX_H
#define NAMEINDEX_H

#include "Std.h"
#include <map>

// Class that maintains a map between a set of string-valued names and
// sequential zero-indexed integers
class NameIndex {
public:
  NameIndex() {}
  int size() const {return m.size();}
  bool has(string s) const {return m.find(s) != m.end();}
  int indexOf(string s) const {
    std::map<string,int>::const_iterator p=m.find(s);
    if (p==m.end())
      return -1;
    else
      return p->second;
  }
  string nameOf(int i) const {
    if (i>=0 && i<m.size())
      return names[i];
    else
      FormatAndThrow<std::runtime_error>() << "NameIndex::nameOf() index out of bounds: " << i;
   return " "; // Should not get here.
  }
  int append(string s) {
    if (has(s)) throw std::runtime_error("NameIndex::append() argument <" + s + "> already exists.");
    names.push_back(s);
    m.insert(std::pair<string,int>(s,names.size()-1));
    return m[s];
  }
private:
  std::map<string, int> m;
  vector<string> names;
};
#endif
