// 2-dimensional array of arbitrary type.  No mathematical intent.
//---------------------------------------------------------------------------
#ifndef ArrayH
#define ArrayH
//---------------------------------------------------------------------------

#include <vector>
#include "Std.h"

template <class T>
class Array {

public:

  explicit Array(uint _n1=0,uint _n2=0) : n1(_n1),n2(_n2),
    data(_n1,myvector<T>(_n2)) {}
  Array(uint _n1,uint _n2,const T& value) : n1(_n1),n2(_n2),
        data(_n1,myvector<T>(_n2,value)) {}
  ~Array() {}
  Array(const Array& rhs) : n1(rhs.n1),n2(rhs.n2),data(rhs.data) {}
  Array<T>& operator=(const Array<T>& rhs)
    { if (this == &rhs) return *this; n1 = rhs.n1; n2 = rhs.n2; 
      data = rhs.data; return *this; }
  std::vector<T>& operator[](uint i) {return data[i]; }
  const std::vector<T>& operator[](uint i) const {return data[i];}

  uint getN1() const {Assert(data.size()==n1); return n1;}
  uint getN2() const {Assert(data.size()==0 || data[0].size()==n2);
        return n2;}
  const T& get(uint i,uint j) const 
    {Assert(i<getN1()); Assert(j<getN2()); return data[i][j]; }
  T& get(uint i,uint j) 
    {Assert(i<getN1()); Assert(j<getN2()); return data[i][j]; }
  const T& operator()(uint i,uint j) const {return get(i,j);}
  T& operator()(uint i,uint j) {return get(i,j);}

  const std::vector<T>& getRow(uint i) const
    {Assert(i<getN1()); return data[i]; }
  std::vector<T> getCol(uint j) const
    {Assert(j<getN2()); 
      std::vector<T> temp(getN1()); 
      for(uint i=0;i<getN1();i++) temp[i] = data[i][j]; 
      return temp; }
  void SetN1(uint newn1);
  void SetN2(uint newn2);
  void SetAllValues(const T& value)
    { for(uint i=0;i<getN1();i++) for(uint j=0;j<getN2();j++) 
        data[i][j] = value; }
  void Set(uint i,uint j,const T& value) {data[i][j] = value;}
  void SetRow(uint i,const std::vector<T>& row) 
    { Assert(i<getN1());
      Assert(row.size()==getN2());
      for(uint j=0;j<getN2();j++) data[i][j] = row[j]; }
  void SetCol(uint j,const std::vector<T>& col)
    { Assert(j<getN2());
      Assert(col.size()==getN1());
      for(uint i=0;i<getN1();i++) data[i][j] = col[i]; }
  void SetRow(uint i,const std::vector<T>& row,uint start,uint end) 
    { Assert(i<getN1());
      Assert(start <= end);
      Assert(end < row.size());
      Assert(end-start+1==getN2());
      for(uint j=0;j<getN2();j++) data[i][j] = row[j+start]; }
  void SetCol(uint j,const std::vector<T>& col,uint start,uint end)
    { Assert(j<getN2());
      Assert(start <= end);
      Assert(end < col.size());
      Assert(end-start+1==getN1());
      for(uint i=0;i<getN1();i++) data[i][j] = col[i+start]; }
  void AddRow(std::vector<T> row) 
    { Assert(row.size()==getN2());
      data.push_back(row); n1++; }
  void AddCol(std::vector<T> col) 
    { Assert(col.size()==getN1());
      for(uint i=0;i<getN1();i++) data[i].push_back(col[i]); n2++; }
  void Reserve(uint n1r,uint n2r);

protected :

  uint n1,n2;
  std::vector<std::vector<T> > data;
};

template <class T>
inline void Array<T>::SetN1(uint newn1)
{
  if (n1 == 0) data = std::vector<std::vector<T> >(newn1,std::vector<T>(n2));
  else if (newn1 > n1) data.insert(data.end(),newn1-n1,std::vector<T>(n2));
  else data.erase(data.begin()+newn1,data.end());
  n1 = newn1;
}

template <class T>
inline void Array<T>::SetN2(uint newn2)
{
  if (n2 == 0) for(uint i=0;i<n1;i++) data[i] = std::vector<T>(newn2);
  else if (newn2 > n2) for(uint i=0;i<n1;i++)
      data[i].insert(data[i].end(),newn2-n2,0.);
  else for(uint i=0;i<n1;i++)
      data[i].erase(data[i].begin()+newn2,data[i].end());
  n2 = newn2;
}

template <class T>
inline void Array<T>::Reserve(uint n1r,uint n2r)
{
  if (n1r > n1) data.reserve(n1r);
  if (n2r > n2) for(uint i=0;i<n1;i++) data[i].reserve(n2r);
}

#endif
