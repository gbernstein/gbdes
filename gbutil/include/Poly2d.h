// 2-dimensional polynomials, using matrix algebra.
#ifndef POLY2D_H
#define POLY2D_H
#include "LinearAlgebra.h"
#include "Std.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif

namespace poly2d {

  class Poly2dError: public std::runtime_error {
  public:
    Poly2dError(const std::string &m=""): std::runtime_error("Poly2d Error: " +m) {}
  };

  class Poly2d {
  public:
    // This enum determines whether polynomial goes to max
    // order in each x and y or only in sum of x & y orders
    enum OrderType {Sum, Each};
    // Construct zero-coeff polynomial with max total summed order
    Poly2d(int order_): 
      orderX(order_), orderY(order_), otype(Sum), cm(order_+1, order_+1, 0.) {}
    // Construct zero-coeff polynomial up to given orders in each x and y:
    Poly2d(int orderX_, int orderY_): 
      orderX(orderX_), orderY(orderY_), otype(Each), cm(orderX_+1, orderY_+1, 0.) {}
    // Construct using coefficients from input matrix
    // Matrix dims give poly orders.  Must be square if otype==Sum.
    Poly2d(linalg::DMatrix m_, OrderType otype_=Sum);
    ~Poly2d() {}
    double evaluate(double x, double y) const {return powers(x,orderX).transpose() * cm * powers(y,orderY);}
    double operator() (double x, double y) const {return evaluate(x,y);}
    double derivx(double x, double y) const {return derivs(x,orderX).transpose() * cm * powers(y,orderY);}
    double derivy(double x, double y) const {return powers(x,orderX).transpose() * cm * derivs(y,orderY);}
    // Access coefficients as vector:
    int nCoeffs() const;
    linalg::DVector getC() const {return vectorFromMatrix(cm);}
    linalg::DMatrix getM() const {return cm;}
    void setC(const linalg::DVector& cv) {fillFromVector(cv);}
    // derivatives w.r.t. coefficient vector:
    linalg::DVector derivC(double x, double y) const {
      return vectorFromMatrix(powers(x,orderX).outer(powers(y,orderY)));} 
    OrderType getOrderType() const {return otype;}
    int getOrderX() const {return orderX;}
    int getOrderY() const {return orderY;}

    // Serialize polynomial to a string:
    // Number of digits to write:
    static const int DEFAULT_PRECISION=8;
    void write(std::ostream& os, int precision=DEFAULT_PRECISION) const;  
    static Poly2d* create(std::istream& is);  // Build polynomial from serialized string

#ifdef USE_YAML
    // Serialize to/from YAML
    void write(YAML::Emitter& os) const;
    static Poly2d* create(const YAML::Node& node,
			  bool& defaulted);
#endif
    
    // If you want to know the orders for given coeffient index or vice-versa:
    // Negative indices/powers returned if inputs are negative.  No checks on upper bounds.
    int vectorIndex(int i, int j) const;  // vector index of coefficient of x^i y^j
    void powersOfIndex(int n, int& i, int& j) const;  // return powers i & j of nth coefficient
  private:
    const int orderX;
    const int orderY;
    const OrderType otype;
    linalg::DMatrix cm;  // coefficients
    linalg::DVector powers(double z, int order) const;
    linalg::DVector derivs(double z, int order) const;
    // Indexing conventions embodied in these:
    linalg::DVector vectorFromMatrix(const linalg::DMatrix& m) const;
    void fillFromVector(const linalg::DVector& v);
  };

} // namespace poly2d

#endif // POLY2D_H
