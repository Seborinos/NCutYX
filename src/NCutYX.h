#ifndef _NCutYX_H
#define _NCutYX_H

#include <RcppEigen.h>//This is new

using namespace Rcpp;
using namespace Eigen;

double NCutY3V1(const NumericMatrix C, const NumericMatrix Wy,
                       const NumericMatrix Wx) ;

double NCutLayer3V1(const NumericMatrix C, const NumericMatrix Wz,
                           const NumericMatrix Wy,const NumericMatrix Wx,
                           const NumericMatrix Wzyx) ;
#endif
