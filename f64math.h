#ifndef __f64math_h
#define __f64math_h

#include <cstdint>


double f64SumNoOpt(uint64_t n, const double * x) __attribute__((optimize("no-tree-vectorize"))) ;

double f64Sum(uint64_t n, const double * x);
double f64Sumi128(uint64_t n, const double * x);

#endif

