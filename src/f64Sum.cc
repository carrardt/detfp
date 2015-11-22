#include <cstdint>
#include "f64math.h"

double f64SumNoOpt(uint64_t n, const double * x) 
{
  double r=0.0;
  for(int i=0; i<n; i++) { r += x[i]; }
  return r;
}

double f64Sum(uint64_t n, const double * x)
{
  double r = 0.0;
#pragma omp parallel for simd reduction(+:r)
  for(int i=0; i<n; i++) { r += x[i]; }
  return r;
}

double f64Sumi128(uint64_t n, const double * x)
{
#ifdef __clang__
  long double r = 0.0;
#else
  __float128 r = 0.0;
#endif

#pragma omp parallel for simd reduction(+:r)
  for(int i=0; i<n; i++) { r += x[i]; }
  return r;
}

