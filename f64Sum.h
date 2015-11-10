#ifndef __f64sum_h
#define __f64sum_h

template<typename AccumT>
double f64sum(int n, const double * x)
{
  AccumT r = 0.0;
#pragma omp simd reduction(+:r)
  for(int i=(n-1); i>=0; i--) { r += x[i]; }
  return r;
}

#endif

