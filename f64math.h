#ifndef __f64math_h
#define __f64math_h

template<typename AccumT>
double f64sum(int n, const double * x)
{
  AccumT r = 0.0;

#ifdef _OPENMP
#pragma omp simd reduction(+:r)
#endif
  for(int i=(n-1); i>=0; i--) { r += x[i]; }

  return r;
}

#endif

