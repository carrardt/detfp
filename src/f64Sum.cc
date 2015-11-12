#include <cstdint>

template<typename AccumT>
double f64SumT(int n, const double * x)
{
  AccumT r = 0.0;
#pragma omp parallel for simd reduction(+:r)
  for(int i=0; i<n; i++) { r += x[i]; }
  return r;
}

double f64SumNoOpt(uint64_t n, const double * x)
{
  double r=0.0;
  for(int i=0; i<n; i++) { r += x[i]; }
  return r;
}


double f64Sum(uint64_t n, const double * x)
{
	return f64SumT<double>(n,x);
}

double f64Sumi128(uint64_t n, const double * x)
{
	return f64SumT<__float128>(n,x);
}

