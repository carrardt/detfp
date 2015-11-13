#ifndef __if64Sum_h
#define __if64Sum_h

#include <cstdint>
#include "IFloat64.h"

double if64Sum(uint64_t n, const double * x);
double if64Sum(uint64_t n, const double * x, IFloat64 & if64 );

#endif

