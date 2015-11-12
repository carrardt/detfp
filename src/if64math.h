#ifndef __if64Sum_h
#define __if64Sum_h

#include <cstdint>
#include "IFloat64.h"

static inline double if64Sum(uint64_t n, const double * x, IFloat64 & if64)
{
    if64.addValues(n,x);
    if64.removeCarries();
    return if64.sumMantissas();
}

static inline double if64Sum(uint64_t n, const double * x)
{
    IFloat64 if64;
    return if64Sum(n,x,if64);
}

#endif

