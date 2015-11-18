#include "IFloat64.h"

double if64Sum(uint64_t n, const double * x, IFloat64 & if64)
{
    if64.addValues(n,x);
    if64.removeCarries();
    return if64.toDouble();
}

double if64Sum(uint64_t n, const double * x)
{
    IFloat64 if64;
    if64.addValues(n,x);
    if64.removeCarries();
    return if64.toDouble();
}

