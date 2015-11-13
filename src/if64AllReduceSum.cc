#include <mpi.h>
#include <assert.h>
#include "IFloat64.h"

double if64AllReduceSum(uint64_t n, const double * x, MPI_Comm comm)
{
    assert( sizeof(int64_t) == sizeof(long long int) );

    IFloat64 tmp;
    tmp.addValues(n,x);
    tmp.removeCarries();

    MPI_Allreduce( MPI_IN_PLACE, tmp.msum , n , MPI_LONG_LONG_INT, MPI_SUM, comm );

    tmp.removeCarries();
    return tmp.sumMantissas();
}

