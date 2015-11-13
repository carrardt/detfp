#include <mpi.h>
#include "IFloat64.h"

double if64AllReduceSum(uint64_t n, const double * x, MPI_Comm comm)
{
    IFloat64 radd;
    radd.addValues( end-start , x+start );
    radd.removeCarries();

    MPI_Allreduce( MPI_IN_PLACE, radd.msum , n , MPI_LONG, MPI_SUM, comm );

    // mantissa normalization
    radd.computeCarriesFromMantissas();
    radd.removeCarries();

    // somme des mantisses
    return radd.sumMantissas();
}

