#include <mpi.h>
#include "IFloat64.h"

#ifdef _OPENMP

#include <omp.h>

double if64AllReduceSum(uint64_t n, const double * x, MPI_Comm comm)
{
    IFloat64 all_radd;

#   pragma omp parallel
    {
        IFloat64 radd;

        const int tid = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        const uint64_t start = (n*tid)/nthreads;
        const uint64_t end = (n*(tid+1))/nthreads;

        radd.insertValues( end-start , x+start );
        radd.removeCarries();

        // add thread mantissas to global bins
#       pragma omp critical
        {
            all_radd.addMantissasNoCarryUpdate( radd );
#           pragma omp flush(all_radd)
        }
	    //printExponentBins(bmin,bmax,msum,mcarry);
    }

    MPI_Allreduce( MPI_IN_PLACE, all_radd.msum , n , MPI_LONG, MPI_SUM, comm );

    // mantissa normalization
    all_radd.computeCarriesFromMantissas();
    all_radd.removeCarries();

    return all_radd.sumMantissas();
}

#else

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

#endif

