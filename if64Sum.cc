#include "IFloat64.h"

#ifdef _OPENMP

#include <omp.h>

double if64Sum(uint64_t n, const double * x)
{
    IFloat64 all_radd;

#   pragma omp parallel
    {
        IFloat64 radd;

        const int tid = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        const uint64_t start = (n*tid)/nthreads;
        const uint64_t end = (n*(tid+1))/nthreads;

        radd.addValues( end-start , x+start );
        radd.removeCarries();

        // add thread mantissas to global bins
#       pragma omp critical
        {
            all_radd.addMantissasNoCarryUpdate( radd );
#           pragma omp flush(all_radd)
        }
	    //printExponentBins(bmin,bmax,msum,mcarry);
    }

    // mantissa normalization
    all_radd.computeCarriesFromMantissas();
    all_radd.removeCarries();

	//printf("All after %d iterations: bmin=%d, bmax=%d, range=%d\n",nIteration,all_bmin,all_bmax,all_bmax-all_bmin);
	//printExponentBins(bmin,bmax,msum,mcarry);

    // somme des mantisses
    return all_radd.sumMantissas();
}

#else

double if64Sum(uint64_t n, const double * x)
{
    F64RAdd radd;
    radd.addValues( n, x );
    radd.removeCarries();
    return radd.sumMantissas();
}

#endif

