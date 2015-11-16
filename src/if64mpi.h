#ifndef __if64mpi_h
#define __if64mpi_h

#include <mpi.h>

void if64MpiReduceSum(uint64_t n, double * x, MPI_Comm comm);
double if64SumMpiReduceSum(uint64_t n, const double * x, MPI_Comm comm);

#endif
