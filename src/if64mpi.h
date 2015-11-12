#ifndef __if64mpi_h
#define __if64mpi_h

#include <mpi.h>
double if64AllReduceSum(uint64_t n, const double * x, MPI_Comm comm);

#endif
