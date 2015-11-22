#include <mpi.h>
#include <assert.h>
#include "IFloat64.h"

double if64SumMpiReduceSum(uint64_t n, const double * x, MPI_Comm comm)
{
    assert( sizeof(int64_t) == sizeof(long long int) );

    IFloat64 tmp;
    tmp.addValues(n,x);
    int64_t buf[IFloat64::EXPSLOTS+1];
    for(int i=0;i<IFloat64::EXPSLOTS;i++) buf[i]=tmp.msum[i];
    buf[IFloat64::EXPSLOTS] = tmp.flags;

    MPI_Allreduce( MPI_IN_PLACE, buf , IFloat64::EXPSLOTS+1 , MPI_LONG_LONG_INT, MPI_SUM, comm );
    for(int i=0;i<IFloat64::EXPSLOTS;i++) tmp.msum[i]=buf[i];
    tmp.flags = buf[IFloat64::EXPSLOTS];
    tmp.normalize();

    return tmp.toDouble();
}

void if64MpiReduceSum(uint64_t n, double * buf, MPI_Comm comm)
{
    assert( sizeof(int64_t) == sizeof(long long int) );

    int rank=0, nproc=1;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );

    uint64_t localStart = (n*rank) / nproc;
    uint64_t localEnd = (n*(rank+1)) / nproc;

    MPI_Request recvRequest[nproc-1];
    double * localValues = 0;
    if( localEnd > localStart )
    {
    	localValues = new double[ (localEnd-localStart)*(nproc-1) ];
    	for(int j=0;j<(nproc-1);j++)
    	{
		int i = (rank+1+j)%nproc;
    		MPI_Irecv( localValues+(localEnd-localStart)*j, localEnd-localStart, MPI_DOUBLE, i, 0, comm, recvRequest+j );
    	}
    }

    MPI_Request sendRequest[nproc-1];
    int sendRequestCount = 0;
    for(int j=0;j<(nproc-1);j++)
    {
	int i = (rank+1+j)%nproc;
    	uint64_t remoteStart = (n*i) / nproc;
    	uint64_t remoteEnd = (n*(i+1)) / nproc;
	if( remoteEnd > remoteStart )
	{
		MPI_Isend( buf+remoteStart, remoteEnd-remoteStart, MPI_DOUBLE, i, 0, comm, sendRequest+sendRequestCount );
		++ sendRequestCount;
	}
    }

    if( localEnd > localStart )
    {
    	MPI_Waitall(nproc-1, recvRequest, MPI_STATUSES_IGNORE);
	double d[nproc];
	for( uint64_t i=localStart; i<localEnd; i++)
	{
		for(int j=0;j<(nproc-1);j++) { d[j] = localValues[ (localEnd-localStart)*j + i-localStart ]; }
		d[nproc-1] = buf[i];
		IFloat64 tmp;
		tmp.addValues(nproc,d);
		buf[i] = tmp.toDouble();
	}
    }

    MPI_Waitall(sendRequestCount,sendRequest,MPI_STATUSES_IGNORE);

    MPI_Request finalRecvRequest[nproc-1];
    int finalRecvRequestCount = 0;
    for(int j=0;j<(nproc-1);j++)
    {
        int i = (rank+1+j)%nproc;
        uint64_t remoteStart = (n*i) / nproc;
        uint64_t remoteEnd = (n*(i+1)) / nproc;
        if( remoteEnd > remoteStart )
        {
		MPI_Irecv( buf+remoteStart, remoteEnd-remoteStart, MPI_DOUBLE, i, 1, comm, finalRecvRequest+finalRecvRequestCount );
		++ finalRecvRequestCount;
	}
    }

    if( localEnd > localStart )
    {
	MPI_Request finalSendRequest[nproc-1];
   	for(int j=0;j<(nproc-1);j++)
    	{
        	int i = (rank+1+j)%nproc;
		MPI_Isend( buf+localStart, localEnd-localStart, MPI_DOUBLE, i, 1, comm, finalSendRequest+j );
	}
	MPI_Waitall(nproc-1,finalSendRequest,MPI_STATUSES_IGNORE);
    }

    MPI_Waitall(finalRecvRequestCount,finalRecvRequest,MPI_STATUSES_IGNORE);
}

