#include <mpi.h>
#include <assert.h>
#include "IFloat64.h"

double if64SumMpiReduceSum(uint64_t n, const double * x, MPI_Comm comm)
{
    assert( sizeof(int64_t) == sizeof(long long int) );

    IFloat64 tmp;
    tmp.addValues(n,x);
    tmp.removeCarries();

    MPI_Allreduce( MPI_IN_PLACE, tmp.msum , IFloat64::EXPRANGESIZE , MPI_LONG_LONG_INT, MPI_SUM, comm );
    tmp.bmin=IFloat64::EXPRANGESIZE-1;
    tmp.bmax=0;
    int i;
    for( i=0 ; i<IFloat64::EXPRANGESIZE && tmp.msum[i]==0 ; ++i )
    if( i < IFloat64::EXPRANGESIZE ) tmp.bmin = i;

    for( i=IFloat64::EXPRANGESIZE-1 ; i>=0 && tmp.msum[i]==0 ; --i )
    if( i >=0 ) tmp.bmax = i;

    tmp.removeCarries();
    return tmp.sumMantissas();
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
	IFloat64 tmp;
	double d[nproc];
	for( uint64_t i=localStart; i<localEnd; i++)
	{
		for(int j=0;j<(nproc-1);j++) { d[j] = localValues[ (localEnd-localStart)*j + i-localStart ]; }
		d[nproc-1] = buf[i];
		tmp.reset();
		tmp.addValues(nproc,d);
		tmp.removeCarries();
		buf[i] = tmp.sumMantissas();
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

