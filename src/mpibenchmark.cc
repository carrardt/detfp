#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include "f64math.h"
#include "if64math.h"
#include "if64mpi.h"

#include <algorithm>
#include <iostream>

static inline double wallclock()
{
  struct timeval timer;
  gettimeofday(&timer, NULL);
  double time = timer.tv_sec + timer.tv_usec * 1.0E-6;
  return time;
}

static inline void decodeFloat64(const double * px, int64_t& s, int16_t& e, int64_t& m)
{
	uint64_t x = * reinterpret_cast<const uint64_t*>(px);
	s = x & (1ULL<<63) ; s = s >> 63; 
	e = ( (x>>52) & ((1ULL<<11)-1) ); 
	m = x & ((1ULL<<52)-1ULL); 
	if(e!=0 || m!=0) { m = m | (1ULL<<52) ; e-=1023; } 
	m = (m^s) + ( s & 1ULL ) ; 
}

template<typename TestFuncT,typename AssessFuncT>
static inline bool runTest( uint64_t n, double* x, const char* methodName, double &Tref, TestFuncT func, AssessFuncT assess )
{
	double t0 = wallclock();
	double r = func( n, x );
	double t1= wallclock();
	int64_t s=0,m=0; int16_t e=0;
	decodeFloat64(&r,s,e,m);
	double t = t1-t0;
	if(Tref==0.0) { Tref = t; }
	bool resultOk = assess( r );
	printf("%-20s : time=%03.4lf, sum=%20.20lf, sign=%ld, exp=%d, mantissa=%017ld %c \tspeedup=%.2lg\n",methodName,t,r,s,e,m, resultOk?' ':'X', Tref/t );
	return resultOk;
}

int main(int argc, char* argv[])
{
	uint64_t N = atol(argv[1]);
	long seed = atol(argv[2]);
	if( N < 1 ) { return 0; }


	MPI_Init(&argc,&argv);
	int rank = 0;
	int nproc = 1;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &nproc );
	std::cout<<"MPI: P#"<<rank<<std::endl;

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
	{
		std::cout<<"N="<<N<<", seed="<<seed<<", NbProc="<<nproc<<std::endl;
	}

	double* x = (double*)malloc(N*sizeof(double));
	if( seed >= 0 )
	{
		seed += rank;
		srand48(seed);
		for(uint64_t i=0;i<N;i++)
		{
			x[i] = (drand48()-0.5) * exp2( static_cast<int>(drand48()*40.0-20.0) );
		}
		x[ static_cast<uint64_t>(drand48()*(N-1)) ] = 0.0;
	}
	else
	{
		for(uint64_t i=0;i<N;i++)
		{
			double E = ( ((double)i) * 40.0 / (double)N ) - 20.0;
			double I = i;
			double S = seed;
			x[i] = S + I * exp2( E );
		}		
	}
	

	// openmp warm-up
	{  double r = f64Sum( N, x ); if(r==0.0) std::cout<<"\n"; }

	double Tref=0.0, sumNoOpt=0.0, sumOpt=0.0, sumi128=0.0, sumif=0.0;
	//IFloat64 sumData1, sumData2;

	auto mpiSumNoOpt = [] (uint64_t n,const double* x)
		{
			double r = f64SumNoOpt(n,x);
			MPI_Allreduce( MPI_IN_PLACE, &r , 1 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
			return r;
		};

	auto mpiSum = [] (uint64_t n,const double* x)
		{
			double r = f64Sum(n,x);
			MPI_Allreduce( MPI_IN_PLACE, &r , 1 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
			return r;
		};

	auto mpiSumI128 = [] (uint64_t n,const double* x)
		{
			double r = f64Sumi128(n,x);
			MPI_Allreduce( MPI_IN_PLACE, &r , 1 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			return r;
		};


	auto mpiSumIF = [] (uint64_t n,const double* x)
		{
			double r = if64AllReduceSum(n,x,MPI_COMM_WORLD);
			return r;
		};

	runTest(N,x,"SumNoOpt",Tref, mpiSumNoOpt , [&sumNoOpt](double r) {sumNoOpt=r; return true;} );
	runTest(N,x,"Sum",Tref, mpiSum, [&sumOpt](double r) {sumOpt=r; return true;} );
	runTest(N,x,"SumI128",Tref, mpiSumI128 , [&sumi128](double r) {sumi128=r; return true;} );
	//runTest(N,x,"SumIF",Tref, [&sumData1](uint64_t n,const double* x) { return if64Sum(n,x,sumData1); } , [&sumif](double r) {sumif=r; return true;} );
	runTest(N,x,"SumIF",Tref, mpiSumIF, [&sumif](double r) {sumif=r; return true;} );
	runTest(N,x,"sort+Sum",Tref,
		[mpiSum](uint64_t N, double* x) -> double
		{
			std::sort( x, x+N, [](double a, double b) { return fabs(a)<fabs(b); } );
			return mpiSum(N,x);
		}
		, [](double) {return true;}
	);

	printf("---- after sort ----\n");

	runTest(N,x,"SumNoOpt",Tref, mpiSumNoOpt , [sumNoOpt](double r) { return r==sumNoOpt; } );
	runTest(N,x,"Sum",Tref, mpiSum , [sumOpt](double r) { return r==sumOpt; } );
	runTest(N,x,"SumI128",Tref, mpiSumI128, [sumi128](double r) { return r==sumi128; } );
	//bool invresult = runTest(N,x,"SumIF",Tref, [&sumData2](uint64_t n,const double* x) { return if64Sum(n,x,sumData2); }, [sumif](double r) { return r==sumif; } );
	bool invresult = runTest(N,x,"SumIF",Tref, mpiSumIF, [sumif](double r) { return r==sumif; } );
/*
	if( !invresult )
	{
		printSumDataDiff(sumData1, sumData2);
	}
*/
	std::cout<<"\n";

	return invresult ? 0 : 1;
}

