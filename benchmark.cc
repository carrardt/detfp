#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include "f64math.h"
#include "if64math.h"

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
static inline void runTest( uint64_t n, double* x, const char* methodName, double &Tref, TestFuncT func, AssessFuncT assess )
{
	double t0 = wallclock();
	double r = func( n, x );
	double t1= wallclock();
	int64_t s=0,m=0; int16_t e=0;
	decodeFloat64(&r,s,e,m);
	double t = t1-t0;
	if(Tref==0.0) { Tref = t; }
	bool resultOk = assess( r );
	printf("%-20s : time=%03.3lf, sum=%20.20lf, sign=%ld, exp=%d, mantissa=%017ld %c \tspeedup=%.2lg\n",methodName,t,r,s,e,m, resultOk?' ':'X', Tref/t );
}

int main(int argc, char* argv[])
{
	long N = atol(argv[1]);
	long seed = atol(argv[2]);
	double* x = (double*)malloc(N*sizeof(double));
	if( seed >= 0 )
	{
		srand48(seed);
		for(uint64_t i=0;i<N;i++)
		{
			x[i] = drand48() * exp2( static_cast<int>(drand48()*40.0-20.0) );
		}
	}
	else
	{
		for(uint64_t i=0;i<N;i++)
		{
			x[i] = seed + i * exp2( (i*40)/N - 20 );
		}		
	}
	
	std::cout<<"N="<<N<<", seed="<<seed<<std::endl;

	// openmp warm-up
	{  double r = f64Sum( N, x ); if(r==0.0) std::cout<<"\n"; }

	double Tref=0.0, sumNoOpt=0.0, sumOpt=0.0, sumi128=0.0, sumif=0.0;
	runTest(N,x,"SumNoOpt",Tref,f64SumNoOpt, [&sumNoOpt](double r)->bool {sumNoOpt=r; return true;} );
	runTest(N,x,"Sum",Tref,f64Sum, [&sumOpt](double r)->bool {sumOpt=r; return true;} );
	runTest(N,x,"SumI128",Tref,f64Sumi128, [&sumi128](double r)->bool {sumi128=r; return true;} );
	runTest(N,x,"SumIF",Tref,if64Sum, [&sumif](double r)->bool {sumif=r; return true;} );
	runTest(N,x,"sort+Sum",Tref,
		[](uint64_t N, double* x) -> double
		{
			std::sort( x, x+N, [](double a, double b) -> bool { return fabs(a)<fabs(b); } );
			return f64Sum(N,x);
		}
		, [](double)->bool{return true;}
	);

	printf("---- after sort ----\n");

	runTest(N,x,"SumNoOpt",Tref,f64SumNoOpt, [sumNoOpt](double r)->bool { return r==sumNoOpt; } );
	runTest(N,x,"Sum",Tref,f64Sum, [sumOpt](double r)->bool { return r==sumOpt; } );
	runTest(N,x,"SumI128",Tref,f64Sumi128, [sumi128](double r)->bool { return r==sumi128; } );
	runTest(N,x,"if64Sum",Tref,if64Sum, [sumif](double r)->bool { return r==sumif; } );

	return 0;
}
