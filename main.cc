#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include <algorithm>
#include "f64math.h"
#include "if64math.h"

static inline double wallclock()
{
  struct timeval timer;
  gettimeofday(&timer, NULL);
  double time = timer.tv_sec + timer.tv_usec * 1.0E-6;
  return time;
}

int main(int argc, char* argv[])
{
	int N = atoi(argv[1]);
	long seed = atoi(argv[2]);
	double* x=(double*)malloc(N*sizeof(double));
	if( seed >= 0 )
	{
		srand48(seed);
		for(int i=0;i<N;i++)
		{
			x[i] = drand48() * exp2( static_cast<int>(drand48()*40.0-20.0) );
		}
	}
	else
	{
		for(int i=0;i<N;i++)
		{
			x[i] = seed + i * exp2( (i*40)/N - 20 );
		}		
	}

	printf("N=%d, seed=%lg\n",N,seed);

	double Tref=0.0;
	{
	  double t0= wallclock();
	  double r = f64Sum<double>( N, x );
	  Tref = wallclock() - t0;
	  printf("f64Sum<double>: time=%lg, result=%.20lg\n",Tref,r);
	}

	{
	  double t0= wallclock();
	  std::sort( x, x+N, [](double a, double b) -> bool { return fabs(a)<fabs(b); } );
	  double r = f64Sum<__float128>( N, x );
	  double t1= wallclock();
	  printf("f64Sum<float128>: time=%lg (x%lg), result=%.20lg\n",t1-t0,(t1-t0)/Tref,r);
	}

	//for(int i=0;i<2;i++)
	{
	  double t0= wallclock();
	  double r = rAddF64( N, x );
	  double t1= wallclock();
	  printf("if64Sum: time=%lg (x%lg), result=%.20lg\n",t1-t0,(t1-t0)/Tref,r);
	}

	{
	  double t0= wallclock();
	  std::sort( x, x+N, [](double a, double b) -> bool { return fabs(a)<fabs(b); } );
	  double r = f64Sum<double>( N, x );
	  double t1= wallclock();
	  printf("sorted f64Sum<double>: time=%lg (x%lg), result=%.20lg\n",t1-t0,(t1-t0)/Tref,r);
	}

	return 0;
}
