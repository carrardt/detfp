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
	
	std::cout<<"N="<<N<<", seed="<<seed<<std::endl;

	// openmp warm-up
	{  double r = f64Sum( N, x ); if(r==0.0) std::cout<<"\n"; }

	double Tref=0.0;
	{
	  double t0= wallclock();
	  double r = f64SumNoOpt( N, x );
	  Tref = wallclock() - t0;
	  std::cout<<"f64SumNoOpt: time="<<Tref<<", result="<<r<<"\n";
	}

	{
	  double t0= wallclock();
	  double r = f64Sum( N, x );
	  double t1= wallclock();
	  std::cout<<"f64Sum: time="<<t1-t0<<"(x"<<(t1-t0)/Tref<<"), result="<<r<<"\n";
	}


	{
	  double t0= wallclock();
	  double r = f64Sumi128( N, x );
	  double t1= wallclock();
	  std::cout<<"f64Sumi128: time="<<t1-t0<<"(x"<<(t1-t0)/Tref<<"), result="<<r<<"\n";
	}

	//for(int i=0;i<2;i++)
	{
	  double t0= wallclock();
	  double r = if64Sum( N, x );
	  double t1= wallclock();
	  std::cout<<"if64Sum: time="<<t1-t0<<"(x"<<(t1-t0)/Tref<<"), result="<<r<<"\n";
	}

	{
	  double t0= wallclock();
	  std::sort( x, x+N, [](double a, double b) -> bool { return fabs(a)<fabs(b); } );
	  double r = f64Sum( N, x );
	  double t1= wallclock();
	  std::cout<<"sorted f64Sum: time="<<t1-t0<<"(x"<<(t1-t0)/Tref<<"), result="<<r<<"\n";
	}

	return 0;
}
