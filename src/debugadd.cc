#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <unistd.h>

#include "IFloat64.h"

#include <algorithm>
#include <iostream>

int main(int argc, char* argv[])
{
	int64_t N = atol(argv[1]);

	std::cout<<"initialization\n";
	IFloat64 if64;
	if64.print( std::cout );

	double* x = 0;
	if( N > 0 )
	{
		long seed = atol(argv[2]);
		std::cout<<"N="<<N<<", seed="<<seed<<std::endl;
		x = (double*)malloc(N*sizeof(double));
		srand48(seed);
		for(uint64_t i=0;i<N;i++)
		{
			x[i] = (drand48()-0.5) * exp2( static_cast<int>(drand48()*40.0-20.0) );
		}
		x[ static_cast<uint64_t>(drand48()*(N-1)) ] = 0.0;
		if64.addValues( N , x );
	}
	else
	{
		N = -N;
		x = (double*)malloc(N*sizeof(double));
		for(uint64_t i=0;i<N;i++)
		{
			std::cin >> x[i];
			std::cout<<"add "<<x[i]<<std::endl;
			if64.addValuesI64( 1, reinterpret_cast<int64_t*>(x+i) );
			if64.print( std::cout );
		}
	}

	double r=0.0;
	for(uint64_t i=0;i<N;i++){ r += x[i]; }
	printf("result = %20.20lf / %20.20lf\n",if64.toDouble(),r);
//	std::cout<<"result = "<<if64.sumMantissas()<<" / "<<r<<"\n";

	return 0;
}

