#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include "IFloat64.h"

#include <algorithm>
#include <iostream>

int main(int argc, char* argv[])
{
	uint64_t N = atol(argv[1]);
	long seed = atol(argv[2]);
	if( N < 1 ) { return 0; }

	double* x = (double*)malloc(N*sizeof(double));
	if( seed >= 0 )
	{
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
	
	std::cout<<"N="<<N<<", seed="<<seed<<std::endl;

	std::cout<<"initialization\n";
	IFloat64 if64;
	if64.print( std::cout );

	std::cout<<"add values\n";
	if64.addValues( N, x );
	if64.print( std::cout );

	int nIterations = 0;
        while( if64.nzCarries() > 0 )
        {
	    ++ nIterations;
	    std::cout<<"non-zero carries = "<< if64.nzCarries()<<"\n";
	    std::cout<<"distributeCarries\n";	
            if64.distributeCarries();
	    if64.print( std::cout );
	    std::cout<<"distributeMantissas\n";	
            if64.distributeMantissas();
	    if64.print( std::cout );
        }
	std::cout<<"iterations = "<<nIterations<<"\n";
	double r=0.0;
	for(uint64_t i=0;i<N;i++){ r += x[i]; }
	std::cout<<"result = "<<if64.sumMantissas()<<" / "<<r<<"\n";

	return 0;
}

