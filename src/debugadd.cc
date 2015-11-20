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
	int64_t N = atoll(argv[1]);
	int64_t mode = atol(argv[2]);
	int64_t flushPeriod = atoll(argv[3]);
	if(flushPeriod<=0) flushPeriod=N;

	std::cout<<"N="<<N<<", mode="<<mode<<", period="<<flushPeriod<<"\n";
	IFloat64 if64;
	//if64.print( std::cout );

	double* x = (double*)malloc(N*sizeof(double));

	if( mode < 0 )
	{
		long seed = -mode-1;
		std::cout<<"seed = "<<seed<<std::endl;
		srand48(seed);
		for(uint64_t i=0;i<N;i++)
		{
			x[i] = (drand48()-0.5) * exp2( static_cast<int>(drand48()*40.0-20.0) );
		}
		x[ static_cast<uint64_t>(drand48()*(N-1)) ] = 0.0;
	}
	else if( mode==0 )
	{
		for(uint64_t i=0;i<N;i++)
		{
			std::cin >> x[i];
		}
	}
	else if( mode==1 )
	{
		int64_t start = atoll(argv[4]);
		int64_t inc = atoll(argv[5]);
		std::cout<<"start = "<<start<<", inc="<<inc<<"\n";
		for(int64_t i=0;i<N;i++)
		{
			x[i] = start;
			start += inc;
		}
	}
	else { return 1; }

	int64_t i=0;
	while( i < N )
	{
		int64_t j = i + flushPeriod;
		if( j>N ) j = N;
		if(flushPeriod==1) std::cout<<x[i]<<"\n";
		if64.addValues( j-i , x+i );
		if64.print( std::cout );
		i += flushPeriod;
	}

	double r=0.0;
	for(uint64_t i=0;i<N;i++){ r += x[i]; }
	printf("if64.toDouble() = %20.20lf , double sum = %20.20lf\n",if64.toDouble(),r);
//	std::cout<<"result = "<<if64.sumMantissas()<<" / "<<r<<"\n";

	return 0;
}

