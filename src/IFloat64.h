#ifndef __IFloat64Common_h
#define __IFloat64Common_h

#include <cstdint>
#include <assert.h>
#include <math.h>
#include <limits>

//#define DBG_ASSERT(x) assert(x)
#define DBG_ASSERT(x) if(!(x))__builtin_unreachable()
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <immintrin.h>

template<int32_t _EXPMIN, int32_t _EXPSLOTS>
struct IFloat64T
{
    static constexpr int32_t EXPSLOTS = _EXPSLOTS;
    static constexpr int32_t SLOTBITS = 32;
    static constexpr int32_t EXPMIN = _EXPMIN;
    static constexpr int32_t EXPMAX = _EXPMIN + EXPSLOTS*SLOTBITS;

    int64_t msum[EXPSLOTS] __attribute__((aligned(64))) ;
    int64_t msum1[EXPSLOTS-1] __attribute__((aligned(64))) ; // msum shifted by one
    int64_t msum2[EXPSLOTS-2] __attribute__((aligned(64))) ; // msum shifted by two
    int32_t emax;
    uint32_t flags;

    inline IFloat64T()
    {
	flags = 0;
	emax = EXPSLOTS-1;
        for(int32_t i=0;i<EXPSLOTS;i++) { msum[i]=0; }
        for(int32_t i=0;i<EXPSLOTS-1;i++) { msum1[i]=0; }
        for(int32_t i=0;i<EXPSLOTS-2;i++) { msum2[i]=0; }
    }

    inline bool isNormalized() const
    {
	int64_t mask = 0;
	for(int32_t i=0;i<(EXPSLOTS-1);i++)
	{
		mask |= msum[i];
	}
	return ( (mask>>32) != 0 );
    }

    static inline uint64_t extractBits(uint64_t x, unsigned int b, unsigned int n)
    {
	
	unsigned int endBit = b+n;
	if( n==0 ) return 0;
	if( n > 64 ) return 0;
	if( (b+n) > 64 ) return 0;
	return ( x << (64-endBit) ) >> (64-n);
    }

    static inline double mantissaAsDouble(int64_t m)
    {
	return exp2(-52) * ((double)m);
    }

    inline void addValuesI64(uint64_t n, const int64_t * x)
    {
      	for(uint64_t i=0;i<n;i++)
   	{
	    int64_t X = x[i];
	    int64_t s = X >> 63; 
	    int32_t e = extractBits( X , 52 , 11 ); 
	    uint64_t m = extractBits( X , 0 , 52 );

	    // detect NaN and infity
	    if(e==0x7FF)
	    {
		if(m) flags |= s ? 4 : 8 ; // +/- NaN
		else flags |= s ? 1 : 2; // +/- Inf
 	    }

	    // denormalized case
	    if(e!=0) { m += (1ULL<<52); e-=1023; }

	    // round to zero
	    if(e<EXPMIN) { e=0; m=0; } 

	    uint32_t E = e - EXPMIN;
	    uint32_t Ebin = E / 32;

	    // round to infinity (exponent overflow)
	    if( Ebin >= (EXPSLOTS-2)  )
	    {
		flags |= s ? 1 : 2;
		m = 0;
		Ebin = 0;
	    }

	    //DBG_ASSERT( Ebin < (EXPSLOTS-2) );
	    unsigned int hbc = E % 32;
	    unsigned int lbc = 32 - hbc;
	    // std::cout<<"E="<<E<<", m="<<m<<", e="<<e<<", s="<<s<<", Ebin="<<Ebin<<", hbc="<<hbc<<", mbc="<<mbc<<", lbc="<<lbc<<"\n";
	
	    int64_t hp = ( extractBits( m, 64-hbc, hbc ) ^ s ) - s;
	    int64_t mp = ( extractBits( m , lbc, 32 ) ^ s ) - s;
	    int64_t lp = ( extractBits( m , 0, lbc ) ^ s ) - s;
	    
	    /* SSE2 exemple to add lp and mp
	    __m128i a = _mm_set_epi64x( lp , mp );
	    __m128i b = _mm_loadu_si128( (__m128i*) (msum+Ebin) );
	    __m128i c = _mm_add_epi64( a , b );
	    _mm_storeu_si128( (__m128i*) (msum+Ebin) , c );
	    */

	    // std::cout<<"Ebin="<<Ebin<<", m="<<mantissaAsDouble(m)<<"("<<m<<")" <<", lp="<<lp<<", mp="<<mp<<", hp="<<hp<<"\n";
	    msum[Ebin] += lp;
	    msum1[Ebin] += mp;
	    msum2[Ebin] += hp;
	}
    }

    inline void addValuesSeq(uint64_t n, const double * dx)
    {
	// we want to avoid carry overflow,
	// so no more than R values are added before a normalization
       	const int64_t * x = reinterpret_cast<const int64_t*>( dx );
	uint64_t rounds = n >> 30;
	for(uint64_t j=0;j<rounds;j++)
	{
		uint64_t offset = j << 30;
		addValuesI64( 1ULL<<30 , x+offset );
		normalize();
	}
	addValuesI64( n & ((1ULL<<30)-1) , x + (rounds<<30) );
	normalize();
    }

    inline void addValues(uint64_t n, const double * dx)
    {
#ifdef _OPENMP
	const uint64_t maxThreads = omp_get_max_threads();
	IFloat64T sharedBuf[maxThreads];
#   	pragma omp parallel
    	{
		const uint64_t tid = omp_get_thread_num();
		uint64_t start = ( n * tid ) / maxThreads;
		uint64_t end = ( n * (tid+1) ) / maxThreads;
		sharedBuf[tid].addValuesSeq( end - start , dx + start );
	}
	for(uint64_t t=0;t<maxThreads;t++)
	{
		for(int32_t i=0;i<EXPSLOTS;i++)
		{
			msum[i] += sharedBuf[t].msum[i];
		}
		flags |= sharedBuf[t].flags;
	}
	normalize();
#else
	addValuesSeq(n,dx);
#endif
    }

    template<typename StreamT>
    inline void print(StreamT& os) const
    {
	int32_t s = 0;
	while( s<EXPSLOTS && msum[s]==0 ) ++s;
	int32_t e = EXPSLOTS-1;
	while( e>=0 && msum[e]==0 ) --e;
	os << "N="<<isNormalized()<<" S="<<s<<" E="<<e<<" :";
        for(int32_t i=s;i<=e;i++)
	    {
		os << " " << (void*)(msum[i]) ;
//		os << " " << exp2(-32)*((double)msum[i]);
	    }
	os <<" = "<<toDouble()<< "\n";
    }

    inline void normalize()
    {
	// flatten shifted arrays
	msum[1] += msum1[0];
	for(uint32_t i=2;i<EXPSLOTS;i++)
	{
		msum[i] += msum1[i-1] + msum2[i-2];
	}
	for(uint32_t i=0;i<EXPSLOTS-1;i++) { msum1[i]=0; }
	for(uint32_t i=0;i<EXPSLOTS-2;i++) { msum2[i]=0; }

	// search for max exponent, to avoid pushing negative sign to the highest exponent (precision concern)
	emax = EXPSLOTS-1;
	while( emax>0 && msum[emax]==0) --emax;

	// propagate carries up to emax
    	int64_t carry = 0;
	for(uint32_t i=0;i<emax;i++)
	{
		int64_t m = msum[i];
		m += carry;
		carry = m >> 32;
		m = m & ((1LL<<32)-1LL);
		msum[i] = m;
	}
	msum[emax] += carry;
    }

    inline double toDouble() const
    {
	if( flags != 0 )
	{
		if( flags==1 ) return - std::numeric_limits<double>::infinity();
		if( flags==2 ) return std::numeric_limits<double>::infinity();
		if( flags==4 ) return - std::numeric_limits<double>::quiet_NaN();
		return std::numeric_limits<double>::quiet_NaN();
	}
	double Sum = 0.0;
        for(int32_t i=0; i<=emax; ++i)
	{
	    Sum += exp2(i*32+EXPMIN-52) * ((double)msum[i]);
	}
        return Sum;
    }

    inline bool operator == (const IFloat64T& rhs) const
    {
	for(uint32_t i=0;i<EXPSLOTS;i++)
	{
		if( msum[i]!=rhs.msum[i] ) return false;
	}
	if( flags!=rhs.flags ) return false;
	return true;
    }

} __attribute__((aligned(64)));

//using IFloat64 = IFloat64T<-256,16>;
//using IFloat64 = IFloat64T<-512,32>;
using IFloat64 = IFloat64T<-1024,64>; // Full DP precision

#endif // __IFloat64Common_h

