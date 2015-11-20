#ifndef __IFloat64Common_h
#define __IFloat64Common_h

#include <cstdint>
#include <assert.h>
#include <math.h>

//#define DBG_ASSERT(x) assert(x)
#define DBG_ASSERT(x) do{}while(0)
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

template<int16_t _EXPMIN=-256, uint16_t _EXPSLOTS=16>
struct IFloat64T
{
    static constexpr int32_t EXPSLOTS = _EXPSLOTS;
    static constexpr int32_t SLOTBITS = 32;
    static constexpr int32_t EXPMIN = _EXPMIN;
    static constexpr int32_t EXPMAX = _EXPMIN + EXPSLOTS*SLOTBITS;

    int64_t msum[EXPSLOTS];
    int32_t emax;

    inline IFloat64T()
    {
	emax = EXPSLOTS-1;
        for(int32_t i=0;i<EXPSLOTS;i++) { msum[i]=0; }
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
	    int32_t e = ( (X >> 52 ) & ((1ULL<<11)-1) ); 
	    uint64_t m = X & ((1ULL<<52)-1ULL);
	    if(e!=0) { m |= (1ULL<<52); e-=1023; } // if not denormalized
	    uint32_t E = e - EXPMIN;

	    DBG_ASSERT( E >= 0 );

	    uint32_t Ebin = E >> 5;

	    DBG_ASSERT( Ebin < (EXPSLOTS-2) );

	    uint32_t hbc = E & 31U;
	    uint32_t lbc = 32 - hbc;
	    // std::cout<<"E="<<E<<", m="<<m<<", e="<<e<<", s="<<s<<", Ebin="<<Ebin<<", hbc="<<hbc<<", mbc="<<mbc<<", lbc="<<lbc<<"\n";
	
	    int64_t hp = extractBits( m, 64-hbc, hbc );
	    int64_t mp = extractBits( m , lbc, 32 );
	    int64_t lp = extractBits( m , 0, lbc );

	    // std::cout<<"Ebin="<<Ebin<<", m="<<mantissaAsDouble(m)<<"("<<m<<")" <<", lp="<<lp<<", mp="<<mp<<", hp="<<hp<<"\n";
	    msum[Ebin] += (lp^s)-s;
	    msum[Ebin+1] += (mp^s)-s;
	    msum[Ebin+2] += (hp^s)-s;

	    //print(std::cout);
	}
    }

    inline void addValues(uint64_t n, const double * dx)
    {
	static constexpr uint64_t R = 1024*1024*1024;
       	const int64_t * x = reinterpret_cast<const int64_t*>( dx );
	uint64_t rounds = n / R;
	for(uint64_t j=0;j<rounds;j++)
	{
		uint64_t offset = j * R;
		addValuesI64( R , x+offset );
		normalize();
	}
	addValuesI64( n % R , x + (rounds*R) );
	normalize();
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
	emax = EXPSLOTS-1;
	while( emax>0 && msum[emax]==0) --emax;
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
	return true;
    }

} __attribute__((aligned(64)));

using IFloat64 = IFloat64T<>;

#endif // __IFloat64Common_h

