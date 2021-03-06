#ifndef __IFloat64Common_h
#define __IFloat64Common_h

#include <cstdint>
#include <assert.h>
#include <math.h>
#include <x86intrin.h>

#define DBG_ASSERT(x) assert(x)
//#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

static inline int log2ui(int64_t x)
{
	return __lzcnt64( (x<0) ? -x : x );
}

template<int16_t _EXPMIN=-128, int16_t _EXPMAX=127>
struct IFloat64T
{
    static constexpr int16_t EXPMIN = _EXPMIN;
    static constexpr int16_t EXPMAX = _EXPMAX;
    static constexpr uint16_t EXPRANGESIZE = (EXPMAX-EXPMIN+1);

    int64_t msum[EXPRANGESIZE];
    int32_t mcarry[EXPRANGESIZE];
    uint16_t bmin; // lowest exponent bin used
    uint16_t bmax;

    inline IFloat64T() : bmin(EXPRANGESIZE-1), bmax(0)
    {
        for(int i=0;i<EXPRANGESIZE;i++)
        {
            msum[i]=0;
            mcarry[i]=0;
        }
    }

    inline void reset()
    {
	for(int i=bmin;i<=bmax;i++)
        {
            msum[i]=0;
            mcarry[i]=0;
        }
 	bmin = EXPRANGESIZE-1;
	bmax = 0;
    }

    inline uint16_t nzCarries() const
    {
	uint16_t n = 0;
	for(int i=bmin;i<=bmax;i++)
	{
		if(mcarry[i]!=0) ++n;
	}
	return n;
    }

    inline void addValues(uint64_t n, const double * dx)
    {
       	const int64_t * x = reinterpret_cast<const int64_t*>( dx );
#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	IFloat64T sharedBuf[maxThreads];
#   	pragma omp parallel
    	{
		IFloat64T privBuf;
		const int tid = omp_get_thread_num(); 
#		pragma omp for nowait
       		for(uint64_t i=0;i<n;i++)
		{
		    int64_t s = x[i] >> 63; 
		    int16_t e = ( (x[i]>>52) & ((1ULL<<11)-1) ); 
		    int64_t m = x[i] & ((1ULL<<52)-1ULL);
		    if(e!=0) { m |= (1ULL<<52); e-=1023; } // if not denormalized
		    m = (m^s) - s ; // make signed normalized mantissa (-1)^s * 1,m

		    // TODO: detect and process nan/inf
		    DBG_ASSERT( e>=EXPMIN && e<=EXPMAX );
		    uint16_t ebin = e - EXPMIN;

		    // add signed mantissa
		    m += privBuf.msum[ebin];

		    // update carry, may be negative.
		    privBuf.mcarry[ebin] += m >> 53;

		    // set remaining mantissa
		    privBuf.msum[ebin] = m & ((1ULL<<53)-1) ;

		    // update exponent range
		    if(ebin<privBuf.bmin) { privBuf.bmin=ebin; }
		    if(ebin>privBuf.bmax) { privBuf.bmax=ebin; }
		}
		for(unsigned int i=privBuf.bmin; i<=privBuf.bmax; i++)
		{
			sharedBuf[tid].msum[i] = privBuf.msum[i];
			sharedBuf[tid].mcarry[i] = privBuf.mcarry[i];
		}
		sharedBuf[tid].bmin = privBuf.bmin;
		sharedBuf[tid].bmax = privBuf.bmax;
//#		pragma omp flush
#		pragma omp barrier
		for(unsigned int t=0;t<maxThreads;t++)
		{
#			pragma omp for
			for(unsigned int i=sharedBuf[t].bmin; i<=sharedBuf[t].bmax; i++)
			{
				int64_t m = msum[i] + sharedBuf[t].msum[i];
				mcarry[i] += sharedBuf[t].mcarry[i];
				mcarry[i] += m >> 53;
				msum[i] = m & ((1ULL<<53)-1) ;
			}
#			pragma omp single
			{
				if( sharedBuf[t].bmin < bmin ) bmin = sharedBuf[t].bmin;
				if( sharedBuf[t].bmax > bmax ) bmax = sharedBuf[t].bmax;
			}
		}
	}
#else
       	for(uint64_t i=0;i<n;i++)
	{
	    int64_t s = x[i] >> 63; 
	    int16_t e = ( (x[i]>>52) & ((1ULL<<11)-1) ); 
	    int64_t m = x[i] & ((1ULL<<52)-1ULL);
	    if(e!=0) { m |= (1ULL<<52); e-=1023; } // if not denormalized
	    m = (m^s) - s ; // make signed normalized mantissa (-1)^s * 1,m

	    // TODO: detect and process nan/inf
	    DBG_ASSERT( e>=EXPMIN && e<=EXPMAX );
	    uint16_t ebin = e - EXPMIN;

	    // add signed mantissa
	    m += msum[ebin];

            // update carry, may be negative.
	    mcarry[ebin] += m >> 53;

            // set remaining mantissa
	    msum[ebin] = m & ((1ULL<<53)-1) ;

            // update exponent range
	    if(ebin<bmin) { bmin=ebin; }
	    if(ebin>bmax) { bmax=ebin; }
	}
#endif
    }
   
    template<typename StreamT>
    inline void print(StreamT& os)
    {
	os << "--- bmin="<<bmin<<", bmax="<<bmax<<" ---\n";
        for(int i=bmin;i<=bmax;i++)
	    {
	        int64_t m = msum[i];
	        int32_t c = mcarry[i];
	        uint16_t cre = (c==0) ? 0 : ( log2ui(c) + 1 );
	        int16_t mre = (m==0) ? 0 : ( log2ui(m) - 52 ) ;
		os << "bin "<<i<<" : cre="<<cre<<", mre="<<mre<<", carry="<<c<<" : sum="<<( (double)m/((double)(1ULL<<52)) ) << "\n";
	    }
    }

    inline void distributeCarries()
    {
    	for(int i=bmin;i<bmax;i++)
    	{
    		int64_t c = mcarry[i];
    		if( c!=0 )
    		{
    			mcarry[i] = 0;
	    		uint16_t cre = log2ui(c) + 1;
    			c = c << (53-cre);
    			uint16_t ebin = i+cre; 
    			DBG_ASSERT( ebin < EXPRANGESIZE );
    			int64_t m = msum[ebin] + c;
	            	mcarry[ebin] += m >> 53;
    			msum[ebin] = m & ((1ULL<<53)-1);
    			if(ebin>bmax) bmax = ebin;
    		}
    	}

	// process highest exponent carry
	{
		int i = bmax;
    		int64_t c = mcarry[i];
    		if( c!=0 )
    		{
    			mcarry[i] = 0;
	    		uint16_t cre = log2ui(c) + 1;
    			c = c << (53-cre);
    			uint16_t ebin = i+cre; 
    			DBG_ASSERT( ebin < EXPRANGESIZE );
    			int64_t m = msum[ebin] + c;
    			msum[ebin] += c;
    			if(ebin>bmax) bmax = ebin;
    		}
	}

    }

    inline void distributeMantissas()
    {
	for(int i=bmin;i<=bmax;i++)
	{
		int64_t m = msum[i];
		int mre = 52 - log2ui(m) ;
		if( m!=0 && mre!=0 )
		{
			msum[i] = 0;
			DBG_ASSERT( i >= mre );
			uint16_t ebin = i - mre;
			m = msum[ebin] + (m << mre);
			mcarry[ebin] += m >> 53;
			msum[ebin] = m & ((1ULL<<53)-1);
			if(ebin<bmin) bmin = ebin ;
		}
	}
	while( bmax>bmin && msum[bmax]==0 && mcarry[bmax]==0 ) --bmax;
    }

    inline void removeCarries()
    {
        int nIteration = 0;
	// print(std::cout);
        while( nzCarries() > 0 )
        {
            ++ nIteration;
            distributeCarries();
            distributeMantissas();
        }
    }

    inline double sumMantissas()
    {
	int64_t mantissaSum = 0;
        for(int i=bmin;i<=bmax;i++)
        {
	    mantissaSum = (mantissaSum>>1) + msum[i];
        }
        return exp2(bmax+EXPMIN-52) * ((double)mantissaSum);
    }

    inline bool operator == (const IFloat64T& rhs)
    {
	if( bmin != rhs.bmin ) return false;
	if( bmax != rhs.bmax ) return false;
	for(int i=bmin;i<=bmax;i++)
	{
		if( msum[i]!=rhs.msum[i] || mcarry[i]!=rhs.mcarry[i] ) return false;
	}
	return true;
    }

} __attribute__((aligned(64)));

using IFloat64 = IFloat64T<>;

#endif // __IFloat64Common_h

