#ifndef __IFloat64Common_h
#define __IFloat64Common_h

#include <cstdint>
#include <assert.h>
#include <math.h>
#include <x86intrin.h>

#define DBG_ASSERT(x) assert(x)
//#include <iostream>

static inline int log2ui(int64_t x)
{
	return __lzcnt64( (x<0) ? -x : x );
}

template<int16_t EXPMIN=-128, int16_t EXPMAX=127>
struct IFloat64T
{
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

    inline void zeroMSum()
    {
        for(int i=0;i<EXPRANGESIZE;i++) msum[i]=0;
    }

    inline void zeroMCarry()
    {
        for(int i=0;i<EXPRANGESIZE;i++) mcarry[i]=0;
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

	    // here after, we take care to keep negative mantissas to prevent excessive carry propagation

            // update carry, may be negative.
	    mcarry[ebin] += ( m >> 53 ) - s;

            // set remaining mantissa
	    msum[ebin] = ( m & ((1ULL<<53)-1) ) - ( s & (1LL<<53) );

            // update exponent range
	    if(ebin<bmin) { bmin=ebin; }
	    if(ebin>bmax) { bmax=ebin; }
	}
    }

    inline void addMantissasNoCarryUpdate( const IFloat64T& radd )
    {
        for(uint16_t i=radd.bmin; i<=radd.bmax; i++)
        {
            msum[i] += radd.msum[i];
        }
        if( radd.bmin < bmin ) bmin = radd.bmin;
        if( radd.bmax > bmax ) bmax = radd.bmax;
    }

    // overwrites carry content
    inline void computeCarriesFromMantissas()
    {
        for(int i=bmin;i<=bmax;i++)
        {
            int64_t m = msum[i];
            mcarry[i] = (m >> 53);
            msum[i] = m & ((1ULL<<53)-1);
        }
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
    	for(int i=bmin;i<=bmax;i++)
    	{
    		int64_t c = mcarry[i];
    		if( c != 0 )
    		{
    			mcarry[i] = 0;
	    		uint16_t cre = log2ui(c) + 1;
    			c = c << (53-cre);
    			uint16_t ebin = i+cre; 
                	//if( ebin >= EXPRANGESIZE ) { printf("i=%d, cre=%d, ebin=%d\n",i,cre,ebin); }
    			DBG_ASSERT( ebin < EXPRANGESIZE );
    			int64_t m = msum[ebin] + c;
	            	mcarry[ebin] += m >> 53;
    			msum[ebin] = m & ((1ULL<<53)-1);

			if( mcarry[ebin] < 0 )
			{
				mcarry[ebin] += 1;
				msum[ebin] -= (1LL<<53);
			}
    			bmax = (ebin>bmax) ? ebin : bmax ;
    		}
    	}
    }

    inline void distributeMantissas()
    {
	    for(int i=bmin;i<=bmax;i++)
	    {
	        int64_t m = msum[i];
		if( m != 0 )
		{
	        	msum[i] = 0;
	        	//uint16_t mre = (m==0) ? 0 : ( 52 - log2ui(m) ) ;
			int mre = log2ui(m) ;
			DBG_ASSERT( mre >= 0 );
			mre = 52 - mre;
	        	DBG_ASSERT( i >= mre );
	        	uint16_t ebin = i - mre;
	        	m = msum[ebin] + (m << mre);
	        	mcarry[ebin] += m >> 53;
	        	msum[ebin] = m & ((1ULL<<53)-1);
                	if( mcarry[ebin] < 0 )
                	{
                		mcarry[ebin] += 1;
                        	msum[ebin] -= (1LL<<53);
                	}
	        	bmin = (ebin<bmin) ? ebin : bmin ;
		}
	    }
    }

    inline void removeCarries()
    {
        int nIteration = 0;
	// print(std::cout);
        while( nzCarries() > 0 )
        {
            ++ nIteration;
            distributeCarries();
	    // print(std::cout);
            distributeMantissas();
	    // print(std::cout);
        }
    }

    inline double sumMantissas()
    {
        double r = 0.0;
        for(int i=bmin;i<=bmax;i++)
        {
            r += exp2(i+EXPMIN-52) *  ((double)msum[i]);
        }
        return r;
    }

} __attribute__((aligned(64)));

using IFloat64 = IFloat64T<>;

#endif // __IFloat64Common_h

