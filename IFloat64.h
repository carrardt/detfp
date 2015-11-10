W#ifndef __IFloat64Common_h
#define __IFloat64Common_h

#include <cstdint>
#include <assert.h>
#include <math.h>
#include <x86intrin.h>

#define DBG_ASSERT(x) assert(x)
//#define DBG_ASSERT(x) (void)0
#define EXPMIN (-128)
#define EXPMAX (127)
#define EXPRANGESIZE (EXPMAX-EXPMIN+1)

static inline int log2ui(int64_t x)
{
	return __lzcnt64( (x<0) ? -x : x );
}

struct IFloat64
{
    int64_t msum[EXPRANGESIZE];
    int32_t mcarry[EXPRANGESIZE];
    uint16_t bmin; // lowest exponent bin used
    uint16_t bmax;
    uint16_t nzCarries;

    inline IFloat64() : bmin(EXPRANGESIZE-1), bmax(0), nzCarries(0)
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

    inline void addValues(uint64_t n, const double * dx)
    {
       const uint64_t * x = reinterpret_cast<const uint64_t*>( dx );
       for(uint64_t i=0;i<n;i++)
	   {
	        int64_t s = x[i] & (1ULL<<63) ; s = s >> 63; 
	        int16_t e = ( (x[i]>>52) & ((1ULL<<11)-1) ); 
	        int64_t m = x[i] & ((1ULL<<52)-1ULL); 
	        if(e!=0 || m!=0) { m = m | (1ULL<<52) ; e-=1023; } 
	        m = (m^s) + ( s & 1ULL ) ; 

	        DBG_ASSERT( e>=EXPMIN && e<=EXPMAX);
	        uint16_t ebin = e - EXPMIN;

	        m += msum[ebin];

            // update carry
            if( mcarry[ebin] != 0 ) { -- nzCarries; }
	        mcarry[ebin] += m >> 53;
            if( mcarry[ebin] != 0 ) { ++ nzCarries; }

            // set mantissa
	        msum[ebin] = m & ((1ULL<<53)-1);	      

            // update exponent range
	        if(ebin<bmin) { bmin=ebin; }
	        if(ebin>bmax) { bmax=ebin; }
	    }
    }

    inline void addMantissasNoCarryUpdate( const IFloat64& radd )
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
        nzCarries = 0;
        for(int i=bmin;i<=bmax;i++)
        {
            int64_t m = msum[i];
            mcarry[i] = (m >> 53);
            if( mcarry[i] != 0 ) { ++ nzCarries; }
            msum[i] = m & ((1ULL<<53)-1);
        }
    }
    
    template<typename StreamT>
    inline void print(StreamT& os)
    {
        for(int i=bmin;i<=bmax;i++)
	    {
	        int64_t m = msum[i];
	        int32_t c = mcarry[i];
	        uint16_t cre = (c==0) ? 0 : ( log2ui(c) + 1 );
	        int16_t mre = (m==0) ? 0 : ( log2ui(m) - 52 ) ;
		os << "bin "<<i<<" : cre="<<cre<<", mre="<<mre,<<", carry="<<c<<" : sum="<<( (double)m/((double)(1ULL<<52)) );
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
                -- nzCarries; 

	    	  	uint16_t cre = log2ui(c) + 1;
    			c = c << (53-cre);
    			uint16_t ebin = i+cre; 
                // if( ebin >= EXPRANGESIZE ) { printf("i=%d, cre=%d, ebin=%d\n",i,cre,ebin); }
    			DBG_ASSERT( ebin < EXPRANGESIZE );
    			int64_t m = msum[ebin] + c;
                if( mcarry[ebin] != 0 ) { -- nzCarries; }
	            mcarry[ebin] += m >> 53;
                if( mcarry[ebin] != 0 ) { ++ nzCarries; }
    			msum[ebin] = m & ((1ULL<<53)-1);
    			bmax = (ebin>bmax) ? ebin : bmax ;
    		}
    	}
    }

    inline void distributeMantissas()
    {
	    for(int i=bmin;i<=bmax;i++)
	    {
	        int64_t m = msum[i];
	        msum[i] = 0;
	        uint16_t mre = (m==0) ? 0 : ( 52 - log2ui(m) ) ;
	        DBG_ASSERT( i >= mre );
	        uint16_t ebin = i - mre;
	        m = msum[ebin] + (m << mre);
            if( mcarry[ebin] != 0 ) { -- nzCarries; }
	        mcarry[ebin] += m >> 53;
            if( mcarry[ebin] != 0 ) { ++ nzCarries; }
	        msum[ebin] = m & ((1ULL<<53)-1);
	        bmin = (ebin<bmin) ? ebin : bmin ;
	    }
    }

    inline void removeCarries()
    {
        int nIteration = 0;
        while( nzCarries != 0 )
        {
            ++ nIteration;
            distributeCarries();
            distributeMantissas();
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


#endif // __IFloat64Common_h

