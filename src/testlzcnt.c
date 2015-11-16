#include <stdio.h>
#include <x86intrin.h>
#include <stdint.h>

static inline int log2ui(int64_t x)
{
        return __lzcnt64( (x<0) ? -x : x );
}

int main(int argc,char* argv[])
{
    int64_t x = atoll(argv[1]);
    printf("%08llX => %d\n",x,log2ui(x));
    return 0;
}

