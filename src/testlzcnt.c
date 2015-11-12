#include <stdio.h>
#include <x86intrin.h>
#include <stdint.h>

static inline int log2ui(int64_t x)
{
        return __lzcnt64( (x<0) ? -x : x );
}

int main(int argc,char* argv[])
{
    long x = atoi(argv[1]);
    printf("%d\n",log2ui(x));
    return 0;
}

