#include <stdio.h>
#include <x86intrin.h>

int main(int argc,char* argv[])
{
    long x = atoi(argv[1]);
    printf("%d\n",__lzcnt64(x));
    return 0;
}

