#define _POSIX_C_SOURCE 200809L
#include "subseqhash.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC,&t);
    return (double)t.tv_sec+(double)t.tv_nsec*1e-9;
}
static int parse_size(const char *text,size_t *value) {
    char *end; unsigned long long x;
    errno=0; x=strtoull(text,&end,10);
    if (errno || !*text || *end || x>(unsigned long long)SIZE_MAX) return -1;
    *value=(size_t)x; return 0;
}
int main(int argc,char **argv) {
    size_t n=50,k=40,d=5120,repetitions=1,i,rep,psi1=0,psi2=0;
    uint8_t *s,*seed1,*seed2,*seed3; double omega1=0,omega2=0,omega3=0,t0,t1,t2,t3;
    size_t psi3=0;
    SubseqHashTables tables;
    if ((argc>1&&parse_size(argv[1],&n)) || (argc>2&&parse_size(argv[2],&k)) ||
        (argc>3&&parse_size(argv[3],&d)) || (argc>4&&parse_size(argv[4],&repetitions)) ||
        argc>5 || !d || !repetitions || k>n) {
        fprintf(stderr,"usage: %s [N [K [D [REPETITIONS]]]]\n",argv[0]);
        return EXIT_FAILURE;
    }
    s=malloc(n?n:1); seed1=malloc(k?k:1); seed2=malloc(k?k:1); seed3=malloc(k?k:1);
    if (!s || !seed1 || !seed2 || !seed3 ||
        subseqhash_tables_init(&tables,k,d,UINT64_C(123456789))) {
        fputs("allocation failed\n",stderr); free(seed3); free(seed2); free(seed1); free(s);
        return EXIT_FAILURE;
    }
    for (i=0;i<n;++i) s[i]=(uint8_t)((i*29+i*i*7+3)%4);
    t0=seconds();
    for (rep=0;rep<repetitions;++rep)
        if (subseqhash_compute(s,n,&tables,seed1,NULL,&psi1,&omega1)) return EXIT_FAILURE;
    t1=seconds();
    for (rep=0;rep<repetitions;++rep)
        if (subseqhash_compute_optimized(s,n,&tables,seed2,NULL,&psi2,&omega2)) return EXIT_FAILURE;
    t2=seconds();
    for (rep=0;rep<repetitions;++rep)
        if (subseqhash_compute_leftmost(s,n,&tables,seed3,NULL,&psi3,&omega3)) return EXIT_FAILURE;
    t3=seconds();
    if (psi1!=psi2 || omega1!=omega2 || memcmp(seed1,seed2,k)) {
        fputs("implementations produced different results\n",stderr); return EXIT_FAILURE;
    }
    if (psi3!=psi1) {
        fputs("leftmost implementation produced a different minimum psi\n",stderr);
        return EXIT_FAILURE;
    }
    printf("n=%zu k=%zu D=%zu repetitions=%zu\n",n,k,d,repetitions);
    printf("original:  %.6f s total, %.6f s/call\n",t1-t0,(t1-t0)/repetitions);
    printf("optimized: %.6f s total, %.6f s/call\n",t2-t1,(t2-t1)/repetitions);
    printf("speedup:   %.2fx\n",(t1-t0)/(t2-t1));
    printf("leftmost:  %.6f s total, %.6f s/call\n",t3-t2,(t3-t2)/repetitions);
    printf("leftmost speedup vs original: %.2fx\n",(t1-t0)/(t3-t2));
    subseqhash_tables_destroy(&tables); free(seed3); free(seed2); free(seed1); free(s);
    return EXIT_SUCCESS;
}
