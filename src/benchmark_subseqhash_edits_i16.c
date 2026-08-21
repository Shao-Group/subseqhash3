#define _POSIX_C_SOURCE 200809L
#include "subseqhash_edits_i16.h"

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
    size_t n=50,edits=10,d=31,repeats=100,rep,i,split;
    uint8_t *s; double t0,t1,t2;
    SubseqHashI16Tables compact;
    SubseqHashTables reference;
    if ((argc>1&&parse_size(argv[1],&n)) ||
        (argc>2&&parse_size(argv[2],&edits)) ||
        (argc>3&&parse_size(argv[3],&d)) ||
        (argc>4&&parse_size(argv[4],&repeats)) ||
        argc>5 || !d || !repeats || edits>n || edits>=32) {
        fprintf(stderr,"usage: %s [n [edits [D [repeats]]]]\n",argv[0]);
        return EXIT_FAILURE;
    }
    s=malloc(n?n:1);
    if (!s || subseqhash_i16_tables_init(
            &compact,n,d,UINT64_C(0x123456789abcdef0)) ||
        subseqhash_i16_tables_make_double(&compact,&reference)) {
        fputs("allocation failed\n",stderr); return EXIT_FAILURE;
    }
    for (i=0;i<n;++i) s[i]=(uint8_t)((i*13+i*i*7+1)%4);
    t0=seconds();
    for (rep=0;rep<repeats;++rep) {
        SubseqHashEditResult result[32];
        if (subseqhash_compute_all_edit_splits(
                s,n,edits,&reference,result)) return EXIT_FAILURE;
        subseqhash_edit_results_destroy(result,edits+1);
    }
    t1=seconds();
    for (rep=0;rep<repeats;++rep) {
        SubseqHashEditResult result[32];
        if (subseqhash_compute_all_edit_splits_i16(
                s,n,edits,&compact,result)) return EXIT_FAILURE;
        subseqhash_edit_results_destroy(result,edits+1);
    }
    t2=seconds();
    {
        SubseqHashEditResult a[32],b[32];
        if (subseqhash_compute_all_edit_splits(s,n,edits,&reference,a) ||
            subseqhash_compute_all_edit_splits_i16(s,n,edits,&compact,b))
            return EXIT_FAILURE;
        for (split=0;split<=edits;++split)
            if (a[split].psi!=b[split].psi ||
                a[split].omega!=b[split].omega ||
                memcmp(a[split].seed,b[split].seed,a[split].length)) {
                fputs("result mismatch\n",stderr); return EXIT_FAILURE;
            }
        subseqhash_edit_results_destroy(b,edits+1);
        subseqhash_edit_results_destroy(a,edits+1);
    }
    printf("n=%zu edits=%zu D=%zu repeats=%zu backend=%s\n",
           n,edits,d,repeats,subseqhash_edits_i16_backend());
    printf("double reference: %.6f s, %.6f s/call\n",
           t1-t0,(t1-t0)/repeats);
    printf("i16/i32:          %.6f s, %.6f s/call\n",
           t2-t1,(t2-t1)/repeats);
    printf("speedup:          %.2fx\n",(t1-t0)/(t2-t1));
    subseqhash_tables_destroy(&reference);
    subseqhash_i16_tables_destroy(&compact); free(s);
    return EXIT_SUCCESS;
}
