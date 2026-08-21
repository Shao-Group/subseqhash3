#define _POSIX_C_SOURCE 200809L
#include "subseqhash_edits_i8.h"
#include "subseqhash_edits_max_i8.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double seconds(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t);
    return (double)t.tv_sec+(double)t.tv_nsec*1e-9;
}
static int parse_size(const char *text,size_t *value) {
    char *end; unsigned long long x;
    errno=0; x=strtoull(text,&end,10);
    if (errno || !*text || *end || x>(unsigned long long)SIZE_MAX) return -1;
    *value=(size_t)x; return 0;
}
int main(int argc,char **argv) {
    size_t n=50,edits=10,d=256,repeats=100,rep,i;
    uint8_t *s; double t0,t1,t2;
    SubseqHashI8Tables old_tables;
    SubseqHashMaxI8Tables max_tables;
    if ((argc>1&&parse_size(argv[1],&n)) ||
        (argc>2&&parse_size(argv[2],&edits)) ||
        (argc>3&&parse_size(argv[3],&d)) ||
        (argc>4&&parse_size(argv[4],&repeats)) || argc>5 || !d ||
        !repeats || edits>n || edits>=32 ||
        n>SUBSEQHASH_MAX_I8_MAX_LENGTH) {
        fprintf(stderr,"usage: %s [n [edits [D [repeats]]]]\n",argv[0]);
        return EXIT_FAILURE;
    }
    s=malloc(n?n:1);
    if (!s || subseqhash_i8_tables_init(&old_tables,n,d,
            UINT64_C(0x123456789abcdef0)) ||
        subseqhash_max_i8_tables_init(&max_tables,n,d,
            UINT64_C(0x123456789abcdef0))) {
        fputs("allocation failed\n",stderr); return EXIT_FAILURE;
    }
    for (i=0;i<n;++i) s[i]=(uint8_t)((i*13+i*i*7+1)%4);
    t0=seconds();
    for (rep=0;rep<repeats;++rep) {
        SubseqHashEditResult result[32];
        if (subseqhash_compute_all_edit_splits_i8(
                s,n,edits,&old_tables,result)) return EXIT_FAILURE;
        subseqhash_edit_results_destroy(result,edits+1);
    }
    t1=seconds();
    for (rep=0;rep<repeats;++rep) {
        SubseqHashEditResult result[32];
        if (subseqhash_compute_all_edit_splits_max_i8(
                s,n,edits,&max_tables,result)) return EXIT_FAILURE;
        subseqhash_edit_results_destroy(result,edits+1);
    }
    t2=seconds();
    printf("n=%zu edits=%zu D=%zu repeats=%zu\n",n,edits,d,repeats);
    printf("old backend: %s\n",subseqhash_edits_i8_backend());
    printf("new backend: %s\n",subseqhash_edits_max_i8_backend());
    printf("old i8 min/max: %.6f s, %.6f s/call\n",
           t1-t0,(t1-t0)/repeats);
    printf("new max-only:   %.6f s, %.6f s/call\n",
           t2-t1,(t2-t1)/repeats);
    printf("speedup:        %.2fx\n",(t1-t0)/(t2-t1));
    subseqhash_max_i8_tables_destroy(&max_tables);
    subseqhash_i8_tables_destroy(&old_tables); free(s);
    return EXIT_SUCCESS;
}
