#define _POSIX_C_SOURCE 200809L
#include "subseqhash_edits_i8.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC,&t);
    return (double)t.tv_sec+(double)t.tv_nsec*1e-9;
}

static int parse_size(const char *text,size_t *value) {
    char *end;
    unsigned long long parsed;
    errno=0;
    parsed=strtoull(text,&end,10);
    if (errno || !*text || *end || parsed>(unsigned long long)SIZE_MAX)
        return -1;
    *value=(size_t)parsed;
    return 0;
}

int main(int argc,char **argv) {
    size_t n=50,edits=10,d=32,repeats=100,rep,i;
    uint8_t *sequence=NULL;
    SubseqHashEditResult *results=NULL;
    SubseqHashI8Tables tables={0};
    double begin,end;
    int status=EXIT_FAILURE;

    if ((argc>1&&parse_size(argv[1],&n)) ||
        (argc>2&&parse_size(argv[2],&edits)) ||
        (argc>3&&parse_size(argv[3],&d)) ||
        (argc>4&&parse_size(argv[4],&repeats)) ||
        argc>5 || !d || !repeats || edits>n || edits==SIZE_MAX) {
        fprintf(stderr,"usage: %s [n [edits [D [repeats]]]]\n",argv[0]);
        return EXIT_FAILURE;
    }
    sequence=malloc(n?n:1);
    results=calloc(edits+1,sizeof(*results));
    if (!sequence || !results || subseqhash_i8_tables_init(
            &tables,n,d,UINT64_C(0x123456789abcdef0))) {
        fputs("allocation failed\n",stderr);
        goto done;
    }
    for (i=0;i<n;++i)
        sequence[i]=(uint8_t)((i*13+i*i*7+1)%4);

    if (subseqhash_compute_all_edit_splits_i8(
            sequence,n,edits,&tables,results)) goto done;
    subseqhash_edit_results_destroy(results,edits+1);

    begin=seconds();
    for (rep=0;rep<repeats;++rep) {
        if (subseqhash_compute_all_edit_splits_i8(
                sequence,n,edits,&tables,results)) goto done;
        subseqhash_edit_results_destroy(results,edits+1);
    }
    end=seconds();

    printf("n=%zu edits=%zu D=%zu repeats=%zu backend=%s\n",
           n,edits,d,repeats,subseqhash_edits_i8_backend());
    printf("elapsed: %.6f s, %.9f s/call\n",
           end-begin,(end-begin)/(double)repeats);
    status=EXIT_SUCCESS;
done:
    subseqhash_edit_results_destroy(results,edits==SIZE_MAX?0:edits+1);
    subseqhash_i8_tables_destroy(&tables);
    free(results);
    free(sequence);
    return status;
}
