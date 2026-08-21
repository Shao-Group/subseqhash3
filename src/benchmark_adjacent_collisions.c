#define _POSIX_C_SOURCE 200809L
#include "subseqhash_adjacent.h"

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
    char *end; unsigned long long parsed;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end || parsed>(unsigned long long)SIZE_MAX)
        return -1;
    *value=(size_t)parsed; return 0;
}
static void trim(char *line,size_t *length) {
    while (*length && (line[*length-1]=='\n' || line[*length-1]=='\r'))
        line[--*length]='\0';
}
static int encode(const char *text,size_t n,uint8_t *s) {
    size_t i;
    for (i=0;i<n;++i) {
        switch (text[i]) {
            case 'A': s[i]=0; break;
            case 'C': s[i]=1; break;
            case 'G': s[i]=2; break;
            case 'T': s[i]=3; break;
            default: return -1;
        }
    }
    return 0;
}

static int benchmark_file(const char *filename,size_t n,size_t k,size_t edits,
                          const SubseqHashAdjacentTables *tables) {
    FILE *input=fopen(filename,"r");
    char *line1=NULL,*line2=NULL;
    size_t cap1=0,cap2=0,pairs=0,split;
    size_t baseline_collisions=0,any_collisions=0,all_collisions=0;
    size_t split_collisions[32]={0};
    double baseline_time=0,edit_time=0;
    uint8_t *s1=NULL,*s2=NULL,*seed1=NULL,*seed2=NULL;
    int status=-1;
    if (!input) { perror(filename); return -1; }
    s1=malloc(n?n:1); s2=malloc(n?n:1);
    seed1=malloc(k?k:1); seed2=malloc(k?k:1);
    if (!s1 || !s2 || !seed1 || !seed2) goto done;
    for (;;) {
        ssize_t read1=getline(&line1,&cap1,input),read2;
        size_t length1,length2,psi; double omega,t0;
        SubseqHashEditResult left[32],right[32];
        int any=0,all=1;
        if (read1<0) { status=ferror(input)?-1:0; break; }
        read2=getline(&line2,&cap2,input);
        if (read2<0) {
            fprintf(stderr,"%s has an odd line count\n",filename); break;
        }
        length1=(size_t)read1; length2=(size_t)read2;
        trim(line1,&length1); trim(line2,&length2);
        if (length1!=n || length2!=n ||
            encode(line1,n,s1) || encode(line2,n,s2)) {
            fprintf(stderr,"invalid pair %zu in %s\n",pairs+1,filename);
            break;
        }

        t0=seconds();
        if (subseqhash_adjacent_compute_edits(
                s1,n,n-k,0,tables,seed1,NULL,NULL,&psi,&omega) ||
            subseqhash_adjacent_compute_edits(
                s2,n,n-k,0,tables,seed2,NULL,NULL,&psi,&omega))
            break;
        baseline_time+=seconds()-t0;
        baseline_collisions+=(size_t)!memcmp(seed1,seed2,k);

        t0=seconds();
        if (subseqhash_adjacent_compute_all_edit_splits(
                s1,n,edits,tables,left) ||
            subseqhash_adjacent_compute_all_edit_splits(
                s2,n,edits,tables,right)) {
            subseqhash_edit_results_destroy(left,edits+1);
            break;
        }
        for (split=0;split<=edits;++split) {
            const int equal=left[split].length==right[split].length &&
                !memcmp(left[split].seed,right[split].seed,
                        left[split].length);
            any|=equal; all&=equal;
            split_collisions[split]+=(size_t)equal;
        }
        subseqhash_edit_results_destroy(right,edits+1);
        subseqhash_edit_results_destroy(left,edits+1);
        edit_time+=seconds()-t0;
        any_collisions+=(size_t)any;
        all_collisions+=(size_t)all;
        ++pairs;
    }
done:
    if (!status) {
        printf("%s,%zu,%.6f,%.3f,%zu,%.6f,%.3f,%zu,%zu",
               filename,pairs,baseline_time,
               baseline_time?(double)pairs/baseline_time:0.0,
               baseline_collisions,edit_time,
               edit_time?(double)pairs/edit_time:0.0,
               any_collisions,all_collisions);
        for (split=0;split<=edits;++split)
            printf(",%zu",split_collisions[split]);
        putchar('\n');
    }
    free(seed2); free(seed1); free(s2); free(s1);
    free(line2); free(line1); fclose(input);
    return status;
}

int main(int argc,char **argv) {
    size_t n,k,d,edits,i;
    SubseqHashAdjacentTables tables={0};
    int status=0;
    if (argc<6 || parse_size(argv[1],&n) || parse_size(argv[2],&k) ||
        parse_size(argv[3],&d) || parse_size(argv[4],&edits) ||
        !d || k>n || edits>n || edits>=32) {
        fprintf(stderr,"usage: %s n k D edits FILE...\n",argv[0]);
        return EXIT_FAILURE;
    }
    if (subseqhash_adjacent_tables_init(
            &tables,n,d,UINT64_C(0x6a09e667f3bcc909))) {
        fputs("adjacent ABC table allocation failed\n",stderr);
        return EXIT_FAILURE;
    }
    printf("file,pairs,adjacent_k%zu_seconds,adjacent_k%zu_pairs_per_second,"
           "adjacent_k%zu_collisions,adjacent_edits_seconds,"
           "adjacent_edits_pairs_per_second,edit_any_seed_collisions,"
           "edit_all_seeds_collisions",k,k,k);
    for (i=0;i<=edits;++i)
        printf(",edit_Es%zu_Ed%zu_collisions",i,edits-i);
    putchar('\n');
    for (i=5;i<(size_t)argc;++i)
        if (benchmark_file(argv[i],n,k,edits,&tables)) {
            status=1; break;
        }
    subseqhash_adjacent_tables_destroy(&tables);
    return status?EXIT_FAILURE:EXIT_SUCCESS;
}
