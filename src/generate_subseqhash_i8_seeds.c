#define _POSIX_C_SOURCE 200809L
#include "subseqhash_edits_i8.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int parse_size(const char *text,size_t *value) {
    char *end; unsigned long long parsed;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end || parsed>(unsigned long long)SIZE_MAX)
        return -1;
    *value=(size_t)parsed; return 0;
}

static int encode_sequence(const char *line,size_t length,uint8_t **encoded) {
    uint8_t *s=malloc(length?length:1); size_t i;
    if (!s) return -1;
    for (i=0;i<length;++i) {
        switch (line[i]) {
        case 'A': case 'a': s[i]=0; break;
        case 'C': case 'c': s[i]=1; break;
        case 'G': case 'g': s[i]=2; break;
        case 'T': case 't': s[i]=3; break;
        default: free(s); return -1;
        }
    }
    *encoded=s; return 0;
}

static int write_result(FILE *out,size_t sequence_index,size_t window_start,
                        const SubseqHashEditResult *r) {
    static const char alphabet[]="ACGT"; size_t i;
    if (fprintf(out,"%zu\t%zu\t%zu\t%zu\t%zu\t",
                sequence_index,window_start,r->substitutions,r->deletions,
                r->length)<0) return -1;
    for (i=0;i<r->length;++i)
        if (fputc(alphabet[r->seed[i]],out)==EOF) return -1;
    return fprintf(out,"\t%zu\t%.0f\n",r->psi,r->omega)<0?-1:0;
}

static int process_sequence(FILE *out,const char *line,size_t length,
                            size_t sequence_index,size_t W,size_t E,
                            const SubseqHashI8Tables *tables) {
    uint8_t *encoded=NULL; SubseqHashEditResult *results=NULL;
    size_t windows,per,total,w,split,index=0; int status=-1;
    if (length<W) {
        fprintf(stderr,"sequence %zu has length %zu, shorter than W=%zu\n",
                sequence_index,length,W);
        return -1;
    }
    if (encode_sequence(line,length,&encoded)) {
        fprintf(stderr,"sequence %zu contains a character outside A/C/G/T\n",
                sequence_index);
        return -1;
    }
    windows=length-W+1; per=E+1;
    if (!per || windows>SIZE_MAX/per ||
        windows*per>SIZE_MAX/sizeof(*results)) goto done;
    total=windows*per;
    results=calloc(total,sizeof(*results));
    if (!results) goto done;
    if (subseqhash_compute_long_edits_i8_fast(
            encoded,length,W,E,tables,results)) goto done;
    for (w=0;w<windows;++w)
        for (split=0;split<per;++split,++index)
            if (write_result(out,sequence_index,w,results+index)) goto done;
    status=0;
done:
    if (results) {
        subseqhash_edit_results_destroy(results,total);
        free(results);
    }
    free(encoded); return status;
}

int main(int argc,char **argv) {
    FILE *input=NULL,*output=NULL; char *line=NULL; size_t capacity=0;
    struct timespec start_time,end_time;
    ssize_t got; size_t W,E,D,sequence_index=0; int status=EXIT_FAILURE;
    SubseqHashI8Tables tables;
    memset(&tables,0,sizeof(tables));
    if (argc!=6 || parse_size(argv[2],&W) || parse_size(argv[3],&E) ||
        parse_size(argv[4],&D) || !W || !D || E>W) {
        fprintf(stderr,"usage: %s INPUT W E D OUTPUT\n",argv[0]);
        return EXIT_FAILURE;
    }
    input=fopen(argv[1],"r");
    if (!input) { perror(argv[1]); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&start_time)) {
        perror("clock_gettime"); goto done;
    }
    output=fopen(argv[5],"w");
    if (!output) { perror(argv[5]); goto done; }
    if (subseqhash_i8_tables_init(
            &tables,W,D,UINT64_C(0x7375627365716831))) {
        fputs("failed to initialize SubseqHash i8 tables\n",stderr); goto done;
    }
    if (fputs("sequence\twindow_start\tsubstitutions\tdeletions"
              "\tseed_length\tseed\tpsi\tomega\n",output)==EOF) goto done;
    while ((got=getline(&line,&capacity,input))>=0) {
        size_t length=(size_t)got;
        while (length && (line[length-1]=='\n'||line[length-1]=='\r')) --length;
        if (!length) {
            fprintf(stderr,"sequence %zu is empty\n",sequence_index); goto done;
        }
        if (process_sequence(output,line,length,sequence_index,W,E,&tables))
            goto done;
        ++sequence_index;
    }
    if (ferror(input)) { perror(argv[1]); goto done; }
    if (fflush(output)) { perror(argv[5]); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&end_time)) {
        perror("clock_gettime"); goto done;
    }
    fprintf(stderr,"elapsed time: %.6f seconds\n",
            (double)(end_time.tv_sec-start_time.tv_sec)+
            (double)(end_time.tv_nsec-start_time.tv_nsec)/1000000000.0);
    status=EXIT_SUCCESS;
done:
    free(line); subseqhash_i8_tables_destroy(&tables);
    if (input && fclose(input)) status=EXIT_FAILURE;
    if (output && fclose(output)) status=EXIT_FAILURE;
    return status;
}
