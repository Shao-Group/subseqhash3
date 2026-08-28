#define _POSIX_C_SOURCE 200809L
#include "table_file.h"

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

static void help(FILE *stream) {
    fputs("Generate edit-aware seeds for sequences and sliding windows.\n\n"
          "Usage:\n  subseqhash3 seeds INPUT W E D R OUTPUT [TABLE_FILE]\n\n"
          "Arguments:\n"
          "  INPUT       DNA sequences, one per line.\n"
          "  W           Sliding-window length.\n"
          "  E           Total edit budget.\n"
          "  D           Number of hash residues.\n"
          "  R           Number of table sets and seed repetitions.\n"
          "  OUTPUT      Destination CSV file.\n"
          "  TABLE_FILE  Optional reusable table file matching D and R.\n\n"
          "If TABLE_FILE is omitted, new tables are written to OUTPUT.abc.\n",
          stream);
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

static int write_result(FILE *out,size_t table_set,size_t sequence_index,size_t window_start,
                        const SubseqHash3EditResult *r) {
    static const char alphabet[]="ACGT"; size_t i;
    if (fprintf(out,"%zu,%zu,%zu,%zu,%zu,",
                table_set,sequence_index,window_start,r->substitutions,
                r->deletions)<0) return -1;
    for (i=0;i<r->length;++i)
        if (fputc(alphabet[r->seed[i]],out)==EOF) return -1;
    return fprintf(out,",%zu,%.0f\n",r->psi,r->omega)<0?-1:0;
}

static int process_sequence(FILE *out,const char *line,size_t length,
                            size_t sequence_index,size_t W,size_t E,
                            const SubseqHash3TableFile *table_file) {
    uint8_t *encoded=NULL; SubseqHash3EditResult *results=NULL;
    size_t windows,per,total,w,split,set; int status=-1;
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
    for (set=0;set<table_file->count;++set) {
        size_t index=0;
        if (subseqhash3_compute_long_edit_splits(
                encoded,length,W,E,table_file->tables+set,results)) goto done;
        for (w=0;w<windows;++w)
            for (split=0;split<per;++split,++index)
                if (write_result(out,set,sequence_index,w,results+index)) goto done;
        subseqhash3_edit_results_destroy(results,total);
    }
    status=0;
done:
    if (results) {
        subseqhash3_edit_results_destroy(results,total);
        free(results);
    }
    free(encoded); return status;
}

int ssh3_seeds_main(int argc,char **argv) {
    FILE *input=NULL,*output=NULL; char *line=NULL; size_t capacity=0;
    struct timespec wall_start,wall_end,cpu_start,cpu_end;
    ssize_t got; size_t W,E,D,R,sequence_index=0; int status=EXIT_FAILURE;
    SubseqHash3TableFile tables={0};
    char *generated_path=NULL;
    uint64_t seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) {
        help(stdout); return EXIT_SUCCESS;
    }
    if ((argc!=7 && argc!=8) || parse_size(argv[2],&W) ||
        parse_size(argv[3],&E) || parse_size(argv[4],&D) ||
        parse_size(argv[5],&R) || !W || !D || !R || E>W) {
        help(stderr);
        return EXIT_FAILURE;
    }
    input=fopen(argv[1],"r");
    if (!input) { perror(argv[1]); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&wall_start) ||
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&cpu_start)) {
        perror("clock_gettime"); goto done;
    }
    output=fopen(argv[6],"w");
    if (!output) { perror(argv[6]); goto done; }
    if (argc==7) {
        size_t length=strlen(argv[6]);
        generated_path=malloc(length+5);
        if (!generated_path) goto done;
        memcpy(generated_path,argv[6],length);
        memcpy(generated_path+length,".abc",5);
        if (subseqhash3_table_file_create(generated_path,W,D,R,seed)) goto done;
    }
    if (subseqhash3_table_file_open(&tables,argc==8?argv[7]:generated_path)) goto done;
    if (tables.n<W || tables.d!=D || tables.count!=R) {
        fputs("table file dimensions do not match W, D, and R\n",stderr); goto done;
    }
    if (fputs("table_set,sequence,window_start,substitutions,deletions,"
              "seed,psi,omega\n",output)==EOF) goto done;
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
    if (fflush(output)) { perror(argv[6]); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&wall_end) ||
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&cpu_end)) {
        perror("clock_gettime"); goto done;
    }
    fprintf(stderr,"elapsed wall time: %.6f seconds, CPU time: %.6f seconds\n",
            (double)(wall_end.tv_sec-wall_start.tv_sec)+
            (double)(wall_end.tv_nsec-wall_start.tv_nsec)/1000000000.0,
            (double)(cpu_end.tv_sec-cpu_start.tv_sec)+
            (double)(cpu_end.tv_nsec-cpu_start.tv_nsec)/1000000000.0);
    status=EXIT_SUCCESS;
done:
    free(generated_path); free(line); subseqhash3_table_file_close(&tables);
    if (input && fclose(input)) status=EXIT_FAILURE;
    if (output && fclose(output)) status=EXIT_FAILURE;
    return status;
}
