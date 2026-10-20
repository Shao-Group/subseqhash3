#define _POSIX_C_SOURCE 200809L
#include "pivot_table_file.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

static int parse_size(const char *text,size_t *value) {
    char *end;
    unsigned long long parsed;
    if (!*text || *text<'0' || *text>'9') return -1;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || *end || parsed>(unsigned long long)SIZE_MAX) return -1;
    *value=(size_t)parsed;
    return 0;
}

static void help(FILE *stream) {
    fputs("Generate single-threaded pivot seeds: one deletion and E substitutions.\n\n"
          "Usage:\n  subseqhash3 repeats INPUT W E D R OUTPUT [TABLE_FILE]\n\n"
          "Arguments:\n"
          "  INPUT       DNA sequences, one per line (A/C/G/T, either case).\n"
          "  W           Sliding-window length, 2..257.\n"
          "  E           Exactly E substitutions of retained bases, 0..W-1.\n"
          "  D           Number of hash residues.\n"
          "  R           Exactly R seed repetitions per window.\n"
          "  OUTPUT      Destination CSV file.\n"
          "  TABLE_FILE  Optional pivot tables, made with 'tables --pivots'.\n\n"
          "Without TABLE_FILE, new pivot tables are saved to OUTPUT.piv.\n"
          "Each family supplies W-1 pivot orders; more repetitions use additional\n"
          "independent families. Legacy seeds/.abc tables are not compatible.\n"
          "CSV indices (repeat, family, pivot, sequence, window_start) are zero-based.\n",
          stream);
}

static int same_file(const char *a,const char *b) {
    struct stat x,y;
    return !strcmp(a,b) || (!stat(a,&x) && !stat(b,&y) &&
        x.st_dev==y.st_dev && x.st_ino==y.st_ino);
}

static int write_result(FILE *out,size_t repeat,size_t k,size_t sequence,
    size_t start,const SubseqHash3EditResult *result) {
    static const char alphabet[]="ACGT";
    size_t i;
    if (fprintf(out,"%zu,%zu,%zu,%zu,%zu,%zu,%zu,",repeat,repeat/k,repeat%k,
            sequence,start,result->substitutions,result->deletions)<0) return -1;
    for (i=0;i<result->length;++i)
        if (fputc(alphabet[result->seed[i]],out)==EOF) return -1;
    return fprintf(out,",%zu,%.0f\n",result->psi,result->omega)<0?-1:0;
}

static int process_sequence(FILE *out,const char *line,size_t length,size_t index,
    size_t w,size_t e,size_t repeats,const SubseqHash3PivotTableFile *tables,
    SubseqHash3PivotWorkspace *workspace,SubseqHash3EditResult *results) {
    uint8_t *encoded;
    size_t i,start;
    if (length<w) {
        fprintf(stderr,"sequence %zu has length %zu, shorter than W=%zu\n",index,length,w);
        return -1;
    }
    encoded=malloc(length);
    if (!encoded) return -1;
    for (i=0;i<length;++i) {
        switch (line[i]) {
        case 'A': case 'a': encoded[i]=0; break;
        case 'C': case 'c': encoded[i]=1; break;
        case 'G': case 'g': encoded[i]=2; break;
        case 'T': case 't': encoded[i]=3; break;
        default:
            fprintf(stderr,"sequence %zu contains a character outside A/C/G/T\n",index);
            free(encoded); return -1;
        }
    }
    for (start=0;start<=length-w;++start) {
        if (subseqhash3_compute_pivot_repeats_with_workspace(encoded+start,w,e,
                tables->tables,tables->count,repeats,results,workspace)) {
            fprintf(stderr,"pivot seeding failed for sequence %zu, window %zu\n",index,start);
            free(encoded); return -1;
        }
        for (i=0;i<repeats;++i)
            if (write_result(out,i,w-1,index,start,results+i)) {
                perror("writing repeats CSV");
                free(encoded); return -1;
            }
        subseqhash3_edit_results_destroy(results,repeats);
    }
    free(encoded);
    return 0;
}

int ssh3_repeats_main(int argc,char **argv) {
    size_t w,e,d,repeats,needed,index=0,capacity=0;
    char *line=NULL,*generated=NULL;
    const char *table_path;
    ssize_t got;
    FILE *input=NULL,*output=NULL;
    SubseqHash3PivotTableFile tables={0};
    SubseqHash3PivotWorkspace *workspace=NULL;
    SubseqHash3EditResult *results=NULL;
    struct timespec wall_start,wall_end,cpu_start,cpu_end;
    int status=EXIT_FAILURE;
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) {
        help(stdout); return EXIT_SUCCESS;
    }
    if ((argc!=7 && argc!=8) || parse_size(argv[2],&w) ||
        parse_size(argv[3],&e) || parse_size(argv[4],&d) ||
        parse_size(argv[5],&repeats) || w<2 || w>257 || !d || !repeats ||
        e>w-1 || repeats>SIZE_MAX/sizeof(*results)) {
        help(stderr); return EXIT_FAILURE;
    }
    needed=1+(repeats-1)/(w-1);
    if (same_file(argv[1],argv[6]) || (argc==8 && same_file(argv[7],argv[6]))) {
        fputs("INPUT, OUTPUT, and TABLE_FILE must not overwrite each other\n",stderr);
        goto done;
    }
    input=fopen(argv[1],"r");
    if (!input) { perror(argv[1]); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&wall_start) ||
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&cpu_start)) {
        perror("clock_gettime"); goto done;
    }
    if (argc==7) {
        const size_t length=strlen(argv[6]);
        const uint64_t seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
        if (length>SIZE_MAX-5) goto done;
        generated=malloc(length+5);
        if (!generated) goto done;
        memcpy(generated,argv[6],length); memcpy(generated+length,".piv",5);
        if (same_file(argv[1],generated) || same_file(argv[6],generated)) {
            fputs("generated pivot table path aliases INPUT or OUTPUT\n",stderr); goto done;
        }
        if (subseqhash3_pivot_table_file_create(generated,w,d,needed,seed)) goto done;
    }
    table_path=argc==8?argv[7]:generated;
    if (subseqhash3_pivot_table_file_open(&tables,table_path)) goto done;
    if (tables.n<w || tables.d!=d || tables.count<needed) {
        fputs("pivot table dimensions do not match W, D, or required family count\n",stderr);
        goto done;
    }
    results=calloc(repeats,sizeof(*results));
    workspace=subseqhash3_pivot_workspace_create();
    if (!results || !workspace) goto done;
    output=fopen(argv[6],"w");
    if (!output) { perror(argv[6]); goto done; }
    if (fputs("repeat,family,pivot,sequence,window_start,substitutions,deletions,"
              "seed,psi,omega\n",output)==EOF) { perror(argv[6]); goto done; }
    while ((got=getline(&line,&capacity,input))>=0) {
        size_t length=(size_t)got;
        while (length && (line[length-1]=='\n' || line[length-1]=='\r')) --length;
        if (!length) {
            fprintf(stderr,"sequence %zu is empty\n",index); goto done;
        }
        if (process_sequence(output,line,length,index,w,e,repeats,&tables,workspace,results))
            goto done;
        ++index;
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
    if (results) subseqhash3_edit_results_destroy(results,repeats);
    free(results); free(generated); free(line);
    subseqhash3_pivot_workspace_destroy(workspace);
    subseqhash3_pivot_table_file_close(&tables);
    if (input && fclose(input)) status=EXIT_FAILURE;
    if (output && fclose(output)) status=EXIT_FAILURE;
    return status;
}
