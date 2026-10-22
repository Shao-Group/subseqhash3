#define _POSIX_C_SOURCE 200809L
#include "ssh3_deletions.h"
#include "pivot_table_file.h"
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

static void help(FILE *f) {
    fputs("Generate up to k=W-E deletion-only pivot seeds per window (outward by default).\n\n"
          "Usage:\n  subseqhash3 deletions INPUT W E D R OUTPUT [TABLE_FILE] [--seed SEED]\n"
          "                       [--subseqhash2-tables TEXT_FILE] [--inward]\n\n"
          "  --inward    Use shared inward DPs (different hash order and seeds).\n"
          "  INPUT       DNA sequences, one per line, A/C/G/T in either case.\n"
          "  W           Window length, 1..256; windows advance one base.\n"
          "  E           Exactly E deletions, 0..W-1; no substitutions.\n"
          "  D           Number of residues, positive.\n"
          "  R           Requested seeds, 1..W-E, using pivots 0..R-1.\n"
          "  OUTPUT      CSV destination, with exact integer scores.\n"
          "  TABLE_FILE  Shared binary pivot tables, made by tables --pivots.\n"
          "              Requires n-1>=W-E and matching D; uses the first family.\n"
          "  --seed SEED Reproducible table generation when TABLE_FILE is omitted.\n\n"
          "  --subseqhash2-tables TEXT_FILE  Import original SubseqHash2 text tables.\n"
          "              Use instead of TABLE_FILE or --seed; k=W-E, matching D.\n\n"
          "Without TABLE_FILE, shared tables are created in OUTPUT.piv and reused.\n"
          "The generation seed is printed to stderr. This command is single-threaded.\n"
          "Legacy .abc tables are not compatible. Text-table permutations are validated;\n"
          "output pivots always use numeric order.\n",f);
}
static int number(const char *text,uint64_t *value) {
    char *end; unsigned long long x;
    if (*text<'0' || *text>'9') return -1;
    errno=0; x=strtoull(text,&end,10);
    if (errno || *end) return -1;
    *value=(uint64_t)x; return 0;
}
static int alias(const char *a,const char *b) {
    struct stat x,y;
    return !strcmp(a,b) || (!stat(a,&x) && !stat(b,&y) && x.st_dev==y.st_dev && x.st_ino==y.st_ino);
}
typedef struct {
    FILE *out;
    size_t sequence;
} DeletionsCSV;
static int write_window(size_t start,const SubseqHash3DeletionResult *results,
    size_t count,void *context) {
    DeletionsCSV *csv=context;
    for (size_t a=0;a<count;++a) {
        const SubseqHash3DeletionResult *r=results+a;
        if (fprintf(csv->out,"%zu,%zu,%zu,%zu,",a,csv->sequence,start,r->deletions)<0) return -1;
        for (size_t j=0;j<r->length;++j)
            if (fputc("ACGT"[r->seed[j]],csv->out)==EOF) return -1;
        if (fprintf(csv->out,",%zu,%" PRId64 "\n",r->psi,r->omega)<0) return -1;
    }
    return 0;
}
int ssh3_deletions_main(int argc,char **argv) {
    const char *args[7],*original=NULL;
    uint64_t values[4],seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    size_t w,e,d,r,seq=0,capacity=0,j;
    int i,count=0,seed_set=0,inward=0,status=EXIT_FAILURE;
    FILE *in=NULL,*out=NULL;
    char *line=NULL,*generated=NULL;
    uint8_t *encoded=NULL;
    ssize_t got;
    SubseqHash3DeletionTables tables={0};
    SubseqHash3PivotTableFile shared={0};
    struct timespec before,after;
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) { help(stdout); return 0; }
    for (i=1;i<argc;++i) {
        if (!strcmp(argv[i],"--inward")) {
            if (inward) goto invalid;
            inward=1;
        } else if (!strcmp(argv[i],"--seed")) {
            if (seed_set || i+1>=argc || number(argv[++i],&seed)) goto invalid;
            seed_set=1;
        } else if (!strcmp(argv[i],"--subseqhash2-tables")) {
            if (original || i+1>=argc) goto invalid;
            original=argv[++i];
        } else {
            if (count==7) goto invalid;
            args[count++]=argv[i];
        }
    }
    if (count<6 || (count==7 && (seed_set || original)) || (original && seed_set)) goto invalid;
    for (i=0;i<4;++i) if (number(args[i+1],values+i) || values[i]>SIZE_MAX) goto invalid;
    w=(size_t)values[0]; e=(size_t)values[1]; d=(size_t)values[2]; r=(size_t)values[3];
    if (!w || w>256 || e>=w || !d || !r || r>w-e) goto invalid;
    if (alias(args[0],args[5]) || (count==7 && alias(args[6],args[5])) ||
        (original && alias(original,args[5]))) {
        fputs("OUTPUT must not alias INPUT or TABLE_FILE\n",stderr); goto done;
    }
    in=fopen(args[0],"r"); if (!in) { perror(args[0]); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&before)) { perror("clock_gettime"); goto done; }
    if (original) {
        if (subseqhash3_deletion_tables_import(&tables,original,w-e,d)) {
            fputs("invalid SubseqHash2 text table or unsupported dimensions/magnitudes\n",stderr); goto done;
        }
    } else {
        const char *path=count==7?args[6]:NULL;
        if (!path) {
            const size_t length=strlen(args[5]);
            /* Match repeats' window capacity, with one extra row for E=0. */
            const size_t maximum=w+(e==0);
            if (length>SIZE_MAX-5) goto done;
            generated=malloc(length+5); if (!generated) goto done;
            memcpy(generated,args[5],length); memcpy(generated+length,".piv",5);
            if (alias(args[0],generated) || alias(args[5],generated)) {
                fputs("generated pivot table path aliases INPUT or OUTPUT\n",stderr); goto done;
            }
            if (subseqhash3_pivot_table_file_create(generated,maximum,d,1,seed)) goto done;
            path=generated;
            fprintf(stderr,"pivot table seed: %" PRIu64 "\n",seed);
        }
        /* Creation can make a previously dangling OUTPUT symlink resolve. */
        if (alias(args[5],path)) {
            fputs("OUTPUT aliases the active pivot table file\n",stderr); goto done;
        }
        if (subseqhash3_pivot_table_file_open(&shared,path)) {
            fputs("expected shared tables from 'tables --pivots'; use --subseqhash2-tables for original text tables\n",stderr);
            goto done;
        }
        if (shared.n-1<w-e || shared.d!=d) {
            fputs("pivot table dimensions must satisfy n-1>=W-E and matching D\n",stderr); goto done;
        }
        if (subseqhash3_deletion_tables_from_pivots(&tables,shared.tables)) goto done;
        subseqhash3_pivot_table_file_close(&shared);
    }
    out=fopen(args[5],"w"); if (!out) { perror(args[5]); goto done; }
    if (fputs("pivot,sequence,window_start,deletions,seed,psi,omega\n",out)==EOF) goto done;
    while ((got=getline(&line,&capacity,in))>=0) {
        size_t length=(size_t)got;
        while (length && (line[length-1]=='\n' || line[length-1]=='\r')) --length;
        if (length<w) { fprintf(stderr,"sequence %zu is empty or shorter than W\n",seq); goto done; }
        encoded=malloc(length); if (!encoded) goto done;
        for (j=0;j<length;++j) {
            switch (line[j]) {
            case 'A': case 'a': encoded[j]=0; break;
            case 'C': case 'c': encoded[j]=1; break;
            case 'G': case 'g': encoded[j]=2; break;
            case 'T': case 't': encoded[j]=3; break;
            default: fprintf(stderr,"sequence %zu contains a character outside A/C/G/T\n",seq); goto done;
            }
        }
        DeletionsCSV csv={out,seq};
        if (subseqhash3_seed_deletion_windows(encoded,length,w,e,&tables,r,
                inward?SUBSEQHASH3_INWARD:SUBSEQHASH3_OUTWARD,write_window,&csv)) goto done;
        free(encoded); encoded=NULL; ++seq;
    }
    if (ferror(in) || fflush(out)) { fputs("deletion input/output error\n",stderr); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&after)) goto done;
    fprintf(stderr,"elapsed wall time: %.6f seconds\n",(double)(after.tv_sec-before.tv_sec)+(double)(after.tv_nsec-before.tv_nsec)/1e9);
    status=EXIT_SUCCESS;
done:
    if (in && fclose(in)) status=EXIT_FAILURE;
    if (out && fclose(out)) status=EXIT_FAILURE;
    subseqhash3_deletion_tables_destroy(&tables);
    subseqhash3_pivot_table_file_close(&shared);
    free(encoded); free(line); free(generated); return status;
invalid:
    help(stderr); return EXIT_FAILURE;
}
