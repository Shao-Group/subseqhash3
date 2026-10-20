#define _POSIX_C_SOURCE 200809L
#include "ssh3_deletions.h"
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

static void help(FILE *f) {
    fputs("Generate up to k=W-E deletion-only SubseqHash2 pivot seeds per window.\n\n"
          "Usage:\n  subseqhash3 deletions INPUT W E D R OUTPUT [TABLE_FILE] [--seed SEED]\n\n"
          "  INPUT       DNA sequences, one per line, A/C/G/T in either case.\n"
          "  W           Window length, 1..256; windows advance one base.\n"
          "  E           Exactly E deletions, 0..W-1; no substitutions.\n"
          "  D           Number of residues, positive.\n"
          "  R           Requested seeds, 1..W-E, using pivots 0..R-1.\n"
          "  OUTPUT      CSV destination, with exact integer scores.\n"
          "  TABLE_FILE  Original SubseqHash2 text table for k=W-E and D.\n"
          "              Its final permutation is validated; pivots use numeric order.\n"
          "  --seed SEED Reproducible table generation when TABLE_FILE is omitted.\n\n"
          "Default tables are randomized once and reused across every input window.\n"
          "The generation seed is printed to stderr. This command is single-threaded.\n"
          "Pivot/legacy .piv/.abc binary files are not accepted here.\n",f);
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
int ssh3_deletions_main(int argc,char **argv) {
    const char *args[7];
    uint64_t values[4],seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    size_t w,e,d,r,seq=0,capacity=0,start,j,a;
    int i,count=0,seed_set=0,status=EXIT_FAILURE;
    FILE *in=NULL,*out=NULL;
    char *line=NULL;
    uint8_t *encoded=NULL;
    ssize_t got;
    SubseqHash3DeletionTables tables={0};
    SubseqHash3DeletionResult *results=NULL;
    SubseqHash3DeletionWorkspace *workspace=NULL;
    struct timespec before,after;
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) { help(stdout); return 0; }
    for (i=1;i<argc;++i) {
        if (!strcmp(argv[i],"--seed")) {
            if (seed_set || i+1>=argc || number(argv[++i],&seed)) goto invalid;
            seed_set=1;
        } else {
            if (count==7) goto invalid;
            args[count++]=argv[i];
        }
    }
    if (count<6 || (count==7 && seed_set)) goto invalid;
    for (i=0;i<4;++i) if (number(args[i+1],values+i) || values[i]>SIZE_MAX) goto invalid;
    w=(size_t)values[0]; e=(size_t)values[1]; d=(size_t)values[2]; r=(size_t)values[3];
    if (!w || w>256 || e>=w || !d || !r || r>w-e) goto invalid;
    if (alias(args[0],args[5]) || (count==7 && alias(args[6],args[5]))) {
        fputs("OUTPUT must not alias INPUT or TABLE_FILE\n",stderr); goto done;
    }
    in=fopen(args[0],"r"); if (!in) { perror(args[0]); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&before)) { perror("clock_gettime"); goto done; }
    if (count==7) {
        if (subseqhash3_deletion_tables_import(&tables,args[6],w-e,d)) {
            fputs("invalid SubseqHash2 text table or unsupported dimensions/magnitudes\n",stderr); goto done;
        }
    } else {
        if (subseqhash3_deletion_tables_init(&tables,w-e,d,seed)) goto done;
        fprintf(stderr,"deletion table seed: %" PRIu64 "\n",seed);
    }
    workspace=subseqhash3_deletion_workspace_create(); results=calloc(r,sizeof(*results));
    if (!workspace || !results) goto done;
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
        for (start=0;start<=length-w;++start) {
            if (subseqhash3_compute_deletion_seeds_with_workspace(encoded+start,w,e,&tables,NULL,r,results,workspace)) goto done;
            for (a=0;a<r;++a) {
                if (fprintf(out,"%zu,%zu,%zu,%zu,",a,seq,start,e)<0) goto done;
                for (j=0;j<w-e;++j) if (fputc("ACGT"[results[a].seed[j]],out)==EOF) goto done;
                if (fprintf(out,",%zu,%" PRId64 "\n",results[a].psi,results[a].omega)<0) goto done;
            }
            subseqhash3_deletion_results_destroy(results,r);
        }
        free(encoded); encoded=NULL; ++seq;
    }
    if (ferror(in) || fflush(out)) { fputs("deletion input/output error\n",stderr); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&after)) goto done;
    fprintf(stderr,"elapsed wall time: %.6f seconds\n",(double)(after.tv_sec-before.tv_sec)+(double)(after.tv_nsec-before.tv_nsec)/1e9);
    status=EXIT_SUCCESS;
done:
    if (in && fclose(in)) status=EXIT_FAILURE;
    if (out && fclose(out)) status=EXIT_FAILURE;
    subseqhash3_deletion_results_destroy(results,results?r:0); free(results);
    subseqhash3_deletion_workspace_destroy(workspace);
    subseqhash3_deletion_tables_destroy(&tables);
    free(encoded); free(line); return status;
invalid:
    help(stderr); return EXIT_FAILURE;
}
