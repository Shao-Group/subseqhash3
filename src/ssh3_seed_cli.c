#define _POSIX_C_SOURCE 200809L
#include "ssh3_commands.h"
#include "ssh3_seeding.h"
#include "ssh3_cli.h"
#include <inttypes.h>
#include <stdio.h>
#include <time.h>

static void help(FILE *out,int ind) {
    fprintf(out,"Generate pivot seeds for %s-heavy sequences (outward by default).\n\n"
        "Usage:\n  subseqhash3 seeds-%s INPUT W E D R S OUTPUT [TABLE_FILE] [--inward] [--seed SEED]\n\n"
        "  INPUT       DNA sequences, one per line, A/C/G/T in either case.\n"
        "  W           Sliding-window length, %s; advance one base per window.\n"
        "  E           Exactly %s.\n"
        "  D           Positive number of hash residues.\n"
        "  R           Seeds per table set, 1..%s; use pivots 0..R-1.\n"
        "  S           Positive number of independent table sets to use.\n"
        "              Each window produces R*S seeds, in set/pivot order.\n"
        "  OUTPUT      Destination CSV; set/pivot/sequence/window_start are zero-based.\n"
        "  TABLE_FILE  Optional shared pivot tables; matching D, at least S sets.\n"
        "  --inward    Alternative inward traversal (changes scores and seeds).\n"
        "  --seed SEED Reproducible table creation when TABLE_FILE is omitted.\n\n"
        "Without TABLE_FILE, create tables in OUTPUT.piv. Generate reusable tables\n"
        "with 'subseqhash3 tables W D S OUTPUT [--seed SEED]'.\n",
        ind?"indel":"substitution",ind?"ind":"sub",ind?"1..256":"2..257",
        ind?"E deletions, no substitutions; 0<=E<W":"one deletion and E substitutions; 0<=E<W",
        ind?"W-E":"W-1");
}

typedef struct { FILE *out; size_t sequence,r; } CSV;
static int write_sub(size_t start,const SubseqHash3EditResult *a,size_t count,void *ctx) {
    CSV *c=ctx;
    for (size_t i=0;i<count;++i) {
        char row[512];
        int used=snprintf(row,sizeof(row),"%zu,%zu,%zu,%zu,%zu,%zu,",i/c->r,i%c->r,c->sequence,start,
            a[i].substitutions,a[i].deletions);
        if (used<0 || (size_t)used+a[i].length>=sizeof(row)) return -1;
        for (size_t j=0;j<a[i].length;++j) row[used++]="ACGT"[a[i].seed[j]];
        const int tail=snprintf(row+used,sizeof(row)-(size_t)used,",%zu,%.0f\n",a[i].psi,a[i].omega);
        if (tail<0 || (size_t)tail>=sizeof(row)-(size_t)used ||
            fwrite(row,1,(size_t)(used+tail),c->out)!=(size_t)(used+tail)) return -1;
    }
    return 0;
}
static int write_ind(size_t start,const SubseqHash3DeletionResult *a,size_t count,void *ctx) {
    CSV *c=ctx;
    for (size_t i=0;i<count;++i) {
        char row[512];
        int used=snprintf(row,sizeof(row),"%zu,%zu,%zu,%zu,0,%zu,",i/c->r,i%c->r,c->sequence,start,a[i].deletions);
        if (used<0 || (size_t)used+a[i].length>=sizeof(row)) return -1;
        for (size_t j=0;j<a[i].length;++j) row[used++]="ACGT"[a[i].seed[j]];
        const int tail=snprintf(row+used,sizeof(row)-(size_t)used,",%zu,%" PRId64 "\n",a[i].psi,a[i].omega);
        if (tail<0 || (size_t)tail>=sizeof(row)-(size_t)used ||
            fwrite(row,1,(size_t)(used+tail),c->out)!=(size_t)(used+tail)) return -1;
    }
    return 0;
}

/* Keep allocations and the converted deletion tables across lines that each
   contain one window. Public single-sequence APIs retain their ownership rules. */
typedef struct {
    SubseqHash3PivotWorkspace *pw;
    SubseqHash3DeletionWorkspace *dw;
    SubseqHash3DeletionTables *tables;
    SubseqHash3EditResult *sub;
    SubseqHash3DeletionResult *ind;
    size_t sets,total;
} SingleWindow;
static void single_destroy(SingleWindow *w) {
    subseqhash3_pivot_workspace_destroy(w->pw);
    subseqhash3_deletion_workspace_destroy(w->dw);
    if (w->tables) for (size_t f=0;f<w->sets;++f) subseqhash3_deletion_tables_destroy(w->tables+f);
    free(w->tables);
    if (w->sub) subseqhash3_edit_results_destroy(w->sub,w->total);
    if (w->ind) subseqhash3_deletion_results_destroy(w->ind,w->total);
    free(w->sub); free(w->ind);
}
static int single_seed(SingleWindow *w,const uint8_t *s,size_t n,size_t e,
    const SubseqHash3PivotTables *tables,size_t r,size_t sets,int ind,int inward,CSV *csv) {
    w->sets=sets; w->total=r*sets;
    if (ind) {
        if (!w->dw) {
            w->dw=subseqhash3_deletion_workspace_create();
            w->tables=calloc(sets,sizeof(*w->tables)); w->ind=calloc(w->total,sizeof(*w->ind));
            if (!w->dw || !w->tables || !w->ind) return -1;
            for (size_t f=0;f<sets;++f)
                if (subseqhash3_deletion_tables_from_pivots(w->tables+f,tables+f)) return -1;
        }
        for (size_t f=0;f<sets;++f)
            if ((inward?subseqhash3_compute_inward_deletion_seeds_with_workspace:
                    subseqhash3_compute_deletion_seeds_with_workspace)(s,n,e,w->tables+f,NULL,r,
                        w->ind+f*r,w->dw)) return -1;
        const int failed=write_ind(0,w->ind,w->total,csv);
        subseqhash3_deletion_results_destroy(w->ind,w->total); return failed;
    }
    if (!w->pw) {
        w->pw=subseqhash3_pivot_workspace_create(); w->sub=calloc(w->total,sizeof(*w->sub));
        if (!w->pw || !w->sub) return -1;
    }
    for (size_t f=0;f<sets;++f)
        if ((inward?subseqhash3_compute_pivot_repeats_with_workspace:
                subseqhash3_compute_outward_pivot_repeats_with_workspace)(s,n,e,tables+f,1,r,
                    w->sub+f*r,w->pw)) return -1;
    const int failed=write_sub(0,w->sub,w->total,csv);
    subseqhash3_edit_results_destroy(w->sub,w->total); return failed;
}

static int seeds_main(int argc,char **argv,int ind) {
    const char *args[8]; size_t argc_pos=0,capacity=0,sequence=0;
    uint64_t values[5],seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    int inward=0,seed_set=0,status=EXIT_FAILURE;
    SubseqHash3PivotTableFile tables={0};
    SingleWindow single={0};
    FILE *input=NULL,*output=NULL;
    char *generated=NULL,*line=NULL; uint8_t *encoded=NULL;
    struct timespec before,after;
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) { help(stdout,ind); return 0; }
    for (int i=1;i<argc;++i) {
        if (!strcmp(argv[i],"--inward")) {
            if (inward) goto invalid;
            inward=1;
        } else if (!strcmp(argv[i],"--seed")) {
            if (seed_set || i+1>=argc || ssh3_cli_number(argv[++i],&seed)) goto invalid;
            seed_set=1;
        } else {
            if (!strncmp(argv[i],"--",2) || argc_pos==8) goto invalid;
            args[argc_pos++]=argv[i];
        }
    }
    if (argc_pos<7 || (argc_pos==8 && seed_set)) goto invalid;
    for (size_t i=0;i<5;++i) if (ssh3_cli_number(args[i+1],values+i) || values[i]>SIZE_MAX) goto invalid;
    const size_t w=(size_t)values[0],e=(size_t)values[1],d=(size_t)values[2],r=(size_t)values[3],sets=(size_t)values[4];
    if (w<(ind?1u:2u) || w>(ind?256u:257u) || e>=w || !d || !r ||
        r>(ind?w-e:w-1) || !sets || sets>SIZE_MAX/r) goto invalid;
    if (ssh3_cli_alias(args[0],args[6]) || (argc_pos==8 && ssh3_cli_alias(args[7],args[6]))) {
        fputs("OUTPUT must not alias INPUT or TABLE_FILE\n",stderr); goto done;
    }
    input=fopen(args[0],"r"); if (!input) { perror(args[0]); goto done; }
    if (clock_gettime(CLOCK_MONOTONIC,&before)) goto done;
    if (ssh3_cli_tables(&tables,argc_pos==8?args[7]:NULL,&generated,args[6],args[0],
        w,d,sets,ind?w-e:w-1,seed)) goto done;
    if (ssh3_cli_alias(args[6],argc_pos==8?args[7]:generated)) {
        fputs("OUTPUT aliases the active table file\n",stderr); goto done;
    }
    output=fopen(args[6],"w"); if (!output) { perror(args[6]); goto done; }
    if (fputs("set,pivot,sequence,window_start,substitutions,deletions,seed,psi,omega\n",output)==EOF) goto done;
    ssize_t got;
    while ((got=getline(&line,&capacity,input))>=0) {
        size_t len=(size_t)got;
        while (len && (line[len-1]=='\n' || line[len-1]=='\r')) --len;
        if (len<w) { fprintf(stderr,"sequence %zu is empty or shorter than W\n",sequence); goto done; }
        encoded=malloc(len); if (!encoded) goto done;
        static const char alphabet[]="ACGTacgt";
        for (size_t j=0;j<len;++j) {
            const char *base=strchr(alphabet,line[j]);
            if (!line[j] || !base) { fprintf(stderr,"invalid DNA in sequence %zu\n",sequence); goto done; }
            encoded[j]=(uint8_t)((base-alphabet)%4);
        }
        CSV csv={output,sequence,r};
        int failed=len==w?single_seed(&single,encoded,w,e,tables.tables,r,sets,ind,inward,&csv):
            ind?subseqhash3_seeds_ind(encoded,len,w,e,tables.tables,r,sets,
                inward?SUBSEQHASH3_INWARD:SUBSEQHASH3_OUTWARD,write_ind,&csv):
            subseqhash3_seeds_sub(encoded,len,w,e,tables.tables,r,sets,
                inward?SUBSEQHASH3_INWARD:SUBSEQHASH3_OUTWARD,write_sub,&csv);
        if (failed) { fputs("pivot seeding failed\n",stderr); goto done; }
        free(encoded); encoded=NULL; ++sequence;
    }
    if (ferror(input) || fflush(output) || clock_gettime(CLOCK_MONOTONIC,&after)) goto done;
    fprintf(stderr,"elapsed wall time: %.6f seconds\n",(double)(after.tv_sec-before.tv_sec)+(double)(after.tv_nsec-before.tv_nsec)/1e9);
    status=EXIT_SUCCESS;
done:
    if (input && fclose(input)) status=EXIT_FAILURE;
    if (output && fclose(output)) status=EXIT_FAILURE;
    single_destroy(&single);
    subseqhash3_pivot_table_file_close(&tables);
    free(generated); free(line); free(encoded); return status;
invalid:
    help(stderr,ind); return EXIT_FAILURE;
}
int ssh3_seeds_sub_main(int argc,char **argv) { return seeds_main(argc,argv,0); }
int ssh3_seeds_ind_main(int argc,char **argv) { return seeds_main(argc,argv,1); }
