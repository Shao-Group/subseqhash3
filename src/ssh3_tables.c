#include "table_file.h"
#include "pivot_table_file.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int parse_u64(const char *text,uint64_t *value) {
    char *end; unsigned long long parsed;
    if (!*text || *text<'0' || *text>'9') return -1;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end) return -1;
    *value=(uint64_t)parsed; return 0;
}

static void help(FILE *stream) {
    fputs("Generate reusable SubseqHash3 or pivot table sets.\n\n"
          "Usage:\n  subseqhash3 tables n D R OUTPUT [--pivots] [--seed SEED]\n\n"
          "Arguments:\n"
          "  n       Legacy: maximum sequence/window length.\n"
          "          With --pivots: n-1 rows/orders per family (2<=n<=257).\n"
          "  D       Number of hash residues.\n"
          "  R       Legacy: number of independent table sets.\n"
          "          With --pivots: number of requested seed repetitions.\n"
          "  OUTPUT  Destination binary table file.\n\n"
          "Options:\n"
          "  --seed SEED  Use a reproducible unsigned random seed.\n"
          "  --pivots     Shared tables for repeats/deletions: ceil(R/(n-1)) families,\n"
          "               each containing n-1 pivot orders; requires 2<=n<=257.\n\n"
          "repeats requires n>=W and enough families for R repetitions.\n"
          "deletions requires n-1>=W-E, uses the first family, and R<=W-E.\n"
          "For E=0 deletions, generate n>=W+1.\n\n"
          "Default legacy files work with seeds/collisions. Pivot files\n"
          "store both side tables and pivot parameters in a distinct format.\n",
          stream);
}

int ssh3_tables_main(int argc,char **argv) {
    uint64_t n,d,r,seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    const char *positional[4];
    int i,count=0,pivots=0,has_seed=0;
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) {
        help(stdout); return EXIT_SUCCESS;
    }
    for (i=1;i<argc;++i) {
        if (!strcmp(argv[i],"--pivots")) {
            if (pivots) goto invalid;
            pivots=1;
        } else if (!strcmp(argv[i],"--seed")) {
            if (has_seed || i+1>=argc || parse_u64(argv[++i],&seed)) goto invalid;
            has_seed=1;
        } else {
            if (count==4) goto invalid;
            positional[count++]=argv[i];
        }
    }
    if (count!=4 || parse_u64(positional[0],&n) ||
        parse_u64(positional[1],&d) || parse_u64(positional[2],&r) ||
        !n || !d || !r || n>SIZE_MAX || d>SIZE_MAX || r>SIZE_MAX) goto invalid;
    if (pivots) {
        size_t families;
        if (n<2 || n>257) goto invalid;
        families=1+((size_t)r-1)/((size_t)n-1);
        return subseqhash3_pivot_table_file_create(positional[3],(size_t)n,
            (size_t)d,families,seed)?EXIT_FAILURE:EXIT_SUCCESS;
    }
    return subseqhash3_table_file_create(positional[3],(size_t)n,(size_t)d,
            (size_t)r,seed)?EXIT_FAILURE:EXIT_SUCCESS;
invalid:
    help(stderr);
    return EXIT_FAILURE;
}
