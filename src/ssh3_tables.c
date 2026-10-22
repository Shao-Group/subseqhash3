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
    fputs("Generate S independent sets of reusable pivot tables.\n\n"
          "Usage:\n  subseqhash3 tables W D S OUTPUT [--seed SEED]\n\n"
          "  W       Target window length, 1..257. Tables reserve min(W,256) rows.\n"
          "  D       Positive number of hash residues.\n"
          "  S       Positive number of independent table sets.\n"
          "  OUTPUT  Shared binary table file for seeds-sub, seeds-ind and collisions.\n"
          "  --seed SEED  Reproducible unsigned random seed.\n\n"
          "Each set supplies multiple pivot orders; R is selected when seeding.\n"
          "The extra row supports seeds-ind with zero deletions.\n",stream);
}

int ssh3_tables_main(int argc,char **argv) {
    uint64_t w,d,sets,seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    const char *args[4]; int count=0,has_seed=0;
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) {
        help(stdout); return EXIT_SUCCESS;
    }
    for (int i=1;i<argc;++i) {
        if (!strcmp(argv[i],"--seed")) {
            if (has_seed || i+1>=argc || parse_u64(argv[++i],&seed)) goto invalid;
            has_seed=1;
        } else {
            if (!strncmp(argv[i],"--",2) || count==4) goto invalid;
            args[count++]=argv[i];
        }
    }
    if (count!=4 || parse_u64(args[0],&w) || parse_u64(args[1],&d) ||
        parse_u64(args[2],&sets) || w<1 || w>257 || !d || !sets ||
        d>SIZE_MAX || sets>SIZE_MAX) goto invalid;
    return subseqhash3_seed_tables_create(args[3],(size_t)w,(size_t)d,
        (size_t)sets,seed)?EXIT_FAILURE:EXIT_SUCCESS;
invalid:
    help(stderr); return EXIT_FAILURE;
}
