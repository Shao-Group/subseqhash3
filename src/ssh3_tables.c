#include "table_file.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int parse_u64(const char *text,uint64_t *value) {
    char *end; unsigned long long parsed;
    errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end) return -1;
    *value=(uint64_t)parsed; return 0;
}

static void help(FILE *stream) {
    fputs("Generate reusable SubseqHash3 table sets in mmap-ready binary format.\n\n"
          "Usage:\n  subseqhash3 tables n D R OUTPUT [--seed SEED]\n\n"
          "Arguments:\n"
          "  n       Maximum supported sequence or window length.\n"
          "  D       Number of hash residues.\n"
          "  R       Number of independently randomized table sets.\n"
          "  OUTPUT  Destination binary table file.\n\n"
          "Options:\n"
          "  --seed SEED  Use a reproducible unsigned random seed.\n",
          stream);
}

int ssh3_tables_main(int argc,char **argv) {
    uint64_t n,d,r,seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32);
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) {
        help(stdout); return EXIT_SUCCESS;
    }
    if ((argc!=5 && argc!=7) || parse_u64(argv[1],&n) ||
        parse_u64(argv[2],&d) || parse_u64(argv[3],&r) || !n || !d || !r ||
        n>SIZE_MAX || d>SIZE_MAX || r>SIZE_MAX ||
        (argc==7 && (strcmp(argv[5],"--seed") || parse_u64(argv[6],&seed)))) {
        help(stderr);
        return EXIT_FAILURE;
    }
    return subseqhash3_table_file_create(argv[4],(size_t)n,(size_t)d,
            (size_t)r,seed)?EXIT_FAILURE:EXIT_SUCCESS;
}
