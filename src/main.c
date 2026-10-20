#include "ssh3_commands.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(FILE *stream,const char *program) {
    fprintf(stream,
            "SubseqHash3: generate and evaluate edit-aware sequence seeds.\n"
            "\n"
            "Usage:\n"
            "  %s <command> [arguments]\n"
            "\n"
            "Commands:\n"
            "  simulate    Generate pairs with exact substitutions, deletions, and insertions.\n"
            "  tables      Generate reusable legacy or --pivots table sets.\n"
            "  seeds       Generate seeds for sequences and sliding windows.\n"
            "  repeats     Generate pivot seeds with one deletion and E substitutions.\n"
            "  deletions   Generate SubseqHash2 pivot seeds with exactly E deletions.\n"
            "  collisions  Estimate seed-collision probabilities from string pairs.\n"
            "  test        Run the built-in correctness tests.\n"
            "\n"
            "Run '%s <command> --help' for command-specific help.\n",
            program,
            program);
}

int main(int argc,char **argv) {
    if (argc<2 || !strcmp(argv[1],"--help") || !strcmp(argv[1],"-h")) {
        usage(argc<2?stderr:stdout,argv[0]);
        return argc<2?EXIT_FAILURE:EXIT_SUCCESS;
    }
    if (!strcmp(argv[1],"simulate")) return ssh3_simulate_main(argc-1,argv+1);
    if (!strcmp(argv[1],"tables")) return ssh3_tables_main(argc-1,argv+1);
    if (!strcmp(argv[1],"seeds")) return ssh3_seeds_main(argc-1,argv+1);
    if (!strcmp(argv[1],"deletions")) return ssh3_deletions_main(argc-1,argv+1);
    if (!strcmp(argv[1],"repeats")) return ssh3_repeats_main(argc-1,argv+1);
    if (!strcmp(argv[1],"collisions")) return ssh3_collisions_main(argc-1,argv+1);
    if (!strcmp(argv[1],"test")) return ssh3_test_main(argc-1,argv+1);
    fprintf(stderr,"subseqhash3: unknown command '%s'\n\n",argv[1]);
    usage(stderr,argv[0]);
    return EXIT_FAILURE;
}
