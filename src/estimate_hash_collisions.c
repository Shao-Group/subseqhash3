#include "hash_collision_estimator.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int parse_size(const char *text,size_t *value) {
    char *end; unsigned long long parsed; errno=0; parsed=strtoull(text,&end,10);
    if (errno || !*text || *end || parsed>(unsigned long long)SIZE_MAX) return -1;
    *value=(size_t)parsed; return 0;
}
int main(int argc,char **argv) {
    size_t k,D; HashCollisionStats stats; int original=0;
    if ((argc!=4 && argc!=5) || parse_size(argv[2],&k) ||
        parse_size(argv[3],&D) || !D ||
        (argc==5 && strcmp(argv[4],"--original"))) {
        fprintf(stderr,"usage: %s FILE K D [--original]\n",argv[0]); return EXIT_FAILURE;
    }
    original=argc==5;
    if ((original?estimate_subseqhash_collisions_original:
                  estimate_subseqhash_collisions)(argv[1],k,D,&stats))
        return EXIT_FAILURE;
    printf("implementation: %s\n",original?"original":"optimized");
    printf("pairs: %zu\nvalid pairs: %zu\ninvalid pairs: %zu\n",stats.pairs,stats.valid_pairs,stats.invalid_pairs);
    printf("hash collisions: %zu\n",stats.collisions);
    printf("collision frequency (all pairs): %.10f\n",stats.frequency);
    printf("collision frequency (valid pairs): %.10f\n",stats.valid_frequency);
    printf("table seed: %llu\n",(unsigned long long)stats.table_seed); return EXIT_SUCCESS;
}
