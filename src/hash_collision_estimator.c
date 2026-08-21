#define _POSIX_C_SOURCE 200809L
#include "hash_collision_estimator.h"
#include "subseqhash.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static void trim_line(char *line, size_t *length) {
    while (*length && (line[*length-1]=='\n' || line[*length-1]=='\r')) line[--*length]='\0';
}
static int encode_dna(const char *text, size_t length, uint8_t *encoded) {
    size_t i;
    for (i=0; i<length; ++i) switch (text[i]) {
        case 'A': case 'a': encoded[i]=0; break;
        case 'C': case 'c': encoded[i]=1; break;
        case 'G': case 'g': encoded[i]=2; break;
        case 'T': case 't': encoded[i]=3; break;
        default: return -1;
    }
    return 0;
}
int estimate_subseqhash_collisions(const char *filename, size_t k, size_t D,
                                   HashCollisionStats *stats) {
    FILE *input=NULL; size_t cap1=0,cap2=0; uint8_t *w1=NULL,*w2=NULL;
    char *line1=NULL,*line2=NULL; SubseqHashTables tables; int status=-1;
    if (!filename || !stats || !D || k>SIZE_MAX/4) return -1;
    memset(stats,0,sizeof(*stats));
    stats->table_seed=(uint64_t)time(NULL)^((uint64_t)clock()<<32)^(uint64_t)(uintptr_t)stats;
    memset(&tables,0,sizeof(tables)); input=fopen(filename,"r");
    w1=malloc(k?k:1); w2=malloc(k?k:1);
    if (!input || !w1 || !w2 || subseqhash_tables_init(&tables,k,D,stats->table_seed)) {
        if (!input) perror(filename); goto cleanup;
    }
    for (;;) {
        ssize_t r1=getline(&line1,&cap1,input),r2; size_t n1,n2,psi; double omega; uint8_t *s1,*s2;
        if (r1<0) { status=ferror(input)?-1:0; break; }
        r2=getline(&line2,&cap2,input);
        if (r2<0) { fprintf(stderr,"%s has an odd number of lines\n",filename); break; }
        n1=(size_t)r1; n2=(size_t)r2; trim_line(line1,&n1); trim_line(line2,&n2); ++stats->pairs;
        if (n1<k || n2<k) { ++stats->invalid_pairs; continue; }
        s1=malloc(n1?n1:1); s2=malloc(n2?n2:1);
        if (!s1 || !s2) { free(s2); free(s1); break; }
        if (encode_dna(line1,n1,s1) || encode_dna(line2,n2,s2)) {
            fprintf(stderr,"invalid DNA character in pair %zu\n",stats->pairs); free(s2); free(s1); break;
        }
        if (subseqhash_compute_optimized(s1,n1,&tables,w1,NULL,&psi,&omega) ||
            subseqhash_compute_optimized(s2,n2,&tables,w2,NULL,&psi,&omega)) {
            fprintf(stderr,"witness computation failed for pair %zu\n",stats->pairs); free(s2); free(s1); break;
        }
        ++stats->valid_pairs; if (!memcmp(w1,w2,k)) ++stats->collisions; free(s2); free(s1);
    }
cleanup:
    if (stats->pairs) stats->frequency=(double)stats->collisions/stats->pairs;
    if (stats->valid_pairs) stats->valid_frequency=(double)stats->collisions/stats->valid_pairs;
    free(line2); free(line1); free(w2); free(w1);
    subseqhash_tables_destroy(&tables); if (input) fclose(input); return status;
}

int estimate_hash_collisions(const char *filename, size_t k, size_t D,
                             HashCollisionStats *stats) {
    return estimate_subseqhash_collisions(filename,k,D,stats);
}
