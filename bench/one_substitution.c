#define _POSIX_C_SOURCE 200809L
#include "subseqhash3.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static uint64_t random_state=1234567;
static unsigned random_base(void) {
    random_state^=random_state<<13;
    random_state^=random_state>>7;
    random_state^=random_state<<17;
    return (unsigned)(random_state&3u);
}

static void compare(const uint8_t *s,size_t n,size_t e,
                    const SubseqHash3Tables *tables) {
    SubseqHash3EditResult fixed={0};
    SubseqHash3EditResult *all=calloc(e+2,sizeof(*all));
    const SubseqHash3EditResult *expected=all+1;
    assert(all);
    assert(!subseqhash3_compute_one_substitution(s,n,e,tables,&fixed));
    assert(!subseqhash3_compute_all_edit_splits(s,n,e+1,tables,all));
    assert(fixed.deletions==e && fixed.substitutions==1);
    assert(fixed.length==n-e && fixed.length==expected->length);
    assert(fixed.psi==expected->psi && fixed.omega==expected->omega);
    assert(!memcmp(fixed.seed,expected->seed,fixed.length));
    assert(!memcmp(fixed.positions,expected->positions,
                   fixed.length*sizeof(*fixed.positions)));
    assert(!memcmp(fixed.substituted,expected->substituted,fixed.length));
    subseqhash3_edit_results_destroy(&fixed,1);
    subseqhash3_edit_results_destroy(all,e+2);
    free(all);
}

static void correctness(void) {
    const size_t residues[]={1,7,31,32,33,47,63,64,96,256,1024,1025};
    size_t trial;
    for (trial=0;trial<360;++trial) {
        const size_t n=1+trial%64;
        const size_t e=(trial/4)%n;
        const size_t d=residues[trial%12];
        SubseqHash3Tables tables={0};
        uint8_t s[64];
        size_t i;
        for (i=0;i<n;++i)
            s[i]=(uint8_t)(trial%3?random_base():trial%4);
        assert(!subseqhash3_tables_init(&tables,n,d,trial+12345));
        compare(s,n,e,&tables);
        subseqhash3_tables_destroy(&tables);
    }
    {
        /* Exercise the largest supported score, absorbing infinity, ties,
           and the search for a reachable residue below D-1. */
        uint8_t s[256]={0};
        SubseqHash3Tables tables={0};
        size_t i;
        assert(!subseqhash3_tables_init(&tables,256,32,123));
        for (i=0;i<256*32*4;++i) {
            tables.A[i]=tables.signed_A[i]=127;
            tables.B1[i]=tables.B2[i]=1;
        }
        for (i=0;i<256*4;++i) {
            tables.B1_masks[i]=UINT32_MAX;
            tables.C[i]=0;
        }
        compare(s,256,0,&tables);
        compare(s,256,8,&tables);
        compare(s,64,63,&tables);
        for (i=0;i<256*32*4;++i) {
            tables.signed_A[i]=-127;
            tables.B2[i]=-1;
        }
        compare(s,256,0,&tables);
        compare(s,256,8,&tables);
        subseqhash3_tables_destroy(&tables);
    }
    {
        uint8_t s[]={0,1,2,3},bad[]={4};
        SubseqHash3Tables tables={0};
        SubseqHash3EditResult out={0};
        assert(!subseqhash3_tables_init(&tables,4,32,1));
        assert(subseqhash3_compute_one_substitution(s,4,4,&tables,&out)==-1);
        assert(subseqhash3_compute_one_substitution(s,0,0,&tables,&out)==-1);
        assert(subseqhash3_compute_one_substitution(bad,1,0,&tables,&out)==-1);
        assert(subseqhash3_compute_one_substitution(NULL,4,1,&tables,&out)==-1);
        assert(subseqhash3_compute_one_substitution(s,4,1,NULL,&out)==-1);
        assert(subseqhash3_compute_one_substitution(s,4,1,&tables,NULL)==-1);
        subseqhash3_tables_destroy(&tables);
    }
    puts("365 exact equivalence cases and invalid-input checks passed");
}

static double seconds(void) {
    struct timespec ts;
    assert(!clock_gettime(CLOCK_MONOTONIC,&ts));
    return ts.tv_sec+ts.tv_nsec*1e-9;
}

static double measure(const uint8_t *sequences,size_t n,size_t e,
                      const SubseqHash3Tables *tables,int fixed) {
    enum { INPUTS=32, LOOPS=4 };
    SubseqHash3EditResult *results=calloc(e+2,sizeof(*results));
    double start;
    size_t i;
    assert(results);
    start=seconds();
    for (i=0;i<INPUTS*LOOPS;++i) {
        const uint8_t *s=sequences+(i%INPUTS)*n;
        if (fixed)
            assert(!subseqhash3_compute_one_substitution(s,n,e,tables,results));
        else
            assert(!subseqhash3_compute_all_edit_splits(s,n,e+1,tables,results));
        subseqhash3_edit_results_destroy(results,fixed?1:e+2);
    }
    start=(seconds()-start)/(INPUTS*LOOPS);
    free(results);
    return start;
}

static double median(double *x) {
    double tmp;
    if (x[0]>x[1]) { tmp=x[0]; x[0]=x[1]; x[1]=tmp; }
    if (x[1]>x[2]) { tmp=x[1]; x[1]=x[2]; x[2]=tmp; }
    return x[0]>x[1]?x[0]:x[1];
}

int main(int argc,char **argv) {
    const size_t cases[][3]={
        {20,5,1024},{30,7,1024},{50,8,1024},
        {50,16,1024},{50,24,1024},{64,32,1024},
        {50,8,32},{50,8,256},{50,8,1025}
    };
    size_t c;
    correctness();
    if (argc==2 && !strcmp(argv[1],"--check-only")) return 0;
    printf("backend: %s\n",subseqhash3_backend());
    puts("n,deletions,D,fixed_us,all_splits_us,speedup");
    for (c=0;c<sizeof(cases)/sizeof(*cases);++c) {
        const size_t n=cases[c][0],e=cases[c][1],d=cases[c][2];
        uint8_t sequences[32*64];
        SubseqHash3Tables tables={0};
        double fixed[3],all[3],f,a;
        size_t i,round;
        assert(!subseqhash3_tables_init(&tables,n,d,c+987654));
        for (i=0;i<32*n;++i) sequences[i]=(uint8_t)random_base();
        for (i=0;i<32;++i) compare(sequences+i*n,n,e,&tables);
        for (round=0;round<3;++round) {
            /* Alternate order to reduce thermal/frequency bias. */
            if (round&1) {
                all[round]=measure(sequences,n,e,&tables,0);
                fixed[round]=measure(sequences,n,e,&tables,1);
            } else {
                fixed[round]=measure(sequences,n,e,&tables,1);
                all[round]=measure(sequences,n,e,&tables,0);
            }
        }
        f=median(fixed); a=median(all);
        printf("%zu,%zu,%zu,%.3f,%.3f,%.2fx\n",n,e,d,f*1e6,a*1e6,a/f);
        fflush(stdout);
        subseqhash3_tables_destroy(&tables);
    }
    return 0;
}
