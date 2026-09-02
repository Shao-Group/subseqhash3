#define _POSIX_C_SOURCE 200809L
#include "subseqhash3.h"
#include "table_file.h"
#include "ssh3_collisions.h"

#include <assert.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct {
    const uint8_t *sequence;
    const SubseqHash3Tables *tables;
    size_t n,k,substitutions;
    int minimum[64],maximum[64];
} BruteForce;

static void brute_symbols(BruteForce *brute,const size_t *positions,
                          size_t row,size_t substitutions,size_t residue,
                          int score) {
    unsigned symbol;
    if (row==brute->k) {
        if (substitutions==brute->substitutions) {
            if (score<brute->minimum[residue])
                brute->minimum[residue]=score;
            if (score>brute->maximum[residue])
                brute->maximum[residue]=score;
        }
        return;
    }
    for (symbol=0;symbol<4;++symbol) {
        const size_t changed=symbol!=brute->sequence[positions[row]];
        size_t next_residue,table_index;
        int next_score;
        if (substitutions+changed>brute->substitutions) continue;
        next_residue=(residue+brute->tables->C[4*row+symbol])%
                     brute->tables->d;
        table_index=(row*4+symbol)*brute->tables->d+next_residue;
        next_score=brute->tables->B1[table_index]*score+
                   brute->tables->A[table_index]*
                   brute->tables->B2[table_index];
        brute_symbols(brute,positions,row+1,substitutions+changed,
                      next_residue,next_score);
    }
}

static void brute_positions(BruteForce *brute,size_t at,size_t used,
                            size_t *positions) {
    size_t position;
    if (used==brute->k) {
        brute_symbols(brute,positions,0,0,0,0);
        return;
    }
    for (position=at;position<brute->n &&
         brute->n-position>=brute->k-used;++position) {
        positions[used]=position;
        brute_positions(brute,position+1,used+1,positions);
    }
}

static void validate_optimum(const uint8_t *sequence,size_t n,
                             const SubseqHash3Tables *tables,
                             const SubseqHash3EditResult *result) {
    BruteForce brute={sequence,tables,n,result->length,
                      result->substitutions,{0},{0}};
    size_t positions[7],residue;
    int expected;
    assert(tables->d<=64);
    for (residue=0;residue<tables->d;++residue) {
        brute.minimum[residue]=INT_MAX;
        brute.maximum[residue]=INT_MIN;
    }
    brute_positions(&brute,0,0,positions);
    {
        float best=0.0f; int have_best=0;
        size_t best_residue=0;
        for (residue=0;residue<tables->d;++residue) {
            size_t magnitude; float transformed;
            if (brute.minimum[residue]==INT_MAX) continue;
            expected=llabs((long long)brute.minimum[residue])>=
                     llabs((long long)brute.maximum[residue])?
                     brute.minimum[residue]:brute.maximum[residue];
            magnitude=(size_t)llabs((long long)expected);
            transformed=tables->F[residue*tables->f_stride+magnitude];
            if (!have_best || transformed>best) {
                have_best=1; best=transformed; best_residue=residue;
            }
        }
        residue=best_residue;
    }
    assert(residue==result->psi);
    expected=llabs((long long)brute.minimum[residue])>=
             llabs((long long)brute.maximum[residue])?
             brute.minimum[residue]:brute.maximum[residue];
    assert((double)expected==result->omega);
}

static void validate_result(const uint8_t *sequence,size_t window_start,
                            size_t window_length,
                            const SubseqHash3Tables *tables,
                            const SubseqHash3EditResult *result) {
    size_t row,residue=0,substitutions=0;
    int score=0;
    assert(result->length+result->deletions==window_length);
    for (row=0;row<result->length;++row) {
        const size_t position=result->positions[row];
        const unsigned symbol=result->seed[row];
        size_t table_index;
        assert(position>=window_start);
        assert(position<window_start+window_length);
        assert(!row || result->positions[row-1]<position);
        assert(symbol<4);
        assert(result->substituted[row]==(symbol!=sequence[position]));
        substitutions+=result->substituted[row];
        residue=(residue+tables->C[4*row+symbol])%tables->d;
        table_index=(row*4+symbol)*tables->d+residue;
        score=tables->B1[table_index]*score+
              tables->A[table_index]*tables->B2[table_index];
    }
    assert(substitutions==result->substitutions);
    assert(residue==result->psi);
    assert((double)score==result->omega);
}

static void test_long_windows(void) {
    const size_t sequence_length=8,window_length=5,edits=2;
    const size_t per_window=edits+1;
    const size_t windows=sequence_length-window_length+1;
    uint8_t sequence[8]={0,1,2,3,1,0,2,3};
    SubseqHash3Tables tables={0};
    SubseqHash3EditResult actual[12]={0},expected[3]={0};
    size_t window,split,row;

    assert(!subseqhash3_tables_init(&tables,window_length,32,777));
    assert(!subseqhash3_compute_long_edit_splits(
        sequence,sequence_length,window_length,edits,&tables,actual));
    for (window=0;window<windows;++window) {
        assert(!subseqhash3_compute_all_edit_splits(
            sequence+window,window_length,edits,&tables,expected));
        for (split=0;split<per_window;++split) {
            const SubseqHash3EditResult *a=actual+window*per_window+split;
            const SubseqHash3EditResult *e=expected+split;
            assert(a->deletions==e->deletions);
            assert(a->substitutions==e->substitutions);
            assert(a->length==e->length);
            assert(a->psi==e->psi);
            assert(a->omega==e->omega);
            assert(!memcmp(a->seed,e->seed,a->length));
            assert(!memcmp(a->substituted,e->substituted,a->length));
            for (row=0;row<a->length;++row)
                assert(a->positions[row]==e->positions[row]+window);
            validate_result(sequence,window,window_length,&tables,a);
        }
        subseqhash3_edit_results_destroy(expected,per_window);
    }
    subseqhash3_edit_results_destroy(actual,windows*per_window);
    subseqhash3_tables_destroy(&tables);
}

static void test_table_file(void) {
    char path[]="/tmp/subseqhash3-tables-XXXXXX";
    const size_t n=7,d=32,count=2,cells=n*d*4;
    const uint64_t seed=UINT64_C(0x123456789abcdef0);
    SubseqHash3TableFile file={0};
    size_t set;
    int descriptor=mkstemp(path);
    assert(descriptor>=0);
    assert(!close(descriptor));
    assert(!subseqhash3_table_file_create(path,n,d,count,seed));
    assert(!subseqhash3_table_file_open(&file,path));
    assert(file.n==n && file.d==d && file.count==count);
    for (set=0;set<count;++set) {
        SubseqHash3Tables expected={0};
        const SubseqHash3Tables *actual=file.tables+set;
        assert(!subseqhash3_tables_init(&expected,n,d,
            seed+set*UINT64_C(0x9e3779b97f4a7c15)));
        assert(!memcmp(actual->A,expected.A,cells));
        assert(!memcmp(actual->B1,expected.B1,cells));
        assert(!memcmp(actual->B2,expected.B2,cells));
        assert(!memcmp(actual->signed_A,expected.signed_A,cells));
        assert(!memcmp(actual->B1_masks,expected.B1_masks,
                       n*4*expected.blocks*sizeof(uint32_t)));
        assert(!memcmp(actual->C,expected.C,n*4*sizeof(size_t)));
        assert(actual->f_stride==expected.f_stride);
        assert(!memcmp(actual->F,expected.F,
                       d*expected.f_stride*sizeof(float)));
        subseqhash3_tables_destroy(&expected);
    }
    subseqhash3_table_file_close(&file);
    assert(!unlink(path));
}

static void test_monotone_function_selection(void) {
    uint8_t sequence[1]={0};
    SubseqHash3Tables tables={0};
    SubseqHash3EditResult results[2]={0};
    size_t target;
    assert(!subseqhash3_tables_init(&tables,1,4,123));
    for (target=0;target<tables.d;++target) {
        size_t x;
        assert(tables.F[target*tables.f_stride]>=0.0f);
        assert(tables.F[target*tables.f_stride]<=1.0f);
        for (x=1;x<tables.f_stride;++x)
            assert(tables.F[target*tables.f_stride+x]>=
                   tables.F[target*tables.f_stride+x-1]);
    }
    /* With one required substitution, symbols 1, 2, and 3 are candidates.
       Give only symbol 2's residue a positive transformed score. */
    target=tables.C[2];
    memset(tables.F,0,tables.d*tables.f_stride*sizeof(*tables.F));
    {
        size_t x;
        for (x=0;x<tables.f_stride;++x)
            tables.F[target*tables.f_stride+x]=1.0f;
    }
    assert(!subseqhash3_compute_all_edit_splits(sequence,1,1,&tables,results));
    assert(results[1].psi==target);
    assert(results[1].seed[0]==2);
    validate_result(sequence,0,1,&tables,results+1);
    subseqhash3_edit_results_destroy(results,2);
    subseqhash3_tables_destroy(&tables);
}

int ssh3_test_main(int argc,char **argv) {
    size_t trial;
    if (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h"))) {
        fputs("Run the built-in SubseqHash3 correctness tests.\n\n"
              "Usage:\n  subseqhash3 test\n",stdout);
        return EXIT_SUCCESS;
    }
    if (argc!=1) {
        fputs("Run the built-in SubseqHash3 correctness tests.\n\n"
              "Usage:\n  subseqhash3 test\n",stderr);
        return EXIT_FAILURE;
    }
    (void)argv;
    assert(ssh3_edit_neighborhood_jaccard("ACGT","ACGT",0,0)==1.0);
    assert(ssh3_edit_neighborhood_jaccard("A","C",0,0)==0.0);
    assert(fabs(ssh3_edit_neighborhood_jaccard("A","C",1,0)-0.5)<1e-12);
    assert(ssh3_edit_neighborhood_jaccard("AC","CA",0,1)==1.0);
    assert(ssh3_edit_neighborhood_jaccard("AC","GT",0,1)==0.0);
    assert(isnan(ssh3_edit_neighborhood_jaccard("AX","AC",0,0)));
    test_monotone_function_selection();
    test_table_file();
    test_long_windows();
    for (trial=0;trial<200;++trial) {
        const size_t n=3+trial%5,edits=trial%3;
        const size_t d=32*(1+trial%2);
        uint8_t sequence[7];
        SubseqHash3Tables tables={0};
        SubseqHash3EditResult results[3]={0};
        size_t i,split;

        for (i=0;i<n;++i)
            sequence[i]=(uint8_t)((trial*11+i*3+i*i)%4);
        assert(!subseqhash3_tables_init(&tables,n,d,12345+trial));
        assert(!subseqhash3_compute_all_edit_splits(
            sequence,n,edits,&tables,results));
        for (split=0;split<=edits;++split) {
            validate_result(sequence,0,n,&tables,results+split);
            validate_optimum(sequence,n,&tables,results+split);
        }
        subseqhash3_edit_results_destroy(results,edits+1);
        subseqhash3_tables_destroy(&tables);
    }
    printf("all SubseqHash3 tests passed (%s)\n",subseqhash3_backend());
    return 0;
}
