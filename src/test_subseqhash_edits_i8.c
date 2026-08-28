#include "subseqhash_edits_i8.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const uint8_t *sequence;
    const SubseqHashI8Tables *tables;
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
                             const SubseqHashI8Tables *tables,
                             const SubseqHashEditResult *result) {
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
    for (residue=0;residue<tables->d &&
         brute.minimum[residue]==INT_MAX;++residue) {}
    assert(residue==result->psi);
    expected=llabs((long long)brute.minimum[residue])>=
             llabs((long long)brute.maximum[residue])?
             brute.minimum[residue]:brute.maximum[residue];
    assert((double)expected==result->omega);
}

static void validate_result(const uint8_t *sequence,size_t window_start,
                            size_t window_length,
                            const SubseqHashI8Tables *tables,
                            const SubseqHashEditResult *result) {
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
    SubseqHashI8Tables tables={0};
    SubseqHashEditResult actual[12]={0},expected[3]={0};
    size_t window,split,row;

    assert(!subseqhash_i8_tables_init(&tables,window_length,32,777));
    assert(!subseqhash_compute_long_edits_i8_fast(
        sequence,sequence_length,window_length,edits,&tables,actual));
    for (window=0;window<windows;++window) {
        assert(!subseqhash_compute_all_edit_splits_i8(
            sequence+window,window_length,edits,&tables,expected));
        for (split=0;split<per_window;++split) {
            const SubseqHashEditResult *a=actual+window*per_window+split;
            const SubseqHashEditResult *e=expected+split;
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
        subseqhash_edit_results_destroy(expected,per_window);
    }
    subseqhash_edit_results_destroy(actual,windows*per_window);
    subseqhash_i8_tables_destroy(&tables);
}

int main(void) {
    size_t trial;
    test_long_windows();
    for (trial=0;trial<200;++trial) {
        const size_t n=3+trial%5,edits=trial%3;
        const size_t d=32*(1+trial%2);
        uint8_t sequence[7];
        SubseqHashI8Tables tables={0};
        SubseqHashEditResult results[3]={0};
        size_t i,split;

        for (i=0;i<n;++i)
            sequence[i]=(uint8_t)((trial*11+i*3+i*i)%4);
        assert(!subseqhash_i8_tables_init(&tables,n,d,12345+trial));
        assert(!subseqhash_compute_all_edit_splits_i8(
            sequence,n,edits,&tables,results));
        for (split=0;split<=edits;++split) {
            validate_result(sequence,0,n,&tables,results+split);
            validate_optimum(sequence,n,&tables,results+split);
        }
        subseqhash_edit_results_destroy(results,edits+1);
        subseqhash_i8_tables_destroy(&tables);
    }
    printf("all i8 edit tests passed (%s)\n",
           subseqhash_edits_i8_backend());
    return 0;
}
