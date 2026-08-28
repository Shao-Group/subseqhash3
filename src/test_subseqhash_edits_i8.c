#include "subseqhash_edits_i8.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    int found;
    uint8_t seed[7];
    size_t positions[7];
} CanonicalBest;

static int candidate_less(const uint8_t *seed,const size_t *positions,
                          const CanonicalBest *best,size_t k) {
    size_t i;
    for (i=0;i<k;++i)
        if (seed[i]!=best->seed[i]) return seed[i]<best->seed[i];
    for (i=0;i<k;++i)
        if (positions[i]!=best->positions[i])
            return positions[i]<best->positions[i];
    return 0;
}

static void brute_symbols(const uint8_t *s,const size_t *positions,size_t k,
                          size_t row,size_t substitutions,size_t target_subs,
                          size_t residue,int omega,size_t target_psi,
                          int target_omega,const SubseqHashI8Tables *tables,
                          uint8_t *seed,CanonicalBest *best) {
    unsigned symbol;
    if (row==k) {
        if (substitutions==target_subs && residue==target_psi &&
            omega==target_omega &&
            (!best->found || candidate_less(seed,positions,best,k))) {
            best->found=1;
            memcpy(best->seed,seed,k);
            memcpy(best->positions,positions,k*sizeof(*positions));
        }
        return;
    }
    for (symbol=0;symbol<4;++symbol) {
        size_t changed=(size_t)(symbol!=s[positions[row]]),next,base;
        int next_omega;
        if (substitutions+changed>target_subs) continue;
        next=(residue+tables->C[4*row+symbol])%tables->d;
        base=(row*4+symbol)*tables->d+next;
        next_omega=tables->B1[base]*omega+
            tables->A[base]*tables->B2[base];
        seed[row]=(uint8_t)symbol;
        brute_symbols(s,positions,k,row+1,substitutions+changed,target_subs,
                      next,next_omega,target_psi,target_omega,tables,seed,best);
    }
}

static void brute_positions(const uint8_t *s,size_t n,size_t k,size_t at,
                            size_t used,size_t target_subs,size_t target_psi,
                            int target_omega,const SubseqHashI8Tables *tables,
                            size_t *positions,uint8_t *seed,
                            CanonicalBest *best) {
    size_t i;
    if (used==k) {
        brute_symbols(s,positions,k,0,0,target_subs,0,0,target_psi,
                      target_omega,tables,seed,best);
        return;
    }
    for (i=at;i<n && n-i>=k-used;++i) {
        positions[used]=i;
        brute_positions(s,n,k,i+1,used+1,target_subs,target_psi,
                        target_omega,tables,positions,seed,best);
    }
}

static void test_long_windows(void) {
    const size_t N=8,n=5,E=2,per=E+1,windows=N-n+1;
    uint8_t s[8]={0,1,2,3,1,0,2,3};
    SubseqHashI8Tables tables;
    SubseqHashEditResult actual[12],expected[3];
    size_t w,split,i;
    assert(!subseqhash_i8_tables_init(&tables,n,32,777));
    assert(!subseqhash_compute_long_edits_i8(
        s,N,n,E,&tables,actual));
    for (w=0;w<windows;++w) {
        assert(!subseqhash_compute_all_edit_splits_i8(
            s+w,n,E,&tables,expected));
        for (split=0;split<per;++split) {
            const SubseqHashEditResult *a=actual+w*per+split;
            const SubseqHashEditResult *e=expected+split;
            assert(a->deletions==e->deletions);
            assert(a->substitutions==e->substitutions);
            assert(a->length==e->length);
            assert(a->psi==e->psi);
            assert(a->omega==e->omega);
            assert(!memcmp(a->seed,e->seed,a->length));
            assert(!memcmp(a->substituted,e->substituted,a->length));
            for (i=0;i<a->length;++i)
                assert(a->positions[i]==e->positions[i]+w);
        }
        subseqhash_edit_results_destroy(expected,per);
    }
    subseqhash_edit_results_destroy(actual,windows*per);
    subseqhash_i8_tables_destroy(&tables);
}
int main(void) {

    size_t trial;
    test_long_windows();
    for (trial=0;trial<100;++trial) {
        const size_t n=3+trial%5,edits=trial%3,d=3+trial%13;
        uint8_t s[7]; size_t i,split;
        SubseqHashI8Tables compact;
        SubseqHashTables reference;
        SubseqHashEditResult actual[3],legacy[3],twopass[3],expected[3];
        for (i=0;i<n;++i) s[i]=(uint8_t)((trial*11+i*3+i*i)%4);
        assert(!subseqhash_i8_tables_init(&compact,n,d,12345+trial));
        assert(!subseqhash_i8_tables_make_double(&compact,&reference));
        assert(!subseqhash_compute_all_edit_splits_i8(
            s,n,edits,&compact,actual));
        assert(!subseqhash_compute_all_edit_splits_i8_legacy(
            s,n,edits,&compact,legacy));
        assert(!subseqhash_compute_all_edit_splits_i8_twopass(
            s,n,edits,&compact,twopass));
        assert(!subseqhash_compute_all_edit_splits(
            s,n,edits,&reference,expected));
        for (split=0;split<=edits;++split) {
            const size_t k=actual[split].length;
            assert(actual[split].deletions==expected[split].deletions);
            assert(actual[split].substitutions==expected[split].substitutions);
            assert(actual[split].psi==expected[split].psi);
            assert(actual[split].omega==expected[split].omega);
            assert(!memcmp(actual[split].seed,expected[split].seed,k));
            assert(!memcmp(actual[split].positions,expected[split].positions,
                           k*sizeof(*actual[split].positions)));
            assert(!memcmp(actual[split].substituted,
                           expected[split].substituted,k));
            assert(actual[split].psi==legacy[split].psi);
            assert(actual[split].omega==legacy[split].omega);
            assert(!memcmp(actual[split].seed,legacy[split].seed,k));
            assert(!memcmp(actual[split].positions,legacy[split].positions,
                           k*sizeof(*actual[split].positions)));
            assert(!memcmp(actual[split].substituted,
                           legacy[split].substituted,k));
            assert(actual[split].psi==twopass[split].psi);
            assert(actual[split].omega==twopass[split].omega);
            {
                size_t positions[7]; uint8_t seed[7];
                CanonicalBest best={0,{0},{0}};
                brute_positions(s,n,k,0,0,split,twopass[split].psi,
                                (int)twopass[split].omega,&compact,
                                positions,seed,&best);
                assert(best.found);
                assert(!memcmp(best.seed,twopass[split].seed,k));
                assert(!memcmp(best.positions,twopass[split].positions,
                               k*sizeof(*best.positions)));
            }
        }
        subseqhash_edit_results_destroy(expected,edits+1);
        subseqhash_edit_results_destroy(twopass,edits+1);
        subseqhash_edit_results_destroy(legacy,edits+1);
        subseqhash_edit_results_destroy(actual,edits+1);
        subseqhash_tables_destroy(&reference);
        subseqhash_i8_tables_destroy(&compact);
    }
    printf("all i8 edit tests passed (%s)\n",
           subseqhash_edits_i8_backend());
    return 0;
}
