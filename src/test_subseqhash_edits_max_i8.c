#include "subseqhash_edits_max_i8.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct { size_t psi; int omega; int found; } Best;

static void visit_symbols(const uint8_t *s,const size_t *positions,size_t k,
                          size_t row,size_t substitutions,size_t target,
                          size_t residue,int omega,
                          const SubseqHashMaxI8Tables *t,Best *best) {
    unsigned a;
    if (row==k) {
        if (substitutions==target &&
            (!best->found || residue<best->psi ||
             (residue==best->psi && omega>best->omega))) {
            best->found=1; best->psi=residue; best->omega=omega;
        }
        return;
    }
    for (a=0;a<4;++a) {
        size_t changed=(size_t)(a!=s[positions[row]]);
        size_t next;
        if (substitutions+changed>target) continue;
        next=(residue+t->C[4*row+a])%t->d;
        visit_symbols(s,positions,k,row+1,substitutions+changed,target,
                      next,omega+t->A[(row*4+a)*t->d+next],t,best);
    }
}

static void visit_positions(const uint8_t *s,size_t n,size_t k,size_t at,
                            size_t used,size_t substitutions,size_t *positions,
                            const SubseqHashMaxI8Tables *t,Best *best) {
    size_t i;
    if (used==k) {
        visit_symbols(s,positions,k,0,0,substitutions,0,0,t,best); return;
    }
    for (i=at;i<n && n-i>=k-used;++i) {
        positions[used]=i;
        visit_positions(s,n,k,i+1,used+1,substitutions,positions,t,best);
    }
}

static void check_witness(const uint8_t *s,size_t n,
                          const SubseqHashMaxI8Tables *t,
                          const SubseqHashEditResult *out) {
    size_t i,residue=0,substitutions=0; int omega=0;
    for (i=0;i<out->length;++i) {
        unsigned a=out->seed[i];
        assert(a<4 && out->positions[i]<n);
        if (i) assert(out->positions[i-1]<out->positions[i]);
        assert(out->substituted[i]==(uint8_t)(a!=s[out->positions[i]]));
        substitutions+=(size_t)(a!=s[out->positions[i]]);
        residue=(residue+t->C[4*i+a])%t->d;
        omega+=t->A[(i*4+a)*t->d+residue];
    }
    assert(residue==out->psi && (double)omega==out->omega);
    assert(substitutions==out->substitutions);
    assert(n-out->length==out->deletions);
}

int main(void) {
    size_t trial;
    for (trial=0;trial<40;++trial) {
        size_t n=3+trial%4,edits=trial%3,d=3+trial%7,i,split;
        uint8_t s[6]; size_t positions[6];
        SubseqHashMaxI8Tables tables;
        SubseqHashEditResult results[3];
        for (i=0;i<n;++i) s[i]=(uint8_t)((trial*7+i*3+i*i)%4);
        assert(!subseqhash_max_i8_tables_init(&tables,n,d,12345+trial));
        assert(!subseqhash_compute_all_edit_splits_max_i8(
            s,n,edits,&tables,results));
        for (split=0;split<=edits;++split) {
            size_t deletions=edits-split,k=n-deletions;
            Best best={0,INT_MIN,0};
            visit_positions(s,n,k,0,0,split,positions,&tables,&best);
            assert(best.found);
            assert(results[split].psi==best.psi);
            assert(results[split].omega==(double)best.omega);
            check_witness(s,n,&tables,results+split);
        }
        subseqhash_edit_results_destroy(results,edits+1);
        subseqhash_max_i8_tables_destroy(&tables);
    }
    {
        SubseqHashMaxI8Tables tables;
        assert(!subseqhash_max_i8_tables_init(
            &tables,SUBSEQHASH_MAX_I8_MAX_LENGTH,1,1));
        subseqhash_max_i8_tables_destroy(&tables);
        assert(subseqhash_max_i8_tables_init(
            &tables,SUBSEQHASH_MAX_I8_MAX_LENGTH+1,1,1));
    }
    printf("all max-only i8 edit tests passed (%s)\n",
           subseqhash_edits_max_i8_backend());
    return 0;
}
