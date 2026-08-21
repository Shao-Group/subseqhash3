#include "subseqhash_adjacent.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct { int found; size_t psi; double omega; } Best;

static void enumerate(const uint8_t *s,size_t n,size_t source,size_t k,
                      size_t chosen,size_t deletions,size_t substitutions,
                      uint8_t *candidate,
                      const SubseqHashAdjacentTables *tables,Best *best) {
    size_t a,psi; double omega;
    if (source==n) {
        if (chosen!=k || deletions || substitutions) return;
        assert(!subseqhash_adjacent_score(candidate,k,tables,&psi,&omega));
        if (!best->found || psi<best->psi ||
            (psi==best->psi && fabs(omega)>fabs(best->omega))) {
            best->found=1; best->psi=psi; best->omega=omega;
        }
        return;
    }
    if (deletions)
        enumerate(s,n,source+1,k,chosen,deletions-1,substitutions,
                  candidate,tables,best);
    if (chosen==k) return;
    for (a=0;a<4;++a) {
        size_t changed=a!=s[source];
        if (changed>substitutions) continue;
        candidate[chosen]=(uint8_t)a;
        enumerate(s,n,source+1,k,chosen+1,deletions,
                  substitutions-changed,candidate,tables,best);
    }
}

int main(void) {
    size_t trial;
    for (trial=0;trial<60;++trial) {
        const size_t n=3+trial%3,total=trial%3,d=3+trial%7;
        uint8_t s[5],candidate[5],fixed_seed[5],fixed_changed[5];
        size_t fixed_positions[5],split,i,fixed_psi,check_psi;
        double fixed_omega,check_omega;
        SubseqHashAdjacentTables tables;
        SubseqHashEditResult results[3];
        for (i=0;i<n;++i) s[i]=(uint8_t)((trial*13+i*3+i*i)%4);
        assert(!subseqhash_adjacent_tables_init(&tables,n,d,7000+trial));
        assert(!subseqhash_adjacent_compute_all_edit_splits(
            s,n,total,&tables,results));
        for (split=0;split<=total;++split) {
            const size_t substitutions=split,deletions=total-split;
            const size_t k=n-deletions;
            size_t changed=0;
            Best best={0,0,0};
            enumerate(s,n,0,k,0,deletions,substitutions,candidate,
                      &tables,&best);
            assert(best.found);
            assert(results[split].psi==best.psi);
            assert(fabs(results[split].omega)==fabs(best.omega));
            assert(!subseqhash_adjacent_score(results[split].seed,k,&tables,
                                               &check_psi,&check_omega));
            assert(check_psi==results[split].psi);
            assert(check_omega==results[split].omega);
            for (i=0;i<k;++i) {
                assert(results[split].positions[i]<n);
                if (i) assert(results[split].positions[i-1]<
                              results[split].positions[i]);
                assert(results[split].substituted[i]==
                       (results[split].seed[i]!=
                        s[results[split].positions[i]]));
                changed+=results[split].substituted[i];
            }
            assert(changed==substitutions);
            assert(!subseqhash_adjacent_compute_edits(
                s,n,deletions,substitutions,&tables,fixed_seed,
                fixed_positions,fixed_changed,&fixed_psi,&fixed_omega));
            assert(fixed_psi==results[split].psi);
            assert(fixed_omega==results[split].omega);
            assert(!memcmp(fixed_seed,results[split].seed,k));
            assert(!memcmp(fixed_positions,results[split].positions,
                           k*sizeof(*fixed_positions)));
            assert(!memcmp(fixed_changed,results[split].substituted,k));
        }
        subseqhash_edit_results_destroy(results,total+1);
        subseqhash_adjacent_tables_destroy(&tables);
    }
    puts("all adjacent subseqhash tests passed");
    return 0;
}
