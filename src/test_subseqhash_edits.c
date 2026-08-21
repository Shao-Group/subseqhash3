#include "subseqhash_edits.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    int found;
    size_t psi;
    double omega;
} EditBest;

static void enumerate_edits(const uint8_t *s,size_t n,size_t source,
                            size_t k,size_t chosen,size_t deletions_left,
                            size_t substitutions_left,uint8_t *candidate,
                            const SubseqHashTables *tables,EditBest *best) {
    size_t a,psi; double omega;
    if (source==n) {
        if (chosen!=k || deletions_left || substitutions_left) return;
        assert(!subseqhash_score(candidate,tables,&psi,&omega));
        if (!best->found || psi<best->psi ||
            (psi==best->psi && fabs(omega)>fabs(best->omega))) {
            best->found=1; best->psi=psi; best->omega=omega;
        }
        return;
    }
    if (deletions_left)
        enumerate_edits(s,n,source+1,k,chosen,deletions_left-1,
                        substitutions_left,candidate,tables,best);
    if (chosen==k) return;
    for (a=0;a<4;++a) {
        const size_t changed=a!=s[source];
        if (changed>substitutions_left) continue;
        candidate[chosen]=(uint8_t)a;
        enumerate_edits(s,n,source+1,k,chosen+1,deletions_left,
                        substitutions_left-changed,candidate,tables,best);
    }
}

int main(void) {
    size_t trial;
    for (trial=0;trial<120;++trial) {
        const size_t n=3+trial%4,deletions=trial%3;
        const size_t k=n-deletions;
        const size_t substitutions=k?(trial/3)%(k<3?k+1:3):0;
        const size_t d=3+trial%8;
        size_t i,psi,check_psi,changed=0;
        double omega,check_omega;
        uint8_t s[6],seed[6],candidate[6],was_substituted[6];
        size_t positions[6];
        SubseqHashTables tables;
        EditBest best={0,0,0};

        for (i=0;i<n;++i) s[i]=(uint8_t)((trial*11+i*3+i*i)%4);
        assert(!subseqhash_tables_init(&tables,k,d,UINT64_C(9001)+trial));
        enumerate_edits(s,n,0,k,0,deletions,substitutions,candidate,
                        &tables,&best);
        assert(best.found);
        assert(!subseqhash_compute_edits(s,n,deletions,substitutions,
                                         &tables,seed,positions,
                                         was_substituted,&psi,&omega));
        assert(psi==best.psi);
        assert(fabs(omega)==fabs(best.omega));
        assert(!subseqhash_score(seed,&tables,&check_psi,&check_omega));
        assert(check_psi==psi); assert(check_omega==omega);
        for (i=0;i<k;++i) {
            assert(positions[i]<n);
            if (i) assert(positions[i-1]<positions[i]);
            assert(was_substituted[i]==(seed[i]!=s[positions[i]]));
            changed+=was_substituted[i];
        }
        assert(changed==substitutions);
        subseqhash_tables_destroy(&tables);
    }
    {
        const uint8_t s[2]={0,1}; uint8_t seed[2];
        size_t psi; double omega; SubseqHashTables tables;
        assert(!subseqhash_tables_init(&tables,1,5,7));
        assert(subseqhash_compute_edits(s,2,0,0,&tables,seed,NULL,NULL,
                                        &psi,&omega));
        assert(subseqhash_compute_edits(s,2,1,2,&tables,seed,NULL,NULL,
                                        &psi,&omega));
        subseqhash_tables_destroy(&tables);
    }
    for (trial=0;trial<60;++trial) {
        const size_t n=3+trial%4,total_edits=trial%3,d=4+trial%7;
        uint8_t s[6],single_seed[6],single_changed[6];
        size_t single_positions[6],i,split,single_psi;
        double single_omega;
        SubseqHashTables tables;
        SubseqHashEditResult results[3];
        for (i=0;i<n;++i) s[i]=(uint8_t)((trial*7+i*5+i*i)%4);
        assert(!subseqhash_tables_init(&tables,n,d,UINT64_C(50000)+trial));
        assert(!subseqhash_compute_all_edit_splits(s,n,total_edits,
                                                   &tables,results));
        for (split=0;split<=total_edits;++split) {
            const size_t es=split,ed=total_edits-split,k=n-ed;
            SubseqHashTables view=tables;
            view.k=k;
            assert(results[split].deletions==ed);
            assert(results[split].substitutions==es);
            assert(results[split].length==k);
            assert(!subseqhash_compute_edits(s,n,ed,es,&view,single_seed,
                                             single_positions,single_changed,
                                             &single_psi,&single_omega));
            assert(results[split].psi==single_psi);
            assert(results[split].omega==single_omega);
            assert(!memcmp(results[split].seed,single_seed,k));
            assert(!memcmp(results[split].positions,single_positions,
                           k*sizeof(*single_positions)));
            assert(!memcmp(results[split].substituted,single_changed,k));
        }
        subseqhash_edit_results_destroy(results,total_edits+1);
        subseqhash_tables_destroy(&tables);
    }
    puts("all subseqhash edit tests passed");
    return 0;
}
