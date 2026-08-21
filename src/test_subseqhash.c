#include "subseqhash.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { int found; size_t psi; double omega; uint8_t *seed; } Best;
static void enumerate(const uint8_t *s,size_t n,const SubseqHashTables *t,size_t start,size_t chosen,uint8_t *candidate,Best *best) {
    size_t p,psi; double omega;
    if (chosen==t->k) { assert(!subseqhash_score(candidate,t,&psi,&omega));
        if (!best->found || psi<best->psi || (psi==best->psi && (fabs(omega)>fabs(best->omega) || (fabs(omega)==fabs(best->omega)&&memcmp(candidate,best->seed,t->k)<0)))) {
            best->found=1; best->psi=psi; best->omega=omega; memcpy(best->seed,candidate,t->k);
        } return;
    }
    for (p=start;p+(t->k-chosen)<=n;++p) { candidate[chosen]=s[p]; enumerate(s,n,t,p+1,chosen+1,candidate,best); }
}
int main(void) {
    size_t trial;
    for (trial=0;trial<100;++trial) {
        size_t n=5+trial%4,k=1+trial%4,d=4+trial%7,i,psi; double omega,check_omega;
        uint8_t s[8],actual[4],fast[4],candidate[4],expected[4]; size_t positions[4],fast_positions[4],check_psi,fast_psi; double fast_omega; Best best={0,0,0,expected}; SubseqHashTables tables;
        for (i=0;i<n;++i) s[i]=(uint8_t)((trial*17+i*5)%4);
        assert(!subseqhash_tables_init(&tables,k,d,UINT64_C(1234)+trial));
        enumerate(s,n,&tables,0,0,candidate,&best);
        assert(!subseqhash_compute(s,n,&tables,actual,positions,&psi,&omega));
        assert(!subseqhash_compute_optimized(s,n,&tables,fast,fast_positions,&fast_psi,&fast_omega));
        assert(psi==best.psi); assert(fabs(omega)==fabs(best.omega)); assert(!memcmp(actual,best.seed,k));
        assert(fast_psi==psi); assert(fast_omega==omega); assert(!memcmp(fast,actual,k));
        assert(!memcmp(fast_positions,positions,k*sizeof(*positions)));
        assert(!subseqhash_score(actual,&tables,&check_psi,&check_omega)); assert(check_psi==psi); assert(check_omega==omega);
        for (i=0;i<k;++i) { assert(positions[i]<n); assert(actual[i]==s[positions[i]]); if(i) assert(positions[i-1]<positions[i]); }
        subseqhash_tables_destroy(&tables);
    }
    puts("all subseqhash tests passed"); return 0;
}
