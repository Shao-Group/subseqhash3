#include "subseqhash2_edits_i8.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int found; size_t psi; int omega; } Best;

static void score_side(const uint8_t *y,size_t begin,size_t end,int reverse,
                       const SubseqHash2I8Tables *t,size_t *residue,int *omega) {
    size_t q=0,r=0; int w=0;
    while (q<end-begin) {
        size_t at=reverse?end-1-q:begin+q,a=y[at],x;
        const int8_t *A=reverse?t->AR:t->AF;
        const int8_t *B1=reverse?t->BR1:t->BF1;
        const int8_t *B2=reverse?t->BR2:t->BF2;
        r=(r+(reverse?t->CR:t->CF)[4*q+a])%t->d;
        x=(q*4+a)*t->d+r; w=B1[x]*w+A[x]*B2[x]; ++q;
    }
    *residue=r; *omega=w;
}

static void consider(const uint8_t *seed,size_t k,size_t order,
                     const SubseqHash2I8Tables *t,Best *best) {
    size_t lr,rr,a=seed[order],r; int lw,rw,w;
    score_side(seed,0,order,1,t,&lr,&lw);
    score_side(seed,order+1,k,0,t,&rr,&rw);
    r=(lr+t->CP[4*order+a]+rr)%t->d;
    w=t->BP1[4*order+a]*lw+t->AP[4*order+a]+
      t->BP2[4*order+a]*rw;
    if (!best->found || r<best->psi || (r==best->psi&&w>best->omega)) {
        best->found=1; best->psi=r; best->omega=w;
    }
}

static void symbols(const uint8_t *x,const size_t *pos,size_t k,size_t row,
                    size_t changed,size_t target,size_t order,
                    const SubseqHash2I8Tables *t,uint8_t *seed,Best *best) {
    unsigned a;
    if (row==k) { if (changed==target) consider(seed,k,order,t,best); return; }
    for (a=0;a<4;++a) {
        size_t delta=(a!=x[pos[row]]);
        if (changed+delta>target) continue;
        seed[row]=(uint8_t)a;
        symbols(x,pos,k,row+1,changed+delta,target,order,t,seed,best);
    }
}

static void positions(const uint8_t *x,size_t start,size_t n,size_t k,
                      size_t at,size_t used,size_t substitutions,size_t order,
                      const SubseqHash2I8Tables *t,size_t *pos,uint8_t *seed,
                      Best *best) {
    size_t i;
    if (used==k) {
        symbols(x,pos,k,0,0,substitutions,order,t,seed,best); return;
    }
    for (i=at;i<start+n && start+n-i>=k-used;++i) {
        pos[used]=i;
        positions(x,start,n,k,i+1,used+1,substitutions,order,t,pos,seed,best);
    }
}

static void validate_result(const uint8_t *x,size_t start,size_t n,
                            const SubseqHash2I8Tables *t,size_t order,
                            const SubseqHashEditResult *r) {
    size_t i,subs=0,psi; int omega;
    for (i=0;i<r->length;++i) {
        assert(r->positions[i]>=start && r->positions[i]<start+n);
        if (i) assert(r->positions[i-1]<r->positions[i]);
        assert(r->substituted[i]==(uint8_t)(r->seed[i]!=x[r->positions[i]]));
        subs+=r->substituted[i];
    }
    assert(subs==r->substitutions);
    {
        size_t lr,rr,a=r->seed[order]; int lw,rw;
        score_side(r->seed,0,order,1,t,&lr,&lw);
        score_side(r->seed,order+1,r->length,0,t,&rr,&rw);
        psi=(lr+t->CP[4*order+a]+rr)%t->d;
        omega=t->BP1[4*order+a]*lw+t->AP[4*order+a]+
              t->BP2[4*order+a]*rw;
    }
    assert(psi==r->psi); assert(omega==(int)r->omega);
}

int main(void) {
    static const size_t residues[]={1,2,3,7,15,31,32,33,63,64,65,127};
    size_t trial;
    for (trial=0;trial<60;++trial) {
        const size_t N=5+trial%2,n=3+trial%2,E=trial%3;
        const size_t d=residues[trial%(sizeof(residues)/sizeof(residues[0]))];
        size_t windows=N-n+1,per=subseqhash2_i8_results_per_window(n,E);
        size_t total=windows*per,w,index=0,s,order,i;
        uint8_t x[6]; SubseqHash2I8Tables tables;
        SubseqHashEditResult *results=calloc(total?total:1,sizeof(*results));
        assert(results);
        for (i=0;i<N;++i) x[i]=(uint8_t)((trial*7+i*3+i*i)%4);
        assert(!subseqhash2_i8_tables_init(&tables,n,d,999+trial));
        assert(!subseqhash2_compute_long_edits_i8(x,N,n,E,&tables,results));
        for (w=0;w<windows;++w) for (s=0;s<=E;++s) {
            size_t k=n-E+s;
            for (order=0;order<k;++order,++index) {
                size_t pos[4]; uint8_t seed[4]; Best best={0,0,0};
                size_t decoded_s,decoded_order;
                assert(!subseqhash2_i8_result_coordinates(
                    n,E,index%per,&decoded_s,&decoded_order));
                assert(decoded_s==s && decoded_order==order);
                positions(x,w,n,k,w,0,s,order,&tables,pos,seed,&best);
                assert(best.found);
                assert(results[index].deletions==E-s);
                assert(results[index].substitutions==s);
                assert(results[index].length==k);
                assert(results[index].psi==best.psi);
                assert((int)results[index].omega==best.omega);
                validate_result(x,w,n,&tables,order,results+index);
            }
        }
        assert(index==total);
        subseqhash_edit_results_destroy(results,total); free(results);
        subseqhash2_i8_tables_destroy(&tables);
    }
    puts("all SubseqHash2 long edit tests passed");
    return 0;
}
