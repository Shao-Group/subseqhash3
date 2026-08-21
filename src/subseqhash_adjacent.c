#include "subseqhash_adjacent.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { PREVIOUS_COUNT = 5, SYMBOL_COUNT = 4 };

static uint64_t next_random(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
static size_t abc_index(const SubseqHashAdjacentTables *t,size_t i,
                        size_t previous,size_t symbol,size_t residue) {
    return (((i*PREVIOUS_COUNT+previous)*SYMBOL_COUNT+symbol)*t->d)+residue;
}
static size_t c_index(size_t i,size_t previous,size_t symbol) {
    return (i*PREVIOUS_COUNT+previous)*SYMBOL_COUNT+symbol;
}
static size_t add_mod(size_t x,size_t y,size_t modulus) {
    return y>=modulus-x?y-(modulus-x):x+y;
}
static size_t subtract_mod(size_t x,size_t y,size_t modulus) {
    return x>=y?x-y:x+(modulus-y);
}

int subseqhash_adjacent_tables_init(SubseqHashAdjacentTables *t,size_t k,
                                    size_t d,uint64_t seed) {
    static const int8_t b1[4]={1,1,-1,-1},b2[4]={1,-1,1,-1};
    size_t contexts,cells,i,p,r,a;
    if (!t || !d || k>SIZE_MAX/PREVIOUS_COUNT ||
        k*PREVIOUS_COUNT>SIZE_MAX/SYMBOL_COUNT ||
        k*PREVIOUS_COUNT*SYMBOL_COUNT>SIZE_MAX/d) return -1;
    contexts=k*PREVIOUS_COUNT*SYMBOL_COUNT;
    cells=contexts*d;
    if (cells>SIZE_MAX/sizeof(*t->A) ||
        contexts>SIZE_MAX/sizeof(*t->C)) return -1;
    memset(t,0,sizeof(*t)); t->k=k; t->d=d;
    t->A=malloc((cells?cells:1)*sizeof(*t->A));
    t->B1=malloc(cells?cells:1); t->B2=malloc(cells?cells:1);
    t->C=malloc((contexts?contexts:1)*sizeof(*t->C));
    if (!t->A || !t->B1 || !t->B2 || !t->C) {
        subseqhash_adjacent_tables_destroy(t); return -1;
    }
    for (i=0;i<cells;++i)
        t->A[i]=1073741824.0+
            (double)(next_random(&seed)&UINT64_C(0xffffffff))/
            4294967296.0*1073741824.0;
    for (i=0;i<k;++i) for (p=0;p<PREVIOUS_COUNT;++p)
        for (r=0;r<d;++r) {
            unsigned permutation[4]={0,1,2,3};
            for (a=4;a>1;--a) {
                size_t q=(size_t)(next_random(&seed)%a);
                unsigned tmp=permutation[a-1];
                permutation[a-1]=permutation[q]; permutation[q]=tmp;
            }
            for (a=0;a<SYMBOL_COUNT;++a) {
                size_t x=abc_index(t,i,p,a,r);
                t->B1[x]=b1[permutation[a]];
                t->B2[x]=b2[permutation[a]];
            }
        }
    for (i=0;i<k;++i) for (p=0;p<PREVIOUS_COUNT;++p)
        for (a=0;a<SYMBOL_COUNT;++a) {
            size_t value;
            do { value=(size_t)(next_random(&seed)%d); }
            while (d>=4 && a &&
                   (value==t->C[c_index(i,p,0)] ||
                    (a>1&&value==t->C[c_index(i,p,1)]) ||
                    (a>2&&value==t->C[c_index(i,p,2)])));
            t->C[c_index(i,p,a)]=value;
        }
    return 0;
}

void subseqhash_adjacent_tables_destroy(SubseqHashAdjacentTables *t) {
    if (!t) return;
    free(t->C); free(t->B2); free(t->B1); free(t->A);
    memset(t,0,sizeof(*t));
}

int subseqhash_adjacent_score(const uint8_t *seed,size_t length,
                              const SubseqHashAdjacentTables *t,
                              size_t *psi,double *omega) {
    size_t i,r=0,p=SUBSEQHASH_ADJACENT_BOS;
    double w=0;
    if (!t || !psi || !omega || !t->d || length>t->k ||
        (length&&!seed)) return -1;
    for (i=0;i<length;++i) {
        size_t a=seed[i],x;
        if (a>=SYMBOL_COUNT) return -1;
        r=add_mod(r,t->C[c_index(i,p,a)],t->d);
        x=abc_index(t,i,p,a,r);
        w=w*(double)t->B1[x]+t->A[x]*(double)t->B2[x];
        p=a;
    }
    *psi=r; *omega=w; return 0;
}

int subseqhash_adjacent_compute_all_edit_splits(
    const uint8_t *s,size_t n,size_t total_edits,
    const SubseqHashAdjacentTables *t,SubseqHashEditResult *results) {
    size_t d,e_count,band,parent_states,l,h,e,q,r,a,p,result_index;
    double *scores=NULL,*old_mn,*old_mx,*new_mn,*new_mx;
    uint8_t *parent=NULL;
    int status=-1;
    if (!t || !results || !t->d || t->k<n || total_edits>n ||
        (n&&!s)) return -1;
    d=t->d; e_count=total_edits+1;
    if (!e_count || e_count>SIZE_MAX/sizeof(*results) ||
        e_count>SIZE_MAX/e_count ||
        e_count*e_count>SIZE_MAX/PREVIOUS_COUNT ||
        e_count*e_count*PREVIOUS_COUNT>SIZE_MAX/d) return -1;
    band=e_count*e_count*PREVIOUS_COUNT*d;
    if (band>SIZE_MAX/(4*sizeof(*scores)) ||
        (n&&band>SIZE_MAX/n)) return -1;
    parent_states=n*band;
    memset(results,0,e_count*sizeof(*results));
    scores=malloc(4*band*sizeof(*scores));
    parent=calloc(parent_states?parent_states:1,1);
    if (!scores || !parent) goto done;
    old_mn=scores; old_mx=old_mn+band;
    new_mn=old_mx+band; new_mx=new_mn+band;
    for (q=0;q<band;++q) {
        old_mn[q]=new_mn[q]=DBL_MAX;
        old_mx[q]=new_mx[q]=-DBL_MAX;
    }
#define ADJ_STATE(H,E,P,R) (((((H)*e_count+(E))*PREVIOUS_COUNT+(P))*d)+(R))
#define ADJ_PARENT(L,H,E,P,R) ((((L)-1)*band)+ADJ_STATE(H,E,P,R))
    old_mn[ADJ_STATE(0,0,SUBSEQHASH_ADJACENT_BOS,0)]=0;
    old_mx[ADJ_STATE(0,0,SUBSEQHASH_ADJACENT_BOS,0)]=0;
    for (l=1;l<=n;++l) {
        const size_t original=s[l-1];
        const size_t previous_max_h=(l-1)<total_edits?(l-1):total_edits;
        if (original>=SYMBOL_COUNT) goto done;
        for (q=0;q<band;++q) { new_mn[q]=DBL_MAX; new_mx[q]=-DBL_MAX; }
        for (h=0;h<=previous_max_h && h<total_edits;++h)
            for (e=0;e+h+1<=total_edits;++e)
                for (p=0;p<PREVIOUS_COUNT;++p) {
                size_t src=ADJ_STATE(h,e,p,0),dst=ADJ_STATE(h+1,e,p,0);
                memcpy(new_mn+dst,old_mn+src,d*sizeof(*new_mn));
                memcpy(new_mx+dst,old_mx+src,d*sizeof(*new_mx));
            }
        for (h=0;h<=previous_max_h;++h) {
            const size_t row=(l-1)-h;
            for (e=0;e+h<=total_edits;++e)
                for (p=0;p<PREVIOUS_COUNT;++p)
                for (a=0;a<SYMBOL_COUNT;++a) {
                    const size_t changed=a!=original;
                    const size_t c=t->C[c_index(row,p,a)];
                    const double *src_mn,*src_mx;
                    double *dst_mn,*dst_mx;
                    uint8_t *decision;
                    if (e+h+changed>total_edits) continue;
                    src_mn=old_mn+ADJ_STATE(h,e,p,0);
                    src_mx=old_mx+ADJ_STATE(h,e,p,0);
                    dst_mn=new_mn+ADJ_STATE(h,e+changed,a,0);
                    dst_mx=new_mx+ADJ_STATE(h,e+changed,a,0);
                    decision=parent+ADJ_PARENT(l,h,e+changed,a,0);
                    for (q=0;q<d;++q) {
                        size_t x; unsigned code,nibble; double lo,hi,add;
                        if (src_mn[q]==DBL_MAX) continue;
                        r=add_mod(q,c,d); x=abc_index(t,row,p,a,r);
                        add=t->A[x]*(double)t->B2[x];
                        nibble=decision[r];
                        if (t->B1[x]>0) {
                            lo=src_mn[q]+add; hi=src_mx[q]+add;
                            if (lo<dst_mn[r]) {
                                dst_mn[r]=lo; code=1u+(unsigned)p*2u;
                                nibble=(nibble&0xf0u)|code;
                            }
                            if (hi>dst_mx[r]) {
                                dst_mx[r]=hi; code=2u+(unsigned)p*2u;
                                nibble=(nibble&0x0fu)|(code<<4);
                            }
                        } else {
                            lo=-src_mx[q]+add; hi=-src_mn[q]+add;
                            if (lo<dst_mn[r]) {
                                dst_mn[r]=lo; code=2u+(unsigned)p*2u;
                                nibble=(nibble&0xf0u)|code;
                            }
                            if (hi>dst_mx[r]) {
                                dst_mx[r]=hi; code=1u+(unsigned)p*2u;
                                nibble=(nibble&0x0fu)|(code<<4);
                            }
                        }
                        decision[r]=(uint8_t)nibble;
                    }
                }
        }
        { double *tmp=old_mn; old_mn=new_mn; new_mn=tmp;
          tmp=old_mx; old_mx=new_mx; new_mx=tmp; }
    }
    for (result_index=0;result_index<e_count;++result_index) {
        SubseqHashEditResult *out=results+result_index;
        const size_t substitutions=result_index;
        const size_t deletions=total_edits-substitutions;
        const size_t k=n-deletions;
        size_t last=SUBSEQHASH_ADJACENT_BOS;
        int type=1,found=0;
        out->deletions=deletions; out->substitutions=substitutions;
        out->length=k;
        out->seed=malloc(k?k:1);
        out->positions=malloc((k?k:1)*sizeof(*out->positions));
        out->substituted=malloc(k?k:1);
        if (!out->seed || !out->positions || !out->substituted) goto done;
        for (r=0;r<d && !found;++r)
            for (p=0;p<PREVIOUS_COUNT;++p)
                if (old_mn[ADJ_STATE(deletions,substitutions,p,r)]!=DBL_MAX)
                    { found=1; break; }
        if (!found) goto done;
        --r;
        for (p=0;p<PREVIOUS_COUNT;++p) {
            const double lo=old_mn[ADJ_STATE(deletions,substitutions,p,r)];
            const double hi=old_mx[ADJ_STATE(deletions,substitutions,p,r)];
            if (lo!=DBL_MAX &&
                (last==SUBSEQHASH_ADJACENT_BOS ||
                 fabs(lo)>fabs(out->omega))) {
                out->omega=lo; last=p; type=1;
            }
            if (hi!=-DBL_MAX &&
                (last==SUBSEQHASH_ADJACENT_BOS ||
                 fabs(hi)>fabs(out->omega))) {
                out->omega=hi; last=p; type=2;
            }
        }
        out->psi=r; l=n; h=deletions; e=substitutions; q=k;
        while (q) {
            const uint8_t packed=parent[ADJ_PARENT(l,h,e,last,r)];
            unsigned code=type==1?(packed&15u):(packed>>4);
            size_t previous;
            if (!code) { --h; --l; continue; }
            --code; previous=code/2u; type=(int)(code%2u)+1;
            out->seed[q-1]=(uint8_t)last;
            out->positions[q-1]=l-1;
            out->substituted[q-1]=(uint8_t)(last!=s[l-1]);
            if (last!=s[l-1]) --e;
            r=subtract_mod(r,t->C[c_index(q-1,previous,last)],d);
            last=previous; --q; --l;
        }
    }
    status=0;
done:
    if (status) subseqhash_edit_results_destroy(results,e_count);
    free(parent); free(scores); return status;
#undef ADJ_PARENT
#undef ADJ_STATE
}

int subseqhash_adjacent_compute_edits(
    const uint8_t *s,size_t n,size_t deletions,size_t substitutions,
    const SubseqHashAdjacentTables *t,uint8_t *seed,size_t *positions,
    uint8_t *substituted,size_t *psi,double *omega) {
    size_t total;
    SubseqHashEditResult *results;
    SubseqHashEditResult *chosen;
    int status=-1;
    if (!t || !psi || !omega || deletions>n ||
        substitutions>n-deletions ||
        substitutions>=SIZE_MAX-deletions ||
        (n-deletions&&!seed)) return -1;
    total=deletions+substitutions;
    results=calloc(total+1,sizeof(*results));
    if (!results) return -1;
    if (subseqhash_adjacent_compute_all_edit_splits(s,n,total,t,results))
        goto done;
    chosen=results+substitutions;
    memcpy(seed,chosen->seed,chosen->length);
    if (positions)
        memcpy(positions,chosen->positions,chosen->length*sizeof(*positions));
    if (substituted)
        memcpy(substituted,chosen->substituted,chosen->length);
    *psi=chosen->psi; *omega=chosen->omega; status=0;
done:
    subseqhash_edit_results_destroy(results,total+1);
    free(results); return status;
}
