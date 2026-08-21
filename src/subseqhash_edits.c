#include "subseqhash_edits.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static size_t add_mod(size_t x,size_t y,size_t modulus) {
    return y>=modulus-x?y-(modulus-x):x+y;
}

static size_t subtract_mod(size_t x,size_t y,size_t modulus) {
    return x>=y?x-y:x+(modulus-y);
}

int subseqhash_compute_edits(const uint8_t *s,size_t n,
                             size_t deletions,size_t substitutions,
                             const SubseqHashTables *t,uint8_t *seed,
                             size_t *positions,uint8_t *substituted,
                             size_t *psi,double *omega) {
    size_t k,d,e_count,states,parent_states,l,j,e,q,r,a;
    double *scores=NULL,*mn,*mx;
    uint8_t *parent=NULL;
    int result=-1;

    if (!t || !psi || !omega || !t->d || deletions>n ||
        t->k!=n-deletions || substitutions>t->k ||
        (n&&!s) || (t->k&&!seed)) return -1;
    k=t->k; d=t->d; e_count=substitutions+1;
    if (!e_count || k+1==0 || k+1>SIZE_MAX/e_count ||
        (k+1)*e_count>SIZE_MAX/d) return -1;
    states=(k+1)*e_count*d;
    if (states>SIZE_MAX/(2*sizeof(*scores))) return -1;
    if (k && (n>SIZE_MAX/k || n*k>SIZE_MAX/e_count ||
              n*k*e_count>SIZE_MAX/d)) return -1;
    parent_states=n*k*e_count*d;

    scores=malloc(2*states*sizeof(*scores));
    parent=calloc(parent_states?parent_states:1,1);
    if (!scores || !parent) goto done;
    mn=scores; mx=mn+states;
    for (q=0;q<states;++q) { mn[q]=DBL_MAX; mx[q]=-DBL_MAX; }
    mn[0]=mx[0]=0.0;

#define STATE(J,E,R) ((((J)*e_count)+(E))*d+(R))
#define PARENT(L,J,E,R) ((((((L)-1)*k+((J)-1))*e_count)+(E))*d+(R))
    for (l=1;l<=n;++l) {
        const size_t original=s[l-1],upper=l<k?l:k;
        if (original>3) goto done;
        for (j=upper;j>0;--j) {
            const size_t max_e=j<substitutions?j:substitutions;
            for (e=max_e+1;e>0;--e) {
                const size_t used=e-1;
                double *dst_mn=mn+STATE(j,used,0);
                double *dst_mx=mx+STATE(j,used,0);
                uint8_t *decision=parent+PARENT(l,j,used,0);
                for (a=0;a<4;++a) {
                    const size_t changed=a!=original;
                    const size_t c=t->C[4*(j-1)+a];
                    const size_t table_base=((j-1)*4+a)*d;
                    const double *src_mn,*src_mx;
                    if (changed>used) continue;
                    src_mn=mn+STATE(j-1,used-changed,0);
                    src_mx=mx+STATE(j-1,used-changed,0);
                    for (q=0;q<d;++q) {
                        size_t x; unsigned code,nibble;
                        double lo,hi,add;
                        if (src_mn[q]==DBL_MAX) continue;
                        r=add_mod(q,c,d); x=table_base+r;
                        add=t->A[x]*(double)t->B2[x];
                        nibble=decision[r];
                        if (t->B1[x]>0) {
                            lo=src_mn[q]+add; hi=src_mx[q]+add;
                            if (lo<dst_mn[r]) {
                                dst_mn[r]=lo; code=1u+(unsigned)a*2u;
                                nibble=(nibble&0xf0u)|code;
                            }
                            if (hi>dst_mx[r]) {
                                dst_mx[r]=hi; code=2u+(unsigned)a*2u;
                                nibble=(nibble&0x0fu)|(code<<4);
                            }
                        } else {
                            lo=-src_mx[q]+add; hi=-src_mn[q]+add;
                            if (lo<dst_mn[r]) {
                                dst_mn[r]=lo; code=2u+(unsigned)a*2u;
                                nibble=(nibble&0xf0u)|code;
                            }
                            if (hi>dst_mx[r]) {
                                dst_mx[r]=hi; code=1u+(unsigned)a*2u;
                                nibble=(nibble&0x0fu)|(code<<4);
                            }
                        }
                        decision[r]=(uint8_t)nibble;
                    }
                }
            }
        }
    }

    for (r=0;r<d && mn[STATE(k,substitutions,r)]==DBL_MAX;++r) {}
    if (r==d) goto done;
    *psi=r;
    {
        int type=fabs(mn[STATE(k,substitutions,r)])>=
                 fabs(mx[STATE(k,substitutions,r)])?1:2;
        *omega=type==1?mn[STATE(k,substitutions,r)]:
                       mx[STATE(k,substitutions,r)];
        l=n; j=k; e=substitutions;
        while (j) {
            const uint8_t packed=parent[PARENT(l,j,e,r)];
            unsigned code=type==1?(packed&15u):(packed>>4);
            if (!code) { --l; continue; }
            --code;
            a=code/2u; type=(int)(code%2u)+1;
            seed[j-1]=(uint8_t)a;
            if (positions) positions[j-1]=l-1;
            if (substituted) substituted[j-1]=(uint8_t)(a!=s[l-1]);
            if (a!=s[l-1]) --e;
            r=subtract_mod(r,t->C[4*(j-1)+a],d);
            --j; --l;
        }
    }
    result=0;
done:
    free(parent); free(scores);
    return result;
#undef PARENT
#undef STATE
}

void subseqhash_edit_results_destroy(SubseqHashEditResult *results,size_t count) {
    size_t i;
    if (!results) return;
    for (i=0;i<count;++i) {
        free(results[i].substituted);
        free(results[i].positions);
        free(results[i].seed);
        memset(results+i,0,sizeof(*results));
    }
}

int subseqhash_compute_all_edit_splits(const uint8_t *s,size_t n,
                                       size_t total_edits,
                                       const SubseqHashTables *t,
                                       SubseqHashEditResult *results) {
    size_t d,e_count,band,parent_states,l,h,e,q,r,a,result_index;
    double *scores=NULL,*old_mn,*old_mx,*new_mn,*new_mx;
    uint8_t *parent=NULL;
    int status=-1;

    if (!t || !results || !t->d || t->k!=n || total_edits>n ||
        (n&&!s)) return -1;
    d=t->d; e_count=total_edits+1;
    if (!e_count || e_count>SIZE_MAX/sizeof(*results) ||
        e_count>SIZE_MAX/e_count || e_count*e_count>SIZE_MAX/d) return -1;
    band=e_count*e_count*d;
    if (band>SIZE_MAX/(4*sizeof(*scores)) ||
        (n && band>SIZE_MAX/n)) return -1;
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
    old_mn[0]=old_mx[0]=0.0;

#define BAND_STATE(H,E,R) ((((H)*e_count)+(E))*d+(R))
#define ALL_PARENT(L,H,E,R) ((((L)-1)*band)+BAND_STATE(H,E,R))
    for (l=1;l<=n;++l) {
        const size_t original=s[l-1];
        const size_t previous_max_h=(l-1)<total_edits?(l-1):total_edits;
        if (original>3) goto done;
        for (q=0;q<band;++q) { new_mn[q]=DBL_MAX; new_mx[q]=-DBL_MAX; }

        /* Delete s[l-1].  This is the unique deletion predecessor. */
        for (h=0;h<previous_max_h+1 && h<total_edits;++h)
            for (e=0;e<e_count;++e) {
                const size_t src=BAND_STATE(h,e,0);
                const size_t dst=BAND_STATE(h+1,e,0);
                memcpy(new_mn+dst,old_mn+src,d*sizeof(*new_mn));
                memcpy(new_mx+dst,old_mx+src,d*sizeof(*new_mx));
            }

        /* Retain s[l-1], either unchanged or changed to one of 3 symbols. */
        for (h=0;h<=previous_max_h;++h) {
            const size_t row=(l-1)-h;
            for (e=0;e<e_count;++e) {
                for (a=0;a<4;++a) {
                    const size_t changed=a!=original;
                    const size_t c=t->C[4*row+a];
                    const size_t table_base=(row*4+a)*d;
                    const double *src_mn,*src_mx;
                    double *dst_mn,*dst_mx;
                    uint8_t *decision;
                    if (e+changed>=e_count) continue;
                    src_mn=old_mn+BAND_STATE(h,e,0);
                    src_mx=old_mx+BAND_STATE(h,e,0);
                    dst_mn=new_mn+BAND_STATE(h,e+changed,0);
                    dst_mx=new_mx+BAND_STATE(h,e+changed,0);
                    decision=parent+ALL_PARENT(l,h,e+changed,0);
                    for (q=0;q<d;++q) {
                        size_t x; unsigned code,nibble;
                        double lo,hi,add;
                        if (src_mn[q]==DBL_MAX) continue;
                        r=add_mod(q,c,d); x=table_base+r;
                        add=t->A[x]*(double)t->B2[x];
                        nibble=decision[r];
                        if (t->B1[x]>0) {
                            lo=src_mn[q]+add; hi=src_mx[q]+add;
                            if (lo<dst_mn[r]) {
                                dst_mn[r]=lo; code=1u+(unsigned)a*2u;
                                nibble=(nibble&0xf0u)|code;
                            }
                            if (hi>dst_mx[r]) {
                                dst_mx[r]=hi; code=2u+(unsigned)a*2u;
                                nibble=(nibble&0x0fu)|(code<<4);
                            }
                        } else {
                            lo=-src_mx[q]+add; hi=-src_mn[q]+add;
                            if (lo<dst_mn[r]) {
                                dst_mn[r]=lo; code=2u+(unsigned)a*2u;
                                nibble=(nibble&0xf0u)|code;
                            }
                            if (hi>dst_mx[r]) {
                                dst_mx[r]=hi; code=1u+(unsigned)a*2u;
                                nibble=(nibble&0x0fu)|(code<<4);
                            }
                        }
                        decision[r]=(uint8_t)nibble;
                    }
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
        int type;
        out->deletions=deletions; out->substitutions=substitutions;
        out->length=k;
        out->seed=malloc(k?k:1);
        out->positions=malloc((k?k:1)*sizeof(*out->positions));
        out->substituted=malloc(k?k:1);
        if (!out->seed || !out->positions || !out->substituted) goto done;
        for (r=0;r<d && old_mn[BAND_STATE(deletions,substitutions,r)]==DBL_MAX;++r) {}
        if (r==d) goto done;
        out->psi=r;
        type=fabs(old_mn[BAND_STATE(deletions,substitutions,r)])>=
             fabs(old_mx[BAND_STATE(deletions,substitutions,r)])?1:2;
        out->omega=type==1?old_mn[BAND_STATE(deletions,substitutions,r)]:
                           old_mx[BAND_STATE(deletions,substitutions,r)];
        l=n; h=deletions; e=substitutions; q=k;
        while (q) {
            const uint8_t packed=parent[ALL_PARENT(l,h,e,r)];
            unsigned code=type==1?(packed&15u):(packed>>4);
            if (!code) { --h; --l; continue; }
            --code; a=code/2u; type=(int)(code%2u)+1;
            out->seed[q-1]=(uint8_t)a;
            out->positions[q-1]=l-1;
            out->substituted[q-1]=(uint8_t)(a!=s[l-1]);
            if (a!=s[l-1]) --e;
            r=subtract_mod(r,t->C[4*(q-1)+a],d);
            --q; --l;
        }
    }
    status=0;
done:
    if (status) subseqhash_edit_results_destroy(results,e_count);
    free(parent); free(scores);
    return status;
#undef ALL_PARENT
#undef BAND_STATE
}
