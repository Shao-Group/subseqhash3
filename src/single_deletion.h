#ifndef SSH3_SINGLE_DELETION_H
#define SSH3_SINGLE_DELETION_H
#include <stddef.h>
#include <stdint.h>
#include <string.h>
/* Exact one-deletion, zero-substitution neighborhoods. Enumerating their n
   members avoids allocating/sweeping D residues. Inward prefixes and suffixes
   are shared across every requested pivot. */
typedef struct {
    const int64_t *wide_A;
    const int8_t *small_A,*B1;
    const size_t *C;
} SingleSide;
typedef struct { size_t psi,deleted; int64_t omega; } SingleChoice;
static size_t single_add(size_t a,size_t b,size_t d) {
    return a>=d-b?a-(d-b):a+b;
}
static int64_t single_A(const SingleSide *side,size_t at) {
    return side->wide_A?side->wide_A[at]:side->small_A[at];
}
static void single_step(const SingleSide *side,size_t d,size_t j,unsigned base,
    size_t *residue,int64_t *score) {
    const size_t row=4*j+base;
    *residue=single_add(*residue,side->C[row],d);
    const size_t at=row*d+*residue;
    *score=side->B1[at]*(*score)+single_A(side,at);
}
static void single_consider(SingleChoice *best,size_t p,size_t del,size_t psi,
    int64_t omega,int inward) {
    int better=psi<best->psi || (psi==best->psi && omega>best->omega);
    if (psi==best->psi && omega==best->omega) {
        const int h=del<=p,oldh=best->deleted<=p;
        if (h<oldh) better=1;
        else if (h==oldh) {
            /* Strict DP ties delete the last-consumed side base first. */
            const int descending=inward?h:!h;
            better=descending?del>best->deleted:del<best->deleted;
        }
    }
    if (better) { best->psi=psi; best->omega=omega; best->deleted=del; }
}
static void single_deletion_choices(const uint8_t *s,size_t n,size_t d,
    const SingleSide *left,const SingleSide *right,
    const int64_t *wide_pivot,const int8_t *small_pivot,
    const int8_t *left_sign,const int8_t *right_sign,const size_t *pivot_C,
    const size_t *pivots,size_t count,int inward,SingleChoice *choices) {
    const size_t k=n-1;
    for (size_t a=0;a<count;++a) { choices[a].psi=d; choices[a].omega=0; choices[a].deleted=n; }
    if (inward) {
        size_t lr[257],rr[257],base_lr[257],base_rr[257];
        int64_t lv[257],rv[257],base_lv[257],base_rv[257];
        size_t maxleft=0,maxright=0;
        for (size_t a=0;a<count;++a) {
            const size_t p=pivots?pivots[a]:a;
            if (p>maxleft) maxleft=p;
            if (k-p-1>maxright) maxright=k-p-1;
        }
        base_lr[0]=base_rr[0]=0; base_lv[0]=base_rv[0]=0;
        for (size_t j=0;j<maxleft;++j) {
            base_lr[j+1]=base_lr[j]; base_lv[j+1]=base_lv[j];
            single_step(left,d,j,s[j],base_lr+j+1,base_lv+j+1);
        }
        for (size_t j=0;j<maxright;++j) {
            base_rr[j+1]=base_rr[j]; base_rv[j+1]=base_rv[j];
            single_step(right,d,j,s[n-1-j],base_rr+j+1,base_rv+j+1);
        }
        for (size_t del=0;del<n;++del) {
            memcpy(lr,base_lr,(maxleft+1)*sizeof(*lr));
            memcpy(lv,base_lv,(maxleft+1)*sizeof(*lv));
            memcpy(rr,base_rr,(maxright+1)*sizeof(*rr));
            memcpy(rv,base_rv,(maxright+1)*sizeof(*rv));
            for (size_t j=del;j<maxleft;++j) {
                lr[j+1]=lr[j]; lv[j+1]=lv[j];
                single_step(left,d,j,s[j+1],lr+j+1,lv+j+1);
            }
            for (size_t j=k-del;j<maxright;++j) {
                rr[j+1]=rr[j]; rv[j+1]=rv[j];
                single_step(right,d,j,s[n-2-j],rr+j+1,rv+j+1);
            }
            for (size_t a=0;a<count;++a) {
                const size_t p=pivots?pivots[a]:a,row=4*p+s[p+(p>=del)],r=k-p-1;
                const size_t psi=single_add(single_add(lr[p],rr[r],d),pivot_C[row],d);
                const int64_t pivot=wide_pivot?wide_pivot[row]:small_pivot[row];
                single_consider(choices+a,p,del,psi,left_sign[row]*lv[p]+right_sign[row]*rv[r]+pivot,1);
            }
        }
    } else {
        /* For an outward side, reuse its unchanged prefix before the deleted
           base. Only the suffix after the deletion needs rescoring. */
        size_t prefix_r[257]; int64_t prefix_v[257];
        for (size_t a=0;a<count;++a) {
            const size_t p=pivots?pivots[a]:a;
            for (size_t h=0;h<=1;++h) {
                const size_t anchor=p+h,length=h?p:k-p-1;
                const SingleSide *variable=h?left:right,*fixed=h?right:left;
                const size_t fixed_length=h?k-p-1:p;
                size_t fixed_r=0; int64_t fixed_v=0;
                for (size_t j=0;j<fixed_length;++j)
                    single_step(fixed,d,j,s[h?anchor+1+j:anchor-1-j],&fixed_r,&fixed_v);
                prefix_r[0]=0; prefix_v[0]=0;
                for (size_t j=0;j<length;++j) {
                    prefix_r[j+1]=prefix_r[j]; prefix_v[j+1]=prefix_v[j];
                    single_step(variable,d,j,s[h?anchor-1-j:anchor+1+j],prefix_r+j+1,prefix_v+j+1);
                }
                const size_t row=4*p+s[anchor];
                const int64_t pivot=wide_pivot?wide_pivot[row]:small_pivot[row];
                /* Update all deletion positions together at each retained rank.
                   The table row is shared, and independent witnesses fill SIMD
                   lanes without changing the seeding algorithm's thread count. */
                for (size_t j=0;j<length;++j) {
                    const size_t tr=4*j+s[h?anchor-2-j:anchor+2+j];
                    const size_t shift=variable->C[tr];
                    size_t q=0;
#if (defined(DEL_SIMD) && DEL_SIMD) || (defined(PIVOT_SIMD) && PIVOT_SIMD)
                    const __m512i shiftv=_mm512_set1_epi64((long long)shift);
                    const __m512i divisor=_mm512_set1_epi64((long long)d);
                    for (;q+8<=j+1;q+=8) {
                        __m512i r=_mm512_add_epi64(_mm512_loadu_si512((const void *)(prefix_r+q)),shiftv);
                        r=_mm512_mask_sub_epi64(r,_mm512_cmp_epu64_mask(r,divisor,_MM_CMPINT_NLT),r,divisor);
                        /* The final retained rank is unused on this side, so the
                           byte-table word gathers have at least one spare row. */
                        const __m256i signs=_mm512_i64gather_epi32(r,variable->B1+tr*d,1);
                        const __mmask8 negative=_mm256_cmp_epi32_mask(_mm256_slli_epi32(signs,24),
                            _mm256_setzero_si256(),_MM_CMPINT_LT);
                        __m512i add;
                        if (variable->wide_A) add=_mm512_i64gather_epi64(r,variable->wide_A+tr*d,8);
                        else {
                            const __m256i bytes=_mm512_i64gather_epi32(r,variable->small_A+tr*d,1);
                            add=_mm512_cvtepi32_epi64(_mm256_srai_epi32(_mm256_slli_epi32(bytes,24),24));
                        }
                        __m512i v=_mm512_loadu_si512((const void *)(prefix_v+q));
                        v=_mm512_mask_sub_epi64(v,negative,_mm512_setzero_si512(),v);
                        v=_mm512_add_epi64(v,add);
                        _mm512_storeu_si512((void *)(prefix_r+q),r);
                        _mm512_storeu_si512((void *)(prefix_v+q),v);
                    }
#endif
                    for (;q<=j;++q) {
                        prefix_r[q]=single_add(prefix_r[q],shift,d);
                        const size_t at=tr*d+prefix_r[q];
                        prefix_v[q]=variable->B1[at]*prefix_v[q]+single_A(variable,at);
                    }
                }
                for (size_t q=0;q<=length;++q) {
                    const size_t vr=prefix_r[q]; const int64_t vv=prefix_v[q];
                    const size_t del=h?anchor-1-q:anchor+1+q;
                    const size_t psi=single_add(single_add(vr,fixed_r,d),pivot_C[row],d);
                    const int64_t omega=(h?left_sign[row]:right_sign[row])*vv+
                        (h?right_sign[row]:left_sign[row])*fixed_v+pivot;
                    single_consider(choices+a,p,del,psi,omega,0);
                }
            }
        }
    }
}
#endif
