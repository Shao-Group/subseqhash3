#include "subseqhash_edits_i8.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if defined(__AVX512F__) && defined(__AVX512BW__)
#include <immintrin.h>
#define HAVE_AVX512_I8_DP 1
#else
#define HAVE_AVX512_I8_DP 0
#endif

static uint64_t next_random(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
static size_t table_index(const SubseqHashI8Tables *t,size_t i,size_t r,
                          size_t a) {
    return (i*4+a)*t->d+r;
}
static size_t subtract_mod(size_t x,size_t y,size_t modulus) {
    return x>=y?x-y:x+(modulus-y);
}

int subseqhash_i8_tables_init(SubseqHashI8Tables *t,size_t k,size_t d,
                               uint64_t seed) {
    static const int8_t b1[4]={1,1,-1,-1},b2[4]={1,-1,1,-1};
    size_t cells,mask_cells,i,r,a,tail,lane;
    if (!t || !d || k>(size_t)((INT16_MAX-1)/127) ||
        k>SIZE_MAX/d || k*d>SIZE_MAX/4 ||
        k*d*4>SIZE_MAX/sizeof(*t->A)) return -1;
    memset(t,0,sizeof(*t)); t->k=k; t->d=d; t->blocks=(d+31)/32;
    cells=k*d*4;
    if (k>SIZE_MAX/4 || k*4>SIZE_MAX/t->blocks) return -1;
    mask_cells=k*4*t->blocks;
    t->A=malloc((cells?cells:1)*sizeof(*t->A));
    t->B1=malloc(cells?cells:1); t->B2=malloc(cells?cells:1);
    t->signed_A=malloc((cells?cells:1)*sizeof(*t->signed_A));
    t->B1_masks=calloc(mask_cells?mask_cells:1,sizeof(*t->B1_masks));
    t->wrap_indices=malloc(32*32*sizeof(*t->wrap_indices));
    t->C=malloc((k?k*4:1)*sizeof(*t->C));
    if (!t->A || !t->B1 || !t->B2 || !t->signed_A ||
        !t->B1_masks || !t->wrap_indices || !t->C) {
        subseqhash_i8_tables_destroy(t); return -1;
    }
    for (tail=0;tail<32;++tail) for (lane=0;lane<32;++lane)
        t->wrap_indices[tail*32+lane]=(uint16_t)(
            lane<tail?lane:32+lane-tail);
    for (i=0;i<cells;++i)
        t->A[i]=(int8_t)(1u+(unsigned)(next_random(&seed)%127u));
    for (i=0;i<k;++i) for (r=0;r<d;++r) {
        unsigned permutation[4]={0,1,2,3};
        for (a=4;a>1;--a) {
            size_t q=(size_t)(next_random(&seed)%a);
            unsigned tmp=permutation[a-1];
            permutation[a-1]=permutation[q]; permutation[q]=tmp;
        }
        for (a=0;a<4;++a) {
            size_t x=table_index(t,i,r,a);
            t->B1[x]=b1[permutation[a]];
            t->B2[x]=b2[permutation[a]];
            t->signed_A[x]=(int8_t)(t->A[x]*t->B2[x]);
            if (t->B1[x]>0)
                t->B1_masks[(i*4+a)*t->blocks+r/32]|=
                    (uint32_t)(UINT32_C(1)<<(r%32));
        }
    }
    for (i=0;i<k;++i) for (a=0;a<4;++a) {
        size_t value;
        do { value=(size_t)(next_random(&seed)%d); }
        while (d>=4 && a &&
               (value==t->C[4*i] || (a>1&&value==t->C[4*i+1]) ||
                (a>2&&value==t->C[4*i+2])));
        t->C[4*i+a]=value;
    }
    return 0;
}

void subseqhash_i8_tables_destroy(SubseqHashI8Tables *t) {
    if (!t) return;
    free(t->C); free(t->wrap_indices); free(t->B1_masks); free(t->signed_A);
    free(t->B2); free(t->B1); free(t->A);
    memset(t,0,sizeof(*t));
}

int subseqhash_i8_tables_make_double(const SubseqHashI8Tables *s,
                                      SubseqHashTables *d) {
    size_t cells,i;
    if (!s || !d || !s->d || s->k>SIZE_MAX/s->d ||
        s->k*s->d>SIZE_MAX/4) return -1;
    cells=s->k*s->d*4;
    memset(d,0,sizeof(*d)); d->k=s->k; d->d=s->d;
    d->A=malloc((cells?cells:1)*sizeof(*d->A));
    d->B1=malloc(cells?cells:1); d->B2=malloc(cells?cells:1);
    d->C=malloc((s->k?s->k*4:1)*sizeof(*d->C));
    if (!d->A || !d->B1 || !d->B2 || !d->C) {
        subseqhash_tables_destroy(d); return -1;
    }
    for (i=0;i<cells;++i) d->A[i]=(double)s->A[i];
    memcpy(d->B1,s->B1,cells); memcpy(d->B2,s->B2,cells);
    memcpy(d->C,s->C,s->k*4*sizeof(*d->C));
    return 0;
}

#if !HAVE_AVX512_I8_DP || !defined(__AVX512VL__)
static void save_parent(uint8_t *decision,size_t r,unsigned a,
                        unsigned min_type,unsigned max_type) {
    unsigned value=decision[r];
    if (min_type) value=(value&0xf0u)|(min_type+(unsigned)a*2u);
    if (max_type) value=(value&0x0fu)|
        ((max_type+(unsigned)a*2u)<<4);
    decision[r]=(uint8_t)value;
}
#endif

#if !HAVE_AVX512_I8_DP || !defined(__AVX512VL__)
static void update_scalar(const int16_t *src_mn,const int16_t *src_mx,
                          int16_t *dst_mn,int16_t *dst_mx,
                          const int8_t *A,const int8_t *B1,
                          const int8_t *B2,uint8_t *decision,
                          size_t q0,size_t r0,size_t count,unsigned symbol) {
    size_t i;
    for (i=0;i<count;++i) {
        size_t q=q0+i,r=r0+i;
        int16_t lo,hi,add;
        unsigned min_type=0,max_type=0;
        if (src_mn[q]==INT16_MAX) continue;
        add=(int16_t)A[r]*(int16_t)B2[r];
        if (B1[r]>0) {
            lo=src_mn[q]+add; hi=src_mx[q]+add;
            if (lo<dst_mn[r]) { dst_mn[r]=lo; min_type=1; }
            if (hi>dst_mx[r]) { dst_mx[r]=hi; max_type=2; }
        } else {
            lo=-src_mx[q]+add; hi=-src_mn[q]+add;
            if (lo<dst_mn[r]) { dst_mn[r]=lo; min_type=2; }
            if (hi>dst_mx[r]) { dst_mx[r]=hi; max_type=1; }
        }
        if (min_type || max_type)
            save_parent(decision,r,symbol,min_type,max_type);
    }
}
#endif

#if HAVE_AVX512_I8_DP
static void update_parent_vector(uint8_t *decision,size_t offset,
                                 __mmask32 active,__mmask32 positive,
                                 __mmask32 minmask,__mmask32 maxmask,
                                 unsigned symbol) {
#if defined(__AVX512VL__)
    const __m256i low_keep=_mm256_set1_epi8((char)0xf0);
    const __m256i high_keep=_mm256_set1_epi8(0x0f);
    const __m256i min_negative=_mm256_set1_epi8(
        (char)(2u+symbol*2u));
    const __m256i min_positive=_mm256_set1_epi8(
        (char)(1u+symbol*2u));
    const __m256i max_negative=min_positive;
    const __m256i max_positive=min_negative;
    __m256i value=_mm256_maskz_loadu_epi8(active,decision+offset);
    __m256i code=_mm256_mask_blend_epi8(
        positive,min_negative,min_positive);
    __m256i updated=_mm256_or_si256(_mm256_and_si256(value,low_keep),code);
    value=_mm256_mask_mov_epi8(value,minmask,updated);
    code=_mm256_mask_blend_epi8(positive,max_negative,max_positive);
    code=_mm256_slli_epi16(code,4);
    updated=_mm256_or_si256(_mm256_and_si256(value,high_keep),code);
    value=_mm256_mask_mov_epi8(value,maxmask,updated);
    _mm256_mask_storeu_epi8(decision+offset,active,value);
#else
    while (minmask || maxmask) {
        unsigned lane=(unsigned)__builtin_ctz(
            (unsigned)(minmask|maxmask));
        unsigned min_type=(minmask&(1u<<lane))?
            ((positive&(1u<<lane))?1u:2u):0u;
        unsigned max_type=(maxmask&(1u<<lane))?
            ((positive&(1u<<lane))?2u:1u):0u;
        save_parent(decision,offset+lane,symbol,min_type,max_type);
        minmask=(__mmask32)(minmask&~(1u<<lane));
        maxmask=(__mmask32)(maxmask&~(1u<<lane));
    }
    (void)active;
#endif
}

#if defined(__AVX512VL__)
static inline __m512i load_wrapped_i16(const int16_t *source,size_t d,
                                       size_t q,const uint16_t *indices) {
    if (q+32<=d)
        return _mm512_loadu_si512((const void *)(source+q));
    {
        const unsigned tail=(unsigned)(d-q);
        const __mmask32 low=(__mmask32)((UINT32_C(1)<<tail)-1u);
        const __m512i end=_mm512_maskz_loadu_epi16(low,source+q);
        const __m512i start=_mm512_loadu_si512((const void *)source);
        const __m512i index=_mm512_loadu_si512(
            (const void *)(indices+tail*32));
        return _mm512_permutex2var_epi16(end,index,start);
    }
}

static void update_avx512_aligned(
    const int16_t *src_mn,const int16_t *src_mx,
    int16_t *dst_mn,int16_t *dst_mx,const int8_t *signed_A,
    const uint32_t *B1_masks,const uint16_t *wrap_indices,
    uint8_t *decision,size_t d,size_t shift,unsigned symbol) {
    const __m512i infinity=_mm512_set1_epi16(INT16_MAX);
    const __m512i zero=_mm512_setzero_si512();
    size_t r;
    for (r=0;r<d;r+=32) {
        const size_t q=r>=shift?r-shift:r+d-shift;
        const __m512i smn=load_wrapped_i16(
            src_mn,d,q,wrap_indices);
        const __m512i smx=load_wrapped_i16(
            src_mx,d,q,wrap_indices);
        const __m512i dmn=_mm512_loadu_si512((const void *)(dst_mn+r));
        const __m512i dmx=_mm512_loadu_si512((const void *)(dst_mx+r));
        const __m256i a8=_mm256_loadu_si256(
            (const __m256i *)(signed_A+r));
        const __m512i add=_mm512_cvtepi8_epi16(a8);
        const __mmask32 positive=(__mmask32)B1_masks[r/32];
        const __mmask32 reachable=
            _mm512_cmpneq_epi16_mask(smn,infinity);
        const __m512i lo=_mm512_add_epi16(
            _mm512_mask_blend_epi16(positive,
                _mm512_sub_epi16(zero,smx),smn),add);
        const __m512i hi=_mm512_add_epi16(
            _mm512_mask_blend_epi16(positive,
                _mm512_sub_epi16(zero,smn),smx),add);
        const __mmask32 minmask=reachable&
            _mm512_cmplt_epi16_mask(lo,dmn);
        const __mmask32 maxmask=reachable&
            _mm512_cmpgt_epi16_mask(hi,dmx);
        _mm512_mask_storeu_epi16(dst_mn+r,minmask,lo);
        _mm512_mask_storeu_epi16(dst_mx+r,maxmask,hi);
        if (minmask || maxmask)
            update_parent_vector(decision,r,(__mmask32)UINT32_MAX,
                                 positive,minmask,maxmask,symbol);
    }
}
#endif

static void update_avx512_block(
    const int16_t *src_mn,const int16_t *src_mx,
    int16_t *dst_mn,int16_t *dst_mx,
    const int8_t *A,const int8_t *B1,const int8_t *B2,
    uint8_t *decision,size_t q,size_t r,size_t count,unsigned symbol) {
    const __mmask32 active=count==32?(__mmask32)UINT32_MAX:
        (__mmask32)((1u<<(unsigned)count)-1u);
    (void)B2;
    const __m512i infinity=_mm512_set1_epi16(INT16_MAX);
    __m512i smn=_mm512_maskz_loadu_epi16(active,src_mn+q);
    __m512i smx=_mm512_maskz_loadu_epi16(active,src_mx+q);
    __m512i dmn=_mm512_maskz_loadu_epi16(active,dst_mn+r);
    __m512i dmx=_mm512_maskz_loadu_epi16(active,dst_mx+r);
    __m256i a16;
    __m256i b1_8;
#if defined(__AVX512VL__)
    a16=_mm256_maskz_loadu_epi8(active,A+r);
    b1_8=_mm256_maskz_loadu_epi8(active,B1+r);
#else
    if (count<32) {
        update_scalar(src_mn,src_mx,dst_mn,dst_mx,A,B1,B2,decision,
                      q,r,count,symbol);
        return;
    }
    a16=_mm256_loadu_si256((const __m256i *)(A+r));
    b1_8=_mm256_loadu_si256((const __m256i *)(B1+r));
#endif
    {
        __m512i add=_mm512_cvtepi8_epi16(a16);
        __m512i b1v=_mm512_cvtepi8_epi16(b1_8);
        __mmask32 reachable=active&_mm512_cmpneq_epi16_mask(smn,infinity);
        __mmask32 positive=active&_mm512_cmpgt_epi16_mask(
            b1v,_mm512_setzero_si512());
        __m512i lo=_mm512_add_epi16(
            _mm512_mask_blend_epi16(positive,
                _mm512_sub_epi16(_mm512_setzero_si512(),smx),smn),add);
        __m512i hi=_mm512_add_epi16(
            _mm512_mask_blend_epi16(positive,
                _mm512_sub_epi16(_mm512_setzero_si512(),smn),smx),add);
        __mmask32 minmask=reachable&_mm512_cmplt_epi16_mask(lo,dmn);
        __mmask32 maxmask=reachable&_mm512_cmpgt_epi16_mask(hi,dmx);
        _mm512_mask_storeu_epi16(dst_mn+r,minmask,lo);
        _mm512_mask_storeu_epi16(dst_mx+r,maxmask,hi);
        if (minmask || maxmask)
            update_parent_vector(decision,r,active,positive,
                                 minmask,maxmask,symbol);
    }
}

static void update_avx512(const int16_t *src_mn,const int16_t *src_mx,
                          int16_t *dst_mn,int16_t *dst_mx,
                          const int8_t *A,const int8_t *B1,
                          const int8_t *B2,uint8_t *decision,
                          size_t q0,size_t r0,size_t count,unsigned symbol) {
    size_t i=0;
    for (;i+32<=count;i+=32)
        update_avx512_block(src_mn,src_mx,dst_mn,dst_mx,A,B1,B2,
                            decision,q0+i,r0+i,32,symbol);
    if (i<count)
        update_avx512_block(src_mn,src_mx,dst_mn,dst_mx,A,B1,B2,
                            decision,q0+i,r0+i,count-i,symbol);
}
#endif

static void update_shifted(const int16_t *src_mn,const int16_t *src_mx,
                           int16_t *dst_mn,int16_t *dst_mx,
                           const int8_t *A,const int8_t *B1,
                           const int8_t *B2,const int8_t *signed_A,
                           const uint32_t *B1_masks,const uint16_t *wrap_indices,
                           uint8_t *decision,
                           size_t d,size_t shift,unsigned symbol,int use_aligned) {
    (void)signed_A;
    (void)B1_masks;
    (void)wrap_indices;
    (void)use_aligned;
#if HAVE_AVX512_I8_DP && defined(__AVX512VL__)
    if (use_aligned && !(d&31u)) {
        update_avx512_aligned(src_mn,src_mx,dst_mn,dst_mx,signed_A,
                              B1_masks,wrap_indices,decision,d,shift,symbol);
        return;
    }
#endif
#if HAVE_AVX512_I8_DP
#define UPDATE_RANGE update_avx512
    A=signed_A;
#else
#define UPDATE_RANGE update_scalar
#endif
    const size_t first=d-shift;
    UPDATE_RANGE(src_mn,src_mx,dst_mn,dst_mx,A,B1,B2,decision,
                 0,shift,first,symbol);
    if (shift)
        UPDATE_RANGE(src_mn,src_mx,dst_mn,dst_mx,A,B1,B2,decision,
                     first,0,shift,symbol);
#undef UPDATE_RANGE
}

static int compute_all_edit_splits_i8_impl(
    const uint8_t *s,size_t n,size_t total_edits,
    const SubseqHashI8Tables *t,SubseqHashEditResult *results,int use_aligned) {
    size_t d,e_count,band,parent_states,l,h,e,q,r,a,result_index;
    int16_t *scores=NULL,*old_mn,*old_mx,*new_mn,*new_mx;
    uint8_t *parent=NULL;
    int status=-1;
    if (!t || !results || !t->d || t->k<n || total_edits>n ||
        (n&&!s)) return -1;
    d=t->d; e_count=total_edits+1;
    if (!e_count || e_count>SIZE_MAX/sizeof(*results) ||
        e_count>SIZE_MAX/e_count || e_count*e_count>SIZE_MAX/d) return -1;
    band=e_count*e_count*d;
    if (band>SIZE_MAX/(4*sizeof(*scores)) || (n&&band>SIZE_MAX/n))
        return -1;
    parent_states=n*band;
    memset(results,0,e_count*sizeof(*results));
    scores=malloc(4*band*sizeof(*scores));
    parent=calloc(parent_states?parent_states:1,1);
    if (!scores || !parent) goto done;
    old_mn=scores; old_mx=old_mn+band;
    new_mn=old_mx+band; new_mx=new_mn+band;
    for (q=0;q<band;++q) {
        old_mn[q]=new_mn[q]=INT16_MAX;
        old_mx[q]=new_mx[q]=INT16_MIN;
    }
#define I8_STATE(H,E,R) ((((H)*e_count+(E))*d)+(R))
#define I8_PARENT(L,H,E,R) ((((L)-1)*band)+I8_STATE(H,E,R))
    old_mn[0]=old_mx[0]=0;
    for (l=1;l<=n;++l) {
        const size_t original=s[l-1];
        const size_t previous_max_h=(l-1)<total_edits?(l-1):total_edits;
        if (original>3) goto done;
        for (q=0;q<band;++q) {
            new_mn[q]=INT16_MAX; new_mx[q]=INT16_MIN;
        }
        for (h=0;h<=previous_max_h && h<total_edits;++h)
            for (e=0;e+h+1<=total_edits;++e) {
                size_t src=I8_STATE(h,e,0),dst=I8_STATE(h+1,e,0);
                memcpy(new_mn+dst,old_mn+src,d*sizeof(*new_mn));
                memcpy(new_mx+dst,old_mx+src,d*sizeof(*new_mx));
            }
        for (h=0;h<=previous_max_h;++h) {
            const size_t row=(l-1)-h;
            for (e=0;e+h<=total_edits;++e) for (a=0;a<4;++a) {
                const size_t changed=a!=original;
                const size_t table_base=(row*4+a)*d;
                if (e+h+changed>total_edits) continue;
                update_shifted(
                    old_mn+I8_STATE(h,e,0),
                    old_mx+I8_STATE(h,e,0),
                    new_mn+I8_STATE(h,e+changed,0),
                    new_mx+I8_STATE(h,e+changed,0),
                    t->A+table_base,t->B1+table_base,t->B2+table_base,
                    t->signed_A+table_base,
                    t->B1_masks+(row*4+a)*t->blocks,
                    t->wrap_indices,
                    parent+I8_PARENT(l,h,e+changed,0),d,
                    t->C[4*row+a],(unsigned)a,use_aligned);
            }
        }
        { int16_t *tmp=old_mn; old_mn=new_mn; new_mn=tmp;
          tmp=old_mx; old_mx=new_mx; new_mx=tmp; }
    }
    for (result_index=0;result_index<e_count;++result_index) {
        SubseqHashEditResult *out=results+result_index;
        size_t substitutions=result_index,deletions=total_edits-substitutions;
        size_t k=n-deletions;
        int type;
        out->deletions=deletions; out->substitutions=substitutions;
        out->length=k; out->seed=malloc(k?k:1);
        out->positions=malloc((k?k:1)*sizeof(*out->positions));
        out->substituted=malloc(k?k:1);
        if (!out->seed || !out->positions || !out->substituted) goto done;
        for (r=0;r<d &&
             old_mn[I8_STATE(deletions,substitutions,r)]==INT16_MAX;++r) {}
        if (r==d) goto done;
        out->psi=r;
        type=llabs((long long)old_mn[I8_STATE(deletions,substitutions,r)])>=
             llabs((long long)old_mx[I8_STATE(deletions,substitutions,r)])
             ?1:2;
        out->omega=(double)(type==1?
            old_mn[I8_STATE(deletions,substitutions,r)]:
            old_mx[I8_STATE(deletions,substitutions,r)]);
        l=n; h=deletions; e=substitutions; q=k;
        while (q) {
            uint8_t packed=parent[I8_PARENT(l,h,e,r)];
            unsigned code=type==1?(packed&15u):(packed>>4);
            if (!code) { --h; --l; continue; }
            --code; a=code/2u; type=(int)(code%2u)+1;
            out->seed[q-1]=(uint8_t)a; out->positions[q-1]=l-1;
            out->substituted[q-1]=(uint8_t)(a!=s[l-1]);
            if (a!=s[l-1]) --e;
            r=subtract_mod(r,t->C[4*(q-1)+a],d);
            --q; --l;
        }
    }
    status=0;
done:
    if (status) subseqhash_edit_results_destroy(results,e_count);
    free(parent); free(scores); return status;
#undef I8_PARENT
#undef I8_STATE
}

int subseqhash_compute_all_edit_splits_i8(
    const uint8_t *s,size_t n,size_t total_edits,
    const SubseqHashI8Tables *t,SubseqHashEditResult *results) {
    return compute_all_edit_splits_i8_impl(s,n,total_edits,t,results,1);
}

int subseqhash_compute_all_edit_splits_i8_legacy(
    const uint8_t *s,size_t n,size_t total_edits,
    const SubseqHashI8Tables *t,SubseqHashEditResult *results) {
    return compute_all_edit_splits_i8_impl(s,n,total_edits,t,results,0);
}

const char *subseqhash_edits_i8_backend(void) {
#if HAVE_AVX512_I8_DP && defined(__AVX512VL__)
    return "AVX-512BW (32 x int16 states; aligned vpermt2w)";
#elif HAVE_AVX512_I8_DP
    return "AVX-512BW (32 x int16 states; scalar parents)";
#else
    return "scalar int16";
#endif
}
