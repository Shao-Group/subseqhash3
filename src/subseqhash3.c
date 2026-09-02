#define _POSIX_C_SOURCE 200112L
#include "subseqhash3.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#if defined(__AVX512F__) && defined(__AVX512BW__)
#include <immintrin.h>
#define HAVE_AVX512_I8_DP 1
#else
#define HAVE_AVX512_I8_DP 0
#endif

#define I8_POSITIVE_INFINITY INT16_C(32639)
#define I8_NEGATIVE_INFINITY INT16_C(-32640)
#define I8_ROW_KERNEL_MAX_D 256

static void *aligned_allocate(size_t alignment,size_t bytes) {
    void *pointer=NULL;
    if (posix_memalign(&pointer,alignment,bytes?bytes:alignment)) return NULL;
    return pointer;
}

static uint64_t next_random(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
static size_t table_index(const SubseqHash3Tables *t,size_t i,size_t r,
                          size_t a) {
    return (i*4+a)*t->d+r;
}
static size_t subtract_mod(size_t x,size_t y,size_t modulus) {
    return x>=y?x-y:x+(modulus-y);
}
static size_t triangular_state(size_t edit_count,size_t deletions,
                               size_t substitutions,size_t d) {
    return (deletions*edit_count-
            deletions*(deletions-1)/2+substitutions)*d;
}

int subseqhash3_tables_init(SubseqHash3Tables *t,size_t k,size_t d,
                               uint64_t seed) {
    static const int8_t b1[4]={1,1,-1,-1},b2[4]={1,-1,1,-1};
    size_t cells,mask_cells,f_cells,i,r,a,tail,lane;
    if (!t || !d || k>(size_t)((I8_POSITIVE_INFINITY-1)/127) ||
        k>SIZE_MAX/d || k*d>SIZE_MAX/4 ||
        k*d*4>SIZE_MAX/sizeof(*t->A)) return -1;
    memset(t,0,sizeof(*t)); t->k=k; t->d=d; t->blocks=(d+31)/32;
    if (k>SIZE_MAX/128 || 128*k==SIZE_MAX) return -1;
    t->f_stride=128*k+1;
    if (d>SIZE_MAX/t->f_stride) return -1;
    f_cells=d*t->f_stride;
    cells=k*d*4;
    if (k>SIZE_MAX/4 || k*4>SIZE_MAX/t->blocks) return -1;
    mask_cells=k*4*t->blocks;
    t->A=malloc((cells?cells:1)*sizeof(*t->A));
    t->B1=malloc(cells?cells:1); t->B2=malloc(cells?cells:1);
    t->signed_A=aligned_allocate(64,(cells?cells:1)*sizeof(*t->signed_A));
    t->B1_masks=calloc(mask_cells?mask_cells:1,sizeof(*t->B1_masks));
    t->wrap_indices=malloc(32*32*sizeof(*t->wrap_indices));
    t->C=malloc((k?k*4:1)*sizeof(*t->C));
    t->F=malloc((f_cells?f_cells:1)*sizeof(*t->F));
    if (!t->A || !t->B1 || !t->B2 || !t->signed_A ||
        !t->B1_masks || !t->wrap_indices || !t->C || !t->F) {
        subseqhash3_tables_destroy(t); return -1;
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
    for (r=0;r<d;++r) {
        float value=0.0f;
        for (i=0;i<t->f_stride;++i) {
            value+=(float)((double)next_random(&seed)/(double)UINT64_MAX);
            t->F[r*t->f_stride+i]=value;
        }
    }
    return 0;
}

void subseqhash3_tables_destroy(SubseqHash3Tables *t) {
    if (!t) return;
    free(t->F); free(t->C); free(t->wrap_indices); free(t->B1_masks);
    free(t->signed_A);
    free(t->B2); free(t->B1); free(t->A);
    memset(t,0,sizeof(*t));
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
        if (src_mn[q]==I8_POSITIVE_INFINITY) continue;
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
        if (decision && (min_type || max_type))
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
    uint8_t *decision,size_t d,size_t shift,unsigned symbol,
    const int16_t *initial_mn,const int16_t *initial_mx,
    int initialize) {
    const __m512i infinity=_mm512_set1_epi16(I8_POSITIVE_INFINITY);
    const __m512i negative_infinity=
        _mm512_set1_epi16(I8_NEGATIVE_INFINITY);
    const __m512i zero=_mm512_setzero_si512();
    size_t r;
    for (r=0;r<d;r+=32) {
        const size_t q=r>=shift?r-shift:r+d-shift;
        const __m512i smn=load_wrapped_i16(
            src_mn,d,q,wrap_indices);
        const __m512i smx=load_wrapped_i16(
            src_mx,d,q,wrap_indices);
        const __m512i dmn=initialize?
            (initial_mn?_mm512_loadu_si512((const void *)(initial_mn+r)):
                        infinity):
            _mm512_loadu_si512((const void *)(dst_mn+r));
        const __m512i dmx=initialize?
            (initial_mx?_mm512_loadu_si512((const void *)(initial_mx+r)):
                        negative_infinity):
            _mm512_loadu_si512((const void *)(dst_mx+r));
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
        if (initialize) {
            _mm512_storeu_si512((void *)(dst_mn+r),
                _mm512_mask_mov_epi16(dmn,minmask,lo));
            _mm512_storeu_si512((void *)(dst_mx+r),
                _mm512_mask_mov_epi16(dmx,maxmask,hi));
        } else {
            _mm512_mask_storeu_epi16(dst_mn+r,minmask,lo);
            _mm512_mask_storeu_epi16(dst_mx+r,maxmask,hi);
        }
        if (decision && (minmask || maxmask))
            update_parent_vector(decision,r,(__mmask32)UINT32_MAX,
                                 positive,minmask,maxmask,symbol);
    }
}

static void update_avx512_substitutions(
    const int16_t *src_mn,const int16_t *src_mx,
    int16_t *dst_mn,int16_t *dst_mx,
    const int8_t *const signed_A[3],const uint32_t *const B1_masks[3],
    const uint16_t *wrap_indices,uint8_t *decision,size_t d,
    const size_t shifts[3],const unsigned symbols[3],
    const int16_t *initial_mn,const int16_t *initial_mx,
    int initialize) {
    const __m512i infinity=_mm512_set1_epi16(I8_POSITIVE_INFINITY);
    const __m512i negative_infinity=
        _mm512_set1_epi16(I8_NEGATIVE_INFINITY);
    const __m512i zero=_mm512_setzero_si512();
    const __m256i low_keep=_mm256_set1_epi8((char)0xf0);
    const __m256i high_keep=_mm256_set1_epi8(0x0f);
    size_t r;
    for (r=0;r<d;r+=32) {
        __m512i dmn=initialize?
            (initial_mn?_mm512_loadu_si512((const void *)(initial_mn+r)):
                        infinity):
            _mm512_loadu_si512((const void *)(dst_mn+r));
        __m512i dmx=initialize?
            (initial_mx?_mm512_loadu_si512((const void *)(initial_mx+r)):
                        negative_infinity):
            _mm512_loadu_si512((const void *)(dst_mx+r));
        __m256i parent=_mm256_loadu_si256(
            (const __m256i *)(decision+r));
        unsigned candidate;
#pragma GCC unroll 3
        for (candidate=0;candidate<3;++candidate) {
            const size_t shift=shifts[candidate];
            const size_t q=r>=shift?r-shift:r+d-shift;
            const __m512i smn=load_wrapped_i16(
                src_mn,d,q,wrap_indices);
            const __m512i smx=load_wrapped_i16(
                src_mx,d,q,wrap_indices);
            const __m256i a8=_mm256_loadu_si256(
                (const __m256i *)(signed_A[candidate]+r));
            const __m512i add=_mm512_cvtepi8_epi16(a8);
            const __mmask32 positive=(__mmask32)
                B1_masks[candidate][r/32];
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
            const unsigned base=symbols[candidate]*2u;
            __m256i code=_mm256_mask_blend_epi8(
                positive,_mm256_set1_epi8((char)(base+2u)),
                _mm256_set1_epi8((char)(base+1u)));
            __m256i updated=_mm256_or_si256(
                _mm256_and_si256(parent,low_keep),code);
            parent=_mm256_mask_mov_epi8(parent,minmask,updated);
            code=_mm256_mask_blend_epi8(
                positive,_mm256_set1_epi8((char)(base+1u)),
                _mm256_set1_epi8((char)(base+2u)));
            code=_mm256_slli_epi16(code,4);
            updated=_mm256_or_si256(
                _mm256_and_si256(parent,high_keep),code);
            parent=_mm256_mask_mov_epi8(parent,maxmask,updated);
            dmn=_mm512_mask_mov_epi16(dmn,minmask,lo);
            dmx=_mm512_mask_mov_epi16(dmx,maxmask,hi);
        }
        _mm512_storeu_si512((void *)(dst_mn+r),dmn);
        _mm512_storeu_si512((void *)(dst_mx+r),dmx);
        _mm256_storeu_si256((__m256i *)(decision+r),parent);
    }
}

static void update_avx512_destination(
    const int16_t *const src_mn[4],const int16_t *const src_mx[4],
    int16_t *dst_mn,int16_t *dst_mx,
    const int8_t *const signed_A[4],const uint32_t *const B1_masks[4],
    uint8_t *decision,size_t d,
    const size_t shifts[4],const unsigned symbols[4],unsigned count,
    const int16_t *initial_mn,const int16_t *initial_mx) {
    const __m512i infinity=_mm512_set1_epi16(I8_POSITIVE_INFINITY);
    const __m512i negative_infinity=
        _mm512_set1_epi16(I8_NEGATIVE_INFINITY);
    const __m256i high_mask=_mm256_set1_epi8((char)0xf0);
    size_t r;
    for (r=0;r<d;r+=32) {
        __m512i dmn=initial_mn?
            _mm512_load_si512((const void *)(initial_mn+r)):infinity;
        __m512i dmx=initial_mx?
            _mm512_load_si512((const void *)(initial_mx+r)):
            negative_infinity;
        __m256i parent_low=_mm256_setzero_si256();
        __m256i parent_high=_mm256_setzero_si256();
        unsigned candidate;
#pragma GCC unroll 4
        for (candidate=0;candidate<count;++candidate) {
            const size_t shift=shifts[candidate];
            const size_t q=r>=shift?r-shift:r+d-shift;
            const __m512i smn=_mm512_loadu_si512(
                (const void *)(src_mn[candidate]+q));
            const __m512i smx=_mm512_loadu_si512(
                (const void *)(src_mx[candidate]+q));
            const __m256i a8=_mm256_load_si256(
                (const __m256i *)(signed_A[candidate]+r));
            const __m512i add=_mm512_cvtepi8_epi16(a8);
            const __mmask32 positive=(__mmask32)
                B1_masks[candidate][r/32];
            const __mmask32 reachable=
                _mm512_cmpneq_epi16_mask(smn,infinity);
            const __m512i lo=_mm512_mask_add_epi16(
                _mm512_sub_epi16(add,smx),positive,smn,add);
            const __m512i hi=_mm512_add_epi16(
                lo,_mm512_sub_epi16(smx,smn));
            const __mmask32 minmask=reachable&
                _mm512_cmplt_epi16_mask(lo,dmn);
            const __mmask32 maxmask=reachable&
                _mm512_cmpgt_epi16_mask(hi,dmx);
            const unsigned base=symbols[candidate]*2u;
            const __m256i code=_mm256_mask_blend_epi8(positive,
                _mm256_set1_epi8((char)(
                    (base+2u)|((base+1u)<<4))),
                _mm256_set1_epi8((char)(
                    (base+1u)|((base+2u)<<4))));
            parent_low=_mm256_mask_mov_epi8(parent_low,minmask,code);
            parent_high=_mm256_mask_mov_epi8(parent_high,maxmask,code);
            dmn=_mm512_mask_mov_epi16(dmn,minmask,lo);
            dmx=_mm512_mask_mov_epi16(dmx,maxmask,hi);
        }
        _mm512_store_si512((void *)(dst_mn+r),dmn);
        _mm512_store_si512((void *)(dst_mx+r),dmx);
        _mm256_store_si256((__m256i *)(decision+r),
            _mm256_ternarylogic_epi32(
                parent_low,parent_high,high_mask,0xd8));
        if (!r) {
            _mm512_store_si512((void *)(dst_mn+d),dmn);
            _mm512_store_si512((void *)(dst_mx+d),dmx);
        }
    }
}

static void update_avx512_row(
    const int16_t *old_mn,const int16_t *old_mx,
    int16_t *new_mn,int16_t *new_mx,uint8_t *parent,
    const SubseqHash3Tables *t,size_t score_stride,size_t state_count,
    size_t h,size_t row,size_t min_e,size_t max_e,unsigned original) {
    const __m512i infinity=_mm512_set1_epi16(I8_POSITIVE_INFINITY);
    const __m512i negative_infinity=
        _mm512_set1_epi16(I8_NEGATIVE_INFINITY);
    const __m256i high_mask=_mm256_set1_epi8((char)0xf0);
    const size_t d=t->d;
    const size_t state_base=triangular_state(state_count,h,0,score_stride);
    const size_t parent_base=triangular_state(state_count,h,0,d);
    const size_t deletion_base=h?
        triangular_state(state_count,h-1,0,score_stride):0;
    const size_t table_row=row*4;
    const size_t original_table=(table_row+original)*d;
    const size_t original_shift=t->C[table_row+original];
    const int8_t *substitution_A[3];
    const uint32_t *substitution_masks[3];
    size_t substitution_shifts[3];
    unsigned substitution_symbols[3];
    unsigned symbol,index=0;
    size_t r,e;
    for (symbol=0;symbol<4;++symbol) if (symbol!=original) {
        const size_t base=(table_row+symbol)*d;
        substitution_A[index]=t->signed_A+base;
        substitution_masks[index]=t->B1_masks+
            (table_row+symbol)*t->blocks;
        substitution_shifts[index]=t->C[table_row+symbol];
        substitution_symbols[index]=symbol;
        ++index;
    }
    for (r=0;r<d;r+=32) {
        __m512i substitution_add[3];
        __mmask32 substitution_positive[3];
        const __m512i original_add=_mm512_cvtepi8_epi16(
            _mm256_load_si256((const __m256i *)(t->signed_A+
                                                original_table+r)));
        const __mmask32 original_positive=(__mmask32)
            t->B1_masks[(table_row+original)*t->blocks+r/32];
        const unsigned original_base=original*2u;
        const __m256i original_code=_mm256_mask_blend_epi8(
            original_positive,
            _mm256_set1_epi8((char)(
                (original_base+2u)|((original_base+1u)<<4))),
            _mm256_set1_epi8((char)(
                (original_base+1u)|((original_base+2u)<<4))));
        unsigned candidate;
#pragma GCC unroll 3
        for (candidate=0;candidate<3;++candidate) {
            substitution_add[candidate]=_mm512_cvtepi8_epi16(
                _mm256_load_si256((const __m256i *)(
                    substitution_A[candidate]+r)));
            substitution_positive[candidate]=(__mmask32)
                substitution_masks[candidate][r/32];
        }
        for (e=min_e;e<=max_e;++e) {
            int16_t *destination_mn=new_mn+state_base+e*score_stride;
            int16_t *destination_mx=new_mx+state_base+e*score_stride;
            uint8_t *destination_parent=parent+parent_base+e*d;
            __m512i dmn=h?_mm512_load_si512((const void *)(
                old_mn+deletion_base+e*score_stride+r)):infinity;
            __m512i dmx=h?_mm512_load_si512((const void *)(
                old_mx+deletion_base+e*score_stride+r)):negative_infinity;
            __m256i parent_low=_mm256_setzero_si256();
            __m256i parent_high=_mm256_setzero_si256();
            if (e) {
                const int16_t *source_mn=old_mn+state_base+
                    (e-1)*score_stride;
                const int16_t *source_mx=old_mx+state_base+
                    (e-1)*score_stride;
#pragma GCC unroll 3
                for (candidate=0;candidate<3;++candidate) {
                    const size_t shift=substitution_shifts[candidate];
                    const size_t q=r>=shift?r-shift:r+d-shift;
                    const __m512i smn=_mm512_loadu_si512(
                        (const void *)(source_mn+q));
                    const __m512i smx=_mm512_loadu_si512(
                        (const void *)(source_mx+q));
                    const __mmask32 reachable=
                        _mm512_cmpneq_epi16_mask(smn,infinity);
                    const __m512i lo=_mm512_mask_add_epi16(
                        _mm512_sub_epi16(substitution_add[candidate],smx),
                        substitution_positive[candidate],smn,
                        substitution_add[candidate]);
                    const __m512i hi=_mm512_add_epi16(
                        lo,_mm512_sub_epi16(smx,smn));
                    const __mmask32 minmask=reachable&
                        _mm512_cmplt_epi16_mask(lo,dmn);
                    const __mmask32 maxmask=reachable&
                        _mm512_cmpgt_epi16_mask(hi,dmx);
                    const unsigned base=substitution_symbols[candidate]*2u;
                    const __m256i code=_mm256_mask_blend_epi8(
                        substitution_positive[candidate],
                        _mm256_set1_epi8((char)(
                            (base+2u)|((base+1u)<<4))),
                        _mm256_set1_epi8((char)(
                            (base+1u)|((base+2u)<<4))));
                    parent_low=_mm256_mask_mov_epi8(
                        parent_low,minmask,code);
                    parent_high=_mm256_mask_mov_epi8(
                        parent_high,maxmask,code);
                    dmn=_mm512_mask_mov_epi16(dmn,minmask,lo);
                    dmx=_mm512_mask_mov_epi16(dmx,maxmask,hi);
                }
            }
            if (e<=row) {
                const int16_t *source_mn=old_mn+state_base+e*score_stride;
                const int16_t *source_mx=old_mx+state_base+e*score_stride;
                const size_t q=r>=original_shift?
                    r-original_shift:r+d-original_shift;
                const __m512i smn=_mm512_loadu_si512(
                    (const void *)(source_mn+q));
                const __m512i smx=_mm512_loadu_si512(
                    (const void *)(source_mx+q));
                const __mmask32 reachable=
                    _mm512_cmpneq_epi16_mask(smn,infinity);
                const __m512i lo=_mm512_mask_add_epi16(
                    _mm512_sub_epi16(original_add,smx),original_positive,
                    smn,original_add);
                const __m512i hi=_mm512_add_epi16(
                    lo,_mm512_sub_epi16(smx,smn));
                const __mmask32 minmask=reachable&
                    _mm512_cmplt_epi16_mask(lo,dmn);
                const __mmask32 maxmask=reachable&
                    _mm512_cmpgt_epi16_mask(hi,dmx);
                parent_low=_mm256_mask_mov_epi8(
                    parent_low,minmask,original_code);
                parent_high=_mm256_mask_mov_epi8(
                    parent_high,maxmask,original_code);
                dmn=_mm512_mask_mov_epi16(dmn,minmask,lo);
                dmx=_mm512_mask_mov_epi16(dmx,maxmask,hi);
            }
            _mm512_store_si512((void *)(destination_mn+r),dmn);
            _mm512_store_si512((void *)(destination_mx+r),dmx);
            _mm256_store_si256((__m256i *)(destination_parent+r),
                _mm256_ternarylogic_epi32(
                    parent_low,parent_high,high_mask,0xd8));
            if (!r) {
                _mm512_store_si512((void *)(destination_mn+d),dmn);
                _mm512_store_si512((void *)(destination_mx+d),dmx);
            }
        }
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
    const __m512i infinity=_mm512_set1_epi16(I8_POSITIVE_INFINITY);
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
        if (decision && (minmask || maxmask))
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
                           size_t d,size_t shift,unsigned symbol) {
    (void)signed_A;
    (void)B1_masks;
    (void)wrap_indices;
#if HAVE_AVX512_I8_DP && defined(__AVX512VL__)
    if (!(d&31u)) {
        update_avx512_aligned(src_mn,src_mx,dst_mn,dst_mx,signed_A,
                              B1_masks,wrap_indices,decision,d,shift,symbol,
                              NULL,NULL,0);
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

typedef struct {
    int16_t *scores;
    uint8_t *parent;
    size_t score_cells,parent_cells;
} I8DPWorkspace;

static int i8_workspace_reserve(I8DPWorkspace *workspace,
                                size_t score_cells,size_t parent_cells) {
    int16_t *scores;
    uint8_t *parent;
    if (score_cells>workspace->score_cells) {
        scores=aligned_allocate(64,score_cells*sizeof(*scores));
        if (!scores) return -1;
        free(workspace->scores); workspace->scores=scores;
        workspace->score_cells=score_cells;
    }
    if (parent_cells>workspace->parent_cells) {
        parent=aligned_allocate(64,parent_cells?parent_cells:1);
        if (!parent) return -1;
        free(workspace->parent); workspace->parent=parent;
        workspace->parent_cells=parent_cells;
    }
    return 0;
}

static void i8_workspace_destroy(I8DPWorkspace *workspace) {
    free(workspace->parent); free(workspace->scores);
    memset(workspace,0,sizeof(*workspace));
}

static int compute_all_edit_splits_impl(
    const uint8_t *s,size_t n,size_t total_edits,
    const SubseqHash3Tables *t,SubseqHash3EditResult *results,
    I8DPWorkspace *workspace) {
    size_t d,e_count,state_count,score_stride,score_band,score_bytes;
    size_t band,parent_states;
    size_t l,h,e,q,r,a,result_index;
    int16_t *scores=NULL,*old_mn,*old_mx,*new_mn,*new_mx;
    uint8_t *parent=NULL;
    int status=-1;
    if (!t || !results || n>SIZE_MAX/128 || !t->d || !t->F ||
        t->f_stride<128*n+1 ||
        t->k<n || total_edits>n ||
        (n&&!s)) return -1;
    d=t->d; e_count=total_edits+1;
    if (!e_count || e_count>SIZE_MAX/sizeof(*results) ||
        e_count==SIZE_MAX || e_count>SIZE_MAX/(e_count+1) ||
        e_count*(e_count+1)/2>SIZE_MAX/d) return -1;
    state_count=e_count*(e_count+1)/2;
    score_stride=d;
#if HAVE_AVX512_I8_DP && defined(__AVX512VL__)
    if (!(d&31u)) {
        if (d>SIZE_MAX-32) return -1;
        score_stride=d+32;
    }
#endif
    if (state_count>SIZE_MAX/score_stride) return -1;
    score_band=state_count*score_stride;
    band=state_count*d;
    if (score_band>SIZE_MAX/(4*sizeof(*scores)) ||
        (n&&band>SIZE_MAX/n)) return -1;
    parent_states=n*band;
    score_bytes=4*score_band*sizeof(*scores);
    if (parent_states>SIZE_MAX-score_bytes) return -1;
    memset(results,0,e_count*sizeof(*results));
    for (l=0;l<n;++l) if (s[l]>3) return -1;
    if (workspace) {
        if (i8_workspace_reserve(
                workspace,4*score_band,parent_states)) goto done;
        scores=workspace->scores; parent=workspace->parent;
    } else {
        scores=aligned_allocate(64,score_bytes+parent_states);
        if (!scores) goto done;
        parent=(uint8_t *)scores+score_bytes;
    }
    if (score_stride==d)
        memset(parent,0,parent_states?parent_states:1);
    old_mn=scores; old_mx=old_mn+score_band;
    new_mn=old_mx+score_band; new_mx=new_mn+score_band;
    /* At layer zero only (h,e)=(0,0) is reachable. Every other state is
       initialized immediately before the first layer that can reach it. */
    memset(old_mn,0x7f,score_stride*sizeof(*old_mn));
    memset(old_mx,0x80,score_stride*sizeof(*old_mx));
#define I8_STATE(H,E,R) (triangular_state(e_count,(H),(E),d)+(R))
#define I8_SCORE_STATE(H,E,R) (triangular_state(e_count,(H),(E),score_stride)+(R))
#define I8_PARENT(L,H,E,R) ((((L)-1)*band)+I8_STATE(H,E,R))
    old_mn[0]=old_mx[0]=0;
    if (score_stride!=d) {
        memcpy(old_mn+d,old_mn,32*sizeof(*old_mn));
        memcpy(old_mx+d,old_mx,32*sizeof(*old_mx));
    }
    for (l=1;l<=n;++l) {
        const size_t original=s[l-1];
        const size_t previous_max_h=(l-1)<total_edits?(l-1):total_edits;
        const size_t deletion_h_count=previous_max_h<total_edits?
            previous_max_h+1:total_edits;
#if HAVE_AVX512_I8_DP && defined(__AVX512VL__)
        if (!(d&31u)) {
            /* Visit each destination once.  Substitution candidates precede
               the original-symbol candidate exactly as in the source-centric
               traversal, preserving strict-comparison tie breaking. */
            const size_t max_h=l<total_edits?l:total_edits;
            for (h=0;h<=max_h;++h) {
                const size_t retained=l-h;
                const size_t max_e=retained<total_edits-h?
                    retained:total_edits-h;
                const size_t remaining=n-l;
                const size_t min_e=total_edits>remaining+h?
                    total_edits-remaining-h:0;
                if (min_e>max_e) continue;
                if (!retained) {
                    const size_t src=I8_SCORE_STATE(h-1,0,0);
                    const size_t dst=I8_SCORE_STATE(h,0,0);
                    memcpy(new_mn+dst,old_mn+src,
                           score_stride*sizeof(*new_mn));
                    memcpy(new_mx+dst,old_mx+src,
                           score_stride*sizeof(*new_mx));
                    memset(parent+I8_PARENT(l,h,0,0),0,d);
                    continue;
                }
                /* Reusing table vectors across edit states wins while a row
                   remains cache-resident. Larger rows retain the existing
                   destination-major traversal for sequential score access. */
                if (d<=I8_ROW_KERNEL_MAX_D) {
                    update_avx512_row(old_mn,old_mx,new_mn,new_mx,
                        parent+(l-1)*band,t,score_stride,e_count,h,
                        retained-1,min_e,max_e,(unsigned)original);
                } else {
                    const size_t row=retained-1;
                    const int8_t *substitution_A[3];
                    const uint32_t *substitution_masks[3];
                    size_t substitution_shifts[3];
                    unsigned substitution_symbols[3],substitution_count=0;
                    const size_t original_base=(row*4+original)*d;
                    for (a=0;a<4;++a) if (a!=original) {
                        const size_t base=(row*4+a)*d;
                        substitution_A[substitution_count]=t->signed_A+base;
                        substitution_masks[substitution_count]=t->B1_masks+
                            (row*4+a)*t->blocks;
                        substitution_shifts[substitution_count]=t->C[4*row+a];
                        substitution_symbols[substitution_count]=(unsigned)a;
                        ++substitution_count;
                    }
                    for (e=min_e;e<=max_e;++e) {
                        const int16_t *candidate_mn[4],*candidate_mx[4];
                        const int8_t *candidate_A[4];
                        const uint32_t *candidate_masks[4];
                        size_t candidate_shifts[4];
                        unsigned candidate_symbols[4],count=0;
                        const int16_t *initial_mn=h?
                            old_mn+I8_SCORE_STATE(h-1,e,0):NULL;
                        const int16_t *initial_mx=h?
                            old_mx+I8_SCORE_STATE(h-1,e,0):NULL;
                        if (e) for (a=0;a<3;++a) {
                            candidate_mn[count]=
                                old_mn+I8_SCORE_STATE(h,e-1,0);
                            candidate_mx[count]=
                                old_mx+I8_SCORE_STATE(h,e-1,0);
                            candidate_A[count]=substitution_A[a];
                            candidate_masks[count]=substitution_masks[a];
                            candidate_shifts[count]=substitution_shifts[a];
                            candidate_symbols[count]=substitution_symbols[a];
                            ++count;
                        }
                        if (e<=row) {
                            candidate_mn[count]=
                                old_mn+I8_SCORE_STATE(h,e,0);
                            candidate_mx[count]=
                                old_mx+I8_SCORE_STATE(h,e,0);
                            candidate_A[count]=t->signed_A+original_base;
                            candidate_masks[count]=t->B1_masks+
                                (row*4+original)*t->blocks;
                            candidate_shifts[count]=t->C[4*row+original];
                            candidate_symbols[count]=(unsigned)original;
                            ++count;
                        }
                        if (count==4)
                            update_avx512_destination(candidate_mn,candidate_mx,
                                new_mn+I8_SCORE_STATE(h,e,0),
                                new_mx+I8_SCORE_STATE(h,e,0),candidate_A,
                                candidate_masks,
                                parent+I8_PARENT(l,h,e,0),d,candidate_shifts,
                                candidate_symbols,4,initial_mn,initial_mx);
                        else if (count==3)
                            update_avx512_destination(candidate_mn,candidate_mx,
                                new_mn+I8_SCORE_STATE(h,e,0),
                                new_mx+I8_SCORE_STATE(h,e,0),candidate_A,
                                candidate_masks,
                                parent+I8_PARENT(l,h,e,0),d,candidate_shifts,
                                candidate_symbols,3,initial_mn,initial_mx);
                        else
                            update_avx512_destination(candidate_mn,candidate_mx,
                                new_mn+I8_SCORE_STATE(h,e,0),
                                new_mx+I8_SCORE_STATE(h,e,0),candidate_A,
                                candidate_masks,
                                parent+I8_PARENT(l,h,e,0),d,candidate_shifts,
                                candidate_symbols,1,initial_mn,initial_mx);
                    }
                }
            }
            goto layer_complete;
        }
#endif
        /* Every reachable h>0 destination is initialized below by copying
           its mandatory deletion predecessor.  Only h=0 lacks such a
           predecessor and needs an infinity fill before its transitions. */
        {
            const size_t reachable_e=l<total_edits?l:total_edits;
            const size_t cells=(reachable_e+1)*d;
            memset(new_mn,0x7f,cells*sizeof(*new_mn));
            memset(new_mx,0x80,cells*sizeof(*new_mx));
        }
        for (h=0;h<deletion_h_count;++h)
            for (e=0;e+h+1<=total_edits && e<=(l-1)-h;++e) {
                size_t src=I8_SCORE_STATE(h,e,0),dst=I8_SCORE_STATE(h+1,e,0);
                memcpy(new_mn+dst,old_mn+src,d*sizeof(*new_mn));
                memcpy(new_mx+dst,old_mx+src,d*sizeof(*new_mx));
            }
        for (h=0;h<=previous_max_h;++h) {
            const size_t row=(l-1)-h;
            const size_t max_e=row<total_edits-h?row:total_edits-h;
            for (e=0;e<=max_e;++e) {
                a=original;
                {
                    const size_t table_base=(row*4+a)*d;
                    update_shifted(old_mn+I8_SCORE_STATE(h,e,0),
                        old_mx+I8_SCORE_STATE(h,e,0),
                        new_mn+I8_SCORE_STATE(h,e,0),new_mx+I8_SCORE_STATE(h,e,0),
                        t->A+table_base,t->B1+table_base,
                        t->B2+table_base,t->signed_A+table_base,
                        t->B1_masks+(row*4+a)*t->blocks,t->wrap_indices,
                        parent+I8_PARENT(l,h,e,0),d,t->C[4*row+a],
                        (unsigned)a);
                }
                if (e+h<total_edits) {
#if HAVE_AVX512_I8_DP && defined(__AVX512VL__)
                    if (!(d&31u)) {
                        const int8_t *substitution_A[3];
                        const uint32_t *substitution_masks[3];
                        size_t substitution_shifts[3];
                        unsigned substitution_symbols[3],index=0;
                        for (a=0;a<4;++a) if (a!=original) {
                            const size_t base=(row*4+a)*d;
                            substitution_A[index]=t->signed_A+base;
                            substitution_masks[index]=t->B1_masks+
                                (row*4+a)*t->blocks;
                            substitution_shifts[index]=t->C[4*row+a];
                            substitution_symbols[index++]=(unsigned)a;
                        }
                        update_avx512_substitutions(
                            old_mn+I8_SCORE_STATE(h,e,0),
                            old_mx+I8_SCORE_STATE(h,e,0),
                            new_mn+I8_SCORE_STATE(h,e+1,0),
                            new_mx+I8_SCORE_STATE(h,e+1,0),substitution_A,
                            substitution_masks,t->wrap_indices,
                            parent+I8_PARENT(l,h,e+1,0),d,
                            substitution_shifts,substitution_symbols,
                            NULL,NULL,0);
                    } else
#endif
                    for (a=0;a<4;++a) if (a!=original) {
                        const size_t table_base=(row*4+a)*d;
                        update_shifted(old_mn+I8_SCORE_STATE(h,e,0),
                            old_mx+I8_SCORE_STATE(h,e,0),
                            new_mn+I8_SCORE_STATE(h,e+1,0),
                            new_mx+I8_SCORE_STATE(h,e+1,0),
                            t->A+table_base,t->B1+table_base,
                            t->B2+table_base,t->signed_A+table_base,
                            t->B1_masks+(row*4+a)*t->blocks,
                            t->wrap_indices,parent+I8_PARENT(l,h,e+1,0),
                            d,t->C[4*row+a],(unsigned)a);
                    }
                }
            }
        }
#if HAVE_AVX512_I8_DP && defined(__AVX512VL__)
layer_complete:
#endif
        { int16_t *tmp=old_mn; old_mn=new_mn; new_mn=tmp;
          tmp=old_mx; old_mx=new_mx; new_mx=tmp; }
    }
    for (result_index=0;result_index<e_count;++result_index) {
        SubseqHash3EditResult *out=results+result_index;
        size_t substitutions=result_index,deletions=total_edits-substitutions;
        size_t k=n-deletions;
        int type=0;
        float best_transformed=0.0f;
        int have_best=0;
        out->deletions=deletions; out->substitutions=substitutions;
        out->length=k; out->seed=malloc(k?k:1);
        out->positions=malloc((k?k:1)*sizeof(*out->positions));
        out->substituted=malloc(k?k:1);
        if (!out->seed || !out->positions || !out->substituted) goto done;
        for (r=0;r<d;++r) {
            int16_t minimum=old_mn[I8_SCORE_STATE(deletions,substitutions,r)];
            int16_t maximum=old_mx[I8_SCORE_STATE(deletions,substitutions,r)];
            size_t magnitude;
            float transformed;
            int candidate_type;
            if (minimum==I8_POSITIVE_INFINITY) continue;
            candidate_type=llabs((long long)minimum)>=llabs((long long)maximum)?1:2;
            magnitude=(size_t)llabs((long long)(candidate_type==1?minimum:maximum));
            if (magnitude>=t->f_stride) goto done;
            transformed=t->F[r*t->f_stride+magnitude];
            if (!have_best || transformed>best_transformed) {
                have_best=1; best_transformed=transformed;
                out->psi=r; type=candidate_type;
            }
        }
        if (!have_best) goto done;
        r=out->psi;
        out->omega=(double)(type==1?
            old_mn[I8_SCORE_STATE(deletions,substitutions,r)]:
            old_mx[I8_SCORE_STATE(deletions,substitutions,r)]);
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
    if (status) subseqhash3_edit_results_destroy(results,e_count);
    if (!workspace) free(scores);
    return status;
#undef I8_PARENT
#undef I8_SCORE_STATE
#undef I8_STATE
}

int subseqhash3_compute_all_edit_splits(
    const uint8_t *s,size_t n,size_t total_edits,
    const SubseqHash3Tables *t,SubseqHash3EditResult *results) {
    return compute_all_edit_splits_impl(s,n,total_edits,t,results,NULL);
}

int subseqhash3_compute_long_edit_splits(
    const uint8_t *s,size_t N,size_t n,size_t total_edits,
    const SubseqHash3Tables *t,SubseqHash3EditResult *results) {
    I8DPWorkspace workspace={0};
    size_t windows,per,total,w,split,i;
    int status=-1;
    if (!t || !results || !n || n>N || total_edits>n || !t->d ||
        t->k<n || (N&&!s) || total_edits==SIZE_MAX) return -1;
    windows=N-n+1; per=total_edits+1;
    if (windows>SIZE_MAX/per) return -1;
    total=windows*per;
    if (total>SIZE_MAX/sizeof(*results)) return -1;
    memset(results,0,total*sizeof(*results));
    for (w=0;w<windows;++w) {
        SubseqHash3EditResult *window_results=results+w*per;
        if (compute_all_edit_splits_impl(
                s+w,n,total_edits,t,window_results,&workspace)) goto done;
        for (split=0;split<per;++split)
            for (i=0;i<window_results[split].length;++i)
                window_results[split].positions[i]+=w;
    }
    status=0;
done:
    if (status) subseqhash3_edit_results_destroy(results,total);
    i8_workspace_destroy(&workspace);
    return status;
}
void subseqhash3_edit_results_destroy(SubseqHash3EditResult *results,
                                      size_t count) {
    size_t i;
    if (!results) return;
    for (i=0;i<count;++i) {
        free(results[i].substituted);
        free(results[i].positions);
        free(results[i].seed);
        memset(results+i,0,sizeof(*results));
    }
}

const char *subseqhash3_backend(void) {
#if HAVE_AVX512_I8_DP && defined(__AVX512VL__)
    return "AVX-512BW (32 x int16 states; padded circular loads)";
#elif HAVE_AVX512_I8_DP
    return "AVX-512BW (32 x int16 states; scalar parents)";
#else
    return "scalar int16";
#endif
}
