#include "subseqhash_edits_max_i8.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#if defined(__AVX512F__) && defined(__AVX512BW__)
#include <immintrin.h>
#define MAX_I8_AVX512 1
#else
#define MAX_I8_AVX512 0
#endif

static uint64_t next_random(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
static size_t subtract_mod(size_t x,size_t y,size_t d) {
    return x>=y?x-y:x+(d-y);
}

int subseqhash_max_i8_tables_init(SubseqHashMaxI8Tables *t,size_t k,
                                  size_t d,uint64_t seed) {
    size_t cells,i,a;
    if (!t || !d || k>SUBSEQHASH_MAX_I8_MAX_LENGTH ||
        k>SIZE_MAX/d || k*d>SIZE_MAX/4) return -1;
    memset(t,0,sizeof(*t)); t->k=k; t->d=d; cells=k*d*4;
    t->A=malloc(cells?cells:1);
    t->C=malloc((k?k*4:1)*sizeof(*t->C));
    if (!t->A || !t->C) { subseqhash_max_i8_tables_destroy(t); return -1; }
    for (i=0;i<cells;++i) {
        unsigned magnitude=1u+(unsigned)(next_random(&seed)%127u);
        t->A[i]=(int8_t)((next_random(&seed)&1u)?
                         (int)magnitude:-(int)magnitude);
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

void subseqhash_max_i8_tables_destroy(SubseqHashMaxI8Tables *t) {
    if (!t) return;
    free(t->C); free(t->A); memset(t,0,sizeof(*t));
}

#if !MAX_I8_AVX512 || !defined(__AVX512VL__)
static void update_scalar(const int16_t *src,int16_t *dst,const int8_t *A,
                          uint8_t *parent,size_t q,size_t r,size_t count,
                          unsigned symbol) {
    size_t i;
    for (i=0;i<count;++i) if (src[q+i]!=INT16_MIN) {
        int16_t candidate=(int16_t)(src[q+i]+A[r+i]);
        if (candidate>dst[r+i]) {
            dst[r+i]=candidate; parent[r+i]=(uint8_t)(symbol+1u);
        }
    }
}
#endif

#if MAX_I8_AVX512
static void update_vector(const int16_t *src,int16_t *dst,const int8_t *A,
                          uint8_t *parent,size_t q0,size_t r0,size_t count,
                          unsigned symbol) {
    const __m512i unreachable=_mm512_set1_epi16(INT16_MIN);
    const __m256i parent_code=_mm256_set1_epi8((char)(symbol+1u));
    size_t i=0;
    for (;i<count;i+=32) {
        size_t lanes=count-i<32?count-i:32;
        __mmask32 active=lanes==32?(__mmask32)UINT32_MAX:
            (__mmask32)((UINT32_C(1)<<lanes)-1u);
        __m512i sv=_mm512_maskz_loadu_epi16(active,src+q0+i);
        __m512i dv=_mm512_maskz_loadu_epi16(active,dst+r0+i);
        __m256i a8;
#if defined(__AVX512VL__)
        a8=_mm256_maskz_loadu_epi8(active,A+r0+i);
#else
        if (lanes<32) {
            update_scalar(src,dst,A,parent,q0+i,r0+i,lanes,symbol);
            continue;
        }
        a8=_mm256_loadu_si256((const __m256i *)(A+r0+i));
#endif
        {
            __m512i candidate=_mm512_add_epi16(
                sv,_mm512_cvtepi8_epi16(a8));
            __mmask32 reachable=active&
                _mm512_cmpneq_epi16_mask(sv,unreachable);
            __mmask32 improve=reachable&
                _mm512_cmpgt_epi16_mask(candidate,dv);
            _mm512_mask_storeu_epi16(dst+r0+i,improve,candidate);
#if defined(__AVX512VL__)
            _mm256_mask_storeu_epi8(parent+r0+i,improve,parent_code);
#else
            while (improve) {
                unsigned lane=(unsigned)__builtin_ctz((unsigned)improve);
                parent[r0+i+lane]=(uint8_t)(symbol+1u);
                improve=(__mmask32)(improve&(improve-1u));
            }
#endif
        }
    }
}
#endif

static void update_shifted(const int16_t *src,int16_t *dst,const int8_t *A,
                           uint8_t *parent,size_t d,size_t shift,
                           unsigned symbol) {
    const size_t first=d-shift;
#if MAX_I8_AVX512
#define MAX_UPDATE update_vector
#else
#define MAX_UPDATE update_scalar
#endif
    MAX_UPDATE(src,dst,A,parent,0,shift,first,symbol);
    if (shift) MAX_UPDATE(src,dst,A,parent,first,0,shift,symbol);
#undef MAX_UPDATE
}

int subseqhash_compute_all_edit_splits_max_i8(
    const uint8_t *s,size_t n,size_t total_edits,
    const SubseqHashMaxI8Tables *t,SubseqHashEditResult *results) {
    size_t d,ec,band,parent_states,l,h,e,q,r,a,result_index;
    int16_t *scores=NULL,*old,*new_scores;
    uint8_t *parent=NULL; int status=-1;
    if (!t || !results || !t->d || t->k<n || n>SUBSEQHASH_MAX_I8_MAX_LENGTH ||
        total_edits>n || (n&&!s)) return -1;
    d=t->d; ec=total_edits+1;
    if (!ec || ec>SIZE_MAX/sizeof(*results) || ec>SIZE_MAX/ec ||
        ec*ec>SIZE_MAX/d) return -1;
    band=ec*ec*d;
    if (band>SIZE_MAX/(2*sizeof(*scores)) || (n&&band>SIZE_MAX/n))
        return -1;
    parent_states=n*band; memset(results,0,ec*sizeof(*results));
    scores=malloc(2*band*sizeof(*scores));
    parent=calloc(parent_states?parent_states:1,1);
    if (!scores || !parent) goto done;
    old=scores; new_scores=old+band;
    for (q=0;q<2*band;++q) scores[q]=INT16_MIN;
#define MAX_STATE(H,E,R) ((((H)*ec+(E))*d)+(R))
#define MAX_PARENT(L,H,E,R) ((((L)-1)*band)+MAX_STATE(H,E,R))
    old[0]=0;
    for (l=1;l<=n;++l) {
        const size_t original=s[l-1];
        const size_t max_h=(l-1)<total_edits?(l-1):total_edits;
        if (original>3) goto done;
        for (q=0;q<band;++q) new_scores[q]=INT16_MIN;
        for (h=0;h<=max_h && h<total_edits;++h)
            for (e=0;e+h+1<=total_edits;++e)
                memcpy(new_scores+MAX_STATE(h+1,e,0),
                       old+MAX_STATE(h,e,0),d*sizeof(*old));
        for (h=0;h<=max_h;++h) {
            const size_t row=(l-1)-h;
            for (e=0;e+h<=total_edits;++e) for (a=0;a<4;++a) {
                const size_t changed=a!=original;
                if (e+h+changed>total_edits) continue;
                update_shifted(old+MAX_STATE(h,e,0),
                    new_scores+MAX_STATE(h,e+changed,0),
                    t->A+(row*4+a)*d,
                    parent+MAX_PARENT(l,h,e+changed,0),d,
                    t->C[4*row+a],(unsigned)a);
            }
        }
        { int16_t *tmp=old; old=new_scores; new_scores=tmp; }
    }
    for (result_index=0;result_index<ec;++result_index) {
        SubseqHashEditResult *out=results+result_index;
        size_t substitutions=result_index;
        size_t deletions=total_edits-substitutions,k=n-deletions;
        out->deletions=deletions; out->substitutions=substitutions;
        out->length=k; out->seed=malloc(k?k:1);
        out->positions=malloc((k?k:1)*sizeof(*out->positions));
        out->substituted=malloc(k?k:1);
        if (!out->seed || !out->positions || !out->substituted) goto done;
        for (r=0;r<d && old[MAX_STATE(deletions,substitutions,r)]==
                         INT16_MIN;++r) {}
        if (r==d) goto done;
        out->psi=r; out->omega=(double)old[MAX_STATE(
            deletions,substitutions,r)];
        l=n; h=deletions; e=substitutions; q=k;
        while (q) {
            unsigned code=parent[MAX_PARENT(l,h,e,r)];
            if (!code) { --h; --l; continue; }
            a=code-1u; out->seed[q-1]=(uint8_t)a;
            out->positions[q-1]=l-1;
            out->substituted[q-1]=(uint8_t)(a!=s[l-1]);
            if (a!=s[l-1]) --e;
            r=subtract_mod(r,t->C[4*(q-1)+a],d);
            --q; --l;
        }
    }
    status=0;
done:
    if (status) subseqhash_edit_results_destroy(results,ec);
    free(parent); free(scores); return status;
#undef MAX_PARENT
#undef MAX_STATE
}

const char *subseqhash_edits_max_i8_backend(void) {
#if MAX_I8_AVX512 && defined(__AVX512VL__)
    return "AVX-512BW (32 x int16 max-only states)";
#elif MAX_I8_AVX512
    return "AVX-512BW (32 x int16 max-only states; scalar tails)";
#else
    return "scalar int16 max-only";
#endif
}
