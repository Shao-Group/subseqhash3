#define _POSIX_C_SOURCE 200112L
#include "ssh3_pivots.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#if defined(__AVX512F__) && defined(__AVX512BW__) && defined(__AVX512VL__)
#include <immintrin.h>
#define PIVOT_SIMD 1
#else
#define PIVOT_SIMD 0
#endif

#define UNREACHABLE UINT16_MAX

typedef struct {
    uint16_t *minimum,*maximum;
    uint8_t *parent;
    size_t stride,states,layer,d,edits;
} PivotDP;

typedef struct {
    size_t psi,input_pivot,left_substitutions,left_residue,right_residue;
    unsigned symbol,cost;
} PivotChoice;

struct SubseqHash3PivotWorkspace {
    PivotDP left,right;
    size_t left_layers,right_layers,scratch_capacity;
    size_t *scratch;
};

SubseqHash3PivotWorkspace *subseqhash3_pivot_workspace_create(void) {
    return calloc(1,sizeof(SubseqHash3PivotWorkspace));
}

static void workspace_clear(SubseqHash3PivotWorkspace *workspace) {
    free(workspace->left.minimum); free(workspace->right.minimum);
    free(workspace->scratch);
    memset(workspace,0,sizeof(*workspace));
}

void subseqhash3_pivot_workspace_destroy(SubseqHash3PivotWorkspace *workspace) {
    if (!workspace) return;
    workspace_clear(workspace);
    free(workspace);
}

static uint64_t pivot_random(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}

static size_t subtract_residue(size_t a,size_t b,size_t d) {
    return a>=b?a-b:a+(d-b);
}

static void *pivot_allocate(size_t bytes) {
    void *p=NULL;
    if (posix_memalign(&p,64,bytes?bytes:64)) return NULL;
    return p;
}

int subseqhash3_pivot_tables_init(SubseqHash3PivotTables *t,
    size_t n,size_t d,uint64_t seed) {
    size_t k,cells,i,a,r;
    if (!t) return -1;
    memset(t,0,sizeof(*t));
    if (n<2 || !d || n-1>256 || d>SIZE_MAX/(4*(n-1))) return -1;
    k=n-1; cells=k*4*d; t->k=k; t->d=d;
    if (subseqhash3_tables_init(&t->left,k,d,pivot_random(&seed)) ||
        subseqhash3_tables_init(&t->right,k,d,pivot_random(&seed)))
        goto failed;
    t->pivot_A=malloc(k*4);
    t->pivot_left_sign=malloc(k*4);
    t->pivot_right_sign=malloc(k*4);
    t->pivot_C=malloc(k*4*sizeof(*t->pivot_C));
    t->right_reversed_A=pivot_allocate(cells);
    t->right_reversed_B1_masks=calloc(k*4*t->right.blocks,
        sizeof(*t->right_reversed_B1_masks));
    if (!t->pivot_A || !t->pivot_left_sign || !t->pivot_right_sign ||
        !t->pivot_C || !t->right_reversed_A || !t->right_reversed_B1_masks)
        goto failed;
    for (i=0;i<k;++i) {
        unsigned permutation[4]={0,1,2,3};
        for (a=4;a>1;--a) {
            size_t q=(size_t)(pivot_random(&seed)%a);
            unsigned tmp=permutation[a-1];
            permutation[a-1]=permutation[q]; permutation[q]=tmp;
        }
        for (a=0;a<4;++a) {
            size_t residue;
            int value=(int)(1+pivot_random(&seed)%127);
            t->pivot_A[4*i+a]=(int8_t)((pivot_random(&seed)&1)?value:-value);
            t->pivot_left_sign[4*i+a]=(permutation[a]&1)?1:-1;
            t->pivot_right_sign[4*i+a]=(permutation[a]&2)?1:-1;
            do { residue=(size_t)(pivot_random(&seed)%d); }
            while (d>=4 && a && (residue==t->pivot_C[4*i] ||
                (a>1 && residue==t->pivot_C[4*i+1]) ||
                (a>2 && residue==t->pivot_C[4*i+2])));
            t->pivot_C[4*i+a]=residue;
        }
    }
    /* Right-side scores are stored at -psi mod D, turning the pivot
       residue join into two contiguous, circular SIMD loads. */
    for (i=0;i<k;++i) for (a=0;a<4;++a) for (r=0;r<d;++r) {
        size_t original=r?d-r:0;
        size_t dst=(4*i+a)*d+r,src=(4*i+a)*d+original;
        t->right_reversed_A[dst]=t->right.signed_A[src];
        if (t->right.B1[src]>0)
            t->right_reversed_B1_masks[(4*i+a)*t->right.blocks+r/32]|=
                UINT32_C(1)<<(r%32);
    }
    return 0;
failed:
    subseqhash3_pivot_tables_destroy(t);
    return -1;
}

void subseqhash3_pivot_tables_destroy(SubseqHash3PivotTables *t) {
    if (!t) return;
    subseqhash3_tables_destroy(&t->left);
    subseqhash3_tables_destroy(&t->right);
    free(t->pivot_A); free(t->pivot_left_sign); free(t->pivot_right_sign);
    free(t->pivot_C); free(t->right_reversed_A);
    free(t->right_reversed_B1_masks);
    memset(t,0,sizeof(*t));
}

static size_t dp_state(const PivotDP *dp,size_t l,size_t h,size_t e) {
    return l*dp->layer+(h*(dp->edits+1)+e)*dp->stride;
}

static size_t dp_parent(const PivotDP *dp,size_t l,size_t h,size_t e) {
    return (l*dp->states+h*(dp->edits+1)+e)*dp->d;
}

static void circular_padding(uint16_t *p,size_t d) {
    size_t done=0;
#if PIVOT_SIMD
    /* A 32-lane load starting at D-1 needs only 31 wrapped elements.
       At D=31 the final stored lane is unused, so one vector copy also
       replaces the old two memcpy calls in that important small-D case. */
    if (d>=31) {
        _mm512_storeu_si512((void *)(p+d),
            _mm512_loadu_si512((const void *)p));
        return;
    }
#endif
    while (done<31) {
        size_t count=d<31-done?d:31-done;
        memcpy(p+d+done,p,count*sizeof(*p));
        done+=count;
    }
}

static int dp_allocate(PivotDP *dp,size_t n,size_t edits,size_t d) {
    size_t cells,bytes;
    memset(dp,0,sizeof(*dp));
    if (d>SIZE_MAX-63 || edits==SIZE_MAX || edits+1>SIZE_MAX/2) return -1;
    dp->d=d; dp->edits=edits;
    dp->stride=((d+31)&~(size_t)31)+32;
    dp->states=2*(edits+1);
    if (dp->states>SIZE_MAX/d || dp->states>SIZE_MAX/dp->stride) return -1;
    dp->layer=dp->states*dp->stride;
    if (n>SIZE_MAX/dp->layer) return -1;
    cells=n*dp->layer;
    if (cells>SIZE_MAX/(2*sizeof(uint16_t))) return -1;
    bytes=2*cells*sizeof(uint16_t);
    /* Small residue spaces favor direct traceback; for large spaces avoid
       recording millions of parents to recover just a few hundred bases. */
    if (d<=64 && cells>SIZE_MAX-bytes) return -1;
    dp->minimum=pivot_allocate(bytes+(d<=64?cells:0));
    if (!dp->minimum) return -1;
    dp->maximum=dp->minimum+cells;
    dp->parent=d<=64?(uint8_t *)(dp->maximum+cells):NULL;
    return 0;
}

static int dp_reserve(PivotDP *dp,size_t *capacity,
    size_t layers,size_t edits,size_t d) {
    if (dp->minimum && dp->d==d && dp->edits==edits && *capacity>=layers)
        return 0;
    free(dp->minimum); dp->minimum=NULL; *capacity=0;
    if (dp_allocate(dp,layers,edits,d)) return -1;
    *capacity=layers;
    return 0;
}

static void empty_state(PivotDP *dp,size_t l,size_t h) {
    const size_t state=dp_state(dp,l,h,0);
    memset(dp->minimum+state,255,dp->stride*sizeof(uint16_t));
    memset(dp->maximum+state,255,dp->stride*sizeof(uint16_t));
    dp->minimum[state]=dp->maximum[state]=0;
    circular_padding(dp->minimum+state,dp->d);
    circular_padding(dp->maximum+state,dp->d);
    if (dp->parent) memset(dp->parent+dp_parent(dp,l,h,0),0,dp->d);
}

static void dp_build(PivotDP *workspace,const uint8_t *s,size_t n,size_t maximum,
    const SubseqHash3PivotTables *tables,int reverse,int record) {
    /* Keep immutable layout metadata out of the score arrays' alias set.
       Explicit alignment also avoids GCC's aligned struct-copy store into
       an insufficiently aligned ASan fake-stack slot. */
    PivotDP layout __attribute__((aligned(64)))=*workspace;
    PivotDP *dp=&layout;
    const SubseqHash3Tables *t=reverse?&tables->right:&tables->left;
    const int8_t *signed_A=reverse?tables->right_reversed_A:t->signed_A;
    const uint32_t *masks=reverse?tables->right_reversed_B1_masks:t->B1_masks;
    const size_t d=dp->d;
    size_t l,h,e,r,a,j;
    empty_state(dp,0,0);
    empty_state(dp,1,1);
    /* Visit both deletion planes at the same retained length. They use the
       same table row; h=1's deletion predecessor is the just-built h=0
       state. The other predecessors have retained length j-1. */
    for (j=1;j<maximum;++j) for (h=0;h<=1;++h) {
        const size_t retained=j,max_e=retained<dp->edits?retained:dp->edits;
        const size_t min_e=dp->edits>n-1-retained?
            dp->edits-(n-1-retained):0;
        l=j+h;
        const unsigned original=s[reverse?n-l:l-1];
        size_t shifts[4];
        for (a=0;a<4;++a) {
            shifts[a]=t->C[4*(retained-1)+a];
            if (reverse && shifts[a]) shifts[a]=d-shifts[a];
        }
#if PIVOT_SIMD
        for (r=0;r<d;r+=32) {
            const size_t lanes=d-r<32?d-r:32;
            const __mmask32 active=lanes==32?(__mmask32)UINT32_MAX:
                (__mmask32)((UINT32_C(1)<<lanes)-1);
            const __m512i infinity=_mm512_set1_epi16(-1);
            const __m512i bias=_mm512_set1_epi16(127);
            const unsigned sym0=original==0?1:0;
            const unsigned sym1=original<=1?2:1;
            const unsigned sym2=original<=2?3:2;
#define PREPARE(ID,SYM) \
            const size_t row##ID=4*(retained-1)+(SYM); \
            const __m512i add##ID=_mm512_cvtepi8_epi16( \
                _mm256_maskz_loadu_epi8(active,signed_A+row##ID*d+r)); \
            const __m512i plus##ID=_mm512_add_epi16(bias,add##ID); \
            const __m512i minus##ID=_mm512_sub_epi16(bias,add##ID); \
            const __mmask32 positive##ID=(__mmask32)masks[row##ID*t->blocks+r/32]; \
            const size_t q##ID=subtract_residue(r,shifts[SYM],d)
            PREPARE(0,sym0); PREPARE(1,sym1); PREPARE(2,sym2);
            PREPARE(3,original);
#undef PREPARE
            for (e=min_e;e<=max_e;++e) {
                const size_t dst=dp_state(dp,l,h,e);
                const size_t deletion=h?dp_state(dp,l-1,h-1,e):0;
                __m512i mn=h?_mm512_load_si512(
                    (const void *)(dp->minimum+deletion+r)):infinity;
                __m512i mx=h?_mm512_load_si512(
                    (const void *)(dp->maximum+deletion+r)):infinity;
                __m256i low=_mm256_setzero_si256(),high=low;
#define PIVOT_UPDATE(ID,SOURCE) do { \
    const size_t src_=(SOURCE)+q##ID; \
    const __m512i smn_=_mm512_loadu_si512((const void *)(dp->minimum+src_)); \
    const __m512i smx_=_mm512_loadu_si512((const void *)(dp->maximum+src_)); \
    const __m512i lo_=_mm512_adds_epu16( \
        _mm512_mask_blend_epi16(positive##ID,smx_,smn_),plus##ID); \
    const __m512i hi_=_mm512_adds_epu16( \
        _mm512_mask_blend_epi16(positive##ID,smn_,smx_),minus##ID); \
    if (record) { \
        const unsigned sym_=(ID)==0?sym0:(ID)==1?sym1:(ID)==2?sym2:original; \
        const unsigned base_=2*sym_; \
        const __m256i code_=_mm256_mask_blend_epi8(positive##ID, \
            _mm256_set1_epi8((char)((base_+2)|((base_+1)<<4))), \
            _mm256_set1_epi8((char)((base_+1)|((base_+2)<<4)))); \
        low=_mm256_mask_mov_epi8(low,_mm512_cmplt_epu16_mask(lo_,mn),code_); \
        high=_mm256_mask_mov_epi8(high,_mm512_cmplt_epu16_mask(hi_,mx),code_); \
    } \
    mn=_mm512_min_epu16(mn,lo_); mx=_mm512_min_epu16(mx,hi_); \
} while (0)
                if (e) {
                    const size_t changed=dp_state(dp,l-1,h,e-1);
                    PIVOT_UPDATE(0,changed); PIVOT_UPDATE(1,changed);
                    PIVOT_UPDATE(2,changed);
                }
                if (e<retained) {
                    const size_t kept=dp_state(dp,l-1,h,e);
                    PIVOT_UPDATE(3,kept);
                }
#undef PIVOT_UPDATE
                _mm512_store_si512((void *)(dp->minimum+dst+r),mn);
                _mm512_store_si512((void *)(dp->maximum+dst+r),mx);
                if (record) _mm256_mask_storeu_epi8(
                    dp->parent+dp_parent(dp,l,h,e)+r,active,
                    _mm256_ternarylogic_epi32(low,high,
                        _mm256_set1_epi8((char)0xf0),0xd8));
            }
        }
#else
        for (e=min_e;e<=max_e;++e) {
            const size_t dst=dp_state(dp,l,h,e);
            for (r=0;r<d;++r) {
                unsigned mn=UNREACHABLE,mx=UNREACHABLE,low=0,high=0;
                unsigned candidate;
                if (h) {
                    const size_t previous=dp_state(dp,l-1,0,e)+r;
                    mn=dp->minimum[previous]; mx=dp->maximum[previous];
                }
                for (candidate=0;candidate<4;++candidate) {
                    const unsigned symbol=candidate<3?
                        candidate+(candidate>=original):original;
                    const size_t changed=symbol!=original;
                    size_t src,table;
                    unsigned lo,hi;
                    int add,positive;
                    if (changed?(!e):(e>=retained)) continue;
                    src=dp_state(dp,l-1,h,e-changed)+
                        subtract_residue(r,shifts[symbol],d);
                    table=(4*(retained-1)+symbol)*d+r;
                    add=signed_A[table];
                    positive=(masks[(4*(retained-1)+symbol)*t->blocks+r/32]>>
                        (r%32))&1u;
                    lo=positive?dp->minimum[src]:dp->maximum[src];
                    hi=positive?dp->maximum[src]:dp->minimum[src];
                    if (lo!=UNREACHABLE && lo+(unsigned)(127+add)<mn) {
                        mn=lo+(unsigned)(127+add);
                        low=2*symbol+(positive?1u:2u);
                    }
                    if (hi!=UNREACHABLE && hi+(unsigned)(127-add)<mx) {
                        mx=hi+(unsigned)(127-add);
                        high=2*symbol+(positive?2u:1u);
                    }
                }
                dp->minimum[dst+r]=(uint16_t)mn;
                dp->maximum[dst+r]=(uint16_t)mx;
                if (record) dp->parent[dp_parent(dp,l,h,e)+r]=
                    (uint8_t)(low|(high<<4));
            }
        }
#endif
        for (e=min_e;e<=max_e;++e) {
            const size_t dst=dp_state(dp,l,h,e);
            circular_padding(dp->minimum+dst,d);
            circular_padding(dp->maximum+dst,d);
        }
    }
}

/* For a fixed output residue, the reflected right-side layout turns the
   residue constraint into a circular shift and pairwise SIMD addition. */
static unsigned join_cost(const uint16_t *left,const uint16_t *right,
    size_t d,size_t shift) {
    size_t r;
#if PIVOT_SIMD
    __m512i minimum=_mm512_set1_epi16(-1);
    for (r=0;r<d;r+=32) {
        const size_t count=d-r<32?d-r:32;
        const __mmask32 active=count==32?(__mmask32)UINT32_MAX:
            (__mmask32)((UINT32_C(1)<<count)-1);
        const size_t q=r>=d-shift?r-(d-shift):r+shift;
        const __m512i sum=_mm512_adds_epu16(
            _mm512_load_si512((const void *)(left+r)),
            _mm512_loadu_si512((const void *)(right+q)));
        minimum=_mm512_mask_min_epu16(minimum,active,minimum,sum);
    }
    {
        const __m512i pairs=_mm512_min_epu16(minimum,
            _mm512_srli_epi32(minimum,16));
        return _mm512_reduce_min_epu32(_mm512_and_si512(
            pairs,_mm512_set1_epi32(65535)));
    }
#else
    unsigned best=UNREACHABLE;
    for (r=0;r<d;++r) {
        const size_t q=r>=d-shift?r-(d-shift):r+shift;
        const unsigned value=(unsigned)left[r]+right[q];
        if (value<best) best=value;
    }
    return best;
#endif
}

/* Resolve the first residue only after the winning score is known. */
static size_t join_witness(const uint16_t *left,const uint16_t *right,
    size_t d,size_t shift,unsigned best) {
    size_t r;
#if PIVOT_SIMD
    const __m512i target=_mm512_set1_epi16((short)best);
    for (r=0;r<d;r+=32) {
        const size_t count=d-r<32?d-r:32;
        const __mmask32 active=count==32?(__mmask32)UINT32_MAX:
            (__mmask32)((UINT32_C(1)<<count)-1);
        const size_t q=r>=d-shift?r-(d-shift):r+shift;
        const __m512i sum=_mm512_adds_epu16(
            _mm512_load_si512((const void *)(left+r)),
            _mm512_loadu_si512((const void *)(right+q)));
        const uint32_t hits=(uint32_t)(_mm512_cmpeq_epi16_mask(sum,target)&active);
        if (hits) return r+(size_t)__builtin_ctz(hits);
    }
#else
    for (r=0;r<d;++r) {
        const size_t q=r>=d-shift?r-(d-shift):r+shift;
        if ((unsigned)left[r]+right[q]==best) return r;
    }
#endif
    return 0; /* The caller has already established a reachable optimum. */
}

/* Find min (left_psi + pivot_psi - reflected_right_psi) mod D in O(D),
   using the previous reachable reflected right residue at every threshold.
   This handles sparse states and does not assume residue zero is reachable. */
static size_t first_residue(const uint16_t *left,const uint16_t *right,
    size_t d,size_t pivot,size_t *previous) {
    size_t r,last=SIZE_MAX,best=d;
    for (r=0;r<d;++r) {
        if (right[r]!=UNREACHABLE) last=r;
        previous[r]=last;
    }
    if (last==SIZE_MAX) return d;
    for (r=0;r<d;++r) if (left[r]!=UNREACHABLE) {
        const size_t sum=r>=d-pivot?r-(d-pivot):r+pivot;
        size_t q=previous[sum];
        size_t value;
        if (q==SIZE_MAX) q=last;
        value=subtract_residue(sum,q,d);
        if (value<best) best=value;
        if (!best) break;
    }
    return best;
}

#if PIVOT_SIMD
static unsigned minimum_word(__m512i values) {
    const __m512i pairs=_mm512_min_epu16(values,_mm512_srli_epi32(values,16));
    return _mm512_reduce_min_epu32(_mm512_and_si512(
        pairs,_mm512_set1_epi32(65535)));
}

/* When the whole residue space fits in one register, keep the winning
   candidate in lanes across ALL edit splits, and reduce just once. */
static int choose_pivot_small(const uint8_t *s,size_t n,size_t edits,size_t pivot,
    const SubseqHash3PivotTables *t,const PivotDP *left,const PivotDP *right,
    PivotChoice *choice) {
    const size_t k=n-1,d=t->d,width=edits+1;
    const __mmask32 active=d==32?(__mmask32)UINT32_MAX:
        (__mmask32)((UINT32_C(1)<<d)-1);
    __m512i minimum=_mm512_set1_epi16(-1),indices=minimum;
    size_t input,e,table,ls,rs,shift,arg,id;
    unsigned symbol,best;
    const uint16_t *lc,*rc;
    for (input=pivot;input<=pivot+1;++input) {
        const size_t left_del=input-pivot,right_del=1-left_del;
        const size_t right_length=k-pivot-1;
        for (symbol=0;symbol<4;++symbol) {
            const size_t changed=symbol!=s[input];
            size_t min_e,max_e;
            __m512i bias;
            if (changed>edits) continue;
            table=4*pivot+symbol;
            bias=_mm512_set1_epi16((short)(127-t->pivot_A[table]));
            min_e=edits-changed>right_length?edits-changed-right_length:0;
            max_e=edits-changed<pivot?edits-changed:pivot;
            for (e=min_e;e<=max_e;++e) {
                const size_t right_e=edits-changed-e;
                __m512i cost;
                __mmask32 better;
                ls=dp_state(left,input,left_del,e);
                rs=dp_state(right,n-input-1,right_del,right_e);
                lc=(t->pivot_left_sign[table]>0?left->maximum:left->minimum)+ls;
                rc=(t->pivot_right_sign[table]>0?right->maximum:right->minimum)+rs;
                cost=_mm512_adds_epu16(_mm512_load_si512((const void *)lc),
                    _mm512_loadu_si512((const void *)(rc+t->pivot_C[table])));
                cost=_mm512_adds_epu16(cost,bias);
                better=active&_mm512_cmplt_epu16_mask(cost,minimum);
                minimum=_mm512_mask_mov_epi16(minimum,better,cost);
                id=((input-pivot)*4+symbol)*width+e;
                indices=_mm512_mask_mov_epi16(indices,better,
                    _mm512_set1_epi16((short)id));
            }
        }
    }
    best=minimum_word(minimum);
    if (best==UNREACHABLE) return -1; /* Use the general sparse-psi search. */
    /* Nonmatching lanes must sort after every finite candidate index. */
    indices=_mm512_mask_mov_epi16(_mm512_set1_epi16(-1),
        _mm512_cmpeq_epi16_mask(minimum,_mm512_set1_epi16((short)best)),indices);
    id=minimum_word(indices);
    choice->input_pivot=pivot+id/(4*width);
    choice->symbol=(unsigned)((id/width)%4);
    choice->left_substitutions=id%width;
    choice->psi=0; choice->cost=best;
    table=4*pivot+choice->symbol;
    ls=dp_state(left,choice->input_pivot,choice->input_pivot-pivot,
        choice->left_substitutions);
    rs=dp_state(right,n-choice->input_pivot-1,
        1-(choice->input_pivot-pivot),
        edits-(choice->symbol!=s[choice->input_pivot])-choice->left_substitutions);
    lc=(t->pivot_left_sign[table]>0?left->maximum:left->minimum)+ls;
    rc=(t->pivot_right_sign[table]>0?right->maximum:right->minimum)+rs;
    shift=t->pivot_C[table];
    arg=join_witness(lc,rc,d,shift,best-(unsigned)(127-t->pivot_A[table]));
    choice->left_residue=arg;
    choice->right_residue=arg>=d-shift?arg-(d-shift):arg+shift;
    return 0;
}
#endif

static int choose_pivot(const uint8_t *s,size_t n,size_t edits,size_t pivot,
    const SubseqHash3PivotTables *t,const PivotDP *left,const PivotDP *right,
    size_t *scratch,PivotChoice *choice) {
    const size_t k=n-1,d=t->d;
    size_t input,e;
    unsigned symbol;
#if PIVOT_SIMD
    if (d<=32 && !choose_pivot_small(s,n,edits,pivot,t,left,right,choice))
        return 0;
#endif
    memset(choice,0,sizeof(*choice)); choice->psi=d; choice->cost=UNREACHABLE;
    for (input=pivot;input<=pivot+1;++input) {
        const size_t left_del=input-pivot,right_del=1-left_del;
        const size_t right_length=k-pivot-1;
        for (symbol=0;symbol<4;++symbol) {
            const size_t changed=symbol!=s[input];
            const size_t table=4*pivot+symbol;
            size_t min_e,max_e;
            if (changed>edits) continue;
            min_e=edits-changed>right_length?edits-changed-right_length:0;
            max_e=edits-changed<pivot?edits-changed:pivot;
            for (e=min_e;e<=max_e;++e) {
                const size_t right_e=edits-changed-e;
                const size_t ls=dp_state(left,input,left_del,e);
                const size_t rs=dp_state(right,n-input-1,right_del,right_e);
                const uint16_t *lc=(t->pivot_left_sign[table]>0?
                    left->maximum:left->minimum)+ls;
                const uint16_t *rc=(t->pivot_right_sign[table]>0?
                    right->maximum:right->minimum)+rs;
                size_t psi=0,shift=t->pivot_C[table];
                unsigned cost=join_cost(lc,rc,d,shift);
                if (cost==UNREACHABLE) {
                    if (!choice->psi) continue;
                    psi=first_residue(lc,rc,d,t->pivot_C[table],scratch);
                    if (psi==d || psi>choice->psi) continue;
                    shift=subtract_residue(t->pivot_C[table],psi,d);
                    cost=join_cost(lc,rc,d,shift);
                    if (cost==UNREACHABLE) continue;
                }
                cost+=(unsigned)(127-t->pivot_A[table]);
                if (psi<choice->psi || (psi==choice->psi && cost<choice->cost)) {
                    choice->psi=psi; choice->cost=cost;
                    choice->input_pivot=input; choice->symbol=symbol;
                    choice->left_substitutions=e;
                }
            }
        }
    }
    if (choice->psi==d) return -1;
    /* Argmin bookkeeping is needed once for the winning join, not for
       every pivot/base/substitution split considered above. */
    {
        const size_t table=4*pivot+choice->symbol;
        const size_t left_del=choice->input_pivot-pivot;
        const size_t changed=choice->symbol!=s[choice->input_pivot];
        const size_t right_e=edits-changed-choice->left_substitutions;
        const size_t ls=dp_state(left,choice->input_pivot,left_del,
            choice->left_substitutions);
        const size_t rs=dp_state(right,n-choice->input_pivot-1,1-left_del,right_e);
        const uint16_t *lc=(t->pivot_left_sign[table]>0?
            left->maximum:left->minimum)+ls;
        const uint16_t *rc=(t->pivot_right_sign[table]>0?
            right->maximum:right->minimum)+rs;
        const size_t shift=subtract_residue(t->pivot_C[table],choice->psi,d);
        const size_t arg=join_witness(lc,rc,d,shift,
            choice->cost-(unsigned)(127-t->pivot_A[table]));
        choice->left_residue=arg;
        choice->right_residue=arg>=d-shift?arg-(d-shift):arg+shift;
    }
    return 0;
}

static void trace_side(const uint8_t *s,size_t n,const PivotDP *dp,
    const SubseqHash3Tables *t,int reverse,size_t l,size_t h,size_t e,
    size_t residue,int type,SubseqHash3EditResult *out) {
    size_t retained=l-h;
    while (retained) {
        const size_t current=dp_state(dp,l,h,e)+residue;
        const unsigned best=type==1?dp->minimum[current]:dp->maximum[current];
        const unsigned original=s[reverse?n-l:l-1];
        unsigned candidate;
        size_t position,output,shift;
        unsigned symbol=original;
        /* Recover only the chosen witnesses from already stored scores.
           Deletion, ascending changed bases, then unchanged reproduce the
           reference's strict-update tie order without storing any parents. */
        if (dp->parent) {
            const unsigned packed=dp->parent[dp_parent(dp,l,h,e)+residue];
            unsigned code=type==1?(packed&15u):(packed>>4);
            if (!code) { --l; --h; continue; }
            --code; symbol=code/2; type=(int)(code%2)+1;
        } else {
            if (h) {
                const size_t previous=dp_state(dp,l-1,h-1,e)+residue;
                const unsigned value=type==1?dp->minimum[previous]:dp->maximum[previous];
                if (value==best) { --l; --h; continue; }
            }
            for (candidate=0;candidate<4;++candidate) {
                const unsigned a=candidate<3?
                    candidate+(candidate>=original):original;
                const size_t changed=a!=original;
                size_t q,src,table;
                int source_type,add;
                unsigned value;
                if (changed?(!e):(e>=retained)) continue;
                shift=t->C[4*(retained-1)+a];
                if (reverse && shift) shift=dp->d-shift;
                q=subtract_residue(residue,shift,dp->d);
                table=(4*(retained-1)+a)*dp->d+
                    (reverse && residue?dp->d-residue:residue);
                source_type=t->B1[table]>0?type:3-type;
                add=127+(type==1?t->signed_A[table]:-t->signed_A[table]);
                src=dp_state(dp,l-1,h,e-changed)+q;
                value=source_type==1?dp->minimum[src]:dp->maximum[src];
                if (value!=UNREACHABLE && value+(unsigned)add==best) {
                    symbol=a; type=source_type;
                    break;
                }
            }
        }
        position=reverse?n-l:l-1;
        output=reverse?out->length-retained:retained-1;
        out->seed[output]=(uint8_t)symbol;
        out->positions[output]=position;
        out->substituted[output]=(uint8_t)(symbol!=s[position]);
        if (symbol!=s[position]) --e;
        shift=t->C[4*(retained-1)+symbol];
        if (reverse && shift) shift=dp->d-shift;
        residue=subtract_residue(residue,shift,dp->d);
        --retained; --l;
    }
}

static int compute_pivot_seeds_impl(const uint8_t *s,size_t n,size_t edits,
    const SubseqHash3PivotTables *t,const size_t *pivots,size_t count,
    SubseqHash3EditResult *results,SubseqHash3PivotWorkspace *workspace) {
    PivotDP *left=&workspace->left,*right=&workspace->right;
    size_t *scratch=NULL,index,k,min_pivot,max_pivot;
    int status=-1;
    if (!results || !count || count>SIZE_MAX/sizeof(*results)) return -1;
    memset(results,0,count*sizeof(*results));
    if (!s || !t || n<2 || n-1>t->k || !t->d || edits>n-1 ||
        (!pivots && count>n-1)) return -1;
    k=n-1;
    for (index=0;index<n;++index) if (s[index]>3) return -1;
    min_pivot=k; max_pivot=0;
    for (index=0;index<count;++index) {
        const size_t pivot=pivots?pivots[index]:index;
        if (pivot>=k) return -1;
        if (pivot<min_pivot) min_pivot=pivot;
        if (pivot>max_pivot) max_pivot=pivot;
    }
    if (t->d>SIZE_MAX/sizeof(*scratch)) return -1;
    if (workspace->scratch_capacity<t->d) {
        scratch=realloc(workspace->scratch,t->d*sizeof(*scratch));
        if (!scratch) goto done;
        workspace->scratch=scratch; workspace->scratch_capacity=t->d;
    }
    scratch=workspace->scratch;
    if (dp_reserve(left,&workspace->left_layers,max_pivot+2,edits,t->d) ||
        dp_reserve(right,&workspace->right_layers,n-min_pivot,edits,t->d))
        goto done;
    if (left->parent) dp_build(left,s,n,max_pivot+1,t,0,1);
    else dp_build(left,s,n,max_pivot+1,t,0,0);
    if (right->parent) dp_build(right,s,n,n-min_pivot-1,t,1,1);
    else dp_build(right,s,n,n-min_pivot-1,t,1,0);
    for (index=0;index<count;++index) {
        const size_t pivot=pivots?pivots[index]:index;
        PivotChoice choice;
        SubseqHash3EditResult *out=results+index;
        size_t table,left_del,changed;
        if (choose_pivot(s,n,edits,pivot,t,left,right,scratch,&choice)) goto done;
        out->seed=malloc(k); out->positions=malloc(k*sizeof(*out->positions));
        out->substituted=malloc(k);
        if (!out->seed || !out->positions || !out->substituted) goto done;
        out->length=k; out->deletions=1; out->substitutions=edits;
        out->psi=choice.psi; out->omega=127.0*k-choice.cost;
        out->seed[pivot]=(uint8_t)choice.symbol;
        out->positions[pivot]=choice.input_pivot;
        out->substituted[pivot]=(uint8_t)(choice.symbol!=s[choice.input_pivot]);
        table=4*pivot+choice.symbol;
        left_del=choice.input_pivot-pivot;
        changed=out->substituted[pivot];
        trace_side(s,n,left,&t->left,0,choice.input_pivot,left_del,
            choice.left_substitutions,choice.left_residue,
            t->pivot_left_sign[table]>0?2:1,out);
        trace_side(s,n,right,&t->right,1,n-choice.input_pivot-1,1-left_del,
            edits-changed-choice.left_substitutions,choice.right_residue,
            t->pivot_right_sign[table]>0?2:1,out);
    }
    status=0;
done:
    if (status) subseqhash3_edit_results_destroy(results,count);
    return status;
}

int subseqhash3_compute_pivot_seeds(const uint8_t *s,size_t n,size_t edits,
    const SubseqHash3PivotTables *t,const size_t *pivots,size_t count,
    SubseqHash3EditResult *results) {
    SubseqHash3PivotWorkspace workspace={0};
    int status=compute_pivot_seeds_impl(s,n,edits,t,pivots,count,results,&workspace);
    workspace_clear(&workspace);
    return status;
}

int subseqhash3_compute_pivot_repeats_with_workspace(const uint8_t *s,size_t n,size_t edits,
    const SubseqHash3PivotTables *tables,size_t table_count,size_t repeats,
    SubseqHash3EditResult *results,SubseqHash3PivotWorkspace *workspace) {
    size_t needed,family,done=0,k;
    if (!workspace || !results || !tables || n<2 || !repeats ||
        repeats>SIZE_MAX/sizeof(*results)) return -1;
    memset(results,0,repeats*sizeof(*results));
    k=n-1; needed=1+(repeats-1)/k;
    if (table_count<needed) return -1;
    for (family=0;family<needed;++family) {
        size_t count=repeats-done<k?repeats-done:k;
        if (compute_pivot_seeds_impl(s,n,edits,tables+family,NULL,
                count,results+done,workspace)) {
            subseqhash3_edit_results_destroy(results,repeats);
            return -1;
        }
        done+=count;
    }
    return 0;
}

int subseqhash3_compute_pivot_repeats(const uint8_t *s,size_t n,size_t edits,
    const SubseqHash3PivotTables *tables,size_t table_count,size_t repeats,
    SubseqHash3EditResult *results) {
    SubseqHash3PivotWorkspace workspace={0};
    int status=subseqhash3_compute_pivot_repeats_with_workspace(s,n,edits,
        tables,table_count,repeats,results,&workspace);
    workspace_clear(&workspace);
    return status;
}
