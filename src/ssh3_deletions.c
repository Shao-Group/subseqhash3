#define _POSIX_C_SOURCE 200112L
#include "ssh3_deletions.h"
#include "ssh3_seeding.h"
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__AVX512F__) && defined(__AVX512BW__) && defined(__AVX512VL__)
#include <immintrin.h>
#define DEL_SIMD 1
#else
#define DEL_SIMD 0
#endif
#include "single_deletion.h"
#define INF (INT64_MAX/2)

static size_t submod(size_t a,size_t b,size_t d) { return a>=b?a-b:d-(b-a); }
static size_t addmod(size_t a,size_t b,size_t d) { return a>=d-b?a-(d-b):a+b; }
static int tables_allocate(SubseqHash3DeletionTables *t,size_t k,size_t d) {
    size_t rows,cells;
    if (!t || !k || k>256 || !d || d>SIZE_MAX/(4*k*sizeof(int64_t))) return -1;
    memset(t,0,sizeof(*t)); t->k=k; t->d=d; rows=4*k; cells=rows*d;
    t->left.A=malloc(cells*sizeof(int64_t)); t->right.A=malloc(cells*sizeof(int64_t));
    t->left.B1=malloc(cells); t->right.B1=malloc(cells);
    t->left.C=malloc(rows*sizeof(size_t)); t->right.C=malloc(rows*sizeof(size_t));
    t->pivot_A=malloc(rows*sizeof(int64_t)); t->pivot_C=malloc(rows*sizeof(size_t));
    t->left_sign=malloc(rows); t->right_sign=malloc(rows);
    if (!t->left.A || !t->right.A || !t->left.B1 || !t->right.B1 ||
        !t->left.C || !t->right.C || !t->pivot_A || !t->pivot_C ||
        !t->left_sign || !t->right_sign) {
        subseqhash3_deletion_tables_destroy(t); return -1;
    }
    return 0;
}
void subseqhash3_deletion_tables_destroy(SubseqHash3DeletionTables *t) {
    if (!t) return;
    free(t->left.A); free(t->left.B1); free(t->left.C);
    free(t->right.A); free(t->right.B1); free(t->right.C);
    free(t->pivot_A); free(t->pivot_C); free(t->left_sign); free(t->right_sign);
    memset(t,0,sizeof(*t));
}
int subseqhash3_deletion_tables_from_pivots(SubseqHash3DeletionTables *t,
    const SubseqHash3PivotTables *source) {
    size_t side,i;
    if (!t || !source || !source->k || source->k>256 || !source->d ||
        source->left.k<source->k || source->right.k<source->k ||
        source->left.d!=source->d || source->right.d!=source->d ||
        !source->left.A || !source->right.A || !source->left.B1 ||
        !source->right.B1 || !source->left.B2 || !source->right.B2 ||
        !source->left.C || !source->right.C || !source->pivot_A ||
        !source->pivot_C || !source->pivot_left_sign || !source->pivot_right_sign)
        return -1;
    if (tables_allocate(t,source->k,source->d)) return -1;
    for (side=0;side<2;++side) {
        const SubseqHash3Tables *s=side?&source->right:&source->left;
        SubseqHash3DeletionSide *dst=side?&t->right:&t->left;
        for (i=0;i<4*t->k*t->d;++i) {
            if (s->A[i]<1 || (s->B1[i]!=1 && s->B1[i]!=-1) ||
                (s->B2[i]!=1 && s->B2[i]!=-1)) goto invalid;
            dst->A[i]=(int64_t)s->A[i]*s->B2[i]; dst->B1[i]=s->B1[i];
        }
        for (i=0;i<4*t->k;++i) {
            if (s->C[i]>=t->d) goto invalid;
            dst->C[i]=s->C[i];
        }
    }
    for (i=0;i<4*t->k;++i) {
        if (!source->pivot_A[i] || source->pivot_A[i]==INT8_MIN ||
            (source->pivot_left_sign[i]!=1 && source->pivot_left_sign[i]!=-1) ||
            (source->pivot_right_sign[i]!=1 && source->pivot_right_sign[i]!=-1) ||
            source->pivot_C[i]>=t->d) goto invalid;
        t->pivot_A[i]=source->pivot_A[i]; t->pivot_C[i]=source->pivot_C[i];
        t->left_sign[i]=source->pivot_left_sign[i];
        t->right_sign[i]=source->pivot_right_sign[i];
    }
    return 0;
invalid:
    subseqhash3_deletion_tables_destroy(t); return -1;
}
int subseqhash3_deletion_tables_init(SubseqHash3DeletionTables *t,size_t k,
    size_t d,uint64_t seed) {
    SubseqHash3PivotTables shared={0};
    int status;
    if (!t || !k || k>256 || !d) return -1;
    if (subseqhash3_pivot_tables_init(&shared,k+1,d,seed)) return -1;
    status=subseqhash3_deletion_tables_from_pivots(t,&shared);
    subseqhash3_pivot_tables_destroy(&shared);
    return status;
}
/* Bound numeric conversion explicitly; scanf integer overflow is undefined. */
static int read_integer(FILE *in,int64_t *value) {
    int c,negative=0;
    uint64_t v=0,limit;
    do { c=fgetc(in); } while (c!=EOF && isspace((unsigned char)c));
    if (c=='-' || c=='+') { negative=c=='-'; c=fgetc(in); }
    if (c<'0' || c>'9') return -1;
    limit=(uint64_t)INT64_MAX+(unsigned)negative;
    do {
        unsigned digit=(unsigned)(c-'0');
        if (v>(limit-digit)/10) return -1;
        v=10*v+digit; c=fgetc(in);
    } while (c>='0' && c<='9');
    if (c!=EOF) ungetc(c,in);
    *value=negative?(v==(uint64_t)INT64_MAX+1?INT64_MIN:-(int64_t)v):(int64_t)v;
    return 0;
}
static int read_pair(FILE *in,int64_t *x,int64_t *y) {
    int c;
    if (read_integer(in,x)) return -1;
    do { c=fgetc(in); } while (c!=EOF && isspace((unsigned char)c));
    return c!=',' || read_integer(in,y)?-1:0;
}
static int read_size(FILE *in,size_t *v) {
    int64_t x;
    if (read_integer(in,&x) || x<0 || (uint64_t)x>SIZE_MAX) return -1;
    *v=(size_t)x; return 0;
}
int subseqhash3_deletion_tables_import(SubseqHash3DeletionTables *t,
    const char *path,size_t k,size_t d) {
    FILE *in;
    size_t i,side,v;
    int ch;
    int64_t value,limit,x,y;
    unsigned char seen[256]={0};
    if (!path || tables_allocate(t,k,d)) return -1;
    in=fopen(path,"r"); if (!in) goto bad;
    limit=(INT64_MAX/4)/(int64_t)k;
    /* Original file order: normal right (A/B/C1), reversed left (revA/B/C2). */
    for (side=0;side<2;++side) {
        SubseqHash3DeletionSide *s=side?&t->left:&t->right;
        for (i=0;i<4*k*d;++i) {
            if (read_integer(in,&value) || value<1 || value>limit) goto invalid;
            s->A[i]=value;
        }
        for (i=0;i<4*k*d;++i) {
            if (read_pair(in,&x,&y) || (x!=1 && x!=-1) || (y!=1 && y!=-1)) goto invalid;
            s->B1[i]=(int8_t)x; s->A[i]*=y;
        }
        for (i=0;i<4*k;++i) {
            if (read_size(in,&v) || v>=d) goto invalid;
            s->C[i]=v;
        }
    }
    for (i=0;i<4*k;++i) {
        if (read_pair(in,&x,&y) || (x!=1 && x!=-1) || (y!=1 && y!=-1)) goto invalid;
        t->right_sign[i]=(int8_t)x; t->left_sign[i]=(int8_t)y;
    }
    for (i=0;i<4*k;++i) {
        if (read_integer(in,&value) || value<1 || value>limit) goto invalid;
        t->pivot_A[i]=value;
    }
    for (i=0;i<4*k;++i) {
        if (read_integer(in,&x) || (x!=1 && x!=-1)) goto invalid;
        t->pivot_A[i]*=x;
    }
    for (i=0;i<4*k;++i) {
        if (read_size(in,&v) || v>=d) goto invalid;
        t->pivot_C[i]=v;
    }
    for (i=0;i<k;++i) {
        if (read_size(in,&v) || v>=k || seen[v]) goto invalid;
        seen[v]=1;
    }
    do { ch=fgetc(in); } while (ch!=EOF && isspace((unsigned char)ch));
    if (ch!=EOF || ferror(in)) goto invalid;
    if (fclose(in)) goto bad;
    return 0;
invalid:
    fclose(in);
bad:
    subseqhash3_deletion_tables_destroy(t); return -1;
}

typedef struct {
    int64_t *low,*high;
    uint8_t *parent;
} SideDP;
typedef struct {
    uint16_t *low,*high;
    size_t capacity,stride;
} InwardDP;
typedef struct {
    InwardDP left,right;
    size_t anchor;
} DeletionWindowCache;
struct SubseqHash3DeletionWorkspace {
    DeletionWindowCache *window_cache;
    size_t window_start,window_slots;
    size_t n,k,e,d,stride,score_capacity,parent_capacity,rolling_capacity;
    SideDP left,right;
    int64_t *rolling;
    size_t *previous,previous_capacity;
    uint16_t *compact;
    uint8_t *packed_A; /* 127+A, exactly representable in one byte. */
    size_t packed_capacity;
    size_t compact_capacity;
    InwardDP inward_left,inward_right;
    uint32_t *inward_masks;
    size_t mask_capacity;
    int packed_parents;
    uint16_t padding[32];
};
static void *aligned_array(size_t bytes) {
    void *p=NULL;
    if (posix_memalign(&p,64,bytes)) return NULL;
    return p;
}
SubseqHash3DeletionWorkspace *subseqhash3_deletion_workspace_create(void) {
    return calloc(1,sizeof(SubseqHash3DeletionWorkspace));
}
void subseqhash3_deletion_workspace_destroy(SubseqHash3DeletionWorkspace *w) {
    if (!w) return;
    free(w->left.low); free(w->right.low);
    free(w->left.parent); free(w->right.parent);
    free(w->rolling); free(w->previous); free(w->compact); free(w->packed_A);
    free(w->inward_left.low); free(w->inward_right.low); free(w->inward_masks); free(w);
}
static int reserve(SubseqHash3DeletionWorkspace *w,size_t n,size_t k,size_t e,size_t d,int mode) {
    size_t scores,parents,rolling,stride;
    if (d>SIZE_MAX-63) return -1;
    stride=((d+7)&~(size_t)7)+8;
    if (n*(e+1)>SIZE_MAX/stride) return -1;
    scores=n*(e+1)*stride;
    if (n*k*(e+1)>SIZE_MAX/d) return -1;
    parents=n*k*(e+1)*(mode==2?((d+31)/32)*8:d);
    rolling=4*(e+1)*stride;
    if (scores>SIZE_MAX/(2*sizeof(int64_t)) || rolling>SIZE_MAX/sizeof(int64_t)) return -1;
    if ((mode==0 || mode==1) && scores>w->score_capacity) {
        free(w->left.low); free(w->right.low);
        w->left.low=aligned_array(2*scores*sizeof(int64_t));
        w->right.low=aligned_array(2*scores*sizeof(int64_t));
        w->score_capacity=0;
        if (!w->left.low || !w->right.low) return -1;
        w->score_capacity=scores;
    }
    if (mode>=1 && parents>w->parent_capacity) {
        free(w->left.parent); free(w->right.parent);
        w->left.parent=malloc(parents); w->right.parent=malloc(parents);
        w->parent_capacity=0;
        if (!w->left.parent || !w->right.parent) return -1;
        w->parent_capacity=parents;
    }
    if ((mode==0 || mode==1) && rolling>w->rolling_capacity) {
        free(w->rolling); w->rolling=aligned_array(rolling*sizeof(int64_t));
        w->rolling_capacity=0;
        if (!w->rolling) return -1;
        w->rolling_capacity=rolling;
    }
    if (d>w->previous_capacity) {
        free(w->previous); w->previous=malloc(d*sizeof(size_t));
        w->previous_capacity=0;
        if (!w->previous) return -1;
        w->previous_capacity=d;
    }
    if (d<31 && w->d!=d) for (size_t r=0;r<32;++r) w->padding[r]=(uint16_t)(r%d);
    w->n=n; w->k=k; w->e=e; w->d=d; w->stride=stride;
    w->left.high=w->left.low?w->left.low+w->score_capacity:NULL;
    w->right.high=w->right.low?w->right.low+w->score_capacity:NULL;
    return 0;
}
static void padding(int64_t *a,size_t d) {
    size_t r;
    if (d>=7) memcpy(a+d,a,7*sizeof(*a));
    else for (r=0;r<7;++r) a[d+r]=a[r%d];
}
static size_t boundary(const SubseqHash3DeletionWorkspace *w,size_t i,size_t h) {
    return (i*(w->e+1)+h)*w->stride;
}
static size_t parent_index(const SubseqHash3DeletionWorkspace *w,size_t i,size_t j,size_t h) {
    return ((i*w->k+j)*(w->e+1)+h)*w->d;
}
/* Build only full-boundary scores. Intermediate scores occupy two retained-
   length planes; traceback stores two choice bits per residue in one byte. */
static void build_side(SubseqHash3DeletionWorkspace *w,const uint8_t *s,
    const SubseqHash3DeletionSide *t,SideDP *out,int left,const uint8_t *anchors,
    const uint8_t *budgets,size_t maximum,int inward) {
    const size_t d=w->d,stride=w->stride,plane=(w->e+1)*stride;
    size_t i,j,h,r;
    for (i=0;i<(inward?1:w->n);++i) if (inward || anchors[i]) {
        const size_t budget=inward?w->e:budgets[i];
        const size_t length=inward?w->n-1:(left?i:w->n-i-1);
        const size_t last=length<maximum?length:maximum;
        int64_t *previous=w->rolling,*previous_high=previous+plane;
        int64_t *current=previous_high+plane,*current_high=current+plane;
        for (h=0;h<=budget && h<=length;++h) {
            int64_t *lo=previous+h*stride,*hi=previous_high+h*stride;
            for (r=0;r<d;++r) { lo[r]=INF; hi[r]=-INF; }
            lo[0]=hi[0]=0; padding(lo,d); padding(hi,d);
            if (inward || h==length) {
                const size_t at=boundary(w,inward?(left?h:w->n-h-1):i,h);
                memcpy(out->low+at,lo,d*sizeof(int64_t));
                memcpy(out->high+at,hi,d*sizeof(int64_t));
                if (inward) { padding(out->low+at,d); padding(out->high+at,d); }
            }
        }
        for (j=1;j<=last;++j) {
            const size_t maxh=budget<length-j?budget:length-j;
            for (h=0;h<=maxh;++h) {
                const size_t consumed=j+h;
                const unsigned a=s[inward?(left?consumed-1:w->n-consumed):
                    (left?i-consumed:i+consumed)];
                const size_t row=4*(j-1)+a,shift=t->C[row];
                int64_t *lo=current+h*stride,*hi=current_high+h*stride;
                const int64_t *pl=previous+h*stride,*ph=previous_high+h*stride;
                uint8_t *parent=inward?NULL:out->parent+parent_index(w,i,j,h);
#if DEL_SIMD
                for (r=0;r<d;r+=8) {
                    const unsigned lanes=(unsigned)(d-r<8?d-r:8);
                    const __mmask8 active=(__mmask8)((1u<<lanes)-1);
                    const size_t q=submod(r,shift,d);
                    const __m512i pinf=_mm512_set1_epi64(INF),ninf=_mm512_set1_epi64(-INF);
                    const __m512i p=_mm512_loadu_si512((const void *)(pl+q));
                    const __m512i v=_mm512_loadu_si512((const void *)(ph+q));
                    const __m128i signs=_mm_maskz_loadu_epi8((__mmask16)active,t->B1+row*d+r);
                    const __mmask8 positive=_mm512_cmp_epi64_mask(_mm512_cvtepi8_epi64(signs),_mm512_setzero_si512(),_MM_CMPINT_GT);
                    const __m512i add=_mm512_maskz_loadu_epi64(active,t->A+row*d+r);
                    const __mmask8 reachable=_mm512_cmp_epi64_mask(p,pinf,_MM_CMPINT_NE);
                    __m512i mn=_mm512_mask_blend_epi64(positive,_mm512_sub_epi64(_mm512_setzero_si512(),v),p);
                    __m512i mx=_mm512_mask_blend_epi64(positive,_mm512_sub_epi64(_mm512_setzero_si512(),p),v);
                    mn=_mm512_mask_mov_epi64(pinf,reachable,_mm512_add_epi64(mn,add));
                    mx=_mm512_mask_mov_epi64(ninf,reachable,_mm512_add_epi64(mx,add));
                    const __m512i dl=h?_mm512_load_si512((const void *)(lo-stride+r)):pinf;
                    const __m512i dh=h?_mm512_load_si512((const void *)(hi-stride+r)):ninf;
                    const __mmask8 lm=_mm512_cmp_epi64_mask(mn,dl,_MM_CMPINT_LT)&active;
                    const __mmask8 hm=_mm512_cmp_epi64_mask(mx,dh,_MM_CMPINT_GT)&active;
                    _mm512_mask_storeu_epi64(lo+r,active,_mm512_min_epi64(mn,dl));
                    _mm512_mask_storeu_epi64(hi+r,active,_mm512_max_epi64(mx,dh));
                    __m128i bits=_mm_or_si128(_mm_maskz_set1_epi8((__mmask16)lm,1),_mm_maskz_set1_epi8((__mmask16)hm,2));
                    if (!inward) _mm_mask_storeu_epi8(parent+r,(__mmask16)active,bits);
                }
#else
                for (r=0;r<d;++r) {
                    const size_t q=submod(r,shift,d),at=row*d+r;
                    int64_t mn=INF,mx=-INF;
                    uint8_t bits=0;
                    lo[r]=h?(lo-stride)[r]:INF; hi[r]=h?(hi-stride)[r]:-INF;
                    if (pl[q]!=INF) {
                        mn=(t->B1[at]>0?pl[q]:-ph[q])+t->A[at];
                        mx=(t->B1[at]>0?ph[q]:-pl[q])+t->A[at];
                    }
                    if (mn<lo[r]) { lo[r]=mn; bits|=1; }
                    if (mx>hi[r]) { hi[r]=mx; bits|=2; }
                    if (!inward) parent[r]=bits;
                }
#endif
                padding(lo,d); padding(hi,d);
                if (inward || consumed==length) {
                    const size_t at=boundary(w,inward?(left?consumed:w->n-consumed-1):i,h);
                    for (r=0;r<d;++r) {
                        const size_t dst=left?r:(r?d-r:0);
                        out->low[at+dst]=lo[r]; out->high[at+dst]=hi[r];
                    }
                    if (inward) { padding(out->low+at,d); padding(out->high+at,d); }
                }
            }
            { int64_t *tmp=previous; previous=current; current=tmp;
              tmp=previous_high; previous_high=current_high; current_high=tmp; }
        }
        /* Pad only the valid full-boundary states used by joins. */
        for (h=0;!inward && h<=budget && h<=length;++h) if (length-h<=last) {
            const size_t at=boundary(w,i,h);
            padding(out->low+at,d); padding(out->high+at,d);
        }
    }
}

/* The shared pivot family has |A| <= 127. Encode both extrema as
   nonnegative costs 127*j +/- omega. Saturating addition propagates the
   unreachable sentinel without reachability masks, and 256*254 < 65535.
   Imported wide tables retain the original exact int64 kernel. */
static void compact_padding(uint16_t *p,size_t d,const uint16_t *indices) {
#if DEL_SIMD
    if (d>=31) _mm512_storeu_si512((void *)(p+d),_mm512_loadu_si512((const void *)p));
    else {
        const __mmask32 active=(UINT32_C(1)<<d)-1;
        const __m512i values=_mm512_maskz_loadu_epi16(active,p);
        _mm512_storeu_si512((void *)(p+d),_mm512_permutexvar_epi16(
            _mm512_loadu_si512((const void *)indices),values));
    }
#else
    (void)indices;
    if (d==1) {
        for (size_t r=1;r<=31;++r) p[r]=p[0];
        return;
    }
    size_t done=0;
    while (done<31) {
        const size_t count=d<31-done?d:31-done;
        memcpy(p+d+done,p,count*sizeof(*p)); done+=count;
    }
#endif
}
/* Reflect residue order for the circular join, keeping residue zero fixed. */
static void reflect_boundary(uint16_t *dst,const uint16_t *src,size_t d) {
    size_t r=1; dst[0]=src[0];
#if DEL_SIMD
    const __m512i reverse=_mm512_set_epi16(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
        16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);
    for (;r+32<=d;r+=32) {
        const __m512i values=_mm512_loadu_si512((const void *)(src+d-r-31));
        _mm512_storeu_si512((void *)(dst+r),_mm512_permutexvar_epi16(reverse,values));
    }
#endif
    for (;r<d;++r) dst[r]=src[d-r];
}

static int build_compact(SubseqHash3DeletionWorkspace *w,const uint8_t *s,
    const SubseqHash3DeletionSide *t,SideDP *out,int left,
    const uint8_t *anchors,const uint8_t *budgets,size_t maximum) {
    const size_t d=w->d,stride=((d+31)&~(size_t)31)+32;
    const size_t plane=(w->e+1)*stride,cells=4*plane,blocks=(d+31)/32;
    InwardDP *bounds=left?&w->inward_left:&w->inward_right;
    const size_t boundary_stride=d==1?1:stride;
    const size_t boundary_cells=w->n*(w->e+1)*boundary_stride;
    if (boundary_cells>bounds->capacity) {
        uint16_t *memory=aligned_array(2*boundary_cells*sizeof(*memory));
        if (!memory) return -1;
        free(bounds->low); bounds->low=memory; bounds->capacity=boundary_cells;
    }
    bounds->high=bounds->low+bounds->capacity; bounds->stride=boundary_stride;
    size_t i,j,h,r;
    if (cells>w->compact_capacity) {
        uint16_t *memory=aligned_array(cells*sizeof(*memory));
        if (!memory) return -1;
        free(w->compact); w->compact=memory; w->compact_capacity=cells;
    }
    if (4*maximum*d>w->packed_capacity) {
        uint8_t *packed=malloc(4*maximum*d*sizeof(*packed));
        if (!packed) return -1;
        free(w->packed_A); w->packed_A=packed; w->packed_capacity=4*maximum*d;
    }
    if (4*maximum*blocks>w->mask_capacity) {
        uint32_t *memory=malloc(4*maximum*blocks*sizeof(*memory));
        if (!memory) return -1;
        free(w->inward_masks); w->inward_masks=memory; w->mask_capacity=4*maximum*blocks;
    }
    int wide=0;
#if DEL_SIMD
    __m512i magnitudes=_mm512_setzero_si512();
#endif
    for (i=0;i<4*maximum;++i) for (r=0;r<d;r+=32) {
        const size_t lanes=d-r<32?d-r:32;
        uint32_t mask=0;
#if DEL_SIMD
        const __mmask32 active=lanes==32?UINT32_MAX:(UINT32_C(1)<<lanes)-1;
        const __m256i signs=_mm256_maskz_loadu_epi8(active,t->B1+i*d+r);
        mask=_mm256_cmpgt_epi8_mask(signs,_mm256_setzero_si256());
        if (lanes==32) {
            const int64_t *a=t->A+i*d+r;
            const __m512i a0=_mm512_loadu_si512((const void *)a);
            const __m512i a1=_mm512_loadu_si512((const void *)(a+8));
            const __m512i a2=_mm512_loadu_si512((const void *)(a+16));
            const __m512i a3=_mm512_loadu_si512((const void *)(a+24));
            magnitudes=_mm512_ternarylogic_epi64(magnitudes,_mm512_abs_epi64(a0),_mm512_abs_epi64(a1),0xfe);
            __m512i abs3=_mm512_abs_epi64(a3);
            magnitudes=_mm512_ternarylogic_epi64(magnitudes,_mm512_abs_epi64(a2),abs3,0xfe);
            __m512i values=_mm512_castsi128_si512(_mm512_cvtepi64_epi16(a0));
            values=_mm512_inserti32x4(values,_mm512_cvtepi64_epi16(a1),1);
            values=_mm512_inserti32x4(values,_mm512_cvtepi64_epi16(a2),2);
            values=_mm512_inserti32x4(values,_mm512_cvtepi64_epi16(a3),3);
            _mm256_storeu_si256((void *)(w->packed_A+i*d+r),_mm512_cvtepi16_epi8(_mm512_add_epi16(values,_mm512_set1_epi16(127))));
        } else
#endif
        for (size_t lane=0;lane<lanes;++lane) {
            wide|=t->A[i*d+r+lane]<-127 || t->A[i*d+r+lane]>127;
            w->packed_A[i*d+r+lane]=(uint8_t)((uint64_t)t->A[i*d+r+lane]+127);
            mask|=(uint32_t)(t->B1[i*d+r+lane]>0)<<lane;
        }
        w->inward_masks[i*blocks+r/32]=mask;
    }
#if DEL_SIMD
    wide|=_mm512_cmp_epu64_mask(magnitudes,_mm512_set1_epi64(127),_MM_CMPINT_GT)!=0;
#endif
    if (wide) return 1;
    /* Parent stores may alias externally owned objects. Keep loop metadata
       in private copies so those stores cannot force repeated pointer loads. */
    SubseqHash3DeletionWorkspace layout __attribute__((aligned(64)))=*w;
    SubseqHash3DeletionSide table=*t;
    SideDP result=*out;
    w=&layout; t=&table; out=&result;
    for (i=0;i<w->n;++i) if (anchors[i]) {
        const size_t budget=budgets[i];
        const size_t length=left?i:w->n-i-1;
        const size_t last=length<maximum?length:maximum;
        size_t unique_residue=0;
        uint16_t *previous=w->compact,*previous_high=previous+plane;
        uint16_t *current=previous_high+plane,*current_high=current+plane;
        for (h=0;h<=budget && h<=length;++h) {
            uint16_t *lo=previous+h*stride,*hi=previous_high+h*stride;
            memset(lo,255,stride*sizeof(*lo)); memset(hi,255,stride*sizeof(*hi));
            lo[0]=hi[0]=0; compact_padding(lo,d,w->padding); compact_padding(hi,d,w->padding);
            if (h==length) {
                const size_t at=(i*(w->e+1)+h)*boundary_stride;
                memcpy(bounds->low+at,lo,d*sizeof(*lo));
                memcpy(bounds->high+at,hi,d*sizeof(*hi));
            }
        }
        for (j=1;j<=last;++j) {
            const size_t maxh=budget<length-j?budget:length-j;
            for (h=0;h<=maxh;++h) {
                const size_t consumed=j+h;
                const unsigned a=s[left?i-consumed:i+consumed];
                const size_t row=4*(j-1)+a,shift=t->C[row];
                uint16_t *lo=current+h*stride,*hi=current_high+h*stride;
                const uint16_t *pl=previous+h*stride,*ph=previous_high+h*stride;
                uint32_t *parent=(uint32_t *)(void *)(out->parent+((i*w->k+j)*(w->e+1)+h)*blocks*8);
                if (!h || j==1) {
                    /* With zero deletions there is exactly one subsequence;
                       with one retained base, the keep predecessor is residue 0.
                       Avoid sweeping unreachable residues in these sparse rows. */
                    const size_t q=h?0:unique_residue;
                    r=addmod(q,shift,d);
                    if (h) {
                        memcpy(lo,lo-stride,d*sizeof(*lo));
                        memcpy(hi,hi-stride,d*sizeof(*hi));
                        memset(parent,0,blocks*8);
                    } else {
                        memset(lo,255,d*sizeof(*lo)); memset(hi,255,d*sizeof(*hi));
                    }
                    const int positive=t->B1[row*d+r]>0,add=(int)t->A[row*d+r];
                    const unsigned mn=(positive?pl[q]:ph[q])+127+add;
                    const unsigned mx=(positive?ph[q]:pl[q])+127-add;
                    if (h) {
                        parent[2*(r/32)]|=(uint32_t)(mn<lo[r])<<(r%32);
                        parent[2*(r/32)+1]|=(uint32_t)(mx<hi[r])<<(r%32);
                        if (mn<lo[r]) lo[r]=(uint16_t)mn;
                        if (mx<hi[r]) hi[r]=(uint16_t)mx;
                    } else {
                        lo[r]=(uint16_t)mn; hi[r]=(uint16_t)mx; unique_residue=r;
                    }
                } else {
#if DEL_SIMD
                const __m512i infinity=_mm512_set1_epi16(-1);
                for (r=0;r<d;r+=32) {
                    const size_t lanes=d-r<32?d-r:32;
                    const __mmask32 active=lanes==32?UINT32_MAX:
                        (UINT32_C(1)<<lanes)-1;
                    const size_t q=submod(r,shift,d);
                    const __m512i add=_mm512_cvtepu8_epi16(_mm256_maskz_loadu_epi8(active,w->packed_A+row*d+r));
                    const __mmask32 positive=w->inward_masks[row*blocks+r/32];
                    const __m512i p=_mm512_loadu_si512((const void *)(pl+q));
                    const __m512i v=_mm512_loadu_si512((const void *)(ph+q));
                    const __m512i selected=_mm512_mask_blend_epi16(positive,v,p);
                    const __m512i other=_mm512_ternarylogic_epi32(p,v,selected,0x96);
                    const __m512i mn=_mm512_adds_epu16(selected,add);
                    const __m512i mx=_mm512_adds_epu16(other,_mm512_sub_epi16(_mm512_set1_epi16(254),add));
                    /* A deletion retains j bases, so its cost has the same bias. */
                    const __m512i dl=h?_mm512_load_si512((const void *)(lo-stride+r)):infinity;
                    const __m512i dh=h?_mm512_load_si512((const void *)(hi-stride+r)):infinity;
                    const __mmask32 lm=_mm512_cmplt_epu16_mask(mn,dl)&active;
                    const __mmask32 hm=_mm512_cmplt_epu16_mask(mx,dh)&active;
                    _mm512_store_si512((void *)(lo+r),_mm512_min_epu16(mn,dl));
                    _mm512_store_si512((void *)(hi+r),_mm512_min_epu16(mx,dh));
                    parent[2*(r/32)]=lm; parent[2*(r/32)+1]=hm;
                }
#else
                for (r=0;r<d;++r) {
                    const size_t q=submod(r,shift,d),at=row*d+r;
                    const int positive=t->B1[at]>0,add=(int)t->A[at];
                    unsigned mn=positive?pl[q]:ph[q],mx=positive?ph[q]:pl[q];
                    const unsigned dl=h?(lo-stride)[r]:UINT16_MAX;
                    const unsigned dh=h?(hi-stride)[r]:UINT16_MAX;
                    if (mn!=UINT16_MAX) mn+=127+add;
                    if (mx!=UINT16_MAX) mx+=127-add;
                    lo[r]=(uint16_t)(mn<dl?mn:dl); hi[r]=(uint16_t)(mx<dh?mx:dh);
                    if (!(r%32)) parent[2*(r/32)]=parent[2*(r/32)+1]=0;
                    parent[2*(r/32)]|=(uint32_t)(mn<dl)<<(r%32);
                    parent[2*(r/32)+1]|=(uint32_t)(mx<dh)<<(r%32);
                }
#endif
                }
                compact_padding(lo,d,w->padding); compact_padding(hi,d,w->padding);
                if (consumed==length) {
                    const size_t at=(i*(w->e+1)+h)*boundary_stride;
                    if (left) {
                        memcpy(bounds->low+at,lo,d*sizeof(*lo));
                        memcpy(bounds->high+at,hi,d*sizeof(*hi));
                    } else {
                        reflect_boundary(bounds->low+at,lo,d);
                        reflect_boundary(bounds->high+at,hi,d);
                    }
                }
            }
            { uint16_t *tmp=previous; previous=current; current=tmp;
              tmp=previous_high; previous_high=current_high; current_high=tmp; }
        }
        for (h=0;h<=budget && h<=length;++h) if (length-h<=last) {
            const size_t at=(i*(w->e+1)+h)*boundary_stride;
            if (d>1) { compact_padding(bounds->low+at,d,w->padding); compact_padding(bounds->high+at,d,w->padding); }
        }
    }
    return 0;
}

/* Compact inward DPs retain every prefix/suffix state for score-based
   traceback. Right residues and tables are reflected for circular joins. */
static size_t inward_state(const InwardDP *dp,size_t l,size_t h,size_t edits) {
    return (l*(edits+1)+h)*dp->stride;
}
/* A one-residue DP has no circular vectors to build or pad. Keep the
   compact scalar states contiguous, retaining the same strict deletion ties. */
static int build_inward_one(SubseqHash3DeletionWorkspace *w,InwardDP *dp,
    const uint8_t *s,const SubseqHash3DeletionSide *t,size_t maximum,int left) {
    const size_t n=w->n,e=w->e,cells=n*(e+1);
    if (cells>dp->capacity) {
        uint16_t *memory=aligned_array(2*cells*sizeof(*memory));
        if (!memory) return -1;
        free(dp->low); dp->low=memory; dp->capacity=cells;
    }
    dp->high=dp->low+dp->capacity; dp->stride=1;
    for (size_t h=0;h<=e;++h) {
        const size_t at=inward_state(dp,h,h,e);
        dp->low[at]=dp->high[at]=0;
    }
    for (size_t j=1;j<=maximum;++j) for (size_t h=0;h<=e && j+h<n;++h) {
        const size_t l=j+h,row=4*(j-1)+s[left?l-1:n-l];
        const int64_t add=t->A[row];
        if (add<-127 || add>127) return 1;
        const int positive=t->B1[row]>0;
        const size_t dst=inward_state(dp,l,h,e),src=inward_state(dp,l-1,h,e);
        const unsigned mn=(positive?dp->low[src]:dp->high[src])+127+(int)add;
        const unsigned mx=(positive?dp->high[src]:dp->low[src])+127-(int)add;
        const size_t del=h?inward_state(dp,l-1,h-1,e):0;
        const unsigned dl=h?dp->low[del]:UINT16_MAX,dh=h?dp->high[del]:UINT16_MAX;
        dp->low[dst]=(uint16_t)(mn<dl?mn:dl); dp->high[dst]=(uint16_t)(mx<dh?mx:dh);
    }
    return 0;
}

/* Preparation is scoped to one call: imported tables can be changed between
   calls using the same public workspace. */
static int prepare_inward(SubseqHash3DeletionWorkspace *w,
    const SubseqHash3DeletionSide *t,size_t maximum,int left) {
    const size_t d=w->d,blocks=(d+31)/32,rows=4*maximum;
    size_t j,r;
    if (rows*d>w->packed_capacity) {
        uint8_t *memory=malloc(rows*d*sizeof(*memory));
        if (!memory) return -1;
        free(w->packed_A); w->packed_A=memory; w->packed_capacity=rows*d;
    }
    if (rows*blocks>w->mask_capacity) {
        uint32_t *memory=malloc(rows*blocks*sizeof(*memory));
        if (!memory) return -1;
        free(w->inward_masks); w->inward_masks=memory; w->mask_capacity=rows*blocks;
    }
    int wide=0;
#if DEL_SIMD
    __m512i magnitudes=_mm512_setzero_si512();
#endif
    for (j=0;j<rows;++j) for (r=0;r<d;r+=32) {
        uint32_t mask=0;
        const size_t lanes=d-r<32?d-r:32;
#if DEL_SIMD
        if (lanes==32) {
            /* Four wide loads narrow to one cost vector. The next unused
               table row supplies one harmless lane at the right wrap. */
            const size_t start=left?r:d-r-31;
            const int64_t *a=t->A+j*d+start;
            const __m512i a0=_mm512_loadu_si512((const void *)a);
            const __m512i a1=_mm512_loadu_si512((const void *)(a+8));
            const __m512i a2=_mm512_loadu_si512((const void *)(a+16));
            const __m512i a3=_mm512_loadu_si512((const void *)(a+24));
            magnitudes=_mm512_ternarylogic_epi64(magnitudes,_mm512_abs_epi64(a0),_mm512_abs_epi64(a1),0xfe);
            __m512i abs3=_mm512_abs_epi64(a3);
            if (!left && !r) {
                abs3=_mm512_maskz_mov_epi64(0x7f,abs3);
                wide|=t->A[j*d]<-127 || t->A[j*d]>127;
            }
            magnitudes=_mm512_ternarylogic_epi64(magnitudes,_mm512_abs_epi64(a2),abs3,0xfe);
            __m512i values=_mm512_castsi128_si512(_mm512_cvtepi64_epi16(a0));
            values=_mm512_inserti32x4(values,_mm512_cvtepi64_epi16(a1),1);
            values=_mm512_inserti32x4(values,_mm512_cvtepi64_epi16(a2),2);
            values=_mm512_inserti32x4(values,_mm512_cvtepi64_epi16(a3),3);
            __m256i signs=_mm256_loadu_si256((const void *)(t->B1+j*d+start));
            if (!left) {
                const __m512i reverse=_mm512_set_epi16(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                    16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);
                values=_mm512_permutexvar_epi16(reverse,values);
                const __m256i bytes=_mm256_setr_epi8(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,
                    15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
                signs=_mm256_permute4x64_epi64(_mm256_shuffle_epi8(signs,bytes),0x4e);
                if (!r) {
                    values=_mm512_mask_set1_epi16(values,1,(short)t->A[j*d]);
                    signs=_mm256_mask_set1_epi8(signs,1,t->B1[j*d]);
                }
            }
            _mm256_storeu_si256((void *)(w->packed_A+j*d+r),_mm512_cvtepi16_epi8(_mm512_add_epi16(values,_mm512_set1_epi16(127))));
            w->inward_masks[j*blocks+r/32]=_mm256_cmpgt_epi8_mask(signs,_mm256_setzero_si256());
            continue;
        }
#endif
        for (size_t lane=0;lane<lanes;++lane) {
            size_t q=r+lane;
            const size_t tr=left?q:(q?d-q:0);
            wide|=t->A[j*d+tr]<-127 || t->A[j*d+tr]>127;
            w->packed_A[j*d+q]=(uint8_t)((uint64_t)t->A[j*d+tr]+127);
            mask|=(uint32_t)(t->B1[j*d+tr]>0)<<lane;
        }
        w->inward_masks[j*blocks+r/32]=mask;
    }
#if DEL_SIMD
    wide|=_mm512_cmp_epu64_mask(magnitudes,_mm512_set1_epi64(127),_MM_CMPINT_GT)!=0;
#endif
    return wide?1:0;
}

static int build_inward(SubseqHash3DeletionWorkspace *w,InwardDP *dp,
    const uint8_t *s,const SubseqHash3DeletionSide *t,size_t maximum,int left,
    const SubseqHash3DeletionWorkspace *prepared) {
    if (w->d==1) return build_inward_one(w,dp,s,t,maximum,left);
    const size_t d=w->d,n=w->n,e=w->e,stride=((d+31)&~(size_t)31)+32;
    const size_t cells=n*(e+1)*stride,blocks=(d+31)/32;
    size_t j,h,r;
    if (cells>SIZE_MAX/(2*sizeof(uint16_t))) return -1;
    if (cells>dp->capacity) {
        uint16_t *memory=aligned_array(2*cells*sizeof(*memory));
        if (!memory) return -1;
        free(dp->low); dp->low=memory; dp->capacity=cells;
    }
    dp->high=dp->low+dp->capacity; dp->stride=stride;
    if (!prepared) {
        int status=prepare_inward(w,t,maximum,left);
        if (status) return status;
        prepared=w;
    }
    uint16_t *low=dp->low,*high=dp->high;
    const uint8_t *packed=prepared->packed_A;
    const uint32_t *masks=prepared->inward_masks;
    const size_t *shifts=t->C;
    for (h=0;h<=e;++h) {
        const size_t at=inward_state(dp,h,h,e);
        memset(low+at,255,stride*sizeof(*low)); memset(high+at,255,stride*sizeof(*high));
        low[at]=high[at]=0; compact_padding(low+at,d,w->padding); compact_padding(high+at,d,w->padding);
    }
    size_t unique_residue=0;
    for (j=1;j<=maximum;++j) for (h=0;h<=e && j+h<n;++h) {
        const size_t l=j+h,row=4*(j-1)+s[left?l-1:n-l];
        const size_t shift=left || !shifts[row]?shifts[row]:d-shifts[row];
        const size_t dst=inward_state(dp,l,h,e),src=inward_state(dp,l-1,h,e);
        uint16_t *lo=low+dst,*hi=high+dst;
        const uint16_t *pl=low+src,*ph=high+src;
        const size_t deletion=h?inward_state(dp,l-1,h-1,e):0;
        if (!h || j==1) {
            const size_t q=h?0:unique_residue;
            r=addmod(q,shift,d);
            if (h) {
                memcpy(lo,low+deletion,d*sizeof(*lo)); memcpy(hi,high+deletion,d*sizeof(*hi));
            } else { memset(lo,255,d*sizeof(*lo)); memset(hi,255,d*sizeof(*hi)); }
            const int positive=(masks[row*blocks+r/32]>>(r%32))&1u;
            const int add=(int)packed[row*d+r]-127;
            const unsigned mn=(positive?pl[q]:ph[q])+127+add;
            const unsigned mx=(positive?ph[q]:pl[q])+127-add;
            if (mn<lo[r]) lo[r]=(uint16_t)mn;
            if (mx<hi[r]) hi[r]=(uint16_t)mx;
            if (!h) unique_residue=r;
        } else {
#if DEL_SIMD
            const __m512i infinity=_mm512_set1_epi16(-1);
            for (r=0;r<d;r+=32) {
                const size_t lanes=d-r<32?d-r:32;
                const __mmask32 active=lanes==32?UINT32_MAX:(UINT32_C(1)<<lanes)-1;
                const __mmask32 positive=masks[row*blocks+r/32];
                const size_t q=submod(r,shift,d);
                const __m512i add=_mm512_cvtepu8_epi16(_mm256_maskz_loadu_epi8(active,packed+row*d+r));
                const __m512i p=_mm512_loadu_si512((const void *)(pl+q));
                const __m512i v=_mm512_loadu_si512((const void *)(ph+q));
                const __m512i selected=_mm512_mask_blend_epi16(positive,v,p);
                const __m512i other=_mm512_ternarylogic_epi32(p,v,selected,0x96);
                const __m512i mn=_mm512_adds_epu16(selected,add);
                const __m512i mx=_mm512_adds_epu16(other,_mm512_sub_epi16(_mm512_set1_epi16(254),add));
                const __m512i dl=h?_mm512_load_si512((const void *)(low+deletion+r)):infinity;
                const __m512i dh=h?_mm512_load_si512((const void *)(high+deletion+r)):infinity;
                _mm512_store_si512((void *)(lo+r),_mm512_min_epu16(mn,dl));
                _mm512_store_si512((void *)(hi+r),_mm512_min_epu16(mx,dh));
            }
#else
            for (r=0;r<d;++r) {
                const size_t q=submod(r,shift,d);
                const int positive=(masks[row*blocks+r/32]>>(r%32))&1u;
                const int add=(int)packed[row*d+r]-127;
                unsigned mn=positive?pl[q]:ph[q],mx=positive?ph[q]:pl[q];
                const unsigned dl=h?low[deletion+r]:UINT16_MAX,dh=h?high[deletion+r]:UINT16_MAX;
                if (mn!=UINT16_MAX) mn+=127+add;
                if (mx!=UINT16_MAX) mx+=127-add;
                lo[r]=(uint16_t)(mn<dl?mn:dl); hi[r]=(uint16_t)(mx<dh?mx:dh);
            }
#endif
        }
        compact_padding(lo,d,w->padding); compact_padding(hi,d,w->padding);
    }
    return 0;
}
static unsigned inward_join(const uint16_t *l,const uint16_t *r,size_t d,size_t shift) {
    if (d==1) {
        const unsigned cost=(unsigned)l[0]+r[0];
        return cost<UINT16_MAX?cost:UINT16_MAX;
    }
    unsigned best=UINT16_MAX;
#if DEL_SIMD
    __m512i minimum=_mm512_set1_epi16(-1);
    for (size_t q=0;q<d;q+=32) {
        const size_t lanes=d-q<32?d-q:32;
        const __mmask32 active=lanes==32?UINT32_MAX:(UINT32_C(1)<<lanes)-1;
        const __m512i cost=_mm512_adds_epu16(_mm512_load_si512((const void *)(l+q)),
            _mm512_loadu_si512((const void *)(r+addmod(q,shift,d))));
        minimum=_mm512_mask_min_epu16(minimum,active,minimum,cost);
    }
    const __m256i m=_mm256_min_epu16(_mm512_castsi512_si256(minimum),
        _mm512_extracti64x4_epi64(minimum,1));
    const __m128i a=_mm_min_epu16(_mm256_castsi256_si128(m),_mm256_extracti128_si256(m,1));
    best=(unsigned)_mm_extract_epi16(_mm_minpos_epu16(a),0);
#else
    for (size_t q=0;q<d;++q) {
        const unsigned cost=(unsigned)l[q]+r[addmod(q,shift,d)];
        if (cost<best) best=cost;
    }
#endif
    return best;
}
static size_t inward_first(SubseqHash3DeletionWorkspace *w,const uint16_t *l,
    const uint16_t *r,size_t c) {
    size_t q,last=w->d,best=w->d;
    for (q=w->d;q>0;--q) if (r[q-1]!=UINT16_MAX) { last=q-1; break; }
    if (last==w->d) return w->d;
    for (q=0;q<w->d;++q) {
        if (r[q]!=UINT16_MAX) last=q;
        w->previous[q]=last;
    }
    for (q=0;q<w->d;++q) if (l[q]!=UINT16_MAX) {
        const size_t x=addmod(q,c,w->d),psi=submod(x,w->previous[x],w->d);
        if (psi<best) best=psi;
    }
    return best;
}
static void inward_trace_compact(const uint8_t *s,size_t n,size_t e,size_t d,
    const SubseqHash3DeletionSide *t,const InwardDP *dp,size_t l,size_t h,
    size_t residue,int high,int left,SubseqHash3DeletionResult *out) {
    size_t j=l-h;
    while (j) {
        if (!h) {
            while (j) {
                const size_t pos=left?j-1:n-j,dst=left?j-1:out->length-j;
                out->seed[dst]=s[pos]; out->positions[dst]=pos; --j;
            }
            break;
        }
        const size_t r=left?residue:(residue?d-residue:0);
        const uint16_t *scores=high?dp->high:dp->low;
        if (h && scores[inward_state(dp,l,h,e)+r]==scores[inward_state(dp,l-1,h-1,e)+r]) {
            --l; --h; continue;
        }
        const size_t pos=left?l-1:n-l,dst=left?j-1:out->length-j;
        const size_t row=4*(j-1)+s[pos];
        out->seed[dst]=s[pos]; out->positions[dst]=pos;
        if (t->B1[row*d+residue]<0) high=!high;
        residue=submod(residue,t->C[row],d); --j; --l;
    }
}
static void trace(const SubseqHash3DeletionWorkspace *,const uint8_t *,
    const SubseqHash3DeletionSide *,const SideDP *,size_t,size_t,size_t,size_t,int,int,size_t,
    SubseqHash3DeletionResult *);
/* Recover an outward witness from the full prefix scores of one cached
   global anchor. Equal scores prefer deletion, matching the parent kernel. */
static void trace_window_side(const uint8_t *s,size_t e,size_t d,
    const SubseqHash3DeletionSide *t,const InwardDP *dp,size_t anchor,
    size_t pivot,size_t l,size_t h,size_t residue,int high,int left,
    SubseqHash3DeletionResult *out) {
    size_t j=l-h;
    while (j) {
        if (!h) {
            while (j) {
                const size_t pos=left?anchor-j:anchor+j;
                const size_t dst=left?pivot-j:pivot+j;
                out->seed[dst]=s[pos]; out->positions[dst]=pos; --j;
            }
            break;
        }
        const size_t r=left?residue:(residue?d-residue:0);
        const uint16_t *scores=high?dp->high:dp->low;
        if (h && scores[inward_state(dp,l,h,e)+r]==
                scores[inward_state(dp,l-1,h-1,e)+r]) { --l; --h; continue; }
        const size_t pos=left?anchor-l:anchor+l,dst=left?pivot-j:pivot+j;
        const size_t row=4*(j-1)+s[pos];
        out->seed[dst]=s[pos]; out->positions[dst]=pos;
        if (t->B1[row*d+residue]<0) high=!high;
        residue=submod(residue,t->C[row],d); --j; --l;
    }
}

static const uint16_t *window_scores(const SubseqHash3DeletionWorkspace *w,
    size_t i,size_t h,int left,int outward,int high) {
    const InwardDP *dp=left?&w->inward_left:&w->inward_right;
    size_t l=left || outward?i:w->n-i-1;
    if (w->window_cache) {
        const DeletionWindowCache *c=w->window_cache+(w->window_start+i)%w->window_slots;
        dp=left?&c->left:&c->right; l=left?i:w->n-i-1;
    }
    return (high?dp->high:dp->low)+inward_state(dp,l,h,w->e);
}

/* Return the first right residue in the original strict witness order. */
static size_t compact_witness(const uint16_t *l,const uint16_t *r,size_t d,
    size_t offset,unsigned cost) {
    if ((unsigned)l[offset]+r[0]==cost) return 0;
#if DEL_SIMD
    const __m512i target=_mm512_set1_epi16((short)cost);
    for (size_t end=d;end>1;) {
        const size_t begin=end>33?end-32:1,lanes=end-begin;
        const __mmask32 active=lanes==32?UINT32_MAX:(UINT32_C(1)<<lanes)-1;
        const __m512i sum=_mm512_adds_epu16(
            _mm512_loadu_si512((const void *)(l+addmod(begin,offset,d))),
            _mm512_loadu_si512((const void *)(r+begin)));
        const uint32_t equal=active&_mm512_cmpeq_epi16_mask(sum,target);
        if (equal) return d-(begin+31-(size_t)__builtin_clz(equal));
        end=begin;
    }
#else
    for (size_t q=1;q<d;++q)
        if ((unsigned)l[submod(offset,q,d)]+r[d-q]==cost) return q;
#endif
    return d;
}

static int compute_compact_seeds(const uint8_t *s,const SubseqHash3DeletionTables *t,
    const size_t *pivots,size_t count,SubseqHash3DeletionResult *out,
    SubseqHash3DeletionWorkspace *w,size_t maxleft,size_t maxright,int outward) {
    const size_t n=w->n,k=w->k,e=w->e,d=w->d;
    if (!outward) {
        int status=build_inward(w,&w->inward_left,s,&t->left,maxleft,1,NULL);
        if (!status) status=build_inward(w,&w->inward_right,s,&t->right,maxright,0,NULL);
        if (status) return status;
    }
    for (size_t a=0;a<count;++a) {
        const size_t p=pivots?pivots[a]:a;
        size_t bestpsi=0,besti=0,bestside=UINT16_MAX,bestq=d;
        int64_t best=-INF;
        /* First try psi=0, then find the exact reachable minimum if needed. */
        for (unsigned pass=0;pass<2;++pass) {
            if (pass) bestpsi=d;
            for (size_t h=0;h<=e;++h) {
                const size_t i=p+h,row=4*p+s[i];
                const uint16_t *l=window_scores(w,i,h,1,outward,t->left_sign[row]>0);
                const uint16_t *r=window_scores(w,i,e-h,0,outward,t->right_sign[row]>0);
                const size_t psi=pass?inward_first(w,l,r,t->pivot_C[row]):0;
                if (psi>bestpsi || psi==d) continue;
                const size_t shift=submod(t->pivot_C[row],psi,d);
                const unsigned cost=p==0?(unsigned)l[0]+r[shift]:
                    p==k-1?(unsigned)l[shift?d-shift:0]+r[0]:inward_join(l,r,d,shift);
                if (cost>=UINT16_MAX) continue;
                const int64_t score=127*(int64_t)(k-1)-(int64_t)cost+t->pivot_A[row];
                if (psi<bestpsi || score>best) {
                    bestpsi=psi; best=score; besti=i; bestside=cost;
                }
            }
            if (best!=-INF) break;
        }
        if (best==-INF) goto failed;
        const size_t h=besti-p,row=4*p+s[besti];
        const uint16_t *l=window_scores(w,besti,h,1,outward,t->left_sign[row]>0);
        const uint16_t *r=window_scores(w,besti,e-h,0,outward,t->right_sign[row]>0);
        bestq=compact_witness(l,r,d,submod(bestpsi,t->pivot_C[row],d),
            (unsigned)bestside);
        if (bestq==d) goto failed;
        out[a].length=k; out[a].deletions=e; out[a].psi=bestpsi; out[a].omega=best;
        out[a].seed=malloc(k); out[a].positions=malloc(k*sizeof(size_t));
        if (!out[a].seed || !out[a].positions) goto failed;
        out[a].seed[p]=s[besti]; out[a].positions[p]=besti;
        if (outward && w->window_cache) {
            DeletionWindowCache *c=w->window_cache+(w->window_start+besti)%w->window_slots;
            trace_window_side(s,e,d,&t->left,&c->left,besti,p,besti,h,
                submod(submod(bestpsi,t->pivot_C[row],d),bestq,d),t->left_sign[row]>0,1,out+a);
            trace_window_side(s,e,d,&t->right,&c->right,besti,p,n-besti-1,e-h,
                bestq,t->right_sign[row]>0,0,out+a);
        } else if (outward) {
            trace(w,s,&t->left,&w->left,besti,p,h,
                submod(submod(bestpsi,t->pivot_C[row],d),bestq,d),t->left_sign[row]>0,1,p,out+a);
            trace(w,s,&t->right,&w->right,besti,k-p-1,e-h,bestq,t->right_sign[row]>0,0,p,out+a);
        } else {
        inward_trace_compact(s,n,e,d,&t->left,&w->inward_left,besti,h,
            submod(submod(bestpsi,t->pivot_C[row],d),bestq,d),t->left_sign[row]>0,1,out+a);
        inward_trace_compact(s,n,e,d,&t->right,&w->inward_right,n-besti-1,e-h,bestq,t->right_sign[row]>0,0,out+a);
        }
    }
    return 0;
failed:
    subseqhash3_deletion_results_destroy(out,count); return -1;
}

static int64_t join(const int64_t *l,const int64_t *r,size_t d,size_t offset,
    int ls,int rs,int64_t pivot) {
    int64_t best=-INF;
    size_t q;
#if DEL_SIMD
    __m512i result=_mm512_set1_epi64(-INF);
    const __m512i li=_mm512_set1_epi64(ls>0?-INF:INF);
    const __m512i ri=_mm512_set1_epi64(rs>0?-INF:INF);
    for (q=0;q<d;q+=8) {
        const __mmask8 active=(__mmask8)((1u<<(d-q<8?d-q:8))-1);
        const __m512i a=_mm512_loadu_si512((const void *)(l+q));
        const __m512i b=_mm512_loadu_si512((const void *)(r+addmod(q,offset,d)));
        const __mmask8 valid=active&_mm512_cmp_epi64_mask(a,li,_MM_CMPINT_NE)&_mm512_cmp_epi64_mask(b,ri,_MM_CMPINT_NE);
        __m512i value=_mm512_add_epi64(ls>0?a:_mm512_sub_epi64(_mm512_setzero_si512(),a),rs>0?b:_mm512_sub_epi64(_mm512_setzero_si512(),b));
        result=_mm512_mask_max_epi64(result,valid,result,value);
    }
    best=_mm512_reduce_max_epi64(result);
    if (best!=-INF) best+=pivot;
#else
    for (q=0;q<d;++q) {
        const int64_t a=l[q],b=r[addmod(q,offset,d)];
        if (a==(ls>0?-INF:INF) || b==(rs>0?-INF:INF)) continue;
        const int64_t v=ls*a+rs*b+pivot;
        if (v>best) best=v;
    }
#endif
    return best;
}
static size_t first_residue(SubseqHash3DeletionWorkspace *w,const int64_t *l,
    const int64_t *r,size_t c,int ls,int rs) {
    size_t q,last=w->d,best=w->d;
    const int64_t li=ls>0?-INF:INF,ri=rs>0?-INF:INF;
    for (q=w->d;q>0;--q) if (r[q-1]!=ri) { last=q-1; break; }
    if (last==w->d) return w->d;
    for (q=0;q<w->d;++q) {
        if (r[q]!=ri) last=q;
        w->previous[q]=last;
    }
    for (q=0;q<w->d;++q) if (l[q]!=li) {
        const size_t x=addmod(q,c,w->d);
        const size_t psi=submod(x,w->previous[x],w->d);
        if (psi<best) best=psi;
    }
    return best;
}
static void trace(const SubseqHash3DeletionWorkspace *w,const uint8_t *s,
    const SubseqHash3DeletionSide *t,const SideDP *dp,size_t i,size_t j,
    size_t h,size_t residue,int high,int left,size_t pivot,
    SubseqHash3DeletionResult *out) {
    while (j) {
        if (!h) {
            while (j) {
                const size_t pos=left?i-j:i+j;
                const size_t dst=left?pivot-j:pivot+j;
                out->seed[dst]=s[pos]; out->positions[dst]=pos; --j;
            }
            break;
        }
        const size_t at=parent_index(w,i,j,h)+residue;
        int keep;
        if (w->packed_parents) {
            const size_t blocks=(w->d+31)/32;
            const uint32_t *bits=(const uint32_t *)(const void *)(dp->parent+((i*w->k+j)*(w->e+1)+h)*blocks*8);
            keep=(int)((bits[2*(residue/32)+(high?1:0)]>>(residue%32))&1u);
        } else keep=!!(dp->parent[at]&(high?2:1));
        if (keep) {
            const size_t pos=left?i-j-h:i+j+h;
            const unsigned a=s[pos];
            const size_t row=4*(j-1)+a;
            const size_t dst=left?pivot-j:pivot+j;
            out->seed[dst]=(uint8_t)a; out->positions[dst]=pos;
            if (t->B1[row*w->d+residue]<0) high=!high;
            residue=submod(residue,t->C[row],w->d); --j;
        } else { --h; }
    }
}
/* Inward states are shared prefixes/suffixes. Recover strict deletion ties
   from their scores, so no per-residue parent storage is needed. */
static void trace_inward(const SubseqHash3DeletionWorkspace *w,const uint8_t *s,
    const SubseqHash3DeletionSide *t,const SideDP *dp,size_t i,size_t j,
    size_t h,size_t residue,int high,int left,size_t pivot,
    SubseqHash3DeletionResult *out) {
    (void)pivot;
    while (j) {
        const size_t r=left?residue:(residue?w->d-residue:0);
        const int64_t *scores=high?dp->high:dp->low;
        const size_t previous=left?i-1:i+1;
        if (h && scores[boundary(w,i,h)+r]==scores[boundary(w,previous,h-1)+r]) {
            i=previous; --h; continue;
        }
        const size_t pos=left?i-1:i+1;
        const unsigned a=s[pos];
        const size_t row=4*(j-1)+a;
        const size_t dst=left?j-1:out->length-j;
        out->seed[dst]=(uint8_t)a; out->positions[dst]=pos;
        if (t->B1[row*w->d+residue]<0) high=!high;
        residue=submod(residue,t->C[row],w->d); --j; i=previous;
    }
}

void subseqhash3_deletion_results_destroy(SubseqHash3DeletionResult *out,size_t count) {
    size_t i;
    if (!out) return;
    for (i=0;i<count;++i) {
        free(out[i].seed); free(out[i].positions); memset(out+i,0,sizeof(*out));
    }
}
/* No combinatorial search is needed for E=0 or k=1. */
static int direct(const uint8_t *s,size_t n,size_t e,
    const SubseqHash3DeletionTables *t,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *out,int inward) {
    const size_t k=n-e,d=t->d;
    size_t a,j;
    for (a=0;a<count;++a) {
        const size_t p=pivots?pivots[a]:a;
        size_t psi=0,chosen=0;
        int64_t omega=0;
        if (k==1) {
            psi=d; omega=-INF;
            for (j=0;j<n;++j) {
                const size_t row=s[j],q=t->pivot_C[row];
                const int64_t v=t->pivot_A[row];
                if (q<psi || (q==psi && v>omega)) {
                    psi=q; omega=v; chosen=j;
                }
            }
        } else {
            size_t lr=0,rr=0;
            int64_t lv=0,rv=0;
            for (j=0;j<p;++j) {
                const size_t row=4*j+s[inward?j:p-1-j];
                lr=addmod(lr,t->left.C[row],d);
                lv=t->left.B1[row*d+lr]*lv+t->left.A[row*d+lr];
            }
            for (j=0;j<k-p-1;++j) {
                const size_t row=4*j+s[inward?k-1-j:p+1+j];
                rr=addmod(rr,t->right.C[row],d);
                rv=t->right.B1[row*d+rr]*rv+t->right.A[row*d+rr];
            }
            const size_t row=4*p+s[p];
            psi=addmod(addmod(lr,rr,d),t->pivot_C[row],d);
            omega=t->left_sign[row]*lv+t->right_sign[row]*rv+t->pivot_A[row];
        }
        out[a].length=k; out[a].deletions=e; out[a].psi=psi; out[a].omega=omega;
        out[a].seed=malloc(k); out[a].positions=malloc(k*sizeof(size_t));
        if (!out[a].seed || !out[a].positions) {
            subseqhash3_deletion_results_destroy(out,count); return -1;
        }
        for (j=0;j<k;++j) { out[a].seed[j]=s[k==1?chosen:j]; out[a].positions[j]=k==1?chosen:j; }
    }
    return 0;
}
static int compute_deletion_seeds(const uint8_t *s,size_t n,
    size_t e,const SubseqHash3DeletionTables *t,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *out,SubseqHash3DeletionWorkspace *w,int inward) {
    size_t k,a,p,h,i,maxleft=0,maxright=0;
    uint8_t anchors[256]={0},left_budget[256]={0},right_budget[256]={0};
    if (!s || !t || !out || !w || !n || n>256 || e>=n || !count) return -1;
    k=n-e;
    if (count>k || t->k<k || !t->d || !t->left.A || !t->right.A ||
        !t->left.B1 || !t->right.B1 || !t->left.C || !t->right.C ||
        !t->pivot_A || !t->pivot_C || !t->left_sign || !t->right_sign) return -1;
    for (i=0;i<n;++i) if (s[i]>3) return -1;
    for (a=0;a<count;++a) {
        p=pivots?pivots[a]:a;
        if (p>=k) return -1;
        if (p>maxleft) maxleft=p;
        if (k-p-1>maxright) maxright=k-p-1;
        for (h=0;h<=e;++h) {
            i=p+h; anchors[i]=1;
            if (h>left_budget[i]) left_budget[i]=(uint8_t)h;
            if (e-h>right_budget[i]) right_budget[i]=(uint8_t)(e-h);
        }
    }
    if (!e || k==1) return direct(s,n,e,t,pivots,count,out,inward);
    if (e==1 && t->d>=128 && t->d>=(inward?4:8)*n) {
        const SingleSide left={t->left.A,NULL,t->left.B1,t->left.C};
        const SingleSide right={t->right.A,NULL,t->right.B1,t->right.C};
        SingleChoice choices[256];
        single_deletion_choices(s,n,t->d,&left,&right,t->pivot_A,NULL,
            t->left_sign,t->right_sign,t->pivot_C,pivots,count,inward,choices);
        for (a=0;a<count;++a) {
            out[a].length=k; out[a].deletions=1;
            out[a].psi=choices[a].psi; out[a].omega=choices[a].omega;
            out[a].seed=malloc(k); out[a].positions=malloc(k*sizeof(size_t));
            if (!out[a].seed || !out[a].positions) goto failed;
            for (i=0;i<k;++i) {
                const size_t pos=i+(i>=choices[a].deleted);
                out[a].seed[i]=s[pos]; out[a].positions[i]=pos;
            }
        }
        return 0;
    }
    /* Range-check while packing the rows actually consumed by the compact
       kernels. A wide row requests the exact int64 fallback before traceback. */
    w->packed_parents=!inward;
    if (reserve(w,n,k,e,t->d,inward?-1:2)) return -1;
    int status;
    if (inward) status=compute_compact_seeds(s,t,pivots,count,out,w,maxleft,maxright,0);
    else {
        status=build_compact(w,s,&t->left,&w->left,1,anchors,left_budget,maxleft);
        if (!status) status=build_compact(w,s,&t->right,&w->right,0,anchors,right_budget,maxright);
        if (!status) return compute_compact_seeds(s,t,pivots,count,out,w,maxleft,maxright,1);
    }
    if (status<=0) return status;
    w->packed_parents=0;
    if (reserve(w,n,k,e,t->d,inward?0:1)) return -1;
    build_side(w,s,&t->left,&w->left,1,anchors,left_budget,maxleft,inward);
    build_side(w,s,&t->right,&w->right,0,anchors,right_budget,maxright,inward);
    for (a=0;a<count;++a) {
        size_t bestpsi=0,besti=0,bestq=0;
        int64_t best=-INF;
        p=pivots?pivots[a]:a;
        /* Most interesting windows reach psi=0. Only calculate the general
           O(D) circular reachability join if NONE of their pivots do. */
        for (h=0;h<=e;++h) {
            i=p+h;
            const size_t row=4*p+s[i];
            const int ls=t->left_sign[row],rs=t->right_sign[row];
            const int64_t *l=(ls>0?w->left.high:w->left.low)+boundary(w,i,h);
            const int64_t *r=(rs>0?w->right.high:w->right.low)+boundary(w,i,e-h);
            const int64_t v=join(l,r,t->d,t->pivot_C[row],ls,rs,t->pivot_A[row]);
            if (v>best) { best=v; besti=i; }
        }
        if (best==-INF) {
            bestpsi=t->d;
            for (h=0;h<=e;++h) {
                i=p+h;
                const size_t row=4*p+s[i];
                const int ls=t->left_sign[row],rs=t->right_sign[row];
                const int64_t *l=(ls>0?w->left.high:w->left.low)+boundary(w,i,h);
                const int64_t *r=(rs>0?w->right.high:w->right.low)+boundary(w,i,e-h);
                const size_t psi=first_residue(w,l,r,t->pivot_C[row],ls,rs);
                if (psi<=bestpsi && psi<t->d) {
                    const size_t off=submod(t->pivot_C[row],psi,t->d);
                    const int64_t v=join(l,r,t->d,off,ls,rs,t->pivot_A[row]);
                    if (psi<bestpsi || v>best) { bestpsi=psi; best=v; besti=i; }
                }
            }
        }
        if (best==-INF) goto failed;
        i=besti; h=i-p;
        const size_t row=4*p+s[i];
        const int ls=t->left_sign[row],rs=t->right_sign[row];
        const int64_t *l=(ls>0?w->left.high:w->left.low)+boundary(w,i,h);
        const int64_t *r=(rs>0?w->right.high:w->right.low)+boundary(w,i,e-h);
        /* Match the reference's tie preference: smallest NORMAL right residue. */
        for (bestq=0;bestq<t->d;++bestq) {
            const size_t lr=submod(submod(bestpsi,t->pivot_C[row],t->d),bestq,t->d);
            const size_t rr=bestq?t->d-bestq:0;
            if (l[lr]!=(ls>0?-INF:INF) && r[rr]!=(rs>0?-INF:INF) &&
                ls*l[lr]+rs*r[rr]+t->pivot_A[row]==best) break;
        }
        if (bestq==t->d) goto failed;
        out[a].length=k; out[a].deletions=e; out[a].psi=bestpsi; out[a].omega=best;
        out[a].seed=malloc(k); out[a].positions=malloc(k*sizeof(size_t));
        if (!out[a].seed || !out[a].positions) goto failed;
        out[a].seed[p]=s[i]; out[a].positions[p]=i;
        (inward?trace_inward:trace)(w,s,&t->left,&w->left,i,p,h,
            submod(submod(bestpsi,t->pivot_C[row],t->d),bestq,t->d),ls>0,1,p,out+a);
        (inward?trace_inward:trace)(w,s,&t->right,&w->right,i,k-p-1,e-h,bestq,rs>0,0,p,out+a);
    }
    return 0;
failed:
    subseqhash3_deletion_results_destroy(out,count); return -1;
}
int subseqhash3_compute_deletion_seeds_with_workspace(const uint8_t *s,size_t n,
    size_t e,const SubseqHash3DeletionTables *t,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *out,SubseqHash3DeletionWorkspace *w) {
    return compute_deletion_seeds(s,n,e,t,pivots,count,out,w,0);
}
int subseqhash3_compute_inward_deletion_seeds_with_workspace(const uint8_t *s,size_t n,
    size_t e,const SubseqHash3DeletionTables *t,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *out,SubseqHash3DeletionWorkspace *w) {
    return compute_deletion_seeds(s,n,e,t,pivots,count,out,w,1);
}
int subseqhash3_compute_inward_deletion_seeds(const uint8_t *s,size_t n,size_t e,
    const SubseqHash3DeletionTables *t,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *out) {
    SubseqHash3DeletionWorkspace *w=subseqhash3_deletion_workspace_create();
    int status;
    if (!w) return -1;
    status=compute_deletion_seeds(s,n,e,t,pivots,count,out,w,1);
    subseqhash3_deletion_workspace_destroy(w); return status;
}
int subseqhash3_compute_deletion_seeds(const uint8_t *s,size_t n,size_t e,
    const SubseqHash3DeletionTables *t,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *out) {
    SubseqHash3DeletionWorkspace *w=subseqhash3_deletion_workspace_create();
    int status;
    if (!w) return -1;
    status=subseqhash3_compute_deletion_seeds_with_workspace(s,n,e,t,pivots,count,out,w);
    subseqhash3_deletion_workspace_destroy(w); return status;
}

/* One cache/workspace per independent set; step all sets before emitting a
   window so callbacks always see a complete set-major R*S batch. */
typedef struct {
    SubseqHash3DeletionWorkspace *w;
    DeletionWindowCache *cache;
    SubseqHash3DeletionWorkspace prepared_left,prepared_right;
    const SubseqHash3DeletionTables *tables;
    SubseqHash3DeletionResult *out;
    size_t slots,count,n,e;
    int shared,inward;
} DeletionWindowEngine;

static void window_engine_destroy(DeletionWindowEngine *v) {
    if (v->cache) for (size_t j=0;j<v->slots;++j) {
        free(v->cache[j].left.low); free(v->cache[j].right.low);
    }
    free(v->cache);
    free(v->prepared_left.packed_A); free(v->prepared_left.inward_masks);
    free(v->prepared_right.packed_A); free(v->prepared_right.inward_masks);
    subseqhash3_deletion_workspace_destroy(v->w);
}

static int window_engine_init(DeletionWindowEngine *v,const uint8_t *s,
    size_t length,size_t n,size_t e,const SubseqHash3DeletionTables *t,
    size_t count,SubseqHash3DeletionResult *out,SubseqHash3Traversal traversal) {
    v->tables=t; v->out=out; v->count=count; v->n=n; v->e=e;
    v->inward=traversal==SUBSEQHASH3_INWARD;
    v->w=subseqhash3_deletion_workspace_create();
    if (!v->w || compute_deletion_seeds(s,n,e,t,NULL,count,out,v->w,v->inward)) return -1;
    v->shared=!v->inward && length>n && e && n-e>1 &&
        !(e==1 && t->d>=128 && t->d>=8*n);
    if (v->shared) for (size_t j=0;j<4*(n-e-1)*t->d;++j)
        if (t->left.A[j]<-127 || t->left.A[j]>127 ||
            t->right.A[j]<-127 || t->right.A[j]>127) { v->shared=0; break; }
    if (v->shared) {
        v->prepared_left.d=v->prepared_right.d=t->d;
        if (t->d>1 && (prepare_inward(&v->prepared_left,&t->left,count-1,1) ||
            prepare_inward(&v->prepared_right,&t->right,n-e-1,0))) return -1;
        v->slots=count+e; v->cache=calloc(v->slots,sizeof(*v->cache));
        if (!v->cache) return -1;
        for (size_t j=0;j<v->slots;++j) v->cache[j].anchor=SIZE_MAX;
        v->w->window_cache=v->cache; v->w->window_slots=v->slots;
    }
    return 0;
}

static int window_engine_step(DeletionWindowEngine *v,const uint8_t *s,
    size_t length,size_t start) {
    SubseqHash3DeletionWorkspace *w=v->w;
    const SubseqHash3DeletionTables *t=v->tables;
    const size_t n=v->n,e=v->e;
    if (!v->shared) return compute_deletion_seeds(s+start,n,e,t,NULL,v->count,v->out,w,v->inward);
    w->window_start=start;
    for (size_t i=0;i<v->slots;++i) {
        const size_t anchor=start+i;
        DeletionWindowCache *c=v->cache+anchor%v->slots;
        if (c->anchor!=anchor) {
            uint8_t input[256]={0};
            for (size_t l=1;l<n && l<=anchor;++l) input[l-1]=s[anchor-l];
            if (build_inward(w,&c->left,input,&t->left,v->count-1,1,&v->prepared_left)) return -1;
            memset(input,0,n);
            for (size_t l=1;l<n && l<length-anchor;++l) input[n-l]=s[anchor+l];
            if (build_inward(w,&c->right,input,&t->right,n-e-1,0,&v->prepared_right)) return -1;
            c->anchor=anchor;
        }
    }
    return compute_compact_seeds(s+start,t,NULL,v->count,v->out,w,0,0,1);
}

/* Process all windows that use an anchor before replacing its score tables.
   Only candidate seeds survive in the ring; the two large DPs are reused. */
static int seed_ind_anchors(const uint8_t *s,size_t length,size_t n,size_t e,
    const SubseqHash3DeletionTables *tables,size_t count,size_t sets,
    SubseqHash3DeletionWindowCallback callback,void *context) {
    SubseqHash3DeletionWorkspace *w=subseqhash3_deletion_workspace_create();
    SubseqHash3DeletionWorkspace *prepared=NULL;
    SubseqHash3DeletionResult *ring=NULL;
    const size_t k=n-e,slots=count+e,total=count*sets,last=length-n;
    int status=-1;
    if (!w || sets>SIZE_MAX/(2*sizeof(*prepared)) || total>SIZE_MAX/slots ||
        total*slots>SIZE_MAX/sizeof(*ring)) goto done;
    prepared=calloc(2*sets,sizeof(*prepared));
    ring=calloc(slots*total,sizeof(*ring));
    if (!prepared || !ring) goto done;
    /* Check compact-score eligibility before emitting any callback. */
    for (size_t f=0;f<sets;++f) {
        SubseqHash3DeletionWorkspace *l=prepared+2*f,*r=l+1;
        l->d=r->d=tables[f].d;
        const int a=prepare_inward(l,&tables[f].left,count-1,1);
        const int b=prepare_inward(r,&tables[f].right,k-1,0);
        if (a<0 || b<0) goto done;
        if (a || b) { status=1; goto done; }
    }
    for (size_t anchor=0;anchor<=last+slots-1;++anchor) {
        for (size_t f=0;f<sets;++f) {
            const SubseqHash3DeletionTables *t=tables+f;
            const size_t d=t->d;
            uint8_t input[256]={0};
            if (reserve(w,n,k,e,d,-1)) goto done;
            for (size_t j=1;j<n && j<=anchor;++j) input[j-1]=s[anchor-j];
            if (build_inward(w,&w->inward_left,input,&t->left,count-1,1,prepared+2*f)) goto done;
            memset(input,0,n);
            for (size_t j=1;j<n && j<length-anchor;++j) input[n-j]=s[anchor+j];
            if (build_inward(w,&w->inward_right,input,&t->right,k-1,0,prepared+2*f+1)) goto done;
            for (size_t p=0;p<count;++p) for (size_t h=0;h<=e;++h) {
                const size_t i=p+h;
                if (anchor<i || anchor-i>last) continue;
                const size_t start=anchor-i,row=4*p+s[anchor];
                SubseqHash3DeletionResult *out=ring+(start%slots)*total+f*count+p;
                const uint16_t *l=(t->left_sign[row]>0?w->inward_left.high:w->inward_left.low)+
                    inward_state(&w->inward_left,i,h,e);
                const uint16_t *r=(t->right_sign[row]>0?w->inward_right.high:w->inward_right.low)+
                    inward_state(&w->inward_right,n-i-1,e-h,e);
                size_t psi=0,shift=t->pivot_C[row];
                unsigned cost=p==0?(unsigned)l[0]+r[shift]:
                    p==k-1?(unsigned)l[shift?d-shift:0]+r[0]:inward_join(l,r,d,shift);
                if (cost>=UINT16_MAX) {
                    if (out->seed && !out->psi) continue;
                    psi=inward_first(w,l,r,shift);
                    if (psi==d) continue;
                    shift=submod(shift,psi,d);
                    cost=p==0?(unsigned)l[0]+r[shift]:
                        p==k-1?(unsigned)l[shift?d-shift:0]+r[0]:inward_join(l,r,d,shift);
                    if (cost>=UINT16_MAX) continue;
                }
                const int64_t score=127*(int64_t)(k-1)-(int64_t)cost+t->pivot_A[row];
                if (out->seed && (psi>out->psi || (psi==out->psi && score<=out->omega))) continue;
                const size_t q=compact_witness(l,r,d,submod(psi,t->pivot_C[row],d),cost);
                if (q==d) goto done;
                if (!out->seed) {
                    out->seed=malloc(k); out->positions=malloc(k*sizeof(size_t));
                    if (!out->seed || !out->positions) goto done;
                }
                out->length=k; out->deletions=e; out->psi=psi; out->omega=score;
                out->seed[p]=s[anchor]; out->positions[p]=i;
                trace_window_side(s+start,e,d,&t->left,&w->inward_left,i,p,i,h,
                    submod(submod(psi,t->pivot_C[row],d),q,d),t->left_sign[row]>0,1,out);
                trace_window_side(s+start,e,d,&t->right,&w->inward_right,i,p,n-i-1,e-h,
                    q,t->right_sign[row]>0,0,out);
            }
        }
        if (anchor>=slots-1) {
            const size_t start=anchor-(slots-1);
            SubseqHash3DeletionResult *out=ring+(start%slots)*total;
            for (size_t a=0;a<total;++a) if (!out[a].seed) goto done;
            if (callback(start,out,total,context)) goto done;
            subseqhash3_deletion_results_destroy(out,total);
        }
    }
    status=0;
done:
    if (prepared) for (size_t f=0;f<2*sets;++f) {
        free(prepared[f].packed_A); free(prepared[f].inward_masks);
    }
    free(prepared);
    subseqhash3_deletion_results_destroy(ring,ring?slots*total:0); free(ring);
    subseqhash3_deletion_workspace_destroy(w); return status;
}

static int seed_ind_windows(const uint8_t *s,size_t length,size_t n,size_t e,
    const SubseqHash3DeletionTables *tables,size_t count,size_t sets,
    SubseqHash3Traversal traversal,SubseqHash3DeletionWindowCallback callback,void *context) {
    DeletionWindowEngine *engines=NULL;
    SubseqHash3DeletionResult *out=NULL;
    int status=-1;
    if (!s || !tables || !callback || !n || n>256 || length<n || e>=n ||
        !count || count>n-e || !sets || sets>SIZE_MAX/count ||
        sets>SIZE_MAX/sizeof(*engines) || count*sets>SIZE_MAX/sizeof(*out) ||
        (traversal!=SUBSEQHASH3_OUTWARD && traversal!=SUBSEQHASH3_INWARD)) return -1;
    for (size_t j=0;j<length;++j) if (s[j]>3) return -1;
    for (size_t f=0;f<sets;++f) {
        const SubseqHash3DeletionTables *t=tables+f;
        if (t->k<n-e || !t->d || !t->left.A || !t->right.A ||
            !t->left.B1 || !t->right.B1 || !t->left.C || !t->right.C ||
            !t->pivot_A || !t->pivot_C || !t->left_sign || !t->right_sign) return -1;
    }
    if (traversal==SUBSEQHASH3_OUTWARD && length-n>=count+e &&
        tables[0].d>=128 && e && n-e>1 &&
        !(e==1 && tables[0].d>=8*n)) {
        const int result=seed_ind_anchors(s,length,n,e,tables,count,sets,callback,context);
        if (result!=1) return result;
    }
    out=calloc(count*sets,sizeof(*out)); engines=calloc(sets,sizeof(*engines));
    if (!out || !engines) goto done;
    for (size_t f=0;f<sets;++f)
        if (window_engine_init(engines+f,s,length,n,e,tables+f,count,out+f*count,traversal)) goto done;
    for (size_t start=0;start<=length-n;++start) {
        if (start) for (size_t f=0;f<sets;++f)
            if (window_engine_step(engines+f,s,length,start)) goto done;
        if (callback(start,out,count*sets,context)) goto done;
        subseqhash3_deletion_results_destroy(out,count*sets);
    }
    status=0;
done:
    if (engines) for (size_t f=0;f<sets;++f) window_engine_destroy(engines+f);
    free(engines); subseqhash3_deletion_results_destroy(out,out?count*sets:0); free(out);
    return status;
}

int subseqhash3_seed_deletion_windows(const uint8_t *s,size_t length,size_t n,
    size_t e,const SubseqHash3DeletionTables *tables,size_t count,
    SubseqHash3Traversal traversal,SubseqHash3DeletionWindowCallback callback,void *context) {
    return seed_ind_windows(s,length,n,e,tables,count,1,traversal,callback,context);
}

int subseqhash3_seeds_ind(const uint8_t *s,size_t length,size_t n,size_t e,
    const SubseqHash3PivotTables *tables,size_t r,size_t sets,
    SubseqHash3Traversal traversal,SubseqHash3IndSeedsCallback callback,void *context) {
    SubseqHash3DeletionTables *converted=NULL;
    int status=-1;
    if (!s || !tables || !callback || !sets || sets>SIZE_MAX/sizeof(*converted) || !r ||
        sets>SIZE_MAX/r || r*sets>SIZE_MAX/sizeof(SubseqHash3DeletionResult) ||
        !n || n>256 || length<n || e>=n || r>n-e ||
        (traversal!=SUBSEQHASH3_OUTWARD && traversal!=SUBSEQHASH3_INWARD)) return -1;
    converted=calloc(sets,sizeof(*converted));
    if (!converted) return -1;
    for (size_t f=0;f<sets;++f) {
        if (tables[f].k<n-e || subseqhash3_deletion_tables_from_pivots(converted+f,tables+f)) goto done;
    }
    status=seed_ind_windows(s,length,n,e,converted,r,sets,traversal,callback,context);
done:
    for (size_t f=0;f<sets;++f) subseqhash3_deletion_tables_destroy(converted+f);
    free(converted); return status;
}
