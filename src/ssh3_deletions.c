#define _POSIX_C_SOURCE 200112L
#include "ssh3_deletions.h"
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
#define INF (INT64_MAX/2)

static uint64_t random64(uint64_t *s) {
    uint64_t z=(*s+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
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
static void permutation4(unsigned *p,uint64_t *seed) {
    size_t a;
    for (a=0;a<4;++a) p[a]=(unsigned)a;
    for (a=4;a>1;--a) {
        const size_t q=(size_t)(random64(seed)%a);
        const unsigned tmp=p[a-1]; p[a-1]=p[q]; p[q]=tmp;
    }
}
static void residue_row(size_t *c,size_t d,uint64_t *seed) {
    size_t a,b;
    if (d<4) {
        unsigned p[4]; permutation4(p,seed);
        for (a=0;a<4;++a) c[a]=p[a]%d;
    } else for (a=0;a<4;++a) {
        do {
            c[a]=(size_t)(random64(seed)%d);
            for (b=0;b<a;++b) if (c[b]==c[a]) break;
        } while (b<a);
    }
}
int subseqhash3_deletion_tables_init(SubseqHash3DeletionTables *t,size_t k,
    size_t d,uint64_t seed) {
    size_t i,side,q,a;
    uint64_t limit;
    if (tables_allocate(t,k,d)) return -1;
    limit=(uint64_t)(INT64_MAX/4)/(uint64_t)k;
    if (limit>(UINT64_C(1)<<55)) limit=UINT64_C(1)<<55;
    for (side=0;side<2;++side) {
        SubseqHash3DeletionSide *s=side?&t->right:&t->left;
        for (i=0;i<k;++i) {
            for (q=0;q<d;++q) {
                unsigned p[4]; permutation4(p,&seed);
                for (a=0;a<4;++a) {
                    const size_t at=(4*i+a)*d+q;
                    const int64_t value=(int64_t)((UINT64_C(1)<<45)+random64(&seed)%(limit-(UINT64_C(1)<<45)+1));
                    s->A[at]=(p[a]&2)?value:-value;
                    s->B1[at]=(p[a]&1)?1:-1;
                }
            }
            residue_row(s->C+4*i,d,&seed);
        }
    }
    for (i=0;i<k;++i) {
        unsigned p[4],signs[4]; permutation4(p,&seed); permutation4(signs,&seed);
        residue_row(t->pivot_C+4*i,d,&seed);
        for (a=0;a<4;++a) {
            const int64_t value=(int64_t)((UINT64_C(1)<<45)+random64(&seed)%(limit-(UINT64_C(1)<<45)+1));
            t->pivot_A[4*i+a]=(signs[a]&1)?value:-value;
            t->left_sign[4*i+a]=(p[a]&2)?1:-1;
            t->right_sign[4*i+a]=(p[a]&1)?1:-1;
        }
    }
    return 0;
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
struct SubseqHash3DeletionWorkspace {
    size_t n,k,e,d,stride,score_capacity,parent_capacity,rolling_capacity;
    SideDP left,right;
    int64_t *rolling;
    size_t *previous,previous_capacity;
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
    free(w->rolling); free(w->previous); free(w);
}
static int reserve(SubseqHash3DeletionWorkspace *w,size_t n,size_t k,size_t e,size_t d) {
    size_t scores,parents,rolling,stride;
    if (d>SIZE_MAX-15) return -1;
    stride=((d+7)&~(size_t)7)+8;
    if (n*(e+1)>SIZE_MAX/stride) return -1;
    scores=n*(e+1)*stride;
    if (n*k*(e+1)>SIZE_MAX/d) return -1;
    parents=n*k*(e+1)*d;
    rolling=4*(e+1)*stride;
    if (scores>SIZE_MAX/(2*sizeof(int64_t)) || rolling>SIZE_MAX/sizeof(int64_t)) return -1;
    if (scores>w->score_capacity) {
        free(w->left.low); free(w->right.low);
        w->left.low=aligned_array(2*scores*sizeof(int64_t));
        w->right.low=aligned_array(2*scores*sizeof(int64_t));
        w->score_capacity=0;
        if (!w->left.low || !w->right.low) return -1;
        w->score_capacity=scores;
    }
    if (parents>w->parent_capacity) {
        free(w->left.parent); free(w->right.parent);
        w->left.parent=malloc(parents); w->right.parent=malloc(parents);
        w->parent_capacity=0;
        if (!w->left.parent || !w->right.parent) return -1;
        w->parent_capacity=parents;
    }
    if (rolling>w->rolling_capacity) {
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
    w->n=n; w->k=k; w->e=e; w->d=d; w->stride=stride;
    w->left.high=w->left.low+w->score_capacity;
    w->right.high=w->right.low+w->score_capacity;
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
    size_t maximum) {
    const size_t d=w->d,stride=w->stride,plane=(w->e+1)*stride;
    size_t i,j,h,r;
    for (i=0;i<w->n;++i) if (anchors[i]) {
        const size_t length=left?i:w->n-i-1;
        const size_t last=length<maximum?length:maximum;
        int64_t *previous=w->rolling,*previous_high=previous+plane;
        int64_t *current=previous_high+plane,*current_high=current+plane;
        for (h=0;h<=w->e && h<=length;++h) {
            int64_t *lo=previous+h*stride,*hi=previous_high+h*stride;
            for (r=0;r<d;++r) { lo[r]=INF; hi[r]=-INF; }
            lo[0]=hi[0]=0; padding(lo,d); padding(hi,d);
            if (h==length) {
                const size_t at=boundary(w,i,h);
                memcpy(out->low+at,lo,d*sizeof(int64_t));
                memcpy(out->high+at,hi,d*sizeof(int64_t));
            }
        }
        for (j=1;j<=last;++j) {
            const size_t maxh=w->e<length-j?w->e:length-j;
            for (h=0;h<=maxh;++h) {
                const size_t consumed=j+h;
                const unsigned a=s[left?i-consumed:i+consumed];
                const size_t row=4*(j-1)+a,shift=t->C[row];
                int64_t *lo=current+h*stride,*hi=current_high+h*stride;
                const int64_t *pl=previous+h*stride,*ph=previous_high+h*stride;
                uint8_t *parent=out->parent+parent_index(w,i,j,h);
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
                    _mm_mask_storeu_epi8(parent+r,(__mmask16)active,bits);
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
                    parent[r]=bits;
                }
#endif
                padding(lo,d); padding(hi,d);
                if (consumed==length) {
                    const size_t at=boundary(w,i,h);
                    for (r=0;r<d;++r) {
                        const size_t dst=left?r:(r?d-r:0);
                        out->low[at+dst]=lo[r]; out->high[at+dst]=hi[r];
                    }
                }
            }
            { int64_t *tmp=previous; previous=current; current=tmp;
              tmp=previous_high; previous_high=current_high; current_high=tmp; }
        }
        /* Pad only the valid full-boundary states used by joins. */
        for (h=0;h<=w->e && h<=length;++h) if (length-h<=last) {
            const size_t at=boundary(w,i,h);
            padding(out->low+at,d); padding(out->high+at,d);
        }
    }
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
        const size_t at=parent_index(w,i,j,h)+residue;
        if (dp->parent[at]&(high?2:1)) {
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
    SubseqHash3DeletionResult *out) {
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
                const size_t row=4*j+s[p-1-j];
                lr=addmod(lr,t->left.C[row],d);
                lv=t->left.B1[row*d+lr]*lv+t->left.A[row*d+lr];
            }
            for (j=0;j<k-p-1;++j) {
                const size_t row=4*j+s[p+1+j];
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
int subseqhash3_compute_deletion_seeds_with_workspace(const uint8_t *s,size_t n,
    size_t e,const SubseqHash3DeletionTables *t,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *out,SubseqHash3DeletionWorkspace *w) {
    size_t k,a,p,h,i,maxleft=0,maxright=0;
    uint8_t anchors[256]={0};
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
        for (h=0;h<=e;++h) anchors[p+h]=1;
    }
    if (!e || k==1) return direct(s,n,e,t,pivots,count,out);
    if (reserve(w,n,k,e,t->d)) return -1;
    build_side(w,s,&t->left,&w->left,1,anchors,maxleft);
    build_side(w,s,&t->right,&w->right,0,anchors,maxright);
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
        trace(w,s,&t->left,&w->left,i,p,h,
            submod(submod(bestpsi,t->pivot_C[row],t->d),bestq,t->d),ls>0,1,p,out+a);
        trace(w,s,&t->right,&w->right,i,k-p-1,e-h,bestq,rs>0,0,p,out+a);
    }
    return 0;
failed:
    subseqhash3_deletion_results_destroy(out,count); return -1;
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
