/* General SubseqHash3 seeding: exactly Ed deletions and exactly Es
   substitutions of retained bases, seed length k=n-Ed. seeds-sub is the
   special case Ed=1 and seeds-ind the special case Es=0; all three share the
   same ABC^k orders (pivot tables converted with deletion_tables_from_pivots).

   Side DP state [len][h][e][j]: len side input bases consumed, h=len-l
   deletions among them, e substitutions among the l retained bases, j residue.
   Scores are compact 16-bit costs (Supplementary Note 3): with M=127,
   Lc = M*l + min omega and Uc = M*l - max omega. Retention with signed
   increment a and sign +1 maps (Lc,Uc) -> (Lc+M+a, Uc+M-a); sign -1 maps
   (Lc,Uc) -> (Uc+M+a, Lc+M-a). All increments are in [0,254] and reachable
   costs are <= 2*M*255 < 65535, so saturating adds keep 0xFFFF unreachable.

   Each state row stores d residues followed by a 32-entry circular tail
   (row[d+i] = row[i mod d]) so every circular shift is one contiguous load.
   The right side is stored at reflected residue v=-j mod d, turning the pivot
   join into a circular vector sum: j1 + j2 + C_P == psi  <=>  v2 = j1+C_P-psi.
   Tracebacks recompute predecessors from exact costs. */
#define _POSIX_C_SOURCE 200112L
#include "ssh3_seeding.h"
#include <stdlib.h>
#include <string.h>
#if defined(__AVX512F__) && defined(__AVX512BW__) && defined(__AVX512VL__)
#include <immintrin.h>
#define G_SIMD 1
#else
#define G_SIMD 0
#endif

#ifdef GPROF_PHASES
#include <x86intrin.h>
#include <stdio.h>
static unsigned long long ph_build,ph_join,ph_trace;
#define PH_START unsigned long long ph0=__rdtsc();
#define PH_END(v) v+=__rdtsc()-ph0;
__attribute__((destructor)) static void ph_print(void){fprintf(stderr,"build %.3g join %.3g trace %.3g Gcyc\n",ph_build/1e9,ph_join/1e9,ph_trace/1e9);}
#else
#define PH_START
#define PH_END(v)
#endif
#define GM 127
#define GSENT 0xFFFFu

typedef struct {
    size_t d,dpad,chunks,rows;
    size_t *C;              /* residue increment of each row (reflected on right) */
    int8_t *A;              /* [row][dpad] signed increment a (B2 folded in) */
    uint32_t *neg;          /* [row][chunks] sign bits, bit set for B1=-1 */
} GTab;

/* State block: Lc row (stride entries) then Uc row. Entries d..d+31 hold the
   circular tail; entry d+32 caches the row minimum (join lower bound), valid
   when entry d+33 is nonzero (computed lazily by the join). */
typedef struct {
    uint16_t *buf;
    size_t cap,H,E1,stride,d,k,Es;
} GSide;

static void *galloc(size_t bytes) {
    void *p=NULL;
    if (!bytes) bytes=64;
    if (posix_memalign(&p,64,(bytes+63)/64*64)) return NULL;
    return p;
}

static void gtab_destroy(GTab *t) {
    free(t->C); free(t->A); free(t->neg);
    memset(t,0,sizeof(*t));
}

/* reflect=1 stores entry j at index (-j mod d) and negates C. */
static int gtab_init(GTab *t,const SubseqHash3DeletionSide *s,size_t k,size_t d,int reflect) {
    memset(t,0,sizeof(*t));
    t->d=d; t->chunks=(d+31)/32; t->dpad=t->chunks*32; t->rows=4*k;
    t->C=malloc(t->rows*sizeof(size_t)); t->A=galloc(t->rows*t->dpad);
    t->neg=calloc(t->rows*t->chunks,sizeof(uint32_t));
    if (!t->C || !t->A || !t->neg) { gtab_destroy(t); return -1; }
    memset(t->A,0,t->rows*t->dpad);
    for (size_t r=0;r<t->rows;++r) {
        const size_t c=s->C[r]%d;
        t->C[r]=reflect?(d-c)%d:c;
        for (size_t v=0;v<d;++v) {
            const size_t j=reflect?(d-v)%d:v;
            const int64_t a=s->A[r*d+j];
            if (a<-GM || a>GM || (s->B1[r*d+j]!=1 && s->B1[r*d+j]!=-1)) { gtab_destroy(t); return -1; }
            t->A[r*t->dpad+v]=(int8_t)a;
            if (s->B1[r*d+j]<0) t->neg[r*t->chunks+v/32]|=UINT32_C(1)<<(v%32);
        }
    }
    return 0;
}

static inline uint16_t *gblk(const GSide *s,size_t len,size_t h,size_t e) {
    return s->buf+((len*s->H+h)*s->E1+e)*2*s->stride;
}

static int gside_reserve(GSide *s,size_t L,size_t H,size_t E1,size_t d,size_t k,size_t Es) {
    const size_t stride=(d+34+31)/32*32,cells=(L+1)*H*E1*2*stride;
    s->H=H; s->E1=E1; s->d=d; s->stride=stride; s->k=k; s->Es=Es;
    if (cells<=s->cap) return 0;
    free(s->buf); s->buf=galloc(cells*sizeof(uint16_t));
    if (!s->buf) { s->cap=0; return -1; }
    s->cap=cells; return 0;
}

static inline uint16_t gadds(uint16_t a,uint16_t b) {
    const uint32_t x=(uint32_t)a+b; return x>GSENT?(uint16_t)GSENT:(uint16_t)x;
}

static inline void gfinish_row(uint16_t *row,size_t d) {
    uint16_t mn=(uint16_t)GSENT;
    for (size_t i=0;i<32;++i) row[d+i]=row[i%d];
    for (size_t j=0;j<d;++j) if (row[j]<mn) mn=row[j];
    row[d+32]=mn; row[d+33]=1;
}

/* Substitutions on a side with l retained bases: e in [emin(l), emax(l)]. */
static inline size_t gemin(const GSide *s,size_t l) { return s->Es>s->k-l?s->Es-(s->k-l):0; }
static inline size_t gemax(const GSide *s,size_t l) { return l<s->Es?l:s->Es; }

#if G_SIMD
typedef struct {
    __m512i plus[4],minus[4];
    __mmask32 ng[4];
    size_t q[4];
} GPrep;

static inline __attribute__((always_inline)) void gprep(GPrep *pr,const GTab *t,
    const size_t *row,const size_t *off,size_t c,size_t d,unsigned need) {
    const __m512i bias=_mm512_set1_epi16(GM);
    const size_t j0=c*32;
    for (unsigned i=0;i<4;++i) {
        if (!((need>>i)&1)) {
            pr->plus[i]=pr->minus[i]=_mm512_setzero_si512(); pr->ng[i]=0; pr->q[i]=0; continue;
        }
        const __m512i av=_mm512_cvtepi8_epi16(_mm256_load_si256((const __m256i *)(t->A+row[i]*t->dpad+j0)));
        pr->plus[i]=_mm512_add_epi16(bias,av); pr->minus[i]=_mm512_sub_epi16(bias,av);
        pr->ng[i]=(__mmask32)t->neg[row[i]*t->chunks+c];
        pr->q[i]=j0+off[i]; if (pr->q[i]>=d) pr->q[i]-=d;
    }
}

/* One residue chunk of all substitution counts of state plane (l,h). */
static inline __attribute__((always_inline)) void gchunk(const GSide *s,const GPrep *pr,
    size_t j0,size_t d,size_t l,size_t h,size_t emin,size_t emax,unsigned base,int subs,int keep) {
    const size_t stride=s->stride,len=l+h;
    const __m512i sent=_mm512_set1_epi16(-1);
    const __mmask32 m=j0+32<=d?(__mmask32)0xFFFFFFFFu:(__mmask32)((UINT32_C(1)<<(d-j0))-1);
    for (size_t e=emin;e<=emax;++e) {
        uint16_t *dst=gblk(s,len,h,e);
        __m512i aL=sent,aU=sent;
        if (h) { const uint16_t *del=gblk(s,len-1,h-1,e);
            aL=_mm512_loadu_si512(del+j0); aU=_mm512_loadu_si512(del+stride+j0); }
#define GUPD(I,SRC) do { \
        const uint16_t *src_=(SRC)+pr->q[I]; \
        const __m512i sL_=_mm512_loadu_si512(src_),sU_=_mm512_loadu_si512(src_+stride); \
        const __m512i selL_=_mm512_mask_blend_epi16(pr->ng[I],sL_,sU_); \
        const __m512i selU_=_mm512_ternarylogic_epi32(sL_,sU_,selL_,0x96); \
        aL=_mm512_min_epu16(aL,_mm512_adds_epu16(selL_,pr->plus[I])); \
        aU=_mm512_min_epu16(aU,_mm512_adds_epu16(selU_,pr->minus[I])); } while (0)
        if (subs && e) {
            const uint16_t *sp=gblk(s,len-1,h,e-1);
            switch (base) { /* constant indices keep the prepared vectors in registers */
            case 0: GUPD(1,sp); GUPD(2,sp); GUPD(3,sp); break;
            case 1: GUPD(0,sp); GUPD(2,sp); GUPD(3,sp); break;
            case 2: GUPD(0,sp); GUPD(1,sp); GUPD(3,sp); break;
            default: GUPD(0,sp); GUPD(1,sp); GUPD(2,sp); break;
            }
        }
        if (keep && e<l) {
            const uint16_t *kp=gblk(s,len-1,h,e);
            switch (base) {
            case 0: GUPD(0,kp); break; case 1: GUPD(1,kp); break;
            case 2: GUPD(2,kp); break; default: GUPD(3,kp); break;
            }
        }
#undef GUPD
        _mm512_mask_storeu_epi16(dst+j0,m,aL); _mm512_mask_storeu_epi16(dst+stride+j0,m,aU);
        if (!j0 && d>=32) { _mm512_storeu_si512(dst+d,aL); _mm512_storeu_si512(dst+stride+d,aU); }
    }
}


#endif

#define GPRUNE_D 128 /* join pruning by row minima pays off for large d */

/* Stream base i is x[i*step]. Builds lengths 0..L with retained length <= maxl
   and substitutions in [emin(l),emax(l)]; other states are never read.
   Iterates retained length l outermost: all states (len=l+h,h,e) of one l use
   the same four table rows, prepared once per residue chunk. State (l,h)
   depends on (l-1,h) and (l,h-1), both already built. */
static void gside_build(GSide *s,const uint8_t *x,ptrdiff_t step,size_t L,
    size_t maxl,size_t Ed,const GTab *t) {
    PH_START
    const size_t d=s->d,stride=s->stride;
    {
        uint16_t *b0=gblk(s,0,0,0);
        for (size_t j=0;j<d;++j) b0[j]=b0[stride+j]=(uint16_t)GSENT;
        b0[0]=b0[stride]=0; gfinish_row(b0,d); gfinish_row(b0+stride,d);
        for (size_t h=1;h<=Ed && h<=L;++h)
            memcpy(gblk(s,h,h,0),b0,2*stride*sizeof(uint16_t));
    }
    const size_t lmax=L<maxl?L:maxl;
    for (size_t l=1;l<=lmax;++l) {
        const size_t emin=gemin(s,l),emax=gemax(s,l);
        if (emin>emax) continue;
        const size_t hmax=L-l<Ed?L-l:Ed;
        size_t row[4],off[4];
        for (unsigned i=0;i<4;++i) { row[i]=4*(l-1)+i; off[i]=(d-t->C[row[i]])%d; }
#if G_SIMD
        const int subs=emax>0,keep=emin<l;
        for (size_t c=0;c<t->chunks;++c) { /* prepare once per chunk, reuse across h */
            GPrep pr; gprep(&pr,t,row,off,c,d,15u);
            for (size_t h=0;h<=hmax;++h)
                gchunk(s,&pr,c*32,d,l,h,emin,emax,x[(ptrdiff_t)(l+h-1)*step],subs,keep);
        }
        for (size_t h=0;h<=hmax;++h) for (size_t e=emin;e<=emax;++e) {
            uint16_t *dst=gblk(s,l+h,h,e);
            if (d<32) { gfinish_row(dst,d); gfinish_row(dst+stride,d); continue; }
            dst[d+33]=dst[stride+d+33]=0; /* minima not yet computed */
        }
#else
        for (size_t h=0;h<=hmax;++h) {
            const size_t len=l+h;
            const unsigned base=x[(ptrdiff_t)(len-1)*step];
            for (size_t e=emin;e<=emax;++e) {
                uint16_t *dst=gblk(s,len,h,e);
                const uint16_t *del=h?gblk(s,len-1,h-1,e):NULL;
                for (size_t j=0;j<d;++j) {
                    uint16_t aL=del?del[j]:(uint16_t)GSENT,aU=del?del[stride+j]:(uint16_t)GSENT;
                    for (unsigned i=0;i<4;++i) {
                        if (i!=base?!e:e>=l) continue;
                        const uint16_t *src=gblk(s,len-1,h,i!=base?e-1:e);
                        const size_t xj=j+off[i]<d?j+off[i]:j+off[i]-d;
                        const int av=t->A[row[i]*t->dpad+j];
                        const int neg=(t->neg[row[i]*t->chunks+j/32]>>(j%32))&1;
                        const uint16_t sL=src[xj],sU=src[stride+xj];
                        const uint16_t cL=gadds(neg?sU:sL,(uint16_t)(GM+av)),cU=gadds(neg?sL:sU,(uint16_t)(GM-av));
                        if (cL<aL) aL=cL;
                        if (cU<aU) aU=cU;
                    }
                    dst[j]=aL; dst[stride+j]=aU;
                }
                gfinish_row(dst,d); gfinish_row(dst+stride,d);
            }
        }
#endif
    }
    PH_END(ph_build)
}

/* Recover one side string with cost c (kindU: 0=Lc, 1=Uc) at residue j. */
static void gside_trace(const GSide *s,const uint8_t *x,ptrdiff_t step,size_t len,
    size_t h,size_t e,size_t j,int kindU,uint16_t c,const GTab *t,
    uint8_t *out_base,size_t *out_off,uint8_t *out_sub) {
    const size_t d=s->d,stride=s->stride;
    while (len) {
        const size_t l=len-h;
        if (h && gblk(s,len-1,h-1,e)[kindU*stride+j]==c) { --len; --h; continue; }
        const unsigned base=x[(ptrdiff_t)(len-1)*step];
        int found=0;
        for (unsigned q=0;q<4 && !found && l;++q) {
            const unsigned sigma=(base+q)&3;
            const size_t delta=sigma!=base;
            if (delta>e || e-delta>l-1) continue;
            const size_t r=4*(l-1)+sigma;
            const int a=t->A[r*t->dpad+j];
            const uint16_t inc=(uint16_t)(kindU?GM-a:GM+a);
            if (c<inc) continue;
            const int neg=(t->neg[r*t->chunks+j/32]>>(j%32))&1;
            const int pk=neg?!kindU:kindU;
            const size_t jp=(j+d-t->C[r])%d;
            if (gblk(s,len-1,h,e-delta)[pk*stride+jp]!=(uint16_t)(c-inc)) continue;
            out_base[l-1]=(uint8_t)sigma; out_off[l-1]=len-1; out_sub[l-1]=(uint8_t)delta;
            --len; e-=delta; j=jp; kindU=pk; c=(uint16_t)(c-inc); found=1;
        }
        if (!found) return; /* unreachable for a consistent DP */
    }
}

/* Row minimum over the d residues, cached in row[d+32] once computed. */
static uint16_t growmin(const uint16_t *crow,size_t d) {
    uint16_t *row=(uint16_t *)crow;
    if (row[d+33]) return row[d+32];
    uint16_t mn=(uint16_t)GSENT;
#if G_SIMD
    const __m512i sent=_mm512_set1_epi16(-1);
    __m512i v=sent;
    for (size_t j0=0;j0<d;j0+=32) {
        const __mmask32 m=j0+32<=d?(__mmask32)0xFFFFFFFFu:(__mmask32)((UINT32_C(1)<<(d-j0))-1);
        v=_mm512_min_epu16(v,_mm512_mask_loadu_epi16(sent,m,row+j0));
    }
    __m256i hh=_mm256_min_epu16(_mm512_castsi512_si256(v),_mm512_extracti64x4_epi64(v,1));
    __m128i qq=_mm_min_epu16(_mm256_castsi256_si128(hh),_mm256_extracti128_si256(hh,1));
    mn=(uint16_t)_mm_cvtsi128_si32(_mm_minpos_epu16(qq));
#else
    for (size_t j=0;j<d;++j) if (row[j]<mn) mn=row[j];
#endif
    row[d+32]=mn; row[d+33]=1;
    return mn;
}

/* min over j<d of a[j] + b[(j+s) mod d] (saturating); b has a circular tail. */
static uint16_t gminsum(const uint16_t *a,const uint16_t *b,size_t s,size_t d) {
#if G_SIMD
    __m512i acc=_mm512_set1_epi16(-1);
    const size_t chunks=(d+31)/32;
    for (size_t c=0;c<chunks;++c) {
        const size_t j0=c*32;
        const __mmask32 m=j0+32<=d?(__mmask32)0xFFFFFFFFu:(__mmask32)((UINT32_C(1)<<(d-j0))-1);
        size_t x=j0+s; if (x>=d) x-=d;
        const __m512i va=_mm512_mask_loadu_epi16(_mm512_set1_epi16(-1),m,a+j0);
        acc=_mm512_min_epu16(acc,_mm512_adds_epu16(va,_mm512_loadu_si512(b+x)));
    }
    __m256i h=_mm256_min_epu16(_mm512_castsi512_si256(acc),_mm512_extracti64x4_epi64(acc,1));
    __m128i q=_mm_min_epu16(_mm256_castsi256_si128(h),_mm256_extracti128_si256(h,1));
    return (uint16_t)_mm_cvtsi128_si32(_mm_minpos_epu16(q));
#else
    uint16_t best=(uint16_t)GSENT;
    for (size_t j=0;j<d;++j) { const uint16_t v=gadds(a[j],b[j+s<d?j+s:j+s-d]); if (v<best) best=v; }
    return best;
#endif
}

/* Smallest psi in 1..limit with a reachable pair, or d. Requires psi=0 failed.
   Pairs need b[(j1+pc-psi) mod d] reachable; prev[] gives the distance back to
   the nearest reachable b entry. scratch holds d entries. */
static size_t gpsi(const uint16_t *a,const uint16_t *b,size_t pc,size_t d,size_t limit,size_t *scratch) {
    size_t last=d,best=d;
    for (size_t v=0;v<d;++v) if (b[v]!=GSENT) last=v;
    if (last==d) return d;
    for (size_t v=0;v<d;++v) {
        if (b[v]!=GSENT) last=v;
        scratch[v]=v>=last?v-last:v+d-last;
    }
    for (size_t j=0;j<d;++j) {
        if (a[j]==GSENT) continue;
        size_t x=j+pc; if (x>=d) x-=d;
        if (scratch[x]<best) best=scratch[x];
    }
    return best<=limit?best:d;
}

typedef struct {
    size_t psi; int64_t omega;
    size_t h,el,s; unsigned sigma; uint16_t mn;  /* winning join, traced later */
} GBest;

typedef struct {
    const GSide *ls,*rs;
    const GTab *lt,*rt;
    const uint8_t *lx,*rx; ptrdiff_t lstep,rstep;
    size_t llen,lh,rlen,rh;
    size_t lpos0,rpos0; ptrdiff_t lposstep,rposstep; /* window position of stream offset */
    int left_reversed,right_reversed; /* side index -> seed index direction */
} GJoin;

/* Window-relative input pivot b=p+h; sides grow outward from the anchor. */
static GJoin gjoin_outward(const GSide *ls,const GSide *rs,const GTab *lt,const GTab *rt,
    const uint8_t *s,size_t a,size_t b,size_t h,size_t n,size_t Ed) {
    GJoin g={ls,rs,lt,rt,a?s+a-1:s,s+a+1,-1,1,b,h,n-1-b,Ed-h,b-1,b+1,-1,1,1,0};
    return g;
}
/* Sides grow inward from the window boundaries x[0] and x[n-1]. */
static GJoin gjoin_inward(const GSide *ls,const GSide *rs,const GTab *lt,const GTab *rt,
    const uint8_t *x,size_t b,size_t h,size_t n,size_t Ed) {
    GJoin g={ls,rs,lt,rt,x,x+n-1,1,-1,b,h,n-1-b,Ed-h,0,n-1,1,-1,0,1};
    return g;
}

/* Score all pivot bases and substitution splits at one input pivot; record
   an improvement in *best without traceback. */
static void gjoin(const GJoin *g,size_t p,size_t k,unsigned xb,size_t h,
    size_t Es,const SubseqHash3DeletionTables *t,GBest *best,size_t *scratch) {
    PH_START
    const size_t d=t->d,lstride=g->ls->stride,rstride=g->rs->stride,nl=p,nr=k-1-p;
    const int prune=d<32 || (Es && d>=GPRUNE_D);
    for (unsigned sigma=0;sigma<4;++sigma) {
        const size_t delta=sigma!=xb;
        if (delta>Es) continue;
        const size_t row=4*p+sigma,pc=t->pivot_C[row];
        const int lu=t->left_sign[row]>0,ru=t->right_sign[row]>0;
        const int64_t pa=t->pivot_A[row];
        for (size_t el=0;el+delta<=Es;++el) {
            const size_t er=Es-delta-el;
            if (el>nl || er>nr) continue;
            const uint16_t *a=gblk(g->ls,g->llen,g->lh,el)+lu*lstride;
            const uint16_t *b=gblk(g->rs,g->rlen,g->rh,er)+ru*rstride;
            if (prune && !best->psi) { /* omega upper bound from the row minima */
                const uint32_t lb=(uint32_t)growmin(a,d)+growmin(b,d);
                if (lb>=GSENT || (int64_t)GM*(int64_t)(k-1)-(int64_t)lb+pa<=best->omega) continue;
            }
            size_t psi=0,s=pc;
            uint16_t mn=gminsum(a,b,s,d);
            if (mn==GSENT) {
                if (!best->psi) continue;
                psi=gpsi(a,b,pc,d,best->psi,scratch);
                if (psi==d) continue;
                s=pc>=psi?pc-psi:pc+d-psi;
                mn=gminsum(a,b,s,d);
            }
            const int64_t omega=(int64_t)GM*(int64_t)(k-1)-(int64_t)mn+pa;
            if (psi<best->psi || (psi==best->psi && omega>best->omega)) {
                best->psi=psi; best->omega=omega;
                best->h=h; best->el=el; best->s=s; best->sigma=sigma; best->mn=mn;
            }
        }
    }
    PH_END(ph_join)
}

/* Trace the recorded winner; g must be the join of the winning anchor/h. */
static void gfinish(const GJoin *g,size_t p,size_t k,size_t pivot_pos,unsigned xb,
    size_t Es,size_t Ed,const SubseqHash3DeletionTables *t,const GBest *best,
    SubseqHash3EditResult *res) {
    PH_START
    const size_t d=t->d,nl=p,nr=k-1-p;
    const unsigned sigma=best->sigma;
    const size_t row=4*p+sigma,el=best->el,er=Es-(sigma!=xb)-el,s=best->s;
    const int lu=t->left_sign[row]>0,ru=t->right_sign[row]>0;
    const uint16_t *a=gblk(g->ls,g->llen,g->lh,el)+lu*g->ls->stride;
    const uint16_t *b=gblk(g->rs,g->rlen,g->rh,er)+ru*g->rs->stride;
    size_t j1=0,v2=0;
    for (;j1<d;++j1) {
        v2=j1+s<d?j1+s:j1+s-d;
        if (gadds(a[j1],b[v2])==best->mn) break;
    }
    uint8_t lb[256],rb[256],lsub[256],rsub[256]; size_t lo[256],ro[256];
    gside_trace(g->ls,g->lx,g->lstep,g->llen,g->lh,el,j1,lu,a[j1],g->lt,lb,lo,lsub);
    gside_trace(g->rs,g->rx,g->rstep,g->rlen,g->rh,er,v2,ru,b[v2],g->rt,rb,ro,rsub);
    for (size_t i=0;i<nl;++i) {
        const size_t z=g->left_reversed?p-1-i:i;
        res->seed[z]=lb[i]; res->substituted[z]=lsub[i];
        res->positions[z]=(size_t)((ptrdiff_t)g->lpos0+g->lposstep*(ptrdiff_t)lo[i]);
    }
    for (size_t i=0;i<nr;++i) {
        const size_t z=g->right_reversed?k-1-i:p+1+i;
        res->seed[z]=rb[i]; res->substituted[z]=rsub[i];
        res->positions[z]=(size_t)((ptrdiff_t)g->rpos0+g->rposstep*(ptrdiff_t)ro[i]);
    }
    res->seed[p]=(uint8_t)sigma; res->substituted[p]=(uint8_t)(sigma!=xb);
    res->positions[p]=pivot_pos;
    res->psi=best->psi; res->omega=(double)best->omega;
    res->length=k; res->deletions=Ed; res->substitutions=Es;
    PH_END(ph_trace)
}

static void slot_reset(GBest *b,size_t count,size_t d) {
    for (size_t i=0;i<count;++i) { b[i].psi=d; b[i].omega=INT64_MIN; }
}

int subseqhash3_seeds(const uint8_t *s,size_t length,size_t n,
    size_t subs,size_t dels,const SubseqHash3PivotTables *tables,
    size_t r,size_t sets,SubseqHash3Traversal traversal,
    SubseqHash3SubSeedsCallback callback,void *context) {
    if (!s || !tables || !callback || !n || n>257 || length<n || dels>=n) return -1;
    const size_t k=n-dels;
    if (k>256 || subs>k || !r || r>k || !sets || sets>SIZE_MAX/r ||
        (traversal!=SUBSEQHASH3_OUTWARD && traversal!=SUBSEQHASH3_INWARD)) return -1;
    for (size_t i=0;i<length;++i) if (s[i]>3) return -1;
    const size_t total=r*sets,d=tables[0].d,H=dels+1,E1=subs+1;
    const int outward=traversal==SUBSEQHASH3_OUTWARD;
    const size_t span=r-1+dels;                   /* last anchor offset of a window */
    const size_t slots=outward?span+1:1;          /* window result slots */
    const size_t ring=outward?dels+1:1;           /* anchors whose DPs are kept */
    SubseqHash3DeletionTables *conv=calloc(sets,sizeof(*conv));
    GTab *lt=calloc(sets,sizeof(*lt)),*rt=calloc(sets,sizeof(*rt));
    SubseqHash3EditResult *res=calloc(slots*total,sizeof(*res));
    GBest *best=calloc(slots*total,sizeof(*best));
    size_t *scratch=malloc((d?d:1)*sizeof(size_t));
    GSide *left=calloc(ring*sets,sizeof(GSide)),*right=calloc(ring*sets,sizeof(GSide));
    int status=-1;
    if (!conv || !lt || !rt || !res || !best || !scratch || !left || !right || !d) goto done;
    for (size_t f=0;f<sets;++f)
        if (tables[f].k<k || tables[f].d!=d ||
            subseqhash3_deletion_tables_from_pivots(conv+f,tables+f) ||
            gtab_init(lt+f,&conv[f].left,k,d,0) || gtab_init(rt+f,&conv[f].right,k,d,1)) goto done;
    for (size_t i=0;i<slots*total;++i) {
        res[i].seed=malloc(k); res[i].positions=malloc(k*sizeof(size_t)); res[i].substituted=malloc(k);
        if (!res[i].seed || !res[i].positions || !res[i].substituted) goto done;
    }
    for (size_t i=0;i<ring*sets;++i)
        if (gside_reserve(left+i,span,H,E1,d,k,subs) || gside_reserve(right+i,n-1,H,E1,d,k,subs)) goto done;
    const size_t windows=length-n+1;
    if (!outward) {
        for (size_t win=0;win<windows;++win) {
            const uint8_t *x=s+win;
            slot_reset(best,total,d);
            for (size_t f=0;f<sets;++f) {
                gside_build(left,x,1,span,r-1,dels,lt+f);
                gside_build(right,x+n-1,-1,n-1,k-1,dels,rt+f);
                for (size_t p=0;p<r;++p) {
                    GBest *bp=best+f*r+p;
                    for (size_t h=0;h<=dels;++h) {
                        const GJoin g=gjoin_inward(left,right,lt+f,rt+f,x,p+h,h,n,dels);
                        gjoin(&g,p,k,x[p+h],h,subs,conv+f,bp,scratch);
                    }
                    const GJoin g=gjoin_inward(left,right,lt+f,rt+f,x,p+bp->h,bp->h,n,dels);
                    gfinish(&g,p,k,p+bp->h,x[p+bp->h],subs,dels,conv+f,bp,res+f*r+p);
                }
            }
            if (callback(win,res,total,context)) goto done;
        }
        status=0; goto done;
    }
    for (size_t a=0;a<length;++a) {
        if (a<windows) slot_reset(best+(a%slots)*total,total,d);
        const size_t wlo=a>span?a-span:0,whi=a<windows-1?a:windows-1;
        if (wlo>whi) continue;
        const size_t bmax=a-wlo,bmin=a-whi,rs=a%ring;
        for (size_t f=0;f<sets;++f) {
            GSide *L=left+rs*sets+f,*R=right+rs*sets+f;
            gside_build(L,a?s+a-1:s,-1,bmax,r-1,dels,lt+f);
            gside_build(R,s+a+1,1,n-1-bmin,k-1,dels,rt+f);
            for (size_t win=wlo;win<=whi;++win) {
                const size_t b=a-win,slot=(win%slots)*total;
                for (size_t h=0;h<=dels && h<=b;++h) {
                    const size_t p=b-h;
                    if (p>=r) continue;
                    const GJoin g=gjoin_outward(L,R,lt+f,rt+f,s,a,b,h,n,dels);
                    gjoin(&g,p,k,s[a],h,subs,conv+f,best+slot+f*r+p,scratch);
                }
            }
            /* Seeds whose last anchor is a: window win, pivot p=a-dels-win. */
            if (a>=dels) for (size_t win=wlo;win<=whi && win<=a-dels;++win) {
                const size_t p=a-dels-win;
                if (p>=r) continue;
                const size_t slot=(win%slots)*total;
                const GBest *bp=best+slot+f*r+p;
                const size_t b=p+bp->h,aw=win+b,ws=aw%ring;
                const GJoin g=gjoin_outward(left+ws*sets+f,right+ws*sets+f,lt+f,rt+f,s,aw,b,bp->h,n,dels);
                gfinish(&g,p,k,b,s[aw],subs,dels,conv+f,bp,res+slot+f*r+p);
            }
        }
        if (a>=span && a-span<windows) {
            const size_t win=a-span;
            if (callback(win,res+(win%slots)*total,total,context)) goto done;
        }
    }
    status=0;
done:
    if (conv) for (size_t f=0;f<sets;++f) subseqhash3_deletion_tables_destroy(conv+f);
    if (lt) for (size_t f=0;f<sets;++f) gtab_destroy(lt+f);
    if (rt) for (size_t f=0;f<sets;++f) gtab_destroy(rt+f);
    if (res) subseqhash3_edit_results_destroy(res,slots*total);
    if (left) for (size_t i=0;i<ring*sets;++i) free(left[i].buf);
    if (right) for (size_t i=0;i<ring*sets;++i) free(right[i].buf);
    free(conv); free(lt); free(rt); free(res); free(best); free(scratch); free(left); free(right);
    return status;
}
