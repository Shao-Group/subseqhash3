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
static unsigned long long ph_lbuild,ph_build,ph_join,ph_trace,ph_gupd,ph_calls,ph_minsum,ph_psi,ph_pruned,ph_psicyc;
#define GCNT (++ph_gupd)
#define GC(x) (++(x))
#define PH_START unsigned long long ph0=__rdtsc();
#define PH_END(v) v+=__rdtsc()-ph0;
__attribute__((destructor)) static void ph_print(void){fprintf(stderr,"lbuild %.3g build %.3g join %.3g trace %.3g Gcyc gupd %.4g calls %.4g minsum %.4g psi %.4g pruned %.4g psicyc %.3g\n",ph_lbuild/1e9,ph_build/1e9,ph_join/1e9,ph_trace/1e9,(double)ph_gupd,(double)ph_calls,(double)ph_minsum,(double)ph_psi,(double)ph_pruned,ph_psicyc/1e9);}
#else
#define GCNT
#define GC(x)
#define PH_START
#define PH_END(v)
#endif
#define GM 127
#define GSENT 0xFFFFu
/* SIMD, d<=16: a block is one vector, Lc in lanes 0..15 and Uc in 16..31
   (row stride 16, entries 0..15 d-periodic). */
#define GPACK(d) (G_SIMD && (d)<=16)

typedef struct {
    size_t d,dpad,chunks,rows;
    size_t *C;              /* residue increment of each row (reflected on right) */
    size_t *off;            /* (d-C) mod d: source offset of each row */
    int8_t *A;              /* [row][dpad] signed increment a (B2 folded in) */
    uint32_t *neg;          /* [row][chunks] sign bits, bit set for B1=-1 */
    uint16_t *prep;         /* [l][chunk][4][32]: M+a of rows 4l..4l+3 */
    uint32_t *negs;         /* [l][chunk][4]: their sign masks */
    uint16_t *tt;           /* [row][d] traceback entry: (M+a) | sign bit << 8 */
    uint16_t *pidx;         /* d<32: [row][Lc,Uc][32] vpermt2w indices: shift + sign swap */
} GTab;

/* State block: Lc row (stride entries) then Uc row. Entries d..d+31 hold the
   circular tail. Row minima (join lower bounds) of the blocks the join reads
   are kept densely in mins[(len,h,e)][Lc/Uc]. */
typedef struct {
    uint16_t *buf,*nil;     /* nil: E1 sentinel blocks (0xFFFF) */
    uint16_t *mins;         /* [len][h][e][2] row minima */
    uint16_t *mm;           /* [len][h][e][2][64] minima per residue mod 32, stored twice */
    uint32_t *qt;           /* [chunk][4] wrapped source offsets of one level */
    size_t cap,H,E1,stride,d,k,Es;
    size_t ready;           /* start blocks (h,h,0), h<ready, are initialized */
    uint16_t per[32];       /* d<32: per[i] = (32+i) mod d (row is d-periodic) */
} GSide;

static void *galloc(size_t bytes) {
    void *p=NULL;
    if (!bytes) bytes=64;
    if (posix_memalign(&p,64,(bytes+63)/64*64)) return NULL;
    return p;
}

static void gtab_destroy(GTab *t) {
    free(t->C); free(t->off); free(t->tt); free(t->pidx); free(t->A); free(t->neg); free(t->prep); free(t->negs);
    memset(t,0,sizeof(*t));
}

/* SIMD kernel layout, interleaved per residue chunk: increments M+a (16-bit;
   the kernel derives M-a=2M-(M+a)) of the four rows of one retained length,
   and their sign masks. Built once per table set; shared by all side DPs. */
#if G_SIMD
static int gtab_prep(GTab *t) {
    const size_t k=t->rows/4,chunks=t->chunks;
    t->prep=galloc(k*chunks*128*sizeof(uint16_t)); t->negs=malloc(k*chunks*4*sizeof(uint32_t));
    if (!t->prep || !t->negs) return -1;
    if (t->d<=16) { /* packed: lane i<16 (Lc) reads source lane (i+off) mod d of
                       Lc, or of Uc (+16) where the target's sign is -1; lane
                       16+i (Uc) the other one. Then increments M+a | M-a. */
        t->pidx=galloc(t->rows*64*sizeof(uint16_t));
        if (!t->pidx) return -1;
        for (size_t r=0;r<t->rows;++r) for (size_t i=0;i<16;++i) {
            const size_t v=i%t->d,src=(t->off[r]+i)%t->d;
            const unsigned ng=(t->neg[r*chunks]>>v)&1u;
            const int a=t->A[r*t->dpad+v];
            t->pidx[r*64+i]=(uint16_t)(src+(ng?16:0));
            t->pidx[r*64+16+i]=(uint16_t)(src+(ng?0:16));
            t->pidx[r*64+32+i]=(uint16_t)(GM+a);
            t->pidx[r*64+48+i]=(uint16_t)(GM-a);
        }
    } else if (t->d<32) { /* lane i of a target reads source lane (i+off) mod d of Lc,
                      or of Uc (index +32) where the target's sign is -1 */
        t->pidx=galloc(t->rows*64*sizeof(uint16_t));
        if (!t->pidx) return -1;
        for (size_t r=0;r<t->rows;++r) for (size_t i=0;i<32;++i) {
            const size_t v=i%t->d,src=(t->off[r]+i)%t->d;
            const unsigned ng=(t->neg[r*chunks]>>v)&1u;
            t->pidx[r*64+i]=(uint16_t)(src+(ng?32:0));
            t->pidx[r*64+32+i]=(uint16_t)(src+(ng?0:32));
        }
    }
    for (size_t l=0;l<k;++l) for (size_t c=0;c<chunks;++c) for (unsigned i=0;i<4;++i) {
        const size_t r=4*l+i;
        uint16_t *pl=t->prep+((l*chunks+c)*4+i)*32;
        if (t->d<32) { /* d-periodic lanes: the kernel then yields d-periodic rows */
            uint32_t m=0;
            for (size_t j=0,v=0;j<32;++j,v=v+1==t->d?0:v+1) {
                pl[j]=(uint16_t)(GM+t->A[r*t->dpad+v]);
                m|=((t->neg[r*chunks]>>v)&1u)<<j;
            }
            t->negs[(l*chunks+c)*4+i]=m;
            continue;
        }
        for (size_t j=0;j<32;++j) pl[j]=(uint16_t)(GM+t->A[r*t->dpad+32*c+j]);
        t->negs[(l*chunks+c)*4+i]=t->neg[r*chunks+c];
    }
    return 0;
}
#endif

/* reflect=1 stores entry j at index (-j mod d) and negates C. */
static int gtab_init(GTab *t,const SubseqHash3DeletionSide *s,size_t k,size_t d,int reflect) {
    memset(t,0,sizeof(*t));
    t->d=d; t->chunks=(d+31)/32; t->dpad=t->chunks*32; t->rows=4*k;
    t->C=malloc(t->rows*sizeof(size_t)); t->off=malloc(t->rows*sizeof(size_t));
    t->tt=malloc(t->rows*d*sizeof(uint16_t)); t->A=galloc(t->rows*t->dpad);
    t->neg=calloc(t->rows*t->chunks,sizeof(uint32_t));
    if (!t->C || !t->off || !t->tt || !t->A || !t->neg) { gtab_destroy(t); return -1; }
    memset(t->A,0,t->rows*t->dpad);
    for (size_t r=0;r<t->rows;++r) {
        const size_t c=s->C[r]%d;
        t->C[r]=reflect?(d-c)%d:c;
        t->off[r]=(d-t->C[r])%d;
        for (size_t v=0;v<d;++v) {
            const size_t j=reflect?(d-v)%d:v;
            const int64_t a=s->A[r*d+j];
            if (a<-GM || a>GM || (s->B1[r*d+j]!=1 && s->B1[r*d+j]!=-1)) { gtab_destroy(t); return -1; }
            t->A[r*t->dpad+v]=(int8_t)a;
            t->tt[r*d+v]=(uint16_t)((GM+a)|(s->B1[r*d+j]<0?256:0));
            if (s->B1[r*d+j]<0) t->neg[r*t->chunks+v/32]|=UINT32_C(1)<<(v%32);
        }
    }
#if G_SIMD
    if (gtab_prep(t)) { gtab_destroy(t); return -1; }
#endif
    return 0;
}

static inline uint16_t *gblk(const GSide *s,size_t len,size_t h,size_t e) {
    return s->buf+((len*s->H+h)*s->E1+e)*2*s->stride;
}

static inline uint16_t *gmins(const GSide *s,size_t len,size_t h,size_t e) {
    return s->mins+((len*s->H+h)*s->E1+e)*2;
}

static inline uint16_t *gmm(const GSide *s,size_t len,size_t h,size_t e) {
    return s->mm+((len*s->H+h)*s->E1+e)*128;
}

#if !G_SIMD
/* Minima of row[0..d) per residue class mod 32 (duplicated to 64 entries so a
   rotation is one unaligned load); meaningful when 32 divides d. */
static void gmm_row(const uint16_t *row,size_t d,uint16_t *out) {
    for (size_t i=0;i<32;++i) {
        uint16_t m=(uint16_t)GSENT;
        for (size_t j=i;j<d;j+=32) if (row[j]<m) m=row[j];
        out[i]=out[i+32]=m;
    }
}
#endif

static int gside_reserve(GSide *s,size_t L,size_t H,size_t E1,size_t d,size_t k,size_t Es) {
    const size_t stride=GPACK(d)?16:(d+32+31)/32*32,cells=(L+1)*H*E1*2*stride;
    s->H=H; s->E1=E1; s->d=d; s->stride=stride; s->k=k; s->Es=Es; s->ready=0;
    for (size_t i=0;i<32;++i) s->per[i]=(uint16_t)((32+i)%d);
    if (cells<=s->cap) return 0;
    free(s->buf); free(s->nil); free(s->mins); free(s->mm);
    s->buf=galloc(cells*sizeof(uint16_t)); s->nil=galloc(E1*2*stride*sizeof(uint16_t));
    s->mins=galloc((L+1)*H*E1*2*sizeof(uint16_t));
    s->mm=galloc((L+1)*H*E1*128*sizeof(uint16_t));
    free(s->qt); s->qt=galloc(4*((d+31)/32)*sizeof(uint32_t));
    if (!s->buf || !s->nil || !s->mins || !s->mm || !s->qt) { s->cap=0; return -1; }
    memset(s->nil,0xFF,E1*2*stride*sizeof(uint16_t));
    s->cap=cells; return 0;
}

static inline uint16_t gadds(uint16_t a,uint16_t b) {
    const uint32_t x=(uint32_t)a+b; return x>GSENT?(uint16_t)GSENT:(uint16_t)x;
}

/* Fill the circular tail; return the row minimum. */
static inline uint16_t gfinish_row(uint16_t *row,size_t d) {
    uint16_t mn=(uint16_t)GSENT;
    for (size_t i=0,m=0;i<32;++i,m=m+1==d?0:m+1) row[d+i]=row[m];
    for (size_t j=0;j<d;++j) if (row[j]<mn) mn=row[j];
    return mn;
}

#if G_SIMD
static inline uint16_t ghmin(__m512i v) {
    __m256i h=_mm256_min_epu16(_mm512_castsi512_si256(v),_mm512_extracti64x4_epi64(v,1));
    __m128i q=_mm_min_epu16(_mm256_castsi256_si128(h),_mm256_extracti128_si256(h,1));
    return (uint16_t)_mm_cvtsi128_si32(_mm_minpos_epu16(q));
}
/* Row minima and per-residue-mod-32 minima (lane-wise over chunks) of a
   target's Lc/Uc rows. */
static inline void gstore_mm(uint16_t *mout,uint16_t *mmo,__m512i mL,__m512i mU) {
    mout[0]=ghmin(mL); mout[1]=ghmin(mU);
    _mm512_storeu_si512(mmo,mL); _mm512_storeu_si512(mmo+32,mL);
    _mm512_storeu_si512(mmo+64,mU); _mm512_storeu_si512(mmo+96,mU);
}
#endif

/* Substitutions on a side with l retained bases: e in [emin(l), emax(l)]. */
static inline size_t gemin(const GSide *s,size_t l) { return s->Es>s->k-l?s->Es-(s->k-l):0; }
static inline size_t gemax(const GSide *s,size_t l) { return l<s->Es?l:s->Es; }

#if G_SIMD
/* Transition with row I: source row pointer bI shifted by this chunk's
   wrapped offset qq[I]; swap Lc/Uc on sign -1 lanes (two non-destructive
   blends), add the increments M+a and M-a=2M-(M+a). */
#define GT(I) do { \
        const uint16_t *r_=b##I+qq[I]; \
        const __m512i sL_=_mm512_loadu_si512(r_),sU_=_mm512_loadu_si512(r_+stride); \
        const __mmask32 n_=(__mmask32)ng[I]; \
        const __m512i pl_=_mm512_load_si512(pp+(I)*32); \
        aL=_mm512_min_epu16(aL,_mm512_adds_epu16(_mm512_mask_blend_epi16(n_,sL_,sU_),pl_)); \
        aU=_mm512_min_epu16(aU,_mm512_adds_epu16(_mm512_mask_blend_epi16(n_,sU_,sL_),_mm512_sub_epi16(c254,pl_))); \
        GCNT; } while (0)

/* One target block (l,h,e) over all residue chunks. SUB: substitution
   transitions from block e-1 (rows i!=B); KEEP: retention from block e (row
   B). qt[4c+i] is row i's source offset (32c+off_i) mod d, shared by every
   target of the level, so the loop has no wrap arithmetic or branches. The
   last chunk is stored unmasked: lanes d..dpad-1 fall into the circular tail,
   rewritten afterwards. MINS: also record the row minima (mout) and the
   minima per residue mod 32 (mmo). Row data lives in L2 at large d, so the
   kernel favours few streams and few uops per chunk over register reuse. */
static inline __attribute__((always_inline)) void gtarget(uint16_t *dst,const uint16_t *del,
    const uint16_t *sp,const uint16_t *kp,const uint16_t *prep,const uint32_t *negs,
    const uint32_t *qt,size_t chunks,size_t d,size_t stride,
    uint16_t *mout,uint16_t *mmo,const unsigned B,const int SUB,const int KEEP,const int MINS) {
    const uint16_t *b0=B==0?kp:sp,*b1=B==1?kp:sp,*b2=B==2?kp:sp,*b3=B==3?kp:sp;
    const __mmask32 mlast=d%32?(__mmask32)((UINT32_C(1)<<(d%32))-1):(__mmask32)0xFFFFFFFFu;
    const __m512i c254=_mm512_set1_epi16(2*GM);
    __m512i mL=_mm512_set1_epi16(-1),mU=mL;
    for (size_t c=0;c<chunks;++c) {
        const size_t j0=32*c;
        const uint16_t *pp=prep+128*c;
        const uint32_t *ng=negs+4*c,*qq=qt+4*c;
        __m512i aL=_mm512_load_si512(del+j0),aU=_mm512_load_si512(del+stride+j0);
        if ((B==0 && KEEP) || (B!=0 && SUB)) GT(0);
        if ((B==1 && KEEP) || (B!=1 && SUB)) GT(1);
        if ((B==2 && KEEP) || (B!=2 && SUB)) GT(2);
        if ((B==3 && KEEP) || (B!=3 && SUB)) GT(3);
        _mm512_store_si512(dst+j0,aL); _mm512_store_si512(dst+stride+j0,aU);
        if (MINS) { /* lanes >= d of the last chunk are not residues */
            const __mmask32 mk=c+1<chunks?(__mmask32)0xFFFFFFFFu:mlast;
            mL=_mm512_mask_min_epu16(mL,mk,mL,aL); mU=_mm512_mask_min_epu16(mU,mk,mU,aU);
        }
    }
    if (MINS) { gstore_mm(mout,mmo,mL,mU); }
}
#undef GT

/* d<=16: the block is one vector [Lc|Uc]; one permute per row applies the
   shift and the sign swap, then the packed increments M+a | M-a. */
static inline __attribute__((always_inline)) void gtarget_packed(uint16_t *dst,const uint16_t *del,
    const uint16_t *sp,const uint16_t *kp,const uint16_t *px,size_t d,uint16_t *mout,uint16_t *mmo,
    const unsigned B,const int SUB,const int KEEP,const int MINS) {
    const __m512i k=_mm512_load_si512(kp),u=_mm512_load_si512(sp);
    __m512i acc=_mm512_load_si512(del);
#define GTP(I) do { const uint16_t *q_=px+(I)*64; \
        acc=_mm512_min_epu16(acc,_mm512_adds_epu16(_mm512_permutexvar_epi16(_mm512_load_si512(q_),B==(I)?k:u),_mm512_load_si512(q_+32))); \
        GCNT; } while (0)
    if ((B==0 && KEEP) || (B!=0 && SUB)) GTP(0);
    if ((B==1 && KEEP) || (B!=1 && SUB)) GTP(1);
    if ((B==2 && KEEP) || (B!=2 && SUB)) GTP(2);
    if ((B==3 && KEEP) || (B!=3 && SUB)) GTP(3);
#undef GTP
    _mm512_store_si512(dst,acc);
    if (MINS) {
        const __mmask32 mk=(__mmask32)((UINT32_C(1)<<d)-1);
        const __m512i sent=_mm512_set1_epi16(-1);
        const __m512i hi=_mm512_permutexvar_epi16(_mm512_add_epi16(_mm512_set_epi16(
            31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0),
            _mm512_set1_epi16(16)),acc);
        gstore_mm(mout,mmo,_mm512_mask_mov_epi16(sent,mk,acc),_mm512_mask_mov_epi16(sent,mk,hi));
    }
}

/* d<32: every row is one vector (entries 0..31 d-periodic). Sources are
   loaded aligned (store forwarding from the previous level always works) and
   one two-source permute per row applies both the shift and the sign swap. */
static inline __attribute__((always_inline)) void gtarget_small(uint16_t *dst,const uint16_t *del,
    const uint16_t *sp,const uint16_t *kp,const uint16_t *prep,const uint16_t *pidx,
    size_t d,size_t stride,uint16_t *mout,uint16_t *mmo,
    const unsigned B,const int SUB,const int KEEP,const int MINS) {
    const __m512i c254=_mm512_set1_epi16(2*GM);
    const __m512i kL=_mm512_load_si512(kp),kU=_mm512_load_si512(kp+stride);
    const __m512i uL=_mm512_load_si512(sp),uU=_mm512_load_si512(sp+stride);
    __m512i aL=_mm512_load_si512(del),aU=_mm512_load_si512(del+stride);
#define GTS(I) do { \
        const __m512i sL_=B==(I)?kL:uL,sU_=B==(I)?kU:uU; \
        const __m512i pl_=_mm512_load_si512(prep+(I)*32); \
        const __m512i tL_=_mm512_permutex2var_epi16(sL_,_mm512_load_si512(pidx+(I)*64),sU_); \
        const __m512i tU_=_mm512_permutex2var_epi16(sL_,_mm512_load_si512(pidx+(I)*64+32),sU_); \
        aL=_mm512_min_epu16(aL,_mm512_adds_epu16(tL_,pl_)); \
        aU=_mm512_min_epu16(aU,_mm512_adds_epu16(tU_,_mm512_sub_epi16(c254,pl_))); \
        GCNT; } while (0)
    if ((B==0 && KEEP) || (B!=0 && SUB)) GTS(0);
    if ((B==1 && KEEP) || (B!=1 && SUB)) GTS(1);
    if ((B==2 && KEEP) || (B!=2 && SUB)) GTS(2);
    if ((B==3 && KEEP) || (B!=3 && SUB)) GTS(3);
#undef GTS
    _mm512_store_si512(dst,aL); _mm512_store_si512(dst+stride,aU);
    if (MINS) {
        const __mmask32 mk=(__mmask32)((UINT32_C(1)<<d)-1);
        const __m512i sent=_mm512_set1_epi16(-1);
        gstore_mm(mout,mmo,_mm512_mask_mov_epi16(sent,mk,aL),_mm512_mask_mov_epi16(sent,mk,aU));
    }
}

static inline __attribute__((always_inline)) void gplane(const GSide *s,const GTab *t,
    const uint16_t *prep,const uint32_t *negs,const uint32_t *qt,size_t l,size_t h,size_t emin,size_t emax,int mins,const unsigned B) {
    const size_t d=s->d,stride=s->stride,len=l+h,chunks=t->chunks,BS=2*stride;
    uint16_t *dst0=gblk(s,len,h,0);
    const uint16_t *src0=gblk(s,len-1,h,0),*del0=h?gblk(s,len-1,h-1,0):s->nil;
    for (size_t e=emin;e<=emax;++e) {
        uint16_t *dst=dst0+e*BS,*mo=gmins(s,len,h,e);
        const uint16_t *del=del0+e*BS,*sp=src0+(e?e-1:0)*BS,*kp=src0+e*BS;
#define GTGT(SUB,KEEP) do { if (mins) gtarget(dst,del,sp,kp,prep,negs,qt,chunks,d,stride,mo,gmm(s,len,h,e),B,SUB,KEEP,1); \
            else gtarget(dst,del,sp,kp,prep,negs,qt,chunks,d,stride,mo,NULL,B,SUB,KEEP,0); } while (0)
        if (GPACK(d)) {
            const uint16_t *px=t->pidx+(l-1)*4*64;
#define GTPK(SUB,KEEP) do { if (mins) gtarget_packed(dst,del,sp,kp,px,d,mo,gmm(s,len,h,e),B,SUB,KEEP,1); \
                else gtarget_packed(dst,del,sp,kp,px,d,mo,NULL,B,SUB,KEEP,0); } while (0)
            if (!e) GTPK(0,1); else if (e>=l) GTPK(1,0); else GTPK(1,1);
#undef GTPK
            continue;
        }
        if (d<32) {
            const uint16_t *px=t->pidx+(l-1)*4*64;
#define GTSM(SUB,KEEP) do { if (mins) gtarget_small(dst,del,sp,kp,prep,px,d,stride,mo,gmm(s,len,h,e),B,SUB,KEEP,1); \
                else gtarget_small(dst,del,sp,kp,prep,px,d,stride,mo,NULL,B,SUB,KEEP,0); } while (0)
            if (!e) GTSM(0,1); else if (e>=l) GTSM(1,0); else GTSM(1,1);
#undef GTSM
            continue;
        }
        if (!e) GTGT(0,1); else if (e>=l) GTGT(1,0); else GTGT(1,1);
#undef GTGT
    }
    if (GPACK(d)) return;
    for (size_t e=emin;e<=emax;++e) {
        uint16_t *dst=dst0+e*BS;
        if (d<32) { /* rows of 64 entries, d-periodic: entries 0..31 already are
                       (periodic tables); only the join reads entries 32..63 */
            if (!mins) continue;
            const __m512i p1=_mm512_loadu_si512(s->per);
            _mm512_store_si512(dst+32,_mm512_permutexvar_epi16(p1,_mm512_load_si512(dst)));
            _mm512_store_si512(dst+stride+32,_mm512_permutexvar_epi16(p1,_mm512_load_si512(dst+stride)));
            continue; }
        _mm512_storeu_si512(dst+d,_mm512_load_si512(dst));
        _mm512_storeu_si512(dst+stride+d,_mm512_load_si512(dst+stride));
    }
}
#endif

#define GPRUNE_D 128 /* join pruning by row minima pays off for large d */

/* Stream base i is x[i*step]. Builds lengths 0..L with retained length <= maxl
   and substitutions in [emin(l),emax(l)]; other states are never read. Row
   minima are recorded for retained lengths >= minl (the join's final states).
   Iterates retained length l outermost: all states (len=l+h,h,e) of one l use
   the same four table rows and the same wrapped offsets qt. State (l,h)
   depends on (l-1,h) and (l,h-1), both already built. */
static void gside_build(GSide *s,const uint8_t *x,ptrdiff_t step,size_t L,
    size_t maxl,size_t minl,size_t Ed,const GTab *t) {
    PH_START
    const size_t d=s->d,stride=s->stride;
    if (s->ready<=Ed) { /* start blocks (h,h,0) are never overwritten: set once */
        uint16_t *b0=gblk(s,0,0,0);
        for (size_t j=0;j<d;++j) b0[j]=b0[stride+j]=(uint16_t)GSENT;
        b0[0]=b0[stride]=0;
        if (GPACK(d)) { for (size_t j=d;j<16;++j) b0[j]=b0[16+j]=b0[j%d]; }
        else { gfinish_row(b0,d); gfinish_row(b0+stride,d); }
        for (size_t h=0;h<=Ed && h<s->H;++h) {
            if (h) memcpy(gblk(s,h,h,0),b0,2*stride*sizeof(uint16_t));
            gmins(s,h,h,0)[0]=gmins(s,h,h,0)[1]=0;
            uint16_t *m=gmm(s,h,h,0); /* row = [0, sentinel...]: class 0 has 0 */
            for (size_t i=0;i<128;++i) m[i]=(uint16_t)GSENT;
            m[0]=m[32]=m[64]=m[96]=0;
        }
        s->ready=Ed+1;
    }
    const size_t lmax=L<maxl?L:maxl;
    for (size_t l=1;l<=lmax;++l) {
        const size_t emin=gemin(s,l),emax=gemax(s,l);
        if (emin>emax) continue;
        const size_t hmax=L-l<Ed?L-l:Ed;
        size_t row[4],off[4];
        for (unsigned i=0;i<4;++i) { row[i]=4*(l-1)+i; off[i]=t->off[row[i]]; }
#if G_SIMD
        const uint16_t *prep=t->prep+(l-1)*t->chunks*128;
        const uint32_t *negs=t->negs+(l-1)*t->chunks*4;
        for (size_t c=0;c<t->chunks;++c) for (unsigned i=0;i<4;++i) {
            const size_t q=32*c+off[i]; s->qt[4*c+i]=(uint32_t)(q>=d?q-d:q);
        }
        for (size_t h=0;h<=hmax;++h) {
            const unsigned base=x[(ptrdiff_t)(l+h-1)*step];
            switch (base) {
            case 0: gplane(s,t,prep,negs,s->qt,l,h,emin,emax,l>=minl,0); break;
            case 1: gplane(s,t,prep,negs,s->qt,l,h,emin,emax,l>=minl,1); break;
            case 2: gplane(s,t,prep,negs,s->qt,l,h,emin,emax,l>=minl,2); break;
            default: gplane(s,t,prep,negs,s->qt,l,h,emin,emax,l>=minl,3); break;
            }
        }
#else
        (void)minl;
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
                gmins(s,len,h,e)[0]=gfinish_row(dst,d); gmins(s,len,h,e)[1]=gfinish_row(dst+stride,d);
                if (l>=minl && d%32==0) { gmm_row(dst,d,gmm(s,len,h,e)); gmm_row(dst+stride,d,gmm(s,len,h,e)+64); }
            }
        }
#endif
    }
    if (maxl<8) { PH_END(ph_lbuild) } else { PH_END(ph_build) }
}

/* Recover one side string with cost c (kindU: 0=Lc, 1=Uc) at residue j.
   Tries the deletion, then retention of the input base, then the three
   substitutions; the common case (plain retention) is predicted well. The
   block offset of (len-1,h,e) is walked incrementally. */
static void gside_trace(const GSide *s,const uint8_t *x,ptrdiff_t step,size_t len,
    size_t h,size_t e,size_t j,int kindU,uint16_t c,const GTab *t,
    uint8_t *out_base,size_t *out_off,uint8_t *out_sub) {
    const size_t d=s->d,stride=s->stride,BS=2*stride,SH=s->E1*BS,SL=s->H*SH;
    const uint16_t *buf=s->buf,*tt=t->tt; const size_t *off=t->off;
    unsigned ku=(unsigned)kindU;
    size_t po=(len-1)*SL+h*SH+e*BS;   /* gblk(len-1,h,e) */
    while (len) {
        const size_t l=len-h;
        if (h && buf[po-SH+ku*stride+j]==c) { --len; --h; po-=SL+SH; continue; }
        if (!l) return; /* unreachable for a consistent DP */
        const unsigned base=x[(ptrdiff_t)(len-1)*step];
        const size_t r0=4*(l-1);
        unsigned q=0;
        for (;q<4;++q) {
            const size_t delta=q!=0;
            if (delta>e || e-delta>l-1) continue;
            const size_t r=r0+((base+q)&3);
            const unsigned w=tt[r*d+j],u=w&255u;
            const uint16_t inc=(uint16_t)(ku?2*GM-u:u);
            if (c<inc) continue;
            const unsigned pk=ku^(w>>8);
            size_t jp=j+off[r]; jp=jp>=d?jp-d:jp;
            const size_t pq=po-delta*BS;
            if (buf[pq+pk*stride+jp]!=(uint16_t)(c-inc)) continue;
            out_base[l-1]=(uint8_t)((base+q)&3); out_off[l-1]=len-1; out_sub[l-1]=(uint8_t)delta;
            --len; e-=delta; j=jp; ku=pk; c=(uint16_t)(c-inc); po=pq-SL;
            break;
        }
        if (q==4) return; /* unreachable for a consistent DP */
    }
}

/* min over j<d of a[j] + b[(j+s) mod d] (saturating); b has a circular tail. */
#if G_SIMD
/* Entries (j+s) mod d, j<32, of a packed row (d<=16). */
static inline __m512i grot16(const uint16_t *b,size_t s,size_t d) {
    const __m512i v=_mm512_add_epi16(_mm512_set_epi16(
        31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0),
        _mm512_set1_epi16((short)s));
    const __mmask32 w=_mm512_cmpge_epu16_mask(v,_mm512_set1_epi16((short)d));
    return _mm512_permutexvar_epi16(_mm512_mask_sub_epi16(v,w,v,_mm512_set1_epi16((short)d)),
        _mm512_maskz_loadu_epi16(0xFFFFu,b));
}
#endif

static uint16_t gminsum(const uint16_t *a,const uint16_t *b,size_t s,size_t d) {
#if G_SIMD
    if (GPACK(d)) {
        const __mmask32 m=(__mmask32)((UINT32_C(1)<<d)-1);
        return ghmin(_mm512_mask_adds_epu16(_mm512_set1_epi16(-1),m,_mm512_maskz_loadu_epi16(m,a),grot16(b,s,d)));
    }
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

/* Bitset of reachable entries (!= sentinel) of row[0..d). */
static void gbits(const uint16_t *row,size_t d,uint64_t *w) {
    const size_t W=(d+63)/64;
    memset(w,0,W*sizeof(uint64_t));
#if G_SIMD
    const __m512i sent=_mm512_set1_epi16(-1);
    for (size_t j0=0;j0<d;j0+=32) {
        const __mmask32 m=j0+32<=d?(__mmask32)0xFFFFFFFFu:(__mmask32)((UINT32_C(1)<<(d-j0))-1);
        w[j0/64]|=(uint64_t)_mm512_mask_cmpneq_epu16_mask(m,_mm512_maskz_loadu_epi16(m,row+j0),sent)<<(j0%64);
    }
#else
    for (size_t j=0;j<d;++j) if (row[j]!=GSENT) w[j/64]|=UINT64_C(1)<<(j%64);
#endif
}

/* Smallest delta<=maxd with bit (t-delta) mod d set, or SIZE_MAX. */
static size_t gpred(const uint64_t *w,size_t d,size_t t,size_t maxd) {
    size_t dist=0,pos=t;
    while (dist<=maxd) {
        const size_t wi=pos>>6,bi=pos&63;
        const uint64_t m=w[wi]&(bi==63?~UINT64_C(0):(UINT64_C(2)<<bi)-1);
        if (m) { const size_t v=wi*64+63-(size_t)__builtin_clzll(m),r=dist+(pos-v); return r<=maxd?r:SIZE_MAX; }
        dist+=bi+1; pos=wi?wi*64-1:d-1;
    }
    return SIZE_MAX;
}

/* Smallest delta<=maxd with bit (u+delta) mod d set, or SIZE_MAX. */
static size_t gsucc(const uint64_t *w,size_t d,size_t u,size_t maxd) {
    size_t dist=0,pos=u;
    while (dist<=maxd) {
        const size_t wi=pos>>6,bi=pos&63;
        const uint64_t m=w[wi]&(~UINT64_C(0)<<bi);
        if (m) { const size_t v=wi*64+(size_t)__builtin_ctzll(m),r=dist+(v-pos); return r<=maxd?r:SIZE_MAX; }
        const size_t end=wi*64+64<d?wi*64+64:d;
        dist+=end-pos; pos=end==d?0:end;
    }
    return SIZE_MAX;
}

/* Reachability bitset of a block (shared by its Lc and Uc rows), computed at
   most once per join call; returns its popcount. */
static size_t gbits_cached(const uint16_t *row,size_t d,uint64_t *w,size_t *cnt) {
    if (*cnt!=SIZE_MAX) return *cnt;
    gbits(row,d,w);
    size_t n=0;
    for (size_t i=0;i<(d+63)/64;++i) n+=(size_t)__builtin_popcountll(w[i]);
    return *cnt=n;
}

/* Smallest psi in 1..limit with a reachable pair, or d. Requires psi=0 failed.
   A pair (j1,v2) gives psi=(j1+pc-v2) mod d. Iterates the sparser of the two
   reachable bitsets and searches the other for the nearest partner within the
   current bound. scratch holds 2*ceil(d/64) words. */
static size_t gpsi(const uint64_t *wa,size_t na,const uint64_t *wb,size_t nb,size_t pc,size_t d,size_t limit) {
    const size_t W=(d+63)/64;
    if (!na || !nb) return d;
    size_t best=limit<d?limit:d-1; int found=0;
    if (na<=nb) {
        for (size_t i=0;i<W && best;++i) for (uint64_t m=wa[i];m;m&=m-1) {
            const size_t j1=i*64+(size_t)__builtin_ctzll(m);
            size_t t=j1+pc; if (t>=d) t-=d;
            const size_t r=gpred(wb,d,t,best);
            if (r!=SIZE_MAX) { best=r; found=1; if (best<=1) return best; }
        }
    } else {
        for (size_t i=0;i<W && best;++i) for (uint64_t m=wb[i];m;m&=m-1) {
            const size_t v2=i*64+(size_t)__builtin_ctzll(m);
            const size_t u=v2>=pc?v2-pc:v2+d-pc;
            const size_t r=gsucc(wa,d,u,best);
            if (r!=SIZE_MAX) { best=r; found=1; if (best<=1) return best; }
        }
    }
    return found?best:d;
}

typedef struct {
    size_t psi; int64_t omega;
    size_t h,el,s; unsigned sigma; uint16_t mn;  /* winning join (for its traceback) */
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

typedef struct { int64_t ub; size_t idx; } GCand;

/* Lower bound of min_j a[j]+b[(j+s) mod d] for 32 | d: residue classes mod 32
   pair as (i, i+s mod 32); ma, mb hold the class minima (mb duplicated). */
static uint16_t gmmbound(const uint16_t *ma,const uint16_t *mb,size_t s32) {
#if G_SIMD
    return ghmin(_mm512_adds_epu16(_mm512_loadu_si512(ma),_mm512_loadu_si512(mb+s32)));
#else
    uint16_t m=(uint16_t)GSENT;
    for (size_t i=0;i<32;++i) { const uint16_t v=gadds(ma[i],mb[i+s32]); if (v<m) m=v; }
    return m;
#endif
}

/* Score all pivot bases and substitution splits at one input pivot; record
   an improvement in *best without traceback. Once some window/pivot already
   has psi=0, candidates are visited in decreasing order of their omega upper
   bound (row minima) so the bound prunes most of them; ties keep the first
   candidate in (sigma,el) order, exactly as the plain sequential scan.
   Returns nonzero when *best changed. */
static int gjoin(const GJoin *g,size_t p,size_t k,unsigned xb,size_t h,
    size_t Es,const SubseqHash3DeletionTables *t,GBest *best,uint64_t *scratch,GCand *cand) {
    PH_START
    const size_t d=t->d,lstride=g->ls->stride,rstride=g->rs->stride,nl=p,nr=k-1-p,E1=Es+1;
    const int prune=d<32 || (Es && d>=GPRUNE_D);
    const int64_t base=(int64_t)GM*(int64_t)(k-1);
    if (prune && !best->psi) {
        size_t n=0;
        for (unsigned sigma=0;sigma<4;++sigma) {
            const size_t delta=sigma!=xb;
            if (delta>Es) continue;
            const size_t row=4*p+sigma;
            const int lu=t->left_sign[row]>0,ru=t->right_sign[row]>0;
            for (size_t el=0;el+delta<=Es;++el) {
                const size_t er=Es-delta-el;
                if (el>nl || er>nr) continue;
                uint32_t lb=(uint32_t)gmins(g->ls,g->llen,g->lh,el)[lu]+gmins(g->rs,g->rlen,g->rh,er)[ru];
                if (lb<GSENT && d%32==0) /* pair minima of residues mod 32 at shift pc */
                    lb=gmmbound(gmm(g->ls,g->llen,g->lh,el)+64*lu,gmm(g->rs,g->rlen,g->rh,er)+64*ru,t->pivot_C[row]%32);
                const int64_t ub=base-(int64_t)lb+t->pivot_A[row];
                if (lb>=GSENT || ub<=best->omega) { GC(ph_pruned); continue; }
                size_t m=n++; /* insertion sort: ub descending, idx ascending */
                while (m && cand[m-1].ub<ub) { cand[m]=cand[m-1]; --m; }
                cand[m].ub=ub; cand[m].idx=sigma*E1+el;
            }
        }
        int here=0; size_t bidx=0; (void)scratch;
        for (size_t c=0;c<n;++c) {
            const int64_t ub=cand[c].ub;
            const size_t idx=cand[c].idx;
            if (ub<best->omega) break;
            if (ub==best->omega && !(here && idx<bidx)) { GC(ph_pruned); continue; }
            const unsigned sigma=(unsigned)(idx/E1);
            const size_t el=idx%E1,er=Es-(sigma!=xb)-el,row=4*p+sigma,pc=t->pivot_C[row];
            const int lu=t->left_sign[row]>0,ru=t->right_sign[row]>0;
            const uint16_t *a=gblk(g->ls,g->llen,g->lh,el)+lu*lstride;
            const uint16_t *b=gblk(g->rs,g->rlen,g->rh,er)+ru*rstride;
            GC(ph_minsum); const uint16_t mn=gminsum(a,b,pc,d);
            if (mn==GSENT) continue; /* psi>0 cannot improve */
            const int64_t omega=base-(int64_t)mn+t->pivot_A[row];
            if (omega>best->omega || (omega==best->omega && here && idx<bidx)) {
                best->omega=omega; best->h=h; best->el=el; best->s=pc; best->sigma=sigma; best->mn=mn;
                here=1; bidx=idx;
            }
        }
        PH_END(ph_join)
        return here;
    }
    int changed=0;
    /* scratch: left then right block bitsets for e=0..Es; counts unknown yet */
    const size_t W=(d+63)/64;
    size_t lcnt[257],rcnt[257];
    if (best->psi) for (size_t e=0;e<=Es;++e) lcnt[e]=rcnt[e]=SIZE_MAX;
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
                const uint32_t lb=(uint32_t)gmins(g->ls,g->llen,g->lh,el)[lu]+gmins(g->rs,g->rlen,g->rh,er)[ru];
                if (lb>=GSENT || base-(int64_t)lb+pa<=best->omega) { GC(ph_pruned); continue; }
            }
            size_t psi=0,s=pc;
            GC(ph_minsum); uint16_t mn=gminsum(a,b,s,d);
            if (mn==GSENT) { GC(ph_psi);
                if (!best->psi) continue;
#ifdef GPROF_PHASES
                unsigned long long t0_=__rdtsc();
#endif
                uint64_t *wl=scratch+el*W,*wr=scratch+(E1+er)*W;
                const size_t na=gbits_cached(a,d,wl,lcnt+el),nb=gbits_cached(b,d,wr,rcnt+er);
                psi=gpsi(wl,na,wr,nb,pc,d,best->psi);
#ifdef GPROF_PHASES
                ph_psicyc+=__rdtsc()-t0_;
#endif
                if (psi==d) continue;
                s=pc>=psi?pc-psi:pc+d-psi;
                mn=gminsum(a,b,s,d);
            }
            const int64_t omega=base-(int64_t)mn+pa;
            if (psi<best->psi || (psi==best->psi && omega>best->omega)) {
                best->psi=psi; best->omega=omega;
                best->h=h; best->el=el; best->s=s; best->sigma=sigma; best->mn=mn;
                changed=1;
            }
        }
    }
    PH_END(ph_join)
    return changed;
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
#if G_SIMD
    if (GPACK(d)) {
        const __mmask32 m=(__mmask32)((UINT32_C(1)<<d)-1);
        const __m512i sum=_mm512_adds_epu16(_mm512_maskz_loadu_epi16(m,a),grot16(b,s,d));
        const __mmask32 eq=_mm512_mask_cmpeq_epu16_mask(m,sum,_mm512_set1_epi16((short)best->mn));
        j1=(size_t)__builtin_ctz((unsigned)eq); v2=j1+s<d?j1+s:j1+s-d;
    } else { /* first residue pair realizing the recorded minimum */
        const __m512i sent=_mm512_set1_epi16(-1),want=_mm512_set1_epi16((short)best->mn);
        for (size_t j0=0;j0<d;j0+=32) {
            const __mmask32 m=j0+32<=d?(__mmask32)0xFFFFFFFFu:(__mmask32)((UINT32_C(1)<<(d-j0))-1);
            size_t x=j0+s; if (x>=d) x-=d;
            const __m512i sum=_mm512_adds_epu16(_mm512_mask_loadu_epi16(sent,m,a+j0),_mm512_loadu_si512(b+x));
            const __mmask32 eq=_mm512_mask_cmpeq_epu16_mask(m,sum,want);
            if (eq) { j1=j0+(size_t)__builtin_ctz((unsigned)eq); break; }
        }
        v2=j1+s<d?j1+s:j1+s-d;
    }
#else
    for (;j1<d;++j1) {
        v2=j1+s<d?j1+s:j1+s-d;
        if (gadds(a[j1],b[v2])==best->mn) break;
    }
#endif
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
    /* Outward: a window's winner is known after its last anchor. Small DPs keep
       the last dels+1 anchors (ring) and trace each winner once; large ones
       trace every improvement immediately while the anchor's single DP is
       hot, avoiding the ring's cache traffic. */
    const size_t stride_=(d+32+31)/32*32,blk=2*stride_*sizeof(uint16_t);
    const size_t dpbytes=((n-1+1)+(span+1))*H*E1*blk;
    const int eager=outward && (dels+1)*dpbytes>(size_t)768*1024;
    const size_t ring=outward && !eager?dels+1:1;
    /* Row minima only serve join pruning (gjoin's prune condition). */
    const int need_mins=d<32 || (subs && d>=GPRUNE_D);
    SubseqHash3DeletionTables *conv=calloc(sets,sizeof(*conv));
    GTab *lt=calloc(sets,sizeof(*lt)),*rt=calloc(sets,sizeof(*rt));
    SubseqHash3EditResult *res=calloc(slots*total,sizeof(*res));
    GBest *best=calloc(slots*total,sizeof(*best));
    uint64_t *scratch=malloc((2*(subs+1)*((d+63)/64)+8)*sizeof(uint64_t));
    GCand *cand=malloc(4*(subs+1)*sizeof(GCand));
    GSide *left=calloc(ring*sets,sizeof(GSide)),*right=calloc(ring*sets,sizeof(GSide));
    int status=-1;
    if (!conv || !lt || !rt || !res || !best || !scratch || !cand || !left || !right || !d) goto done;
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
                gside_build(left,x,1,span,r-1,need_mins?0:SIZE_MAX,dels,lt+f);
                gside_build(right,x+n-1,-1,n-1,k-1,need_mins?k-r:SIZE_MAX,dels,rt+f);
                for (size_t p=0;p<r;++p) {
                    GBest *bp=best+f*r+p;
                    for (size_t h=0;h<=dels;++h) {
                        const GJoin g=gjoin_inward(left,right,lt+f,rt+f,x,p+h,h,n,dels);
                        gjoin(&g,p,k,x[p+h],h,subs,conv+f,bp,scratch,cand);
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
            gside_build(L,a?s+a-1:s,-1,bmax,r-1,need_mins?0:SIZE_MAX,dels,lt+f);
            gside_build(R,s+a+1,1,n-1-bmin,k-1,need_mins?k-r:SIZE_MAX,dels,rt+f);
            for (size_t win=wlo;win<=whi;++win) {
                const size_t b=a-win,slot=(win%slots)*total;
                for (size_t h=0;h<=dels && h<=b;++h) {
                    const size_t p=b-h;
                    if (p>=r) continue;
                    const GJoin g=gjoin_outward(L,R,lt+f,rt+f,s,a,b,h,n,dels);
                    GBest *bp=best+slot+f*r+p;
                    /* Trace improvements while this anchor's DPs are live:
                       only one DP per side is kept, and it stays in cache. */
                    if (gjoin(&g,p,k,s[a],h,subs,conv+f,bp,scratch,cand) && eager)
                        gfinish(&g,p,k,b,s[a],subs,dels,conv+f,bp,res+slot+f*r+p);
                }
            }
            /* Ring mode: seeds whose last anchor is a (window win, pivot
               p=a-dels-win), traced from the winning anchor's kept DPs. */
            if (!eager && a>=dels) for (size_t win=wlo;win<=whi && win<=a-dels;++win) {
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
    if (left) for (size_t i=0;i<ring*sets;++i) { free(left[i].buf); free(left[i].nil); free(left[i].mins); free(left[i].mm); free(left[i].qt); }
    if (right) for (size_t i=0;i<ring*sets;++i) { free(right[i].buf); free(right[i].nil); free(right[i].mins); free(right[i].mm); free(right[i].qt); }
    free(conv); free(lt); free(rt); free(res); free(best); free(scratch); free(cand); free(left); free(right);
    return status;
}
