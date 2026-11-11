/* General SubseqHash3 seeding: exactly Ed deletions and exactly Es
   substitutions of retained bases, seed length k=n-Ed. seeds-sub is the
   special case Ed=1 and seeds-ind the special case Es=0; all three share the
   same ABC^k orders (pivot tables converted with deletion_tables_from_pivots).

   Side DP state: [len][h][e][j] where len is the number of side input bases
   consumed, h=len-l the deletions among them, e the substitutions among the l
   retained bases and j the side residue. Each state keeps the exact int64
   minimum and maximum omega; tracebacks recompute predecessors from scores. */
#include "ssh3_seeding.h"
#include <stdlib.h>
#include <string.h>

#define GINF (INT64_MAX/4)

typedef struct {
    int64_t *mn,*mx;
    size_t cap,H,E1,d;
} GSide;

typedef struct {
    GSide left,right;
} GWork;

static size_t gidx(const GSide *s,size_t len,size_t h,size_t e) {
    return ((len*s->H+h)*s->E1+e)*s->d;
}

static int gside_reserve(GSide *s,size_t L,size_t H,size_t E1,size_t d) {
    const size_t cells=(L+1)*H*E1*d;
    s->H=H; s->E1=E1; s->d=d;
    if (cells<=s->cap) return 0;
    free(s->mn); free(s->mx);
    s->mn=malloc(cells*sizeof(int64_t)); s->mx=malloc(cells*sizeof(int64_t));
    if (!s->mn || !s->mx) { free(s->mn); free(s->mx); s->mn=s->mx=NULL; s->cap=0; return -1; }
    s->cap=cells; return 0;
}

/* Stream base i is x[i*step]. Builds lengths 0..L with retained length <= maxl. */
static void gside_build(GSide *s,const uint8_t *x,ptrdiff_t step,size_t L,
    size_t maxl,size_t Ed,size_t Es,const SubseqHash3DeletionSide *t) {
    const size_t d=s->d;
    for (size_t e=0;e<=Es;++e) {
        int64_t *mn=s->mn+gidx(s,0,0,e),*mx=s->mx+gidx(s,0,0,e);
        for (size_t j=0;j<d;++j) { mn[j]=GINF; mx[j]=-GINF; }
    }
    s->mn[gidx(s,0,0,0)]=0; s->mx[gidx(s,0,0,0)]=0;
    for (size_t len=1;len<=L;++len) {
        const unsigned base=x[(ptrdiff_t)(len-1)*step];
        const size_t hmax=len<Ed?len:Ed;
        for (size_t h=0;h<=hmax;++h) {
            const size_t l=len-h;
            for (size_t e=0;e<=Es;++e) {
                int64_t *mn=s->mn+gidx(s,len,h,e),*mx=s->mx+gidx(s,len,h,e);
                if (h) {
                    memcpy(mn,s->mn+gidx(s,len-1,h-1,e),d*sizeof(int64_t));
                    memcpy(mx,s->mx+gidx(s,len-1,h-1,e),d*sizeof(int64_t));
                } else {
                    for (size_t j=0;j<d;++j) { mn[j]=GINF; mx[j]=-GINF; }
                }
                if (!l || l>maxl || e>l) {
                    if (l>maxl || e>l) for (size_t j=0;j<d;++j) { mn[j]=GINF; mx[j]=-GINF; }
                    continue;
                }
                for (unsigned sigma=0;sigma<4;++sigma) {
                    const size_t delta=sigma!=base;
                    if (delta>e) continue;
                    const int64_t *pmn=s->mn+gidx(s,len-1,h,e-delta),*pmx=s->mx+gidx(s,len-1,h,e-delta);
                    const size_t row=4*(l-1)+sigma,c=t->C[row];
                    const int64_t *A=t->A+row*d; const int8_t *B=t->B1+row*d;
                    for (size_t j=0;j<d;++j) {
                        const size_t jp=j>=c?j-c:j+d-c;
                        if (pmn[jp]==GINF) continue;
                        int64_t a,b;
                        if (B[j]>0) { a=pmn[jp]+A[j]; b=pmx[jp]+A[j]; }
                        else { a=A[j]-pmx[jp]; b=A[j]-pmn[jp]; }
                        if (a<mn[j]) mn[j]=a;
                        if (b>mx[j]) mx[j]=b;
                    }
                }
            }
        }
    }
}

/* Recover one side string reaching value v (kind 0=min, 1=max). Side index i
   receives out_base[i], stream offset out_off[i] and substitution flag. */
static void gside_trace(const GSide *s,const uint8_t *x,ptrdiff_t step,size_t len,
    size_t h,size_t e,size_t j,int kind,int64_t v,const SubseqHash3DeletionSide *t,
    uint8_t *out_base,size_t *out_off,uint8_t *out_sub) {
    const size_t d=s->d;
    while (len) {
        const size_t l=len-h;
        if (h) {
            const int64_t *p=(kind?s->mx:s->mn)+gidx(s,len-1,h-1,e);
            if (p[j]==v) { --len; --h; continue; }
        }
        const unsigned base=x[(ptrdiff_t)(len-1)*step];
        int found=0;
        for (unsigned sigma=0;sigma<4 && !found;++sigma) {
            const size_t delta=sigma!=base;
            if (delta>e) continue;
            const size_t row=4*(l-1)+sigma,c=t->C[row];
            const size_t jp=j>=c?j-c:j+d-c;
            const int64_t A=t->A[row*d+j];
            const int pkind=t->B1[row*d+j]>0?kind:!kind;
            const int64_t pv=t->B1[row*d+j]>0?v-A:A-v;
            const int64_t *p=(pkind?s->mx:s->mn)+gidx(s,len-1,h,e-delta);
            if (p[jp]==GINF || p[jp]==-GINF || p[jp]!=pv) continue;
            out_base[l-1]=(uint8_t)sigma; out_off[l-1]=len-1; out_sub[l-1]=(uint8_t)delta;
            --len; e-=delta; j=jp; kind=pkind; v=pv; found=1;
        }
        if (!found) return; /* unreachable for a consistent DP */
    }
}

typedef struct {
    size_t psi; int64_t omega;
} GBest;

/* Join one input pivot. Sides are built from the pivot outward or from the
   window boundaries inward; lens/hs select states. Updates *res on improvement. */
typedef struct {
    const GSide *ls,*rs;
    const uint8_t *lx,*rx; ptrdiff_t lstep,rstep;
    size_t llen,lh,rlen,rh;
    size_t lpos0,rpos0; ptrdiff_t lposstep,rposstep; /* window position of stream offset */
    int left_reversed,right_reversed; /* side index -> seed index direction */
} GJoin;

static void gjoin(const GJoin *g,size_t p,size_t k,size_t pivot_pos,unsigned xb,
    size_t Es,const SubseqHash3DeletionTables *t,GBest *best,SubseqHash3EditResult *res,
    size_t Ed) {
    const size_t d=t->d;
    for (unsigned sigma=0;sigma<4;++sigma) {
        const size_t delta=sigma!=xb;
        if (delta>Es) continue;
        const size_t row=4*p+sigma,pc=t->pivot_C[row];
        const int ls=t->left_sign[row],rsg=t->right_sign[row];
        const int64_t pa=t->pivot_A[row];
        for (size_t el=0;el+delta<=Es;++el) {
            const size_t er=Es-delta-el;
            if (el>=g->ls->E1 || er>=g->rs->E1) continue;
            const int64_t *Lmn=g->ls->mn+gidx(g->ls,g->llen,g->lh,el),*Lmx=g->ls->mx+gidx(g->ls,g->llen,g->lh,el);
            const int64_t *Rmn=g->rs->mn+gidx(g->rs,g->rlen,g->rh,er),*Rmx=g->rs->mx+gidx(g->rs,g->rlen,g->rh,er);
            /* Smallest reachable combined residue, not exceeding the current best. */
            size_t psi=d;
            for (size_t target=0;target<=best->psi && target<d && psi==d;++target) {
                const size_t base=(target+d-pc)%d; /* j1+j2 == base mod d */
                for (size_t j1=0;j1<d;++j1) {
                    if (Lmn[j1]==GINF) continue;
                    const size_t j2=base>=j1?base-j1:base+d-j1;
                    if (Rmn[j2]!=GINF) { psi=target; break; }
                }
            }
            if (psi==d) continue;
            const size_t base=(psi+d-pc)%d;
            int64_t bo=0; size_t bj1=d;
            for (size_t j1=0;j1<d;++j1) {
                if (Lmn[j1]==GINF) continue;
                const size_t j2=base>=j1?base-j1:base+d-j1;
                if (Rmn[j2]==GINF) continue;
                const int64_t o=(ls>0?Lmx[j1]:-Lmn[j1])+(rsg>0?Rmx[j2]:-Rmn[j2])+pa;
                if (bj1==d || o>bo) { bo=o; bj1=j1; }
            }
            if (psi<best->psi || (psi==best->psi && bo>best->omega)) {
                best->psi=psi; best->omega=bo;
                const size_t j2=base>=bj1?base-bj1:base+d-bj1;
                uint8_t lb[256],rb[256],lsub[256],rsub[256]; size_t lo[256],ro[256];
                const size_t nl=p,nr=k-1-p;
                gside_trace(g->ls,g->lx,g->lstep,g->llen,g->lh,el,bj1,ls>0,
                    ls>0?Lmx[bj1]:Lmn[bj1],&t->left,lb,lo,lsub);
                gside_trace(g->rs,g->rx,g->rstep,g->rlen,g->rh,er,j2,rsg>0,
                    rsg>0?Rmx[j2]:Rmn[j2],&t->right,rb,ro,rsub);
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
                res->seed[p]=(uint8_t)sigma; res->substituted[p]=(uint8_t)delta;
                res->positions[p]=pivot_pos;
                res->psi=psi; res->omega=(double)bo;
                res->length=k; res->deletions=Ed; res->substitutions=Es;
            }
        }
    }
}

static void slot_reset(SubseqHash3EditResult *r,GBest *b,size_t count,size_t d) {
    for (size_t i=0;i<count;++i) { b[i].psi=d; b[i].omega=INT64_MIN; r[i].psi=d; }
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
    const size_t slots=outward?span+1:1;
    SubseqHash3DeletionTables *conv=calloc(sets,sizeof(*conv));
    SubseqHash3EditResult *res=calloc(slots*total,sizeof(*res));
    GBest *best=calloc(slots*total,sizeof(*best));
    GWork w={{0},{0}};
    int status=-1;
    if (!conv || !res || !best) goto done;
    for (size_t f=0;f<sets;++f)
        if (tables[f].k<k || tables[f].d!=d ||
            subseqhash3_deletion_tables_from_pivots(conv+f,tables+f)) goto done;
    for (size_t i=0;i<slots*total;++i) {
        res[i].seed=malloc(k); res[i].positions=malloc(k*sizeof(size_t)); res[i].substituted=malloc(k);
        if (!res[i].seed || !res[i].positions || !res[i].substituted) goto done;
    }
    if (gside_reserve(&w.left,n,H,E1,d) || gside_reserve(&w.right,n,H,E1,d)) goto done;
    const size_t windows=length-n+1;
    if (!outward) {
        for (size_t win=0;win<windows;++win) {
            const uint8_t *x=s+win;
            slot_reset(res,best,total,d);
            for (size_t f=0;f<sets;++f) {
                const SubseqHash3DeletionTables *t=conv+f;
                gside_build(&w.left,x,1,span,r-1,dels,subs,&t->left);
                gside_build(&w.right,x+n-1,-1,n-1,k-1,dels,subs,&t->right);
                for (size_t p=0;p<r;++p) for (size_t h=0;h<=dels;++h) {
                    const size_t b=p+h;
                    GJoin g={&w.left,&w.right,x,x+n-1,1,-1,b,h,n-1-b,dels-h,
                        0,n-1,1,-1,0,1};
                    gjoin(&g,p,k,b,x[b],subs,t,best+f*r+p,res+f*r+p,dels);
                }
            }
            if (callback(win,res,total,context)) goto done;
        }
        status=0; goto done;
    }
    for (size_t a=0;a<length;++a) {
        if (a<windows) slot_reset(res+(a%slots)*total,best+(a%slots)*total,total,d);
        const size_t wlo=a>span?a-span:0,whi=a<windows-1?a:windows-1;
        if (wlo<=whi) {
            const size_t bmax=a-wlo,bmin=a-whi;
            for (size_t f=0;f<sets;++f) {
                const SubseqHash3DeletionTables *t=conv+f;
                gside_build(&w.left,a?s+a-1:s,-1,bmax,r-1,dels,subs,&t->left);
                gside_build(&w.right,s+a+1,1,n-1-bmin,k-1,dels,subs,&t->right);
                for (size_t win=wlo;win<=whi;++win) {
                    const size_t b=a-win,slot=(win%slots)*total;
                    for (size_t h=0;h<=dels && h<=b;++h) {
                        const size_t p=b-h;
                        if (p>=r) continue;
                        GJoin g={&w.left,&w.right,a?s+a-1:s,s+a+1,-1,1,b,h,n-1-b,dels-h,
                            b-1,b+1,-1,1,1,0};
                        gjoin(&g,p,k,b,s[a],subs,t,best+slot+f*r+p,res+slot+f*r+p,dels);
                    }
                }
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
    if (res) subseqhash3_edit_results_destroy(res,slots*total);
    free(conv); free(res); free(best);
    free(w.left.mn); free(w.left.mx); free(w.right.mn); free(w.right.mx);
    return status;
}
