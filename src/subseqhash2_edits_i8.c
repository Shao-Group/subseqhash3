#include "subseqhash2_edits_i8.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define SH2_POS_INF INT16_C(32639)
#define SH2_NEG_INF INT16_C(-32640)

typedef struct {
    size_t N,n,E,d,band;
    size_t *offset;
    int16_t *fmn,*fmx,*rmn,*rmx;
    uint8_t *fp,*rp;
} SharedDP;

static uint64_t next_random(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}

static size_t tri(size_t E,size_t h,size_t e) {
    return h*(E+1)-h*(h-1)/2+e;
}

static size_t table_index(size_t d,size_t row,size_t r,unsigned a) {
    return (row*4+(size_t)a)*d+r;
}

static size_t block_state(const SharedDP *p,size_t u,size_t b,
                          size_t h,size_t e,size_t r) {
    return p->offset[b]+u*p->band+tri(p->E,h,e)*p->d+r;
}

static size_t sub_mod(size_t r,size_t shift,size_t d) {
    return r>=shift?r-shift:r+d-shift;
}

static void save_parent(uint8_t *slot,int which,unsigned symbol,
                        int previous_type) {
    unsigned code=1u+(unsigned)symbol*2u+(unsigned)(previous_type-1);
    if (which==1) *slot=(uint8_t)((*slot&0xf0u)|code);
    else *slot=(uint8_t)((*slot&0x0fu)|(code<<4));
}

static int allocate_tables(SubseqHash2I8Tables *t,size_t k,size_t d) {
    size_t side,pivot;
    if (!d || k>SIZE_MAX/d || k*d>SIZE_MAX/4) return -1;
    side=k*d*4; pivot=k*4;
    if (pivot>SIZE_MAX/sizeof(*t->CF)) return -1;
    t->AF=malloc(side?side:1); t->BF1=malloc(side?side:1);
    t->BF2=malloc(side?side:1); t->AR=malloc(side?side:1);
    t->BR1=malloc(side?side:1); t->BR2=malloc(side?side:1);
    t->AP=malloc(pivot?pivot:1); t->BP1=malloc(pivot?pivot:1);
    t->BP2=malloc(pivot?pivot:1);
    t->CF=malloc((pivot?pivot:1)*sizeof(*t->CF));
    t->CR=malloc((pivot?pivot:1)*sizeof(*t->CR));
    t->CP=malloc((pivot?pivot:1)*sizeof(*t->CP));
    return t->AF&&t->BF1&&t->BF2&&t->AR&&t->BR1&&t->BR2&&
           t->AP&&t->BP1&&t->BP2&&t->CF&&t->CR&&t->CP?0:-1;
}

int subseqhash2_i8_tables_init(SubseqHash2I8Tables *t,size_t k,size_t d,
                               uint64_t seed) {
    static const int8_t b1[4]={1,1,-1,-1},b2[4]={1,-1,1,-1};
    size_t row,r,a,side,pivot;
    if (!t || k>(size_t)((SH2_POS_INF-1)/127)) return -1;
    memset(t,0,sizeof(*t)); t->k=k; t->d=d;
    if (allocate_tables(t,k,d)) { subseqhash2_i8_tables_destroy(t); return -1; }
    side=k*d*4; pivot=k*4;
    for (row=0;row<side;++row) {
        t->AF[row]=(int8_t)(1u+next_random(&seed)%127u);
        t->AR[row]=(int8_t)(1u+next_random(&seed)%127u);
    }
    for (row=0;row<k;++row) for (r=0;r<d;++r) {
        unsigned pf[4]={0,1,2,3},pr[4]={0,1,2,3};
        for (a=4;a>1;--a) {
            size_t q=(size_t)(next_random(&seed)%a); unsigned x=pf[a-1];
            pf[a-1]=pf[q]; pf[q]=x;
            q=(size_t)(next_random(&seed)%a); x=pr[a-1];
            pr[a-1]=pr[q]; pr[q]=x;
        }
        for (a=0;a<4;++a) {
            size_t x=table_index(d,row,r,(unsigned)a);
            t->BF1[x]=b1[pf[a]]; t->BF2[x]=b2[pf[a]];
            t->BR1[x]=b1[pr[a]]; t->BR2[x]=b2[pr[a]];
        }
    }
    for (row=0;row<pivot;++row) {
        t->AP[row]=(int8_t)(1u+next_random(&seed)%127u);
        t->BP1[row]=(next_random(&seed)&1u)?1:-1;
        t->BP2[row]=(next_random(&seed)&1u)?1:-1;
        t->CF[row]=(size_t)(next_random(&seed)%d);
        t->CR[row]=(size_t)(next_random(&seed)%d);
        t->CP[row]=(size_t)(next_random(&seed)%d);
    }
    return 0;
}

void subseqhash2_i8_tables_destroy(SubseqHash2I8Tables *t) {
    if (!t) return;
    free(t->CP); free(t->CR); free(t->CF);
    free(t->BP2); free(t->BP1); free(t->AP);
    free(t->BR2); free(t->BR1); free(t->AR);
    free(t->BF2); free(t->BF1); free(t->AF);
    memset(t,0,sizeof(*t));
}

size_t subseqhash2_i8_results_per_window(size_t n,size_t E) {
    size_t count;
    if (E>n || E==SIZE_MAX) return SIZE_MAX;
    if (n==E) count=0;
    else if (E+1>SIZE_MAX/(n-E)) return SIZE_MAX;
    else count=(E+1)*(n-E);
    if (E && E+1>SIZE_MAX/E) return SIZE_MAX;
    if (count>SIZE_MAX-E*(E+1)/2) return SIZE_MAX;
    return count+E*(E+1)/2;
}

int subseqhash2_i8_result_coordinates(size_t n,size_t E,size_t index,
                                      size_t *subs,size_t *order) {
    size_t s;
    if (!subs || !order || E>n) return -1;
    for (s=0;s<=E;++s) {
        size_t k=n-E+s;
        if (index<k) { *subs=s; *order=index; return 0; }
        index-=k;
    }
    return -1;
}

static int allocate_dp(SharedDP *p,size_t N,size_t n,size_t E,size_t d) {
    size_t b,total=0,tri_count;
    if (E==SIZE_MAX || E+1>SIZE_MAX/(E+2)) return -1;
    tri_count=(E+1)*(E+2)/2;
    if (tri_count>SIZE_MAX/d) return -1;
    p->N=N; p->n=n; p->E=E; p->d=d; p->band=tri_count*d;
    p->offset=malloc((n+2)*sizeof(*p->offset));
    if (!p->offset) return -1;
    for (b=0;b<=n;++b) {
        size_t starts=N-b+1;
        p->offset[b]=total;
        if (starts>SIZE_MAX/p->band || total>SIZE_MAX-starts*p->band)
            return -1;
        total+=starts*p->band;
    }
    p->offset[n+1]=total;
    if (total>SIZE_MAX/sizeof(*p->fmn)) return -1;
    p->fmn=malloc(total*sizeof(*p->fmn)); p->fmx=malloc(total*sizeof(*p->fmx));
    p->rmn=malloc(total*sizeof(*p->rmn)); p->rmx=malloc(total*sizeof(*p->rmx));
    p->fp=calloc(total?total:1,1); p->rp=calloc(total?total:1,1);
    if (!p->fmn||!p->fmx||!p->rmn||!p->rmx||!p->fp||!p->rp) return -1;
    for (b=0;b<total;++b) {
        p->fmn[b]=p->rmn[b]=SH2_POS_INF;
        p->fmx[b]=p->rmx[b]=SH2_NEG_INF;
    }
    return 0;
}

static void free_dp(SharedDP *p) {
    free(p->rp); free(p->fp); free(p->rmx); free(p->rmn);
    free(p->fmx); free(p->fmn); free(p->offset); memset(p,0,sizeof(*p));
}

static void update_interval(int16_t pmn,int16_t pmx,int8_t A,int8_t b1,
                            int8_t b2,int16_t *mn,int16_t *mx,
                            uint8_t *parent,unsigned symbol) {
    int16_t add,lo,hi; int lt,ut;
    if (pmn==SH2_POS_INF) return;
    add=(int16_t)A*(int16_t)b2;
    if (b1>0) { lo=(int16_t)(pmn+add); hi=(int16_t)(pmx+add); lt=1; ut=2; }
    else { lo=(int16_t)(-pmx+add); hi=(int16_t)(-pmn+add); lt=2; ut=1; }
    if (lo<*mn) { *mn=lo; save_parent(parent,1,symbol,lt); }
    if (hi>*mx) { *mx=hi; save_parent(parent,2,symbol,ut); }
}

static void fill_dp(const uint8_t *x,const SubseqHash2I8Tables *t,SharedDP *p) {
    size_t b,u,h,e,r,a;
    for (u=0;u<=p->N;++u) {
        size_t z=block_state(p,u,0,0,0,0);
        p->fmn[z]=p->fmx[z]=p->rmn[z]=p->rmx[z]=0;
    }
    for (b=1;b<=p->n;++b) for (u=0;u+b<=p->N;++u) {
        unsigned fs=x[u+b-1],rs=x[u];
        for (h=0;h<=p->E && h<=b;++h) {
            size_t l=b-h,maxe=p->E-h;
            if (l>t->k) continue;
            if (maxe>l) maxe=l;
            for (e=0;e<=maxe;++e) for (r=0;r<p->d;++r) {
                size_t dst=block_state(p,u,b,h,e,r);
                if (h) {
                    size_t q=block_state(p,u,b-1,h-1,e,r);
                    p->fmn[dst]=p->fmn[q]; p->fmx[dst]=p->fmx[q];
                    q=block_state(p,u+1,b-1,h-1,e,r);
                    p->rmn[dst]=p->rmn[q]; p->rmx[dst]=p->rmx[q];
                }
                if (!l) continue;
                for (a=0;a<4;++a) {
                    size_t df=(a!=fs),dr=(a!=rs),q,base;
                    if (df<=e) {
                        q=sub_mod(r,t->CF[4*(l-1)+a],p->d);
                        base=table_index(p->d,l-1,r,(unsigned)a);
                        q=block_state(p,u,b-1,h,e-df,q);
                        update_interval(p->fmn[q],p->fmx[q],t->AF[base],
                            t->BF1[base],t->BF2[base],p->fmn+dst,p->fmx+dst,
                            p->fp+dst,(unsigned)a);
                    }
                    if (dr<=e) {
                        q=sub_mod(r,t->CR[4*(l-1)+a],p->d);
                        base=table_index(p->d,l-1,r,(unsigned)a);
                        q=block_state(p,u+1,b-1,h,e-dr,q);
                        update_interval(p->rmn[q],p->rmx[q],t->AR[base],
                            t->BR1[base],t->BR2[base],p->rmn+dst,p->rmx+dst,
                            p->rp+dst,(unsigned)a);
                    }
                }
            }
        }
    }
}

static int traceback_forward(const uint8_t *x,const SubseqHash2I8Tables *t,
    const SharedDP *p,size_t u,size_t b,size_t h,size_t e,size_t r,int type,
    uint8_t *seed,size_t *positions,size_t output_offset) {
    size_t q=b-h;
    while (q) {
        size_t state=block_state(p,u,b,h,e,r),pos=u+b-1,row=q-1;
        unsigned code=type==1?(p->fp[state]&15u):(p->fp[state]>>4),a;
        if (!code) { if (!h) return -1; --h; --b; continue; }
        --code; a=code/2u; type=(int)(code%2u)+1;
        seed[output_offset+q-1]=(uint8_t)a;
        positions[output_offset+q-1]=pos;
        if (a!=x[pos]) --e;
        r=sub_mod(r,t->CF[4*row+a],p->d); --q; --b;
    }
    return 0;
}

static int traceback_reverse(const uint8_t *x,const SubseqHash2I8Tables *t,
    const SharedDP *p,size_t u,size_t b,size_t h,size_t e,size_t r,int type,
    uint8_t *seed,size_t *positions) {
    size_t cursor=0,q=b-h;
    while (q) {
        size_t state=block_state(p,u,b,h,e,r),row=q-1;
        unsigned code=type==1?(p->rp[state]&15u):(p->rp[state]>>4),a;
        if (!code) { if (!h) return -1; --h; ++u; --b; continue; }
        --code; a=code/2u; type=(int)(code%2u)+1;
        seed[cursor]=(uint8_t)a; positions[cursor]=u;
        if (a!=x[u]) --e;
        r=sub_mod(r,t->CR[4*row+a],p->d);
        ++cursor; ++u; --q; --b;
    }
    return 0;
}

static int emit_one(const uint8_t *x,size_t w,size_t n,size_t S,size_t order,
                    const SubseqHash2I8Tables *t,const SharedDP *p,
                    SubseqHashEditResult *out,size_t *next_right) {
    size_t H=p->E-S,K=n-H,i=order+1,b,a,eL,rL;
    int have=0; size_t best_r=0,best_b=0,best_a=0,best_eL=0;
    size_t best_rL=0,best_rR=0; int best_lt=0,best_rt=0; int16_t best_w=0;
    for (b=i;b<=n-(K-i);++b) {
        size_t hL=b-i,hR=(n-b)-(K-i),pivot=w+b-1;
        for (a=0;a<4;++a) {
            size_t delta=(a!=x[pivot]);
            if (delta>S) continue;
            for (eL=0;eL+delta<=S;++eL) {
                size_t eR=S-delta-eL,q,min_right;
                int lt=t->BP1[4*order+a]>0?2:1;
                int rt=t->BP2[4*order+a]>0?2:1;
                next_right[p->d]=SIZE_MAX;
                for (q=p->d;q>0;--q) {
                    size_t rs=block_state(p,w+b,n-b,hR,eR,q-1);
                    int16_t rv=rt==1?p->fmn[rs]:p->fmx[rs];
                    int reachable=rt==1?rv!=SH2_POS_INF:rv!=SH2_NEG_INF;
                    next_right[q-1]=reachable?q-1:next_right[q];
                }
                min_right=next_right[0];
                if (min_right==SIZE_MAX) continue;
                for (rL=0;rL<p->d;++rL) {
                    size_t ls=block_state(p,w,b-1,hL,eL,rL);
                    int16_t lv=lt==1?p->rmn[ls]:p->rmx[ls];
                    size_t base,threshold,rR,rs,residue; int16_t rv,score;
                    if ((lt==1&&lv==SH2_POS_INF)||(lt==2&&lv==SH2_NEG_INF))
                        continue;
                    base=rL+t->CP[4*order+a];
                    if (base>=p->d) base-=p->d;
                    threshold=base? p->d-base:0;
                    rR=next_right[threshold];
                    if (rR==SIZE_MAX) rR=min_right;
                    rs=block_state(p,w+b,n-b,hR,eR,rR);
                    rv=rt==1?p->fmn[rs]:p->fmx[rs];
                    residue=base+rR;
                    if (residue>=p->d) residue-=p->d;
                    score=(int16_t)(t->BP1[4*order+a]*lv+t->AP[4*order+a]+
                                    t->BP2[4*order+a]*rv);
                    if (!have || residue<best_r ||
                        (residue==best_r && score>best_w)) {
                        have=1; best_r=residue; best_w=score; best_b=b;
                        best_a=a; best_eL=eL; best_rL=rL; best_rR=rR;
                        best_lt=lt; best_rt=rt;
                    }
                }
            }
        }
    }
    if (!have) return -1;
    out->deletions=H; out->substitutions=S; out->length=K;
    out->psi=best_r; out->omega=best_w;
    out->seed=malloc(K?K:1); out->positions=malloc((K?K:1)*sizeof(*out->positions));
    out->substituted=malloc(K?K:1);
    if (!out->seed||!out->positions||!out->substituted) return -1;
    {
        size_t hL=best_b-i,hR=(n-best_b)-(K-i);
        size_t eR=S-(best_a!=x[w+best_b-1])-best_eL,q;
        if (traceback_reverse(x,t,p,w,best_b-1,hL,best_eL,best_rL,
                              best_lt,out->seed,out->positions) ||
            traceback_forward(x,t,p,w+best_b,n-best_b,hR,eR,best_rR,
                              best_rt,out->seed,out->positions,i)) return -1;
        out->seed[i-1]=(uint8_t)best_a; out->positions[i-1]=w+best_b-1;
        for (q=0;q<K;++q)
            out->substituted[q]=(uint8_t)(out->seed[q]!=x[out->positions[q]]);
    }
    return 0;
}

int subseqhash2_compute_long_edits_i8(const uint8_t *x,size_t N,size_t n,
    size_t E,const SubseqHash2I8Tables *t,SubseqHashEditResult *results) {
    SharedDP p={0}; size_t windows,per,total,w,s,order,index=0; int status=-1;
    size_t *next_right=NULL;
    if (!t||!results||!t->d||!n||n>N||E>n||t->k<n||(N&&!x)) return -1;
    for (w=0;w<N;++w) if (x[w]>3) return -1;
    windows=N-n+1; per=subseqhash2_i8_results_per_window(n,E);
    if (per==SIZE_MAX || windows>SIZE_MAX/per) return -1;
    total=windows*per;
    if (total>SIZE_MAX/sizeof(*results)) return -1;
    memset(results,0,total*sizeof(*results));
    if (allocate_dp(&p,N,n,E,t->d)) goto done;
    if (t->d==SIZE_MAX ||
        !(next_right=malloc((t->d+1)*sizeof(*next_right)))) goto done;
    fill_dp(x,t,&p);
    for (w=0;w<windows;++w) for (s=0;s<=E;++s) {
        size_t K=n-E+s;
        for (order=0;order<K;++order,++index)
            if (emit_one(x,w,n,s,order,t,&p,results+index,next_right)) goto done;
    }
    status=0;
done:
    if (status) subseqhash_edit_results_destroy(results,total);
    free(next_right);
    free_dp(&p); return status;
}
