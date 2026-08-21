#include "subseqhash.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static uint64_t next_random(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
static size_t table_index(const SubseqHashTables *t,size_t i,size_t r,size_t a) {
    return (i*4+a)*t->d+r;
}
int subseqhash_tables_init(SubseqHashTables *t,size_t k,size_t d,uint64_t seed) {
    size_t cells,i,r,a;
    static const int8_t b1[4]={1,1,-1,-1}, b2[4]={1,-1,1,-1};
    if (!t || !d || k>SIZE_MAX/d || k*d>SIZE_MAX/4 ||
        k*d*4>SIZE_MAX/sizeof(*t->A)) return -1;
    memset(t,0,sizeof(*t)); t->k=k; t->d=d; cells=k*d*4;
    t->A=malloc((cells?cells:1)*sizeof(*t->A));
    t->B1=malloc(cells?cells:1); t->B2=malloc(cells?cells:1);
    t->C=malloc((k?k*4:1)*sizeof(*t->C));
    if (!t->A || !t->B1 || !t->B2 || !t->C) { subseqhash_tables_destroy(t); return -1; }
    for (i=0;i<cells;++i)
        t->A[i]=1073741824.0+(double)(next_random(&seed)&UINT64_C(0xffffffff))/4294967296.0*1073741824.0;
    for (i=0;i<k;++i) for (r=0;r<d;++r) {
        unsigned permutation[4]={0,1,2,3};
        for (a=4;a>1;--a) { size_t q=(size_t)(next_random(&seed)%a); unsigned x=permutation[a-1]; permutation[a-1]=permutation[q]; permutation[q]=x; }
        for (a=0;a<4;++a) { size_t x=table_index(t,i,r,a); t->B1[x]=b1[permutation[a]]; t->B2[x]=b2[permutation[a]]; }
    }
    for (i=0;i<k;++i) for (a=0;a<4;++a) {
        size_t value;
        do { value=(size_t)(next_random(&seed)%d); }
        while (d>=4 && a>0 && (value==t->C[4*i] || (a>1&&value==t->C[4*i+1]) || (a>2&&value==t->C[4*i+2])));
        t->C[4*i+a]=value;
    }
    return 0;
}
void subseqhash_tables_destroy(SubseqHashTables *t) {
    if (!t) return; free(t->C); free(t->B2); free(t->B1); free(t->A); memset(t,0,sizeof(*t));
}
int subseqhash_score(const uint8_t *seed,const SubseqHashTables *t,size_t *psi,double *omega) {
    size_t i,r=0; double w=0;
    if (!t || !psi || !omega || (t->k && !seed) || !t->d) return -1;
    for (i=0;i<t->k;++i) { size_t x; if (seed[i]>3) return -1; r=(r+t->C[4*i+seed[i]])%t->d; x=table_index(t,i,r,seed[i]); w=w*t->B1[x]+t->A[x]*t->B2[x]; }
    *psi=r; *omega=w; return 0;
}
int subseqhash_compute(const uint8_t *s,size_t n,const SubseqHashTables *t,uint8_t *seed,size_t *positions,size_t *psi,double *omega) {
    size_t states,l,j,r,k,d; double *mn,*mx; uint8_t *pm,*px; int result=-1;
    if (!t || !psi || !omega || !t->d || (n&&!s) || (t->k&&!seed) || t->k>n) return -1;
    k=t->k; d=t->d;
    if ((n+1)>SIZE_MAX/(k+1) || (n+1)*(k+1)>SIZE_MAX/d) return -1;
    states=(n+1)*(k+1)*d;
    if (states>SIZE_MAX/sizeof(*mn)) return -1;
    mn=malloc(states*sizeof(*mn)); mx=malloc(states*sizeof(*mx)); pm=calloc(states,1); px=calloc(states,1);
    if (!mn || !mx || !pm || !px) goto done;
#define IDX(L,J,R) (((L)*(k+1)+(J))*d+(R))
    for (l=0;l<states;++l) { mn[l]=DBL_MAX; mx[l]=-DBL_MAX; }
    mn[IDX(0,0,0)]=mx[IDX(0,0,0)]=0;
    for (l=1;l<=n;++l) {
        size_t a=s[l-1], upper=l<k?l:k; if (a>3) goto done;
        for (j=0;j<=upper;++j) for (r=0;r<d;++r) {
            size_t cur=IDX(l,j,r),skip=IDX(l-1,j,r);
            if (mn[skip]!=DBL_MAX) { mn[cur]=mn[skip]; mx[cur]=mx[skip]; }
            if (j) {
                size_t c=t->C[4*(j-1)+a],q=(r+d-c)%d,prev=IDX(l-1,j-1,q),x=table_index(t,j-1,r,a);
                if (mn[prev]!=DBL_MAX) {
                    double amin=t->B1[x]>0?mn[prev]:-mx[prev], amax=t->B1[x]>0?mx[prev]:-mn[prev];
                    amin+=t->A[x]*t->B2[x]; amax+=t->A[x]*t->B2[x];
                    if (amin<mn[cur]) { mn[cur]=amin; pm[cur]=(uint8_t)(t->B1[x]>0?1:2); }
                    if (amax>mx[cur]) { mx[cur]=amax; px[cur]=(uint8_t)(t->B1[x]>0?2:1); }
                }
            }
        }
    }
    for (r=0;r<d && mn[IDX(n,k,r)]==DBL_MAX;++r) {}
    if (r==d) goto done;
    *psi=r;
    { int type=fabs(mn[IDX(n,k,r)])>=fabs(mx[IDX(n,k,r)])?1:2; *omega=type==1?mn[IDX(n,k,r)]:mx[IDX(n,k,r)]; l=n; j=k;
      while (j) { size_t cur=IDX(l,j,r); uint8_t parent=type==1?pm[cur]:px[cur]; if (!parent) { --l; continue; }
        seed[j-1]=s[l-1]; if (positions) positions[j-1]=l-1; r=(r+d-t->C[4*(j-1)+s[l-1]])%d; type=parent; --j; --l; }
    }
    result=0;
done:
    free(px); free(pm); free(mx); free(mn); return result;
#undef IDX
}

int subseqhash_compute_optimized(const uint8_t *s,size_t n,const SubseqHashTables *t,
                                 uint8_t *seed,size_t *positions,size_t *psi,double *omega) {
    size_t k,d,plane_cells,parent_cells,parent_stride,l,j,q,r;
    double *storage=NULL,*mn,*mx;
    uint8_t *parent=NULL;
    int result=-1;

    if (!t || !psi || !omega || !t->d || (n&&!s) ||
        (t->k&&!seed) || t->k>n) return -1;
    k=t->k; d=t->d;
    if (k+1==0 || k+1>SIZE_MAX/d) return -1;
    plane_cells=(k+1)*d;
    if (plane_cells>SIZE_MAX/(2*sizeof(*storage))) return -1;
    parent_stride=d/2+d%2;
    if (k && (n>SIZE_MAX/k || n*k>SIZE_MAX/parent_stride)) return -1;
    parent_cells=n*k*parent_stride;

    storage=malloc(2*plane_cells*sizeof(*storage));
    parent=calloc(parent_cells?parent_cells:1,1);
    if (!storage || !parent) goto done;
    mn=storage;
    mx=mn+plane_cells;
    for (q=0;q<plane_cells;++q) {
        mn[q]=DBL_MAX;
        mx[q]=-DBL_MAX;
    }
    mn[0]=mx[0]=0.0;

    for (l=1;l<=n;++l) {
        const size_t a=s[l-1],upper=l<k?l:k;
        if (a>3) goto done;
        for (j=upper;j>0;--j) {
            double *dst_mn=mn+j*d,*dst_mx=mx+j*d;
            const double *src_mn=mn+(j-1)*d,*src_mx=mx+(j-1)*d;
            const size_t c=t->C[4*(j-1)+a];
            const size_t table_base=((j-1)*4+a)*d;
            uint8_t *decisions=parent+((l-1)*k+(j-1))*parent_stride;
            for (q=0;q<d;++q) {
                size_t x,byte_index; unsigned shift,nibble; double lo,hi,add;
                uint8_t *decision_byte;
                if (src_mn[q]==DBL_MAX) continue;
                r=q+c; if (r>=d) r-=d;
                x=table_base+r;
                add=t->A[x]*(double)t->B2[x];
                byte_index=r>>1; shift=(unsigned)(r&1u)*4u;
                decision_byte=decisions+byte_index;
                nibble=(unsigned)(*decision_byte>>shift)&15u;
                if (t->B1[x]>0) {
                    lo=src_mn[q]+add; hi=src_mx[q]+add;
                    if (lo<dst_mn[r]) { dst_mn[r]=lo; nibble=(nibble&12u)|1u; }
                    if (hi>dst_mx[r]) { dst_mx[r]=hi; nibble=(nibble&3u)|8u; }
                } else {
                    lo=-src_mx[q]+add; hi=-src_mn[q]+add;
                    if (lo<dst_mn[r]) { dst_mn[r]=lo; nibble=(nibble&12u)|2u; }
                    if (hi>dst_mx[r]) { dst_mx[r]=hi; nibble=(nibble&3u)|4u; }
                }
                *decision_byte=(uint8_t)((*decision_byte&~(15u<<shift))|(nibble<<shift));
            }
        }
    }

    for (r=0;r<d && mn[k*d+r]==DBL_MAX;++r) {}
    if (r==d) goto done;
    *psi=r;
    {
        int type=fabs(mn[k*d+r])>=fabs(mx[k*d+r])?1:2;
        *omega=type==1?mn[k*d+r]:mx[k*d+r];
        l=n; j=k;
        while (j) {
            const uint8_t packed=parent[((l-1)*k+(j-1))*parent_stride+(r>>1)];
            const uint8_t decision=(uint8_t)((packed>>((r&1u)*4u))&15u);
            const unsigned p=type==1?(decision&3u):((decision>>2)&3u);
            if (!p) { --l; continue; }
            seed[j-1]=s[l-1];
            if (positions) positions[j-1]=l-1;
            q=t->C[4*(j-1)+s[l-1]];
            r=r>=q?r-q:r+d-q;
            type=(int)p; --j; --l;
        }
    }
    result=0;
done:
    free(parent); free(storage); return result;
}

int subseqhash_compute_leftmost(const uint8_t *s,size_t n,const SubseqHashTables *t,
                                uint8_t *seed,size_t *positions,size_t *psi,double *omega) {
    size_t k,d,states,next_cells,i,j,r,a,best;
    size_t *dp=NULL,*next=NULL,*parent_residue=NULL;
    uint8_t *parent_symbol=NULL;
    int result=-1;

    if (!t || !psi || !omega || !t->d || (n&&!s) ||
        (t->k&&!seed) || t->k>n) return -1;
    k=t->k; d=t->d;
    if (k+1==0 || k+1>SIZE_MAX/d || n==SIZE_MAX ||
        n+1>SIZE_MAX/4) return -1;
    states=(k+1)*d; next_cells=(n+1)*4;
    if (states>SIZE_MAX/sizeof(*dp) ||
        states>SIZE_MAX/sizeof(*parent_residue) ||
        next_cells>SIZE_MAX/sizeof(*next)) return -1;
    dp=malloc(states*sizeof(*dp));
    next=malloc(next_cells*sizeof(*next));
    parent_residue=malloc(states*sizeof(*parent_residue));
    parent_symbol=malloc(states?states:1);
    if (!dp || !next || !parent_residue || !parent_symbol) goto done;

    for (a=0;a<4;++a) next[n*4+a]=n;
    for (i=n;i>0;--i) {
        if (s[i-1]>3) goto done;
        memcpy(next+(i-1)*4,next+i*4,4*sizeof(*next));
        next[(i-1)*4+s[i-1]]=i-1;
    }
    for (i=0;i<states;++i) {
        dp[i]=n+1;
        parent_residue[i]=SIZE_MAX;
        parent_symbol[i]=UINT8_MAX;
    }
    dp[0]=0;

    for (j=1;j<=k;++j) {
        size_t *previous=dp+(j-1)*d,*current=dp+j*d;
        for (r=0;r<d;++r) {
            const size_t previous_position=previous[r];
            if (previous_position>n) continue;
            for (a=0;a<4;++a) {
                const size_t found=next[previous_position*4+a];
                size_t nr,candidate,index;
                if (found==n) continue;
                nr=r+t->C[4*(j-1)+a]; if (nr>=d) nr-=d;
                candidate=found+1; index=j*d+nr;
                if (candidate<current[nr] ||
                    (candidate==current[nr] && a<parent_symbol[index])) {
                    current[nr]=candidate;
                    parent_residue[index]=r;
                    parent_symbol[index]=(uint8_t)a;
                }
            }
        }
    }
    for (best=0;best<d && dp[k*d+best]>n;++best) {}
    if (best==d) goto done;
    *psi=best;
    r=best;
    for (j=k;j>0;--j) {
        const size_t index=j*d+r;
        seed[j-1]=parent_symbol[index];
        if (positions) positions[j-1]=dp[index]-1;
        r=parent_residue[index];
    }
    if (subseqhash_score(seed,t,psi,omega)) goto done;
    result=0;
done:
    free(parent_symbol); free(parent_residue); free(next); free(dp);
    return result;
}
