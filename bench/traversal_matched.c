#define _POSIX_C_SOURCE 200809L
#include "ssh3_deletions.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static uint64_t checksum;
static uint64_t random64(uint64_t *s) {
    uint64_t z=(*s+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
static double now(void) {
    struct timespec t; assert(!clock_gettime(CLOCK_MONOTONIC,&t));
    return t.tv_sec+t.tv_nsec*1e-9;
}
static int compare(const void *a,const void *b) {
    double x=*(const double *)a,y=*(const double *)b;
    return (x>y)-(x<y);
}
/* Fixed length pairs: distinct substitutions, or b random deletions followed
   by b random insertions. Insertions can cancel deletions by chance. */
static void pair(uint8_t *x,uint8_t *y,size_t n,unsigned kind,size_t edits,uint64_t *rng) {
    size_t order[256],length=n;
    for (size_t j=0;j<n;++j) { x[j]=(uint8_t)(random64(rng)&3); order[j]=j; }
    memcpy(y,x,n);
    if (kind==2) {
        for (size_t j=0;j<n;++j) y[j]=(uint8_t)(random64(rng)&3);
        return;
    }
    for (size_t j=0;j<edits;++j) {
        size_t q=j+random64(rng)%(n-j),tmp=order[j]; order[j]=order[q]; order[q]=tmp;
    }
    if (kind==0) {
        for (size_t j=0;j<edits;++j) y[order[j]]=(uint8_t)((x[order[j]]+1+random64(rng)%3)&3);
    } else {
        uint8_t deleted[256]={0};
        for (size_t j=0;j<edits;++j) deleted[order[j]]=1;
        length=0;
        for (size_t j=0;j<n;++j) if (!deleted[j]) y[length++]=x[j];
        for (size_t j=0;j<edits;++j) {
            size_t q=random64(rng)%(length+1);
            memmove(y+q+1,y+q,length-q); y[q]=(uint8_t)(random64(rng)&3); ++length;
        }
    }
}

/* Mode 0 is outward, mode 1 inward. Each has its own D and workspace. */
typedef struct {
    size_t n,e,r; int repeats, inward;
    SubseqHash3PivotTables pt;
    SubseqHash3DeletionTables dt;
    SubseqHash3PivotWorkspace *pw;
    SubseqHash3DeletionWorkspace *dw;
    SubseqHash3EditResult pr[256];
    SubseqHash3DeletionResult dr[256];
} Engine;
static void init(Engine *a,int repeats,int inward,size_t n,size_t e,size_t r,size_t d,uint64_t seed) {
    memset(a,0,sizeof(*a));
    a->n=n; a->e=e; a->r=r; a->repeats=repeats; a->inward=inward;
    if (repeats) {
        assert(!subseqhash3_pivot_tables_init(&a->pt,n,d,seed));
        a->pw=subseqhash3_pivot_workspace_create(); assert(a->pw);
    } else {
        assert(!subseqhash3_deletion_tables_init(&a->dt,n-e,d,seed));
        a->dw=subseqhash3_deletion_workspace_create(); assert(a->dw);
    }
}
static void compute(Engine *a,const uint8_t *x) {
    if (a->repeats) {
        if (a->inward) assert(!subseqhash3_compute_pivot_repeats_with_workspace(x,a->n,a->e,&a->pt,1,a->r,a->pr,a->pw));
        else assert(!subseqhash3_compute_outward_pivot_repeats_with_workspace(x,a->n,a->e,&a->pt,1,a->r,a->pr,a->pw));
    } else {
        if (a->inward) assert(!subseqhash3_compute_inward_deletion_seeds_with_workspace(x,a->n,a->e,&a->dt,NULL,a->r,a->dr,a->dw));
        else assert(!subseqhash3_compute_deletion_seeds_with_workspace(x,a->n,a->e,&a->dt,NULL,a->r,a->dr,a->dw));
    }
}
static const uint8_t *seed(Engine *a,size_t j) { return a->repeats?a->pr[j].seed:a->dr[j].seed; }
static void clear(Engine *a) {
    if (a->repeats) subseqhash3_edit_results_destroy(a->pr,a->r);
    else subseqhash3_deletion_results_destroy(a->dr,a->r);
}
static void destroy(Engine *a) {
    if (a->repeats) {
        subseqhash3_pivot_workspace_destroy(a->pw); subseqhash3_pivot_tables_destroy(&a->pt);
    } else {
        subseqhash3_deletion_workspace_destroy(a->dw); subseqhash3_deletion_tables_destroy(&a->dt);
    }
}
static void runtime(int repeats,size_t n,size_t e,size_t r,size_t di,size_t dout) {
    Engine a[2]; uint8_t x[16][256]; uint64_t rng=812733; double times[2][7];
    for (int m=0;m<2;++m) init(&a[m],repeats,m,n,e,r,m?di:dout,43);
    for (size_t q=0;q<16;++q) {
        for (size_t j=0;j<n;++j) x[q][j]=(uint8_t)(random64(&rng)&3);
        for (int m=0;m<2;++m) { compute(&a[m],x[q]); clear(&a[m]); }
    }
    for (size_t round=0;round<7;++round) for (size_t step=0;step<2;++step) {
        size_t m=(round+step)%2; double start=now();
        for (size_t pass=0;pass<3;++pass) for (size_t q=0;q<16;++q) {
            compute(&a[m],x[q]); checksum+=seed(&a[m],0)[0]; clear(&a[m]);
        }
        times[m][round]=(now()-start)*1e6/48;
    }
    for (int m=0;m<2;++m) { qsort(times[m],7,sizeof(double),compare); destroy(&a[m]); }
    puts("inward_us,outward_us"); printf("%.3f,%.3f\n",times[1][3],times[0][3]);
}
static void collisions(int repeats,size_t n,size_t e,size_t r,size_t di,size_t dout,
                       size_t families,size_t batch,unsigned kind,size_t changes) {
    const size_t k=repeats?n-1:n-e;
    puts("family,table_seed,pairs,nonidentical_pairs,outward_same,inward_same,outward_any,inward_any,outward_per,inward_per,inward_only,outward_only");
    for (size_t f=0;f<families;++f) {
        uint64_t ts=UINT64_C(872341)+f*UINT64_C(739391);
        uint64_t rng=UINT64_C(918273)+n*1009+e*101+UINT64_C(1024)*17+kind*9127+changes*7919+f*UINT64_C(893417);
        Engine a[2]; size_t same[2]={0},any[2]={0},per[2]={0},only[2]={0},different=0;
        for (int m=0;m<2;++m) init(&a[m],repeats,m,n,e,r,m?di:dout,ts);
        for (size_t b=0;b<batch;++b) {
            uint8_t x[256],y[256],xs[256][256]; int matched[2];
            pair(x,y,n,kind,changes,&rng); different+=memcmp(x,y,n)!=0;
            for (size_t step=0;step<2;++step) {
                size_t m=(step+b)%2; int h=0,ha=0;
                compute(&a[m],x);
                for (size_t j=0;j<r;++j) memcpy(xs[j],seed(&a[m],j),k);
                clear(&a[m]); compute(&a[m],y);
                for (size_t j=0;j<r;++j) {
                    if (!memcmp(xs[j],seed(&a[m],j),k)) { h=1; ++per[m]; }
                    if (!ha) for (size_t l=0;l<r;++l) if (!memcmp(xs[j],seed(&a[m],l),k)) { ha=1; break; }
                }
                clear(&a[m]); same[m]+=h; any[m]+=ha; matched[m]=h;
            }
            only[0]+=matched[0] && !matched[1]; only[1]+=matched[1] && !matched[0];
        }
        printf("%zu,%llu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu\n",f,(unsigned long long)ts,batch,different,same[0],same[1],any[0],any[1],per[0],per[1],only[1],only[0]);
        fflush(stdout);
        for (int m=0;m<2;++m) destroy(&a[m]);
    }
}
static size_t argument(const char *s) {
    char *end; unsigned long long v=strtoull(s,&end,10);
    assert(*s && *s!='-' && !*end && v<=SIZE_MAX); return (size_t)v;
}
int main(int argc,char **argv) {
    if (argc!=8 && argc!=12) {
        fprintf(stderr,"Usage: %s runtime|collisions deletions|repeats N E R INWARD_D OUTWARD_D [FAMILIES PAIRS KIND CHANGES]\n",argv[0]); return 1;
    }
    assert(!strcmp(argv[2],"repeats") || !strcmp(argv[2],"deletions"));
    int rep=!strcmp(argv[2],"repeats");
    size_t n=argument(argv[3]),e=argument(argv[4]),r=argument(argv[5]),di=argument(argv[6]),dout=argument(argv[7]);
    assert(n>1 && n<=256 && e<n && r && r<=(rep?n-1:n-e) && di && dout);
    if (!strcmp(argv[1],"runtime") && argc==8) runtime(rep,n,e,r,di,dout);
    else {
        assert(!strcmp(argv[1],"collisions") && argc==12);
        size_t f=argument(argv[8]),b=argument(argv[9]),kind=argument(argv[10]),c=argument(argv[11]);
        assert(f>=32 && b && kind<=2 && c<=n);
        collisions(rep,n,e,r,di,dout,f,b,(unsigned)kind,c);
    }
    return 0;
}
