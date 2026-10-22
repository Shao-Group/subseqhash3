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

typedef int (*Compute)(const uint8_t *,size_t,size_t,
    const SubseqHash3DeletionTables *,const size_t *,size_t,
    SubseqHash3DeletionResult *,SubseqHash3DeletionWorkspace *);
static Compute compute[2]={subseqhash3_compute_deletion_seeds_with_workspace,
    subseqhash3_compute_inward_deletion_seeds_with_workspace};
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
static void runtime(size_t n,size_t e,size_t d,size_t count) {
    SubseqHash3PivotTables shared={0};
    SubseqHash3DeletionTables t={0};
    SubseqHash3DeletionWorkspace *w[2]={subseqhash3_deletion_workspace_create(),
        subseqhash3_deletion_workspace_create()};
    SubseqHash3PivotWorkspace *pw=subseqhash3_pivot_workspace_create();
    SubseqHash3DeletionResult out[256]={{0}};
    SubseqHash3EditResult rep[256]={{0}};
    uint8_t inputs[16][256]; uint64_t rng=812733;
    double times[3][7];
    assert(w[0] && w[1] && pw);
    assert(!subseqhash3_pivot_tables_init(&shared,n+1,d,43));
    assert(!subseqhash3_deletion_tables_from_pivots(&t,&shared));
    for (size_t q=0;q<16;++q) {
        for (size_t j=0;j<n;++j) inputs[q][j]=(uint8_t)(random64(&rng)&3);
        for (size_t mode=0;mode<2;++mode) {
            assert(!compute[mode](inputs[q],n,e,&t,NULL,count,out,w[mode]));
            subseqhash3_deletion_results_destroy(out,count);
        }
        assert(!subseqhash3_compute_pivot_repeats_with_workspace(inputs[q],n,e,&shared,1,count,rep,pw));
        subseqhash3_edit_results_destroy(rep,count);
    }
    for (size_t round=0;round<7;++round) for (size_t step=0;step<3;++step) {
        const size_t mode=(round+step)%3;
        const double start=now();
        for (size_t pass=0;pass<3;++pass) for (size_t q=0;q<16;++q) {
            if (mode==2) {
                assert(!subseqhash3_compute_pivot_repeats_with_workspace(inputs[q],n,e,&shared,1,count,rep,pw));
                checksum+=rep[0].psi+(uint64_t)rep[0].omega;
                subseqhash3_edit_results_destroy(rep,count);
            } else {
                assert(!compute[mode](inputs[q],n,e,&t,NULL,count,out,w[mode]));
                checksum+=out[0].psi+(uint64_t)out[0].omega;
                subseqhash3_deletion_results_destroy(out,count);
            }
        }
        times[mode][round]=(now()-start)*1e6/48;
    }
    for (size_t mode=0;mode<3;++mode) qsort(times[mode],7,sizeof(double),compare);
    printf("%zu,%zu,%zu,%zu,%.3f,%.3f,%.2f,%.3f,%.2f\n",n,e,d,count,
        times[0][3],times[1][3],times[0][3]/times[1][3],times[2][3],times[2][3]/times[1][3]);
    fflush(stdout);
    subseqhash3_pivot_workspace_destroy(pw);
    for (size_t m=0;m<2;++m) subseqhash3_deletion_workspace_destroy(w[m]);
    subseqhash3_deletion_tables_destroy(&t); subseqhash3_pivot_tables_destroy(&shared);
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
typedef struct { size_t same,any,per; } Hits;
static Hits hits(const SubseqHash3DeletionResult *x,const SubseqHash3DeletionResult *y,
    size_t count,size_t k) {
    Hits h={0};
    for (size_t i=0;i<count;++i) {
        if (!memcmp(x[i].seed,y[i].seed,k)) { h.same=1; ++h.per; }
        if (!h.any) for (size_t j=0;j<count;++j)
            if (!memcmp(x[i].seed,y[j].seed,k)) { h.any=1; break; }
    }
    return h;
}
/* Student t 0.975 at df=31 is about 2.04; conservative for F >= 32. */
static double halfwidth(double sum,double squares,size_t families) {
    const double mean=sum/families;
    double variance=(squares-families*mean*mean)/(families-1);
    if (variance<0) variance=0;
    return 2.04*sqrt(variance/families);
}
static void collisions(size_t n,size_t e,size_t d,size_t count,unsigned kind,
    size_t changes,size_t families,size_t batch,FILE *raw) {
    const char *names[]={"substitutions","delete_insert","unrelated"};
    const size_t k=n-e,total=families*batch;
    SubseqHash3DeletionWorkspace *w[2]={subseqhash3_deletion_workspace_create(),
        subseqhash3_deletion_workspace_create()};
    SubseqHash3DeletionResult xout[256]={{0}},yout[256]={{0}};
    size_t same[2]={0},any[2]={0},per[2]={0};
    size_t inward_only=0,outward_only=0,nonidentical=0;
    double sum=0,squares=0,anysum=0,anysquares=0;
    assert(w[0] && w[1]);
    for (size_t f=0;f<families;++f) {
        uint64_t table_seed=UINT64_C(872341)+f*UINT64_C(739391);
        uint64_t rng=UINT64_C(918273)+n*1009+e*101+d*17+kind*9127+changes*7919+f*UINT64_C(893417);
        SubseqHash3DeletionTables t={0};
        Hits family[2]={{0},{0}};
        assert(!subseqhash3_deletion_tables_init(&t,k,d,table_seed));
        for (size_t b=0;b<batch;++b) {
            uint8_t x[256],y[256]; Hits h[2];
            pair(x,y,n,kind,changes,&rng); nonidentical+=memcmp(x,y,n)!=0;
            /* Alternate call order; both modes see identical pairs and tables. */
            for (size_t step=0;step<2;++step) {
                size_t mode=(step+b)%2;
                assert(!compute[mode](x,n,e,&t,NULL,count,xout,w[mode]));
                assert(!compute[mode](y,n,e,&t,NULL,count,yout,w[mode]));
                h[mode]=hits(xout,yout,count,k);
                family[mode].same+=h[mode].same; family[mode].any+=h[mode].any;
                family[mode].per+=h[mode].per;
                subseqhash3_deletion_results_destroy(xout,count);
                subseqhash3_deletion_results_destroy(yout,count);
            }
            inward_only+=h[1].same && !h[0].same;
            outward_only+=h[0].same && !h[1].same;
        }
        double delta=((double)family[1].same-family[0].same)/batch;
        double anydelta=((double)family[1].any-family[0].any)/batch;
        sum+=delta; squares+=delta*delta; anysum+=anydelta; anysquares+=anydelta*anydelta;
        for (size_t m=0;m<2;++m) {
            same[m]+=family[m].same; any[m]+=family[m].any; per[m]+=family[m].per;
        }
        if (raw) fprintf(raw,"%zu,%zu,%zu,%zu,%s,%zu,%zu,%llu,%zu,%zu,%zu,%zu,%zu,%zu,%zu\n",
            n,e,d,count,names[kind],changes,f,(unsigned long long)table_seed,batch,
            family[0].same,family[1].same,family[0].any,family[1].any,family[0].per,family[1].per);
        subseqhash3_deletion_tables_destroy(&t);
    }
    const double delta=100*sum/families,hw=100*halfwidth(sum,squares,families);
    const double adelta=100*anysum/families,ahw=100*halfwidth(anysum,anysquares,families);
    printf("%zu,%zu,%zu,%zu,%s,%zu,%zu,%zu,%zu,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%zu,%zu\n",
        n,e,d,count,names[kind],changes,families,total,nonidentical,
        100.0*same[0]/total,100.0*same[1]/total,delta,delta-hw,delta+hw,
        100.0*any[0]/total,100.0*any[1]/total,adelta,adelta-ahw,adelta+ahw,
        100.0*per[0]/(total*count),100.0*per[1]/(total*count),inward_only,outward_only);
    fflush(stdout);
    for (size_t m=0;m<2;++m) subseqhash3_deletion_workspace_destroy(w[m]);
}
int main(int argc,char **argv) {
    if (argc==2 && !strcmp(argv[1],"--runtime")) {
        puts("n,deletions,D,seeds,outward_us,inward_us,outward_over_inward,repeats_us,repeats_over_inward");
        for (size_t d=31;d<=1024;d=d==31?1024:1025) {
            runtime(30,4,d,10); runtime(30,4,d,26);
            runtime(50,8,d,10); runtime(50,8,d,42);
        }
        runtime(30,1,1024,10);
        fprintf(stderr,"runtime checksum=%llu\n",(unsigned long long)checksum);
    } else if ((argc==2 || argc==5) && !strcmp(argv[1],"--collisions")) {
        size_t families=32,batch=64; FILE *raw=NULL;
        if (argc==5) {
            families=(size_t)strtoul(argv[2],NULL,10); batch=(size_t)strtoul(argv[3],NULL,10);
            assert(families>=32 && families<=1024 && batch && batch<=100000);
            raw=fopen(argv[4],"w"); assert(raw);
            fputs("n,deletions,D,seeds,mutation,changes,family,table_seed,pairs,outward_same,inward_same,outward_any,inward_any,outward_per,inward_per\n",raw);
        }
        puts("n,deletions,D,seeds,mutation,changes,families,pairs,nonidentical_pairs,outward_same_pct,inward_same_pct,delta_same_pp,delta_same_low95,delta_same_high95,outward_any_pct,inward_any_pct,delta_any_pp,delta_any_low95,delta_any_high95,outward_per_pct,inward_per_pct,inward_only,outward_only");
        for (size_t config=0;config<4;++config) {
            size_t n=config<2?30:50,e=config<2?4:8;
            size_t count=config%2?n-e:10;
            for (unsigned kind=0;kind<3;++kind) {
                for (size_t scenario=0;scenario<(kind==2?1:2);++scenario) {
                    size_t changes=kind==0?(scenario?e:e/2):kind==1?scenario+1:0;
                    collisions(n,e,1024,count,kind,changes,families,batch,raw);
                }
            }
        }
        if (raw) assert(!fclose(raw));
    } else {
        fprintf(stderr,"Usage: %s --runtime | --collisions [FAMILIES PAIRS_PER_FAMILY RAW_CSV]\n",argv[0]);
        return 1;
    }
}
