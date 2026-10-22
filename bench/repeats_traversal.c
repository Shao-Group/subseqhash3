#define _POSIX_C_SOURCE 200809L
#include "ssh3_pivots.h"
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
    const SubseqHash3PivotTables *,size_t,size_t,
    SubseqHash3EditResult *,SubseqHash3PivotWorkspace *);
static Compute compute[2]={subseqhash3_compute_outward_pivot_repeats_with_workspace,
    subseqhash3_compute_pivot_repeats_with_workspace};
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
    SubseqHash3PivotTables t={0};
    SubseqHash3PivotWorkspace *w[2]={subseqhash3_pivot_workspace_create(),
        subseqhash3_pivot_workspace_create()};
    SubseqHash3EditResult out[256]={{0}};
    uint8_t inputs[16][257]; uint64_t rng=812733;
    double times[2][7];
    assert(w[0] && w[1]);
    assert(!subseqhash3_pivot_tables_init(&t,n,d,43));
    for (size_t q=0;q<16;++q) {
        for (size_t j=0;j<n;++j) inputs[q][j]=(uint8_t)(random64(&rng)&3);
        for (size_t mode=0;mode<2;++mode) {
            assert(!compute[mode](inputs[q],n,e,&t,1,count,out,w[mode]));
            subseqhash3_edit_results_destroy(out,count);
        }
    }
    for (size_t round=0;round<7;++round) for (size_t step=0;step<2;++step) {
        const size_t mode=(round+step)%2;
        const double start=now();
        for (size_t pass=0;pass<3;++pass) for (size_t q=0;q<16;++q) {
            assert(!compute[mode](inputs[q],n,e,&t,1,count,out,w[mode]));
            checksum+=out[0].psi+(uint64_t)out[0].omega;
            subseqhash3_edit_results_destroy(out,count);
        }
        times[mode][round]=(now()-start)*1e6/48;
    }
    for (size_t mode=0;mode<2;++mode) qsort(times[mode],7,sizeof(double),compare);
    printf("%zu,%zu,%zu,%zu,%.3f,%.3f,%.2f\n",n,e,d,count,
        times[0][3],times[1][3],times[0][3]/times[1][3]);
    fflush(stdout);
    for (size_t m=0;m<2;++m) subseqhash3_pivot_workspace_destroy(w[m]);
    subseqhash3_pivot_tables_destroy(&t);
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
static Hits hits(const SubseqHash3EditResult *x,const SubseqHash3EditResult *y,
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
    size_t changes,size_t families,size_t batch,size_t pair_offset,FILE *raw) {
    const char *names[]={"substitutions","delete_insert","unrelated"};
    const size_t k=n-1,total=families*batch;
    SubseqHash3PivotWorkspace *w[2]={subseqhash3_pivot_workspace_create(),
        subseqhash3_pivot_workspace_create()};
    SubseqHash3EditResult xout[256]={{0}},yout[256]={{0}};
    size_t same[2]={0},any[2]={0},per[2]={0};
    size_t inward_only=0,outward_only=0,nonidentical=0;
    double sum=0,squares=0,anysum=0,anysquares=0;
    assert(w[0] && w[1]);
    for (size_t f=0;f<families;++f) {
        uint64_t table_seed=UINT64_C(872341)+f*UINT64_C(739391);
        uint64_t rng=UINT64_C(918273)+n*1009+e*101+UINT64_C(1024)*17+kind*9127+changes*7919+f*UINT64_C(893417);
        SubseqHash3PivotTables t={0};
        Hits family[2]={{0},{0}};
        size_t family_nonidentical=0,family_inward_only=0,family_outward_only=0;
        assert(!subseqhash3_pivot_tables_init(&t,n,d,table_seed));
        /* The canonical pair stream is independent of D. Offsets let workers
           evaluate disjoint chunks without changing the experiment. */
        for (size_t b=0;b<pair_offset;++b) {
            uint8_t x[256],y[256]; pair(x,y,n,kind,changes,&rng);
        }
        for (size_t b=0;b<batch;++b) {
            uint8_t x[256],y[256]; Hits h[2];
            pair(x,y,n,kind,changes,&rng);
            const size_t different=memcmp(x,y,n)!=0;
            nonidentical+=different; family_nonidentical+=different;
            /* Alternate call order; both modes see identical pairs and tables. */
            for (size_t step=0;step<2;++step) {
                size_t mode=(step+b)%2;
                assert(!compute[mode](x,n,e,&t,1,count,xout,w[mode]));
                assert(!compute[mode](y,n,e,&t,1,count,yout,w[mode]));
                h[mode]=hits(xout,yout,count,k);
                family[mode].same+=h[mode].same; family[mode].any+=h[mode].any;
                family[mode].per+=h[mode].per;
                subseqhash3_edit_results_destroy(xout,count);
                subseqhash3_edit_results_destroy(yout,count);
            }
            const size_t in_only=h[1].same && !h[0].same;
            const size_t out_only=h[0].same && !h[1].same;
            inward_only+=in_only; outward_only+=out_only;
            family_inward_only+=in_only; family_outward_only+=out_only;
        }
        double delta=((double)family[1].same-family[0].same)/batch;
        double anydelta=((double)family[1].any-family[0].any)/batch;
        sum+=delta; squares+=delta*delta; anysum+=anydelta; anysquares+=anydelta*anydelta;
        for (size_t m=0;m<2;++m) {
            same[m]+=family[m].same; any[m]+=family[m].any; per[m]+=family[m].per;
        }
        if (raw) fprintf(raw,"%zu,%zu,%zu,%zu,%s,%zu,%zu,%llu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu\n",
            n,e,d,count,names[kind],changes,f,(unsigned long long)table_seed,batch,
            family[0].same,family[1].same,family[0].any,family[1].any,family[0].per,family[1].per,
            family_nonidentical,family_inward_only,family_outward_only,pair_offset);
        subseqhash3_pivot_tables_destroy(&t);
    }
    const double delta=100*sum/families,hw=100*halfwidth(sum,squares,families);
    const double adelta=100*anysum/families,ahw=100*halfwidth(anysum,anysquares,families);
    printf("%zu,%zu,%zu,%zu,%s,%zu,%zu,%zu,%zu,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%zu,%zu\n",
        n,e,d,count,names[kind],changes,families,total,nonidentical,
        100.0*same[0]/total,100.0*same[1]/total,delta,delta-hw,delta+hw,
        100.0*any[0]/total,100.0*any[1]/total,adelta,adelta-ahw,adelta+ahw,
        100.0*per[0]/(total*count),100.0*per[1]/(total*count),inward_only,outward_only);
    fflush(stdout);
    for (size_t m=0;m<2;++m) subseqhash3_pivot_workspace_destroy(w[m]);
}
static size_t argument(const char *value) {
    char *end;
    const unsigned long long number=strtoull(value,&end,10);
    assert(*value && *value!='-' && !*end && number<=SIZE_MAX);
    return (size_t)number;
}
int main(int argc,char **argv) {
    if ((argc==2 || argc==3) && !strcmp(argv[1],"--runtime")) {
        puts("n,substitutions,D,seeds,outward_us,inward_us,outward_over_inward");
        const size_t residues[]={argc==3?argument(argv[2]):31,1024};
        assert(residues[0]);
        for (size_t index=0;index<(argc==3?1:2);++index) {
            const size_t d=residues[index];
            runtime(30,4,d,10); runtime(30,4,d,29);
            runtime(50,8,d,10); runtime(50,8,d,49);
        }
        if (argc==2) runtime(30,0,1024,10);
        fprintf(stderr,"runtime checksum=%llu\n",(unsigned long long)checksum);
    } else if ((argc==2 || (argc>=5 && argc<=8)) && !strcmp(argv[1],"--collisions")) {
        size_t families=32,batch=64; FILE *raw=NULL;
        const size_t d=argc>=6?argument(argv[5]):1024;
        const size_t selected=argc>=7?argument(argv[6]):SIZE_MAX;
        const size_t offset=argc>=8?argument(argv[7]):0;
        assert(d && (selected==SIZE_MAX || selected<4) && offset<=100000);
        if (argc>=5) {
            families=argument(argv[2]); batch=argument(argv[3]);
            assert(families>=32 && families<=1024 && batch && batch<=100000);
            raw=fopen(argv[4],"w"); assert(raw);
            fputs("n,substitutions,D,seeds,mutation,changes,family,table_seed,pairs,outward_same,inward_same,outward_any,inward_any,outward_per,inward_per,nonidentical_pairs,inward_only,outward_only,pair_offset\n",raw);
        }
        puts("n,substitutions,D,seeds,mutation,changes,families,pairs,nonidentical_pairs,outward_same_pct,inward_same_pct,delta_same_pp,delta_same_low95,delta_same_high95,outward_any_pct,inward_any_pct,delta_any_pp,delta_any_low95,delta_any_high95,outward_per_pct,inward_per_pct,inward_only,outward_only");
        for (size_t config=0;config<4;++config) {
            if (selected!=SIZE_MAX && config!=selected) continue;
            size_t n=config<2?30:50,e=config<2?4:8;
            size_t count=config%2?n-1:10;
            for (unsigned kind=0;kind<3;++kind) {
                for (size_t scenario=0;scenario<(kind==2?1:2);++scenario) {
                    size_t changes=kind==0?(scenario?e:e/2):kind==1?scenario+1:0;
                    collisions(n,e,d,count,kind,changes,families,batch,offset,raw);
                }
            }
        }
        if (raw) assert(!fclose(raw));
    } else {
        fprintf(stderr,"Usage: %s --runtime [D] | --collisions [FAMILIES PAIRS_PER_FAMILY RAW_CSV [D [CONFIG [PAIR_OFFSET]]]]\n",argv[0]);
        return 1;
    }
}
