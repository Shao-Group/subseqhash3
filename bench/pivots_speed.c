#define _POSIX_C_SOURCE 200809L
#include "ssh3_pivots.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Link a frozen implementation with renamed exports. Both versions use
   identical tables, inputs, budgets, and output pivot orders. */
#ifndef PIVOT_BASELINE_LABEL
#define PIVOT_BASELINE_LABEL "ca35567"
#endif
extern int baseline_compute_pivot_seeds(const uint8_t *,size_t,size_t,
    const SubseqHash3PivotTables *,const size_t *,size_t,SubseqHash3EditResult *);
extern int baseline_compute_pivot_repeats(const uint8_t *,size_t,size_t,
    const SubseqHash3PivotTables *,size_t,size_t,SubseqHash3EditResult *);

static uint64_t random_state=314159265;
static unsigned random_base(void) {
    random_state^=random_state<<13; random_state^=random_state>>7;
    random_state^=random_state<<17;
    return (unsigned)(random_state&3);
}
static double seconds(void) {
    struct timespec t;
    assert(!clock_gettime(CLOCK_MONOTONIC,&t));
    return t.tv_sec+t.tv_nsec*1e-9;
}
static void same(const SubseqHash3EditResult *a,const SubseqHash3EditResult *b) {
    assert(a->length==b->length && a->psi==b->psi && a->omega==b->omega);
    assert(a->deletions==b->deletions && a->substitutions==b->substitutions);
    assert(!memcmp(a->seed,b->seed,a->length));
    assert(!memcmp(a->positions,b->positions,a->length*sizeof(*a->positions)));
    assert(!memcmp(a->substituted,b->substituted,a->length));
}
static void check_ties(void) {
    const size_t pivots[]={255,0,128},budgets[]={0,1,256};
    size_t d,mode,budget;
    for (d=31;d<=32;++d) for (mode=0;mode<3;++mode) {
        const size_t n=mode==2?9:257,k=n-1,shift=mode==2?1:0;
        const int sign=mode==1?-1:1;
        uint8_t s[257]={0};
        SubseqHash3PivotTables t={0};
        size_t i;
        assert(!subseqhash3_pivot_tables_init(&t,n,d,888));
        for (i=0;i<k*4;++i) {
            t.left.C[i]=t.right.C[i]=t.pivot_C[i]=shift;
            t.pivot_left_sign[i]=t.pivot_right_sign[i]=1;
            t.pivot_A[i]=(int8_t)(sign*127);
        }
        memset(t.left.B1_masks,255,k*4*t.left.blocks*sizeof(uint32_t));
        memset(t.right.B1_masks,255,k*4*t.right.blocks*sizeof(uint32_t));
        memset(t.right_reversed_B1_masks,255,k*4*t.right.blocks*sizeof(uint32_t));
        for (i=0;i<k*4*d;++i) {
            t.left.A[i]=t.right.A[i]=127;
            t.left.B1[i]=t.right.B1[i]=1;
            t.left.B2[i]=t.right.B2[i]=(int8_t)sign;
            t.left.signed_A[i]=t.right.signed_A[i]=(int8_t)(sign*127);
            t.right_reversed_A[i]=(int8_t)(sign*127);
        }
        for (budget=0;budget<(mode==2?1:3);++budget) {
            const size_t e=mode==2?3:budgets[budget];
            const size_t *selected=mode==2?NULL:pivots;
            SubseqHash3EditResult a[3]={{0}},b[3]={{0}};
            assert(!subseqhash3_compute_pivot_seeds(s,n,e,&t,selected,3,a));
            assert(!baseline_compute_pivot_seeds(s,n,e,&t,selected,3,b));
            for (i=0;i<3;++i) {
                same(a+i,b+i);
                assert(a[i].psi==(k*shift)%d);
                assert(a[i].omega==sign*127.0*k);
            }
            subseqhash3_edit_results_destroy(a,3);
            subseqhash3_edit_results_destroy(b,3);
        }
        subseqhash3_pivot_tables_destroy(&t);
    }
}
static void check(void) {
    const size_t residues[]={1,3,7,30,31,32,33,63,64,65,257,1024,1025};
    size_t trial;
    for (trial=0;trial<260;++trial) {
        const size_t n=2+trial%49,k=n-1,e=(trial/3)%(k+1),r=1+trial%k;
        uint8_t s[50];
        SubseqHash3PivotTables t={0};
        SubseqHash3EditResult a[49]={{0}},b[49]={{0}};
        size_t i;
        size_t pivots[49];
        for (i=0;i<n;++i) s[i]=(uint8_t)(trial%5?random_base():trial%4);
        for (i=0;i<r;++i) pivots[i]=(17*i+7*trial)%k;
        assert(!subseqhash3_pivot_tables_init(&t,n,residues[trial%13],trial+678));
        assert(!subseqhash3_compute_pivot_seeds(s,n,e,&t,trial&1?pivots:NULL,r,a));
        assert(!baseline_compute_pivot_seeds(s,n,e,&t,trial&1?pivots:NULL,r,b));
        for (i=0;i<r;++i) same(a+i,b+i);
        subseqhash3_edit_results_destroy(a,r);
        subseqhash3_edit_results_destroy(b,r);
        subseqhash3_pivot_tables_destroy(&t);
    }
    for (trial=0;trial<3;++trial) {
        const size_t budgets[]={0,128,256},pivots[]={0,127,255};
        uint8_t s[257];
        SubseqHash3PivotTables t={0};
        SubseqHash3EditResult a[3]={{0}},b[3]={{0}};
        size_t i;
        for (i=0;i<257;++i) s[i]=(uint8_t)random_base();
        assert(!subseqhash3_pivot_tables_init(&t,257,65,trial+789));
        assert(!subseqhash3_compute_pivot_seeds(s,257,budgets[trial],&t,pivots,3,a));
        assert(!baseline_compute_pivot_seeds(s,257,budgets[trial],&t,pivots,3,b));
        for (i=0;i<3;++i) same(a+i,b+i);
        subseqhash3_edit_results_destroy(a,3);
        subseqhash3_edit_results_destroy(b,3);
        subseqhash3_pivot_tables_destroy(&t);
    }
    check_ties();
    puts("277 exact baseline comparisons passed");
}
static double measure(const uint8_t *s,size_t n,size_t e,size_t r,
    const SubseqHash3PivotTables *t,int old) {
    SubseqHash3EditResult out[49]={{0}};
    size_t i;
    double start=0;
    for (i=0;i<=128;++i) {
        const uint8_t *input=s+((i?i-1:0)%32)*n;
        if (i==1) start=seconds();
        if (old) assert(!baseline_compute_pivot_repeats(input,n,e,t,1,r,out));
        else assert(!subseqhash3_compute_pivot_repeats(input,n,e,t,1,r,out));
        subseqhash3_edit_results_destroy(out,r);
    }
    return (seconds()-start)*1e6/128;
}
static double median(double *a) {
    size_t i;
    for (i=1;i<7;++i) {
        const double value=a[i];
        size_t j=i;
        while (j && a[j-1]>value) { a[j]=a[j-1]; --j; }
        a[j]=value;
    }
    return a[3];
}
int main(int argc,char **argv) {
    const size_t cases[][3]={{30,5,31},{50,8,31},{30,5,1024},{50,8,1024},
        {30,5,65},{50,8,65},{30,5,257},{50,8,257},{30,5,512},{50,8,512}};
    const size_t count=argc==2 && !strcmp(argv[1],"--sweep")?10:4;
    size_t c;
    printf("# baseline=%s, threads=1\n",PIVOT_BASELINE_LABEL);
    check();
    if (argc==2 && !strcmp(argv[1],"--check-only")) return 0;
    puts("n,E,D,r,baseline_us,current_us,further_speedup");
    for (c=0;c<count;++c) {
        const size_t n=cases[c][0],e=cases[c][1],d=cases[c][2],r=10;
        uint8_t s[32*50];
        SubseqHash3PivotTables t={0};
        double old[7],current[7],a,b;
        size_t i,round;
        for (i=0;i<32*n;++i) s[i]=(uint8_t)random_base();
        assert(!subseqhash3_pivot_tables_init(&t,n,d,1234));
        for (i=0;i<32;++i) {
            SubseqHash3EditResult a[10]={{0}},b[10]={{0}};
            size_t p;
            assert(!subseqhash3_compute_pivot_repeats(s+i*n,n,e,&t,1,r,a));
            assert(!baseline_compute_pivot_repeats(s+i*n,n,e,&t,1,r,b));
            for (p=0;p<r;++p) same(a+p,b+p);
            subseqhash3_edit_results_destroy(a,r);
            subseqhash3_edit_results_destroy(b,r);
        }
        for (round=0;round<7;++round) {
            if (round&1) {
                current[round]=measure(s,n,e,r,&t,0);
                old[round]=measure(s,n,e,r,&t,1);
            } else {
                old[round]=measure(s,n,e,r,&t,1);
                current[round]=measure(s,n,e,r,&t,0);
            }
            printf("# sample,%zu,%zu,%zu,10,%zu,%.6f,%.6f\n",
                n,e,d,round+1,old[round],current[round]);
        }
        a=median(old); b=median(current);
        printf("%zu,%zu,%zu,10,%.3f,%.3f,%.3f\n",n,e,d,a,b,a/b);
        fflush(stdout);
        subseqhash3_pivot_tables_destroy(&t);
    }
    return 0;
}
