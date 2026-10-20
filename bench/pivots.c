#define _POSIX_C_SOURCE 200809L
#include "ssh3_pivots.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static uint64_t random_state=UINT64_C(123456789);
static unsigned random_base(void) {
    random_state^=random_state<<13; random_state^=random_state>>7;
    random_state^=random_state<<17;
    return (unsigned)(random_state&3u);
}

/* Literal evaluation of the written order definition, independent of DP
   cost transforms, mirrored residues, joins, and traceback. */
static void evaluate(const uint8_t *seed,size_t k,size_t pivot,
    const SubseqHash3PivotTables *t,size_t *psi,int *omega) {
    size_t i,lr=0,rr=0;
    int lo=0,ro=0;
    const size_t pt=4*pivot+seed[pivot];
    for (i=0;i<pivot;++i) {
        const unsigned a=seed[i];
        size_t at;
        lr=(lr+t->left.C[4*i+a])%t->d;
        at=(4*i+a)*t->d+lr;
        lo=t->left.B1[at]*lo+t->left.A[at]*t->left.B2[at];
    }
    for (i=0;i<k-pivot-1;++i) {
        const unsigned a=seed[k-1-i];
        size_t at;
        rr=(rr+t->right.C[4*i+a])%t->d;
        at=(4*i+a)*t->d+rr;
        ro=t->right.B1[at]*ro+t->right.A[at]*t->right.B2[at];
    }
    *psi=(lr+rr+t->pivot_C[pt])%t->d;
    *omega=t->pivot_left_sign[pt]*lo+t->pivot_right_sign[pt]*ro+t->pivot_A[pt];
}

static void validate(const uint8_t *s,size_t n,size_t edits,size_t pivot,
    const SubseqHash3PivotTables *t,const SubseqHash3EditResult *out) {
    size_t i,changed=0,psi;
    int omega;
    assert(out->length==n-1 && out->deletions==1 && out->substitutions==edits);
    for (i=0;i<out->length;++i) {
        assert(out->positions[i]<n && (!i || out->positions[i]>out->positions[i-1]));
        assert(out->seed[i]<4);
        assert(out->substituted[i]==(out->seed[i]!=s[out->positions[i]]));
        changed+=out->substituted[i];
    }
    assert(changed==edits);
    evaluate(out->seed,n-1,pivot,t,&psi,&omega);
    assert(psi==out->psi && (double)omega==out->omega);
}

typedef struct {
    const uint8_t *s;
    size_t n,edits,deleted;
    const SubseqHash3PivotTables *tables;
    uint8_t seed[8];
    size_t psi[8];
    int omega[8];
} Brute;

static void enumerate(Brute *b,size_t row,size_t changed) {
    unsigned a;
    const size_t k=b->n-1;
    if (row==k) {
        size_t p;
        if (changed!=b->edits) return;
        for (p=0;p<k;++p) {
            size_t psi;
            int omega;
            evaluate(b->seed,k,p,b->tables,&psi,&omega);
            if (psi<b->psi[p] || (psi==b->psi[p] && omega>b->omega[p])) {
                b->psi[p]=psi; b->omega[p]=omega;
            }
        }
        return;
    }
    for (a=0;a<4;++a) {
        const size_t next=changed+(a!=b->s[row+(row>=b->deleted)]);
        if (next>b->edits || next+k-row-1<b->edits) continue;
        b->seed[row]=(uint8_t)a;
        enumerate(b,row+1,next);
    }
}

static void same(const SubseqHash3EditResult *a,const SubseqHash3EditResult *b) {
    assert(a->length==b->length && a->psi==b->psi && a->omega==b->omega);
    assert(!memcmp(a->seed,b->seed,a->length));
    assert(!memcmp(a->positions,b->positions,a->length*sizeof(*a->positions)));
    assert(!memcmp(a->substituted,b->substituted,a->length));
}

static void correctness(void) {
    const size_t ds[]={1,3,7,31,32,33,64,257,1024,1025};
    SubseqHash3PivotWorkspace *workspace=subseqhash3_pivot_workspace_create();
    size_t trial,p;
    assert(workspace);
    for (trial=0;trial<140;++trial) {
        const size_t n=2+trial%6,k=n-1,edits=(trial/3)%(k+1);
        const size_t d=ds[trial%10];
        SubseqHash3PivotTables tables={0};
        SubseqHash3EditResult out[6]={0},reused[6]={0};
        uint8_t s[7];
        Brute b={0};
        size_t i;
        for (i=0;i<n;++i) s[i]=(uint8_t)(trial%4?random_base():0);
        assert(!subseqhash3_pivot_tables_init(&tables,n,d,12345+trial));
        assert(!subseqhash3_compute_pivot_seeds(s,n,edits,&tables,NULL,k,out));
        assert(!subseqhash3_compute_pivot_repeats_with_workspace(
            s,n,edits,&tables,1,k,reused,workspace));
        for (p=0;p<k;++p) same(out+p,reused+p);
        subseqhash3_edit_results_destroy(reused,k);
        b.s=s; b.n=n; b.edits=edits; b.tables=&tables;
        for (p=0;p<k;++p) { b.psi[p]=d; b.omega[p]=INT_MIN; }
        for (i=0;i<n;++i) { b.deleted=i; enumerate(&b,0,0); }
        for (p=0;p<k;++p) {
            SubseqHash3EditResult single={0};
            validate(s,n,edits,p,&tables,out+p);
            assert(out[p].psi==b.psi[p] && out[p].omega==b.omega[p]);
            assert(!subseqhash3_compute_pivot_seeds(s,n,edits,&tables,&p,1,&single));
            same(out+p,&single);
            subseqhash3_edit_results_destroy(&single,1);
        }
        subseqhash3_edit_results_destroy(out,k);
        subseqhash3_pivot_tables_destroy(&tables);
    }
    for (trial=0;trial<30;++trial) {
        const size_t n=8+trial%25,k=n-1,edits=trial%(k+1);
        SubseqHash3PivotTables tables={0};
        SubseqHash3EditResult out[32]={0};
        uint8_t s[33];
        size_t i;
        for (i=0;i<n;++i) s[i]=(uint8_t)random_base();
        assert(!subseqhash3_pivot_tables_init(&tables,n,ds[trial%10],54321+trial));
        assert(!subseqhash3_compute_pivot_seeds(s,n,edits,&tables,NULL,k,out));
        for (p=0;p<k;++p) validate(s,n,edits,p,&tables,out+p);
        subseqhash3_edit_results_destroy(out,k);
        subseqhash3_pivot_tables_destroy(&tables);
    }
    {
        const uint8_t s[]={0,0,1,2,3};
        SubseqHash3PivotTables tables[3]={{0}};
        SubseqHash3EditResult out[10]={{0}},single={0};
        size_t i;
        for (i=0;i<3;++i)
            assert(!subseqhash3_pivot_tables_init(tables+i,5,31,100+i));
        assert(!subseqhash3_compute_pivot_repeats(s,5,2,tables,3,10,out));
        for (i=0;i<10;++i) {
            size_t pivot=i%4;
            assert(!subseqhash3_compute_pivot_seeds(s,5,2,tables+i/4,&pivot,1,&single));
            same(out+i,&single);
            subseqhash3_edit_results_destroy(&single,1);
        }
        subseqhash3_edit_results_destroy(out,10);
        assert(!subseqhash3_compute_pivot_repeats_with_workspace(
            s,5,2,tables,3,10,out,workspace));
        for (i=0;i<10;++i) validate(s,5,2,i%4,tables+i/4,out+i);
        subseqhash3_edit_results_destroy(out,10);
        assert(subseqhash3_compute_pivot_repeats(s,5,2,tables,2,10,out)==-1);
        assert(subseqhash3_compute_pivot_seeds(s,5,5,tables,NULL,1,out)==-1);
        assert(subseqhash3_compute_pivot_seeds(s,1,0,tables,NULL,1,out)==-1);
        assert(subseqhash3_compute_pivot_seeds(s,5,2,tables,NULL,5,out)==-1);
        i=4;
        assert(subseqhash3_compute_pivot_seeds(s,5,2,tables,&i,1,out)==-1);
        {
            const uint8_t bad[]={0,4,1,2,3};
            assert(subseqhash3_compute_pivot_seeds(bad,5,2,tables,NULL,1,out)==-1);
            assert(subseqhash3_compute_pivot_seeds(NULL,5,2,tables,NULL,1,out)==-1);
            assert(subseqhash3_compute_pivot_seeds(s,5,2,NULL,NULL,1,out)==-1);
            assert(subseqhash3_compute_pivot_seeds(s,5,2,tables,NULL,1,NULL)==-1);
            tables[1].k=1;
            assert(subseqhash3_compute_pivot_repeats_with_workspace(
                s,5,2,tables,3,10,out,workspace)==-1);
            for (i=0;i<10;++i) assert(!out[i].seed && !out[i].positions);
            tables[1].k=4;
        }
        for (i=0;i<3;++i) subseqhash3_pivot_tables_destroy(tables+i);
    }
    {
        uint8_t s[257]={0};
        SubseqHash3PivotTables t={0};
        SubseqHash3EditResult out[3]={{0}};
        const size_t pivots[]={255,0,128};
        size_t i,pass;
        assert(!subseqhash3_pivot_tables_init(&t,257,33,888));
        for (i=0;i<t.k*4;++i) {
            t.left.C[i]=t.right.C[i]=t.pivot_C[i]=0;
            t.pivot_left_sign[i]=t.pivot_right_sign[i]=1;
        }
        memset(t.left.B1_masks,255,t.k*4*t.left.blocks*sizeof(uint32_t));
        memset(t.right.B1_masks,255,t.k*4*t.right.blocks*sizeof(uint32_t));
        memset(t.right_reversed_B1_masks,255,t.k*4*t.right.blocks*sizeof(uint32_t));
        for (pass=0;pass<2;++pass) {
            const int sign=pass?-1:1;
            for (i=0;i<t.k*4*33;++i) {
                t.left.A[i]=t.right.A[i]=127;
                t.left.B1[i]=t.right.B1[i]=1;
                t.left.B2[i]=t.right.B2[i]=(int8_t)sign;
                t.left.signed_A[i]=t.right.signed_A[i]=(int8_t)(sign*127);
                t.right_reversed_A[i]=(int8_t)(sign*127);
            }
            memset(t.pivot_A,sign*127,t.k*4);
            assert(!subseqhash3_compute_pivot_seeds(s,257,1,&t,pivots,3,out));
            for (i=0;i<3;++i) {
                validate(s,257,1,pivots[i],&t,out+i);
                assert(out[i].psi==0 && out[i].omega==sign*32512.0);
            }
            subseqhash3_edit_results_destroy(out,3);
            assert(!subseqhash3_compute_pivot_repeats_with_workspace(
                s,5,4,&t,1,3,out,workspace));
            for (i=0;i<3;++i) validate(s,5,4,i,&t,out+i);
            subseqhash3_edit_results_destroy(out,3);
        }
        subseqhash3_pivot_tables_destroy(&t);
    }
    subseqhash3_pivot_workspace_destroy(workspace);
    puts("pivot checks: 140 exhaustive neighborhoods, 30 larger witnesses, score limits, workspace and repeat batching passed");
}

static double seconds(void) {
    struct timespec ts;
    assert(!clock_gettime(CLOCK_MONOTONIC,&ts));
    return ts.tv_sec+ts.tv_nsec*1e-9;
}

static double measure(const uint8_t *inputs,size_t n,size_t edits,size_t r,
    const SubseqHash3PivotTables *tables,size_t families,
    const SubseqHash3Tables *old_tables,int method,size_t input_count,
    size_t loops) {
    SubseqHash3EditResult *out=calloc(r,sizeof(*out));
    SubseqHash3PivotWorkspace *workspace=NULL;
    double start;
    size_t i,j;
    assert(out);
    if (method==3) {
        workspace=subseqhash3_pivot_workspace_create();
        assert(workspace);
    }
    start=seconds();
    /* One untimed warmup per method/sample, followed by identical inputs. */
    for (i=0;i<=input_count*loops;++i) {
        const uint8_t *s=inputs+((i?i-1:0)%input_count)*n;
        if (i==1) start=seconds();
        if (method==3) {
            assert(!subseqhash3_compute_pivot_repeats_with_workspace(
                s,n,edits,tables,families,r,out,workspace));
        } else if (!method) {
            assert(!subseqhash3_compute_pivot_repeats(s,n,edits,tables,families,r,out));
        } else for (j=0;j<r;++j) {
            if (method==1) {
                const size_t pivot=j%(n-1);
                assert(!subseqhash3_compute_pivot_seeds(s,n,edits,tables+j/(n-1),&pivot,1,out+j));
            } else {
                /* Earlier API: ONE substitution and E deletions. Different
                   neighborhood; this is only the requested timing baseline. */
                assert(!subseqhash3_compute_one_substitution(s,n,edits,old_tables+j,out+j));
            }
        }
        subseqhash3_edit_results_destroy(out,r);
    }
    start=(seconds()-start)/(input_count*loops);
    subseqhash3_pivot_workspace_destroy(workspace);
    free(out);
    return start;
}

static double median(double *a,size_t count) {
    size_t i;
    for (i=1;i<count;++i) {
        const double value=a[i];
        size_t j=i;
        while (j && a[j-1]>value) { a[j]=a[j-1]; --j; }
        a[j]=value;
    }
    return a[count/2];
}

int main(int argc,char **argv) {
    const size_t default_cases[][3]={{20,3,31},{30,5,31},{30,5,1024},{50,8,1024}};
    const size_t r10_cases[][3]={{30,5,31},{50,8,31},{30,5,1024},{50,8,1024}};
    const int focused=argc==2 && !strcmp(argv[1],"--r10");
    const size_t (*cases)[3]=focused?r10_cases:default_cases;
    const size_t input_count=focused?32:8,loops=focused?4:2,rounds=focused?7:3;
    size_t c;
    if (argc>2 || (argc==2 && !focused && strcmp(argv[1],"--check-only"))) {
        fputs("usage: bench-pivots [--check-only | --r10]\n",stderr);
        return EXIT_FAILURE;
    }
    correctness();
    if (argc==2 && !strcmp(argv[1],"--check-only")) return 0;
    printf("# inputs=%zu, loops=%zu, samples=%zu, untimed_warmup=1\n",
           input_count,loops,rounds);
    puts("n,E,D,r,shared_us,reuse_us,separate_pivots_us,previous_api_us,speedup_separate,speedup_previous,speedup_previous_reuse");
    for (c=0;c<4;++c) {
        const size_t n=cases[c][0],edits=cases[c][1],d=cases[c][2];
        const size_t repeats[]={1,4,8,n-1,2*(n-1)+3};
        uint8_t inputs[32*50];
        size_t i,ri;
        for (i=0;i<input_count*n;++i) inputs[i]=(uint8_t)random_base();
        for (ri=0;ri<(focused?1:sizeof(repeats)/sizeof(*repeats));++ri) {
            const size_t r=focused?10:repeats[ri],families=1+(r-1)/(n-1);
            SubseqHash3PivotTables *tables=calloc(families,sizeof(*tables));
            SubseqHash3Tables *old_tables=calloc(r,sizeof(*old_tables));
            double measurements[4][7],shared,reuse,separate,old;
            size_t round,method;
            assert(tables && old_tables);
            for (i=0;i<families;++i)
                assert(!subseqhash3_pivot_tables_init(tables+i,n,d,1234+i));
            for (i=0;i<r;++i)
                assert(!subseqhash3_tables_init(old_tables+i,n,d,9876+i));
            for (round=0;round<rounds;++round) for (i=0;i<4;++i) {
                method=(round+i)%4;
                measurements[method][round]=measure(inputs,n,edits,r,
                    tables,families,old_tables,(int)method,input_count,loops);
                if (focused) {
                    static const char *names[]={"shared","separate","simple","reuse"};
                    printf("# sample,%zu,%zu,%zu,%zu,%zu,%s,%.6f\n",
                        n,edits,d,r,round+1,names[method],
                        measurements[method][round]*1e6);
                }
            }
            shared=median(measurements[0],rounds);
            separate=median(measurements[1],rounds);
            old=median(measurements[2],rounds); reuse=median(measurements[3],rounds);
            printf("%zu,%zu,%zu,%zu,%.3f,%.3f,%.3f,%.3f,%.2fx,%.2fx,%.2fx\n",
                n,edits,d,r,shared*1e6,reuse*1e6,separate*1e6,old*1e6,
                separate/shared,old/shared,old/reuse);
            fflush(stdout);
            for (i=0;i<families;++i) subseqhash3_pivot_tables_destroy(tables+i);
            for (i=0;i<r;++i) subseqhash3_tables_destroy(old_tables+i);
            free(tables); free(old_tables);
        }
    }
    return 0;
}
