#define _POSIX_C_SOURCE 200809L
#include "ssh3_deletions.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int baseline_compute_deletion_seeds_with_workspace(const uint8_t *,size_t,size_t,
    const SubseqHash3DeletionTables *,const size_t *,size_t,
    SubseqHash3DeletionResult *,SubseqHash3DeletionWorkspace *);
SubseqHash3DeletionWorkspace *baseline_deletion_workspace_create(void);
void baseline_deletion_workspace_destroy(SubseqHash3DeletionWorkspace *);
static uint64_t rng=12345,checksum;
static double now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t);
    return t.tv_sec+t.tv_nsec*1e-9;
}
static int compare(const void *a,const void *b) {
    const double x=*(const double *)a,y=*(const double *)b;
    return (x>y)-(x<y);
}
static void run(size_t n,size_t e,size_t d,size_t count,int check) {
    SubseqHash3PivotTables shared={0};
    SubseqHash3DeletionTables tables={0};
    SubseqHash3DeletionWorkspace *w=subseqhash3_deletion_workspace_create();
    SubseqHash3DeletionWorkspace *old=baseline_deletion_workspace_create();
    SubseqHash3PivotWorkspace *pw=subseqhash3_pivot_workspace_create();
    SubseqHash3DeletionResult a[256]={{0}},b[256]={{0}};
    SubseqHash3EditResult r[256]={{0}};
    uint8_t inputs[16][256];
    double times[3][5];
    assert(w && old && pw);
    assert(!subseqhash3_pivot_tables_init(&shared,n+1,d,42));
    assert(!subseqhash3_deletion_tables_from_pivots(&tables,&shared));
    for (size_t q=0;q<16;++q) {
        for (size_t i=0;i<n;++i) {
            rng^=rng<<13; rng^=rng>>7; rng^=rng<<17;
            inputs[q][i]=(uint8_t)(rng&3);
        }
        assert(!baseline_compute_deletion_seeds_with_workspace(inputs[q],n,e,&tables,NULL,count,a,old));
        assert(!subseqhash3_compute_deletion_seeds_with_workspace(inputs[q],n,e,&tables,NULL,count,b,w));
        for (size_t p=0;p<count;++p) {
            assert(a[p].psi==b[p].psi && a[p].omega==b[p].omega);
            assert(!memcmp(a[p].seed,b[p].seed,n-e));
            assert(!memcmp(a[p].positions,b[p].positions,(n-e)*sizeof(size_t)));
        }
        subseqhash3_deletion_results_destroy(a,count);
        subseqhash3_deletion_results_destroy(b,count);
    }
    if (!check) {
        for (size_t round=0;round<5;++round) for (size_t step=0;step<3;++step) {
            const size_t mode=(step+round)%3;
            double start=now();
            for (size_t batch=0;batch<3;++batch) for (size_t q=0;q<16;++q) {
                if (mode==2) {
                    assert(!subseqhash3_compute_pivot_repeats_with_workspace(inputs[q],n,e,&shared,1,count,r,pw));
                    checksum+=r[0].psi+(uint64_t)r[0].omega;
                    subseqhash3_edit_results_destroy(r,count);
                } else {
                    assert(!(mode?subseqhash3_compute_deletion_seeds_with_workspace:
                        baseline_compute_deletion_seeds_with_workspace)(inputs[q],n,e,&tables,NULL,count,b,mode?w:old));
                    checksum+=b[0].psi+(uint64_t)b[0].omega;
                    subseqhash3_deletion_results_destroy(b,count);
                }
            }
            times[mode][round]=(now()-start)*1e6/48;
        }
        for (size_t mode=0;mode<3;++mode) qsort(times[mode],5,sizeof(double),compare);
        printf("%zu,%zu,%zu,%zu,%.3f,%.3f,%.2f,%.3f,%.2f\n",n,e,d,count,
            times[0][2],times[1][2],times[0][2]/times[1][2],times[2][2],times[2][2]/times[1][2]);
        fflush(stdout);
    }
    subseqhash3_pivot_workspace_destroy(pw);
    baseline_deletion_workspace_destroy(old);
    subseqhash3_deletion_workspace_destroy(w);
    subseqhash3_deletion_tables_destroy(&tables);
    subseqhash3_pivot_tables_destroy(&shared);
}
int main(int argc,char **argv) {
    const int check=argc==2 && !strcmp(argv[1],"--check-only");
    if (!check) puts("n,edits,D,seeds,baseline_us,deletions_us,speedup,repeats_us,repeats_over_deletions");
    for (size_t d=31;d<=1024;d=d==31?1024:1025) {
        run(30,4,d,10,check); run(30,4,d,26,check);
        run(50,8,d,10,check); run(50,8,d,42,check);
    }
    run(20,1,33,19,check); run(20,18,65,2,check);
    fprintf(stderr,"160 windows match baseline scores, seeds and positions; checksum=%llu\n",(unsigned long long)checksum);
}
