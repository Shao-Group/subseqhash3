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
static uint64_t checksum;
static size_t verified_windows;
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
/* Each implementation owns its tables, workspace and results. */
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
SubseqHash3PivotWorkspace *baseline_subseqhash3_pivot_workspace_create(void);
void baseline_subseqhash3_pivot_workspace_destroy(SubseqHash3PivotWorkspace *workspace);
int baseline_subseqhash3_pivot_tables_init(SubseqHash3PivotTables *tables,
    size_t maximum_input_length,size_t residues,uint64_t seed);
void baseline_subseqhash3_pivot_tables_destroy(SubseqHash3PivotTables *tables);
int baseline_subseqhash3_compute_pivot_seeds(const uint8_t *sequence,size_t n,
    size_t substitutions,const SubseqHash3PivotTables *tables,
    const size_t *pivots,size_t count,SubseqHash3EditResult *results);
int baseline_subseqhash3_compute_pivot_repeats(const uint8_t *sequence,size_t n,
    size_t substitutions,const SubseqHash3PivotTables *tables,
    size_t table_count,size_t repeats,SubseqHash3EditResult *results);
int baseline_subseqhash3_compute_pivot_repeats_with_workspace(const uint8_t *sequence,
    size_t n,size_t substitutions,const SubseqHash3PivotTables *tables,
    size_t table_count,size_t repeats,SubseqHash3EditResult *results,
    SubseqHash3PivotWorkspace *workspace);
int baseline_subseqhash3_compute_outward_pivot_seeds(const uint8_t *,size_t n,
    size_t substitutions,const SubseqHash3PivotTables *,const size_t *pivots,
    size_t count,SubseqHash3EditResult *);
int baseline_subseqhash3_compute_outward_pivot_repeats(const uint8_t *,size_t n,
    size_t substitutions,const SubseqHash3PivotTables *,size_t table_count,
    size_t repeats,SubseqHash3EditResult *);
int baseline_subseqhash3_compute_outward_pivot_repeats_with_workspace(const uint8_t *,size_t n,
    size_t substitutions,const SubseqHash3PivotTables *,size_t table_count,
    size_t repeats,SubseqHash3EditResult *,SubseqHash3PivotWorkspace *);
int baseline_subseqhash3_deletion_tables_init(SubseqHash3DeletionTables *,size_t k,
    size_t residues,uint64_t seed);
int baseline_subseqhash3_deletion_tables_from_pivots(SubseqHash3DeletionTables *,
    const SubseqHash3PivotTables *);
int baseline_subseqhash3_deletion_tables_import(SubseqHash3DeletionTables *,
    const char *path,size_t k,size_t residues);
void baseline_subseqhash3_deletion_tables_destroy(SubseqHash3DeletionTables *);
SubseqHash3DeletionWorkspace *baseline_subseqhash3_deletion_workspace_create(void);
void baseline_subseqhash3_deletion_workspace_destroy(SubseqHash3DeletionWorkspace *);
void baseline_subseqhash3_deletion_results_destroy(SubseqHash3DeletionResult *,size_t count);
int baseline_subseqhash3_compute_deletion_seeds_with_workspace(const uint8_t *,size_t n,
    size_t deletions,const SubseqHash3DeletionTables *,const size_t *pivots,
    size_t count,SubseqHash3DeletionResult *,SubseqHash3DeletionWorkspace *);
int baseline_subseqhash3_compute_deletion_seeds(const uint8_t *,size_t n,size_t deletions,
    const SubseqHash3DeletionTables *,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *);
int baseline_subseqhash3_compute_inward_deletion_seeds_with_workspace(const uint8_t *,size_t n,
    size_t deletions,const SubseqHash3DeletionTables *,const size_t *pivots,
    size_t count,SubseqHash3DeletionResult *,SubseqHash3DeletionWorkspace *);
int baseline_subseqhash3_compute_inward_deletion_seeds(const uint8_t *,size_t n,size_t deletions,
    const SubseqHash3DeletionTables *,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *);
static void baseline_init(Engine *a,int repeats,int inward,size_t n,size_t e,size_t r,size_t d,uint64_t seed) {
    memset(a,0,sizeof(*a));
    a->n=n; a->e=e; a->r=r; a->repeats=repeats; a->inward=inward;
    if (repeats) {
        assert(!baseline_subseqhash3_pivot_tables_init(&a->pt,n,d,seed));
        a->pw=baseline_subseqhash3_pivot_workspace_create(); assert(a->pw);
    } else {
        assert(!baseline_subseqhash3_deletion_tables_init(&a->dt,n-e,d,seed));
        a->dw=baseline_subseqhash3_deletion_workspace_create(); assert(a->dw);
    }
}
static void baseline_compute(Engine *a,const uint8_t *x) {
    if (a->repeats) {
        if (a->inward) assert(!baseline_subseqhash3_compute_pivot_repeats_with_workspace(x,a->n,a->e,&a->pt,1,a->r,a->pr,a->pw));
        else assert(!baseline_subseqhash3_compute_outward_pivot_repeats_with_workspace(x,a->n,a->e,&a->pt,1,a->r,a->pr,a->pw));
    } else {
        if (a->inward) assert(!baseline_subseqhash3_compute_inward_deletion_seeds_with_workspace(x,a->n,a->e,&a->dt,NULL,a->r,a->dr,a->dw));
        else assert(!baseline_subseqhash3_compute_deletion_seeds_with_workspace(x,a->n,a->e,&a->dt,NULL,a->r,a->dr,a->dw));
    }
}

static void baseline_clear(Engine *a) {
    if (a->repeats) subseqhash3_edit_results_destroy(a->pr,a->r);
    else baseline_subseqhash3_deletion_results_destroy(a->dr,a->r);
}
static void baseline_destroy(Engine *a) {
    if (a->repeats) {
        baseline_subseqhash3_pivot_workspace_destroy(a->pw); baseline_subseqhash3_pivot_tables_destroy(&a->pt);
    } else {
        baseline_subseqhash3_deletion_workspace_destroy(a->dw); baseline_subseqhash3_deletion_tables_destroy(&a->dt);
    }
}

static void compare_run(int rep,int inward,size_t n,size_t e,size_t r,size_t d,int check) {
    Engine a[2]; uint8_t x[16][256]; uint64_t rng=812733; double times[2][7];
    baseline_init(&a[0],rep,inward,n,e,r,d,43); init(&a[1],rep,inward,n,e,r,d,43);
    for (size_t q=0;q<19;++q) {
        uint8_t structured[256];
        const uint8_t *input;
        if (q<16) {
            for (size_t j=0;j<n;++j) x[q][j]=(uint8_t)(random64(&rng)&3);
            input=x[q];
        } else {
            for (size_t j=0;j<n;++j) structured[j]=(uint8_t)(q==16?0:q==17?3:j%4);
            input=structured;
        }
        baseline_compute(&a[0],input); compute(&a[1],input);
        for(size_t j=0;j<r;++j) {
            assert(!memcmp(seed(&a[0],j),seed(&a[1],j),rep?n-1:n-e));
            if(rep) {
                assert(a[0].pr[j].psi==a[1].pr[j].psi && a[0].pr[j].omega==a[1].pr[j].omega);
                assert(!memcmp(a[0].pr[j].positions,a[1].pr[j].positions,(n-1)*sizeof(size_t)));
                assert(!memcmp(a[0].pr[j].substituted,a[1].pr[j].substituted,n-1));
            } else {
                assert(a[0].dr[j].psi==a[1].dr[j].psi && a[0].dr[j].omega==a[1].dr[j].omega);
                assert(!memcmp(a[0].dr[j].positions,a[1].dr[j].positions,(n-e)*sizeof(size_t)));
            }
        }
        baseline_clear(&a[0]); clear(&a[1]); ++verified_windows;
    }
    if(!check) {
        for (size_t round=0;round<7;++round) for (size_t step=0;step<2;++step) {
            size_t m=(round+step)%2; double start=now();
            for (size_t pass=0;pass<3;++pass) for (size_t q=0;q<16;++q) {
                if(m) compute(&a[m],x[q]); else baseline_compute(&a[m],x[q]);
                checksum+=seed(&a[m],0)[0];
                if(m) clear(&a[m]); else baseline_clear(&a[m]);
            }
            times[m][round]=(now()-start)*1e6/48;
        }
        for(int m=0;m<2;++m)qsort(times[m],7,sizeof(double),compare);
        printf("%s,%s,%zu,%zu,%zu,%zu,%.3f,%.3f,%.3f\n",rep?"repeats":"deletions",inward?"inward":"outward",n,e,d,r,times[0][3],times[1][3],times[0][3]/times[1][3]);
        fflush(stdout);
    }
    baseline_destroy(&a[0]); destroy(&a[1]);
}
int main(int argc,char **argv) {
    if (argc>3 || (argc>1 && strcmp(argv[1],"--check-only") &&
        strcmp(argv[1],"--single") && strcmp(argv[1],"--residues")) ||
        (argc>1 && !strcmp(argv[1],"--residues") && argc!=3)) {
        fprintf(stderr,"Usage: %s [--check-only | --residues D | --single [D]]\n",argv[0]);
        return 1;
    }
    int check=argc>1 && !strcmp(argv[1],"--check-only");
    size_t residues=1024;
    if(argc==3) {
        char *end;
        unsigned long long value=strtoull(argv[2],&end,10);
        assert(*argv[2] && *argv[2]!='-' && !*end && value && value<=SIZE_MAX);
        residues=(size_t)value;
    }
    puts("algorithm,traversal,n,E,D,seeds,baseline_us,optimized_us,speedup");
    if(argc>1 && !strcmp(argv[1],"--single")) {
        for(int rep=0;rep<2;++rep) for(int inward=0;inward<2;++inward)
            for(size_t n=30;n<=50;n+=20)for(int all=0;all<2;++all)
                compare_run(rep,inward,n,rep?0:1,all?n-1:10,residues,0);
        return 0;
    }
    for(int rep=0;rep<2;++rep) for(int inward=0;inward<2;++inward)
        for(size_t n=30;n<=50;n+=20)for(int all=0;all<2;++all) {
            size_t e=n==30?4:8,r=all?n-(rep?1:e):10;
            compare_run(rep,inward,n,e,r,residues,check);
        }
    if(check) for(int rep=0;rep<2;++rep) for(int inward=0;inward<2;++inward)
        for(size_t d=1;d<=2048;d=d<33?d+16:d<1024?d*2:2049)
            compare_run(rep,inward,20,rep?3:4,rep?19:16,d,1);
    if(check) for(int rep=0;rep<2;++rep) for(int inward=0;inward<2;++inward)
        for(size_t n=8;n<=56;n+=16) for(size_t d=128;d<=1025;d=d==128?129:d==129?1024:d==1024?1025:1026)
            compare_run(rep,inward,n,rep?0:1,n-1,d,1);
    if(check) for(int inward=0;inward<2;++inward) {
        compare_run(1,inward,36,31,3,1,1);
        compare_run(1,inward,36,32,3,1,1);
        compare_run(0,inward,36,31,5,1,1);
    }
    if (check) fprintf(stderr,"%zu configurations, %zu windows match frozen scores, seeds, positions and substitution flags\n",
        verified_windows/19,verified_windows);
    return 0;
}
