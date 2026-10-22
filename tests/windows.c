#include "ssh3_deletions.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const uint8_t *sequence;
    size_t window,edits,count,table_count,seen;
    int inward,stop;
    SubseqHash3PivotTables *pt;
    SubseqHash3DeletionTables *dt;
    SubseqHash3PivotWorkspace *pw;
    SubseqHash3DeletionWorkspace *dw;
} Check;
static int repeats_check(size_t start,const SubseqHash3EditResult *out,size_t count,void *context) {
    Check *c=context;
    assert(start==c->seen++ && count==c->count);
    if (c->stop) return 1;
    SubseqHash3EditResult *ref=calloc(count,sizeof(*ref)); assert(ref);
    assert(!(c->inward?subseqhash3_compute_pivot_repeats_with_workspace:
        subseqhash3_compute_outward_pivot_repeats_with_workspace)(c->sequence+start,
        c->window,c->edits,c->pt,c->table_count,count,ref,c->pw));
    for (size_t a=0;a<count;++a) {
        assert(out[a].length==ref[a].length && out[a].psi==ref[a].psi &&
            out[a].omega==ref[a].omega && out[a].substitutions==ref[a].substitutions &&
            out[a].deletions==ref[a].deletions);
        assert(!memcmp(out[a].seed,ref[a].seed,ref[a].length));
        assert(!memcmp(out[a].positions,ref[a].positions,ref[a].length*sizeof(size_t)));
        assert(!memcmp(out[a].substituted,ref[a].substituted,ref[a].length));
    }
    subseqhash3_edit_results_destroy(ref,count); free(ref); return 0;
}
static int deletions_check(size_t start,const SubseqHash3DeletionResult *out,size_t count,void *context) {
    Check *c=context;
    assert(start==c->seen++ && count==c->count);
    if (c->stop) return 1;
    SubseqHash3DeletionResult ref[256]={{0}};
    assert(!(c->inward?subseqhash3_compute_inward_deletion_seeds_with_workspace:
        subseqhash3_compute_deletion_seeds_with_workspace)(c->sequence+start,c->window,
        c->edits,c->dt,NULL,count,ref,c->dw));
    for (size_t a=0;a<count;++a) {
        assert(out[a].length==ref[a].length && out[a].psi==ref[a].psi &&
            out[a].omega==ref[a].omega && out[a].deletions==ref[a].deletions);
        assert(!memcmp(out[a].seed,ref[a].seed,ref[a].length));
        assert(!memcmp(out[a].positions,ref[a].positions,ref[a].length*sizeof(size_t)));
    }
    subseqhash3_deletion_results_destroy(ref,count); return 0;
}
static size_t configurations,windows;
static void check(size_t n,size_t e,size_t d,int pattern,int wide) {
    const size_t length=3*n+7;
    uint8_t sequence[800];
    unsigned rng=31943;
    for (size_t j=0;j<length;++j) {
        rng=1664525u*rng+1013904223u;
        sequence[j]=(uint8_t)(pattern==0?(rng>>30):pattern==1?0:j%4);
    }
    SubseqHash3PivotTables pt[2]={{0}};
    SubseqHash3DeletionTables dt={0};
    assert(!subseqhash3_pivot_tables_init(pt,n+1,d,43));
    assert(!subseqhash3_pivot_tables_init(pt+1,n+1,d,73));
    assert(!subseqhash3_deletion_tables_from_pivots(&dt,pt));
    if (wide) for (size_t j=0;j<4*dt.k*d;++j) {
        dt.left.A[j]*=INT64_C(10000000000000); dt.right.A[j]*=INT64_C(10000000000000);
    }
    Check c={0}; c.sequence=sequence; c.window=n; c.edits=e;
    c.pt=pt; c.dt=&dt; c.table_count=2;
    c.pw=subseqhash3_pivot_workspace_create(); c.dw=subseqhash3_deletion_workspace_create();
    assert(c.pw && c.dw);
    for (c.inward=0;c.inward<2;++c.inward) for (unsigned partial=0;partial<2;++partial) {
        c.count=partial?1:n-e; c.seen=0;
        assert(!subseqhash3_seed_deletion_windows(sequence,length,n,e,&dt,c.count,
            (SubseqHash3Traversal)c.inward,deletions_check,&c));
        assert(c.seen==length-n+1); windows+=c.seen; ++configurations;
        c.count=partial?1:n+1; c.seen=0;
        assert(!subseqhash3_seed_repeats_windows(sequence,length,n,e,pt,2,c.count,
            (SubseqHash3Traversal)c.inward,repeats_check,&c));
        assert(c.seen==length-n+1); windows+=c.seen; ++configurations;
    }
    c.stop=1; c.seen=0; c.count=1;
    assert(subseqhash3_seed_repeats_windows(sequence,length,n,e,pt,2,1,
        SUBSEQHASH3_OUTWARD,repeats_check,&c)==-1 && c.seen==1);
    c.seen=0;
    assert(subseqhash3_seed_deletion_windows(sequence,length,n,e,&dt,1,
        SUBSEQHASH3_OUTWARD,deletions_check,&c)==-1 && c.seen==1);
    c.seen=0; sequence[length-1]=4;
    assert(subseqhash3_seed_repeats_windows(sequence,length,n,e,pt,2,1,
        SUBSEQHASH3_OUTWARD,repeats_check,&c)==-1 && !c.seen);
    assert(subseqhash3_seed_deletion_windows(sequence,length,n,e,&dt,1,
        SUBSEQHASH3_OUTWARD,deletions_check,&c)==-1 && !c.seen);
    sequence[length-1]=0;
    assert(subseqhash3_seed_repeats_windows(sequence,n-1,n,e,pt,2,1,
        SUBSEQHASH3_OUTWARD,repeats_check,&c)==-1);
    assert(subseqhash3_seed_deletion_windows(sequence,length,n,e,&dt,1,
        (SubseqHash3Traversal)9,deletions_check,&c)==-1);
    subseqhash3_pivot_workspace_destroy(c.pw); subseqhash3_deletion_workspace_destroy(c.dw);
    subseqhash3_deletion_tables_destroy(&dt);
    subseqhash3_pivot_tables_destroy(pt); subseqhash3_pivot_tables_destroy(pt+1);
}
static void endpoints(void) {
    uint8_t sequence[260];
    for (size_t j=0;j<260;++j) sequence[j]=(uint8_t)(j%4);
    SubseqHash3PivotTables pt={0}; SubseqHash3DeletionTables dt={0};
    assert(!subseqhash3_pivot_tables_init(&pt,257,1,43));
    assert(!subseqhash3_deletion_tables_from_pivots(&dt,&pt));
    Check c={0}; c.sequence=sequence; c.pt=&pt; c.dt=&dt; c.count=1; c.table_count=1;
    c.pw=subseqhash3_pivot_workspace_create(); c.dw=subseqhash3_deletion_workspace_create();
    assert(c.pw && c.dw);
    for (c.inward=0;c.inward<2;++c.inward) {
        c.window=257; c.edits=256; c.seen=0;
        assert(!subseqhash3_seed_repeats_windows(sequence,260,257,256,&pt,1,1,
            (SubseqHash3Traversal)c.inward,repeats_check,&c));
        assert(c.seen==4); windows+=c.seen; ++configurations;
        const size_t ns[]={1,256,256},es[]={0,0,255};
        for (size_t j=0;j<3;++j) {
            c.window=ns[j]; c.edits=es[j]; c.seen=0;
            assert(!subseqhash3_seed_deletion_windows(sequence,260,c.window,c.edits,&dt,1,
                (SubseqHash3Traversal)c.inward,deletions_check,&c));
            assert(c.seen==261-c.window); windows+=c.seen; ++configurations;
        }
    }
    subseqhash3_pivot_workspace_destroy(c.pw); subseqhash3_deletion_workspace_destroy(c.dw);
    subseqhash3_pivot_tables_destroy(&pt); subseqhash3_deletion_tables_destroy(&dt);
}

int main(void) {
    endpoints();
    const size_t ds[]={1,2,7,31,32,33,128,129};
    for (size_t d=0;d<sizeof(ds)/sizeof(*ds);++d)
        for (size_t e=0;e<5;++e) for (int pattern=0;pattern<3;++pattern)
            check(5,e,ds[d],pattern,0);
    check(9,3,33,0,1); check(9,1,1024,0,1);
    check(12,4,32,0,0); check(12,8,1,0,0);
    printf("window APIs: %zu configurations, %zu windows match exact per-window scores, seeds and witnesses; callback cancellation and validation passed\n",configurations,windows);
    return 0;
}
