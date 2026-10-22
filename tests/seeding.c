#include "ssh3_seeding.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint8_t input[64];
    SubseqHash3PivotTables tables[3];
    size_t n,e,r,sets,seen;
    int inward,cancel;
} Check;
static int sub_check(size_t start,const SubseqHash3EditResult *out,size_t count,void *ctx) {
    Check *c=ctx; assert(count==c->r*c->sets && start==c->seen++);
    for (size_t f=0;f<c->sets;++f) {
        SubseqHash3EditResult expected[32]={{0}};
        assert(!(c->inward?subseqhash3_compute_pivot_seeds:subseqhash3_compute_outward_pivot_seeds)(
            c->input+start,c->n,c->e,c->tables+f,NULL,c->r,expected));
        for (size_t p=0;p<c->r;++p) {
            const SubseqHash3EditResult *a=out+f*c->r+p,*b=expected+p;
            assert(a->length==b->length && a->psi==b->psi && a->omega==b->omega);
            assert(!memcmp(a->seed,b->seed,a->length));
            assert(!memcmp(a->positions,b->positions,a->length*sizeof(size_t)));
            assert(!memcmp(a->substituted,b->substituted,a->length));
        }
        subseqhash3_edit_results_destroy(expected,c->r);
    }
    return c->cancel;
}
static int ind_check(size_t start,const SubseqHash3DeletionResult *out,size_t count,void *ctx) {
    Check *c=ctx; assert(count==c->r*c->sets && start==c->seen++);
    for (size_t f=0;f<c->sets;++f) {
        SubseqHash3DeletionTables table={0}; SubseqHash3DeletionResult expected[32]={{0}};
        assert(!subseqhash3_deletion_tables_from_pivots(&table,c->tables+f));
        assert(!(c->inward?subseqhash3_compute_inward_deletion_seeds:subseqhash3_compute_deletion_seeds)(
            c->input+start,c->n,c->e,&table,NULL,c->r,expected));
        for (size_t p=0;p<c->r;++p) {
            const SubseqHash3DeletionResult *a=out+f*c->r+p,*b=expected+p;
            assert(a->length==b->length && a->psi==b->psi && a->omega==b->omega);
            assert(!memcmp(a->seed,b->seed,a->length));
            assert(!memcmp(a->positions,b->positions,a->length*sizeof(size_t)));
        }
        subseqhash3_deletion_results_destroy(expected,c->r);
        subseqhash3_deletion_tables_destroy(&table);
    }
    return c->cancel;
}
int main(void) {
    const size_t ds[]={1,7,32,33,129,1024};
    size_t cases=0;
    for (size_t d=0;d<sizeof(ds)/sizeof(*ds);++d) {
        Check c={0}; c.n=7; c.sets=3;
        uint64_t state=812733;
        for (size_t i=0;i<sizeof(c.input);++i) { state^=state<<13;state^=state>>7;state^=state<<17;c.input[i]=(uint8_t)(state&3); }
        for (size_t f=0;f<c.sets;++f) assert(!subseqhash3_pivot_tables_init(c.tables+f,8,ds[d],43+f));
        for (c.inward=0;c.inward<2;++c.inward) for (c.e=0;c.e<c.n;++c.e)
            for (int all=0;all<2;++all) for (int ind=0;ind<2;++ind) {
                c.r=all?(ind?c.n-c.e:c.n-1):1;c.seen=0;
                int status=ind?subseqhash3_seeds_ind(c.input,21,c.n,c.e,c.tables,c.r,c.sets,
                    (SubseqHash3Traversal)c.inward,ind_check,&c):
                    subseqhash3_seeds_sub(c.input,21,c.n,c.e,c.tables,c.r,c.sets,
                    (SubseqHash3Traversal)c.inward,sub_check,&c);
                assert(!status && c.seen==21-c.n+1);++cases;
            }
        c.r=2;c.e=1;c.inward=0;c.cancel=1;c.seen=0;
        assert(subseqhash3_seeds_sub(c.input,21,c.n,c.e,c.tables,c.r,c.sets,SUBSEQHASH3_OUTWARD,sub_check,&c)==-1 && c.seen==1);
        c.seen=0;
        assert(subseqhash3_seeds_ind(c.input,21,c.n,c.e,c.tables,c.r,c.sets,SUBSEQHASH3_OUTWARD,ind_check,&c)==-1 && c.seen==1);
        c.cancel=0;c.seen=0;
        c.input[20]=4;
        assert(subseqhash3_seeds_sub(c.input,21,c.n,c.e,c.tables,c.r,c.sets,SUBSEQHASH3_OUTWARD,sub_check,&c)==-1 && !c.seen);
        assert(subseqhash3_seeds_ind(c.input,21,c.n,c.e,c.tables,c.r,c.sets,SUBSEQHASH3_OUTWARD,ind_check,&c)==-1 && !c.seen);
        c.input[20]=0;
        const size_t saved_d=c.tables[2].d;
        c.tables[2].d=0;
        assert(subseqhash3_seeds_sub(c.input,21,c.n,1,c.tables,2,c.sets,SUBSEQHASH3_OUTWARD,sub_check,&c)==-1 && !c.seen);
        assert(subseqhash3_seeds_ind(c.input,21,c.n,1,c.tables,2,c.sets,SUBSEQHASH3_OUTWARD,ind_check,&c)==-1 && !c.seen);
        c.tables[2].d=saved_d;
        assert(subseqhash3_seeds_sub(c.input,21,c.n,1,c.tables,c.n,c.sets,SUBSEQHASH3_OUTWARD,sub_check,&c)==-1);
        assert(subseqhash3_seeds_ind(c.input,21,c.n,1,c.tables,c.n,c.sets,SUBSEQHASH3_OUTWARD,ind_check,&c)==-1);
        assert(subseqhash3_seeds_sub(c.input,21,c.n,1,c.tables,2,SIZE_MAX,SUBSEQHASH3_OUTWARD,sub_check,&c)==-1);
        assert(subseqhash3_seeds_ind(c.input,21,c.n,1,c.tables,2,SIZE_MAX,SUBSEQHASH3_OUTWARD,ind_check,&c)==-1);
        for (size_t f=0;f<c.sets;++f) subseqhash3_pivot_tables_destroy(c.tables+f);
    }
    printf("R*S seeding: %zu configurations match independent per-set APIs, including exact witnesses and cancellation\n",cases);
    return 0;
}
