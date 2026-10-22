#define _POSIX_C_SOURCE 200809L
#include "ssh3_deletions.h"
#include "ssh3_commands.h"
#include "pivot_table_file.h"
#include "table_file.h"
#include <inttypes.h>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static uint64_t state=917;
static unsigned base(void) {
    state^=state<<13; state^=state>>7; state^=state<<17;
    return (unsigned)(state&3);
}
/* Literal SubseqHash2 order; no DP states, reflection, or cost encoding. */
static void evaluate(const uint8_t *s,size_t k,size_t p,
    const SubseqHash3DeletionTables *t,size_t *psi,int64_t *omega,int inward) {
    size_t j,lr=0,rr=0;
    int64_t lv=0,rv=0;
    for (j=0;j<p;++j) {
        const size_t row=4*j+s[inward?j:p-1-j];
        lr=(lr+t->left.C[row])%t->d;
        lv=t->left.B1[row*t->d+lr]*lv+t->left.A[row*t->d+lr];
    }
    for (j=0;j<k-p-1;++j) {
        const size_t row=4*j+s[inward?k-1-j:p+1+j];
        rr=(rr+t->right.C[row])%t->d;
        rv=t->right.B1[row*t->d+rr]*rv+t->right.A[row*t->d+rr];
    }
    const size_t row=4*p+s[p];
    *psi=(lr+rr+t->pivot_C[row])%t->d;
    *omega=t->left_sign[row]*lv+t->right_sign[row]*rv+t->pivot_A[row];
}
typedef struct {
    const uint8_t *s;
    size_t n,k;
    int inward;
    const SubseqHash3DeletionTables *t;
    uint8_t seed[256];
    size_t psi[256];
    int64_t omega[256];
} Brute;
static void enumerate(Brute *b,size_t input,size_t kept) {
    if (kept==b->k) {
        size_t p;
        for (p=0;p<b->k;++p) {
            size_t psi; int64_t omega;
            evaluate(b->seed,b->k,p,b->t,&psi,&omega,b->inward);
            if (psi<b->psi[p] || (psi==b->psi[p] && omega>b->omega[p])) {
                b->psi[p]=psi; b->omega[p]=omega;
            }
        }
        return;
    }
    for (;input+b->k-kept<=b->n;++input) {
        b->seed[kept]=b->s[input]; enumerate(b,input+1,kept+1);
    }
}
static void validate(const uint8_t *s,size_t n,size_t e,size_t p,
    const SubseqHash3DeletionTables *t,const SubseqHash3DeletionResult *out,int inward) {
    size_t j,psi; int64_t omega;
    assert(out->length==n-e && out->deletions==e);
    for (j=0;j<out->length;++j) {
        assert(out->positions[j]<n && (!j || out->positions[j]>out->positions[j-1]));
        assert(out->seed[j]==s[out->positions[j]]);
    }
    evaluate(out->seed,out->length,p,t,&psi,&omega,inward);
    assert(psi==out->psi && omega==out->omega);
}
static void same(const SubseqHash3DeletionResult *a,const SubseqHash3DeletionResult *b) {
    assert(a->length==b->length && a->psi==b->psi && a->omega==b->omega);
    assert(!memcmp(a->seed,b->seed,a->length));
    assert(!memcmp(a->positions,b->positions,a->length*sizeof(size_t)));
}
static void run_order(size_t n,size_t e,size_t d,unsigned mode,
    SubseqHash3DeletionWorkspace *workspace,int inward) {
    SubseqHash3DeletionTables t={0};
    SubseqHash3DeletionResult out[256]={{0}},subset[4]={{0}};
    uint8_t s[256];
    size_t i,k=n-e,pivots[4]={k-1,0,k/2,k-1};
    Brute b={0};
    for (i=0;i<n;++i) s[i]=(uint8_t)(mode==1?0:base());
    assert(!subseqhash3_deletion_tables_init(&t,k,d,71+mode));
    if (mode==1) { /* Keep coverage beyond double precision after unification. */
        for (i=0;i<4*k*d;++i) {
            t.left.A[i]=t.left.A[i]*(INT64_C(1)<<45)+7;
            t.right.A[i]=t.right.A[i]*(INT64_C(1)<<45)+11;
        }
        for (i=0;i<4*k;++i) t.pivot_A[i]=t.pivot_A[i]*(INT64_C(1)<<45)+3;
    }
    if (mode==4) { /* Compact-cost limits, sign flips and long finite costs. */
        for (i=0;i<4*k*d;++i) {
            t.left.A[i]=(i&1)?127:-127;
            t.right.A[i]=(i&2)?127:-127;
            t.left.B1[i]=(i&4)?1:-1; t.right.B1[i]=(i&8)?1:-1;
        }
    }
    if (mode==5) { /* An unused wide coefficient must not affect selected witnesses. */
        t.left.A[4*k*d-1]=128; t.right.A[4*k*d-1]=-128;
    }
    if (mode==2 || mode==3) {
        for (i=0;i<4*k*d;++i) {
            t.left.A[i]=t.right.A[i]=1;
            t.left.B1[i]=t.right.B1[i]=1;
        }
        for (i=0;i<4*k;++i) {
            t.left.C[i]=t.right.C[i]=0; t.pivot_C[i]=mode==3?d-1:0;
            t.pivot_A[i]=1; t.left_sign[i]=t.right_sign[i]=1;
        }
    }
    assert(!(inward?subseqhash3_compute_inward_deletion_seeds_with_workspace:
        subseqhash3_compute_deletion_seeds_with_workspace)(s,n,e,&t,NULL,k,out,workspace));
    if (n<=10) {
        b.s=s; b.n=n; b.k=k; b.t=&t; b.inward=inward;
        for (i=0;i<k;++i) { b.psi[i]=d; b.omega[i]=INT64_MIN; }
        enumerate(&b,0,0);
    }
    for (i=0;i<k;++i) {
        validate(s,n,e,i,&t,out+i,inward);
        if (n<=10) assert(out[i].psi==b.psi[i] && out[i].omega==b.omega[i]);
    }
    const size_t count=k<4?k:4;
    assert(!(inward?subseqhash3_compute_inward_deletion_seeds:
        subseqhash3_compute_deletion_seeds)(s,n,e,&t,pivots,count,subset));
    for (i=0;i<count;++i) same(out+pivots[i],subset+i);
    subseqhash3_deletion_results_destroy(subset,count);
    subseqhash3_deletion_results_destroy(out,k);
    subseqhash3_deletion_tables_destroy(&t);
}
static void run(size_t n,size_t e,size_t d,unsigned mode,
    SubseqHash3DeletionWorkspace *workspace) {
    run_order(n,e,d,mode,workspace,0);
    run_order(n,e,d,mode,workspace,1);
}
static void import_checks(void) {
    char path[]="/tmp/ssh3-deletion-table-XXXXXX";
    int fd=mkstemp(path);
    SubseqHash3DeletionTables t={0};
    assert(fd>=0); assert(!close(fd));
    for (unsigned bad=0;bad<10;++bad) {
        FILE *f=fopen(path,"w"); assert(f);
        for (size_t side=0;side<2;++side) {
            for (size_t i=0;i<4;++i) {
                if (bad==1 && side==0 && i==0) fputs("9999999999999999999999999999999999\n",f);
                else fprintf(f,"%d\n",bad==2 && side==0 && i==0?0:1);
            }
            for (size_t i=0;i<4;++i) fprintf(f,"%d,1\n",bad==3 && i==0?2:1);
            for (size_t i=0;i<4;++i) fprintf(f,"%d\n",bad==4 && i==0?1:0);
        }
        for (size_t i=0;i<4;++i) fprintf(f,"1,%d\n",bad==5 && i==0?0:-1);
        for (size_t i=0;i<4;++i) fputs("1\n",f);
        for (size_t i=0;i<4;++i) fputs("-1\n",f);
        for (size_t i=0;i<4;++i) fputs("0\n",f);
        if (bad!=6) fprintf(f,"%d\n",bad==7?1:0);
        if (bad==8) fputs("extra\n",f);
        if (bad==9) fputs("999999999999999999999999999999\n",f);
        assert(!fclose(f));
        if (bad) assert(subseqhash3_deletion_tables_import(&t,path,1,1)==-1 && !t.k);
        else {
            assert(!subseqhash3_deletion_tables_import(&t,path,1,1));
            assert(t.k==1 && t.d==1 && t.left.A[0]==1 && t.pivot_A[0]==-1);
            subseqhash3_deletion_tables_destroy(&t);
        }
    }
    assert(!unlink(path));
}
static void compare_parameters(const SubseqHash3PivotTables *p,const SubseqHash3DeletionTables *t) {
    assert(p->k==t->k && p->d==t->d);
    for (size_t i=0;i<4*t->k*t->d;++i) {
        assert(t->left.A[i]==(int64_t)p->left.A[i]*p->left.B2[i]);
        assert(t->right.A[i]==(int64_t)p->right.A[i]*p->right.B2[i]);
        assert(t->left.B1[i]==p->left.B1[i] && t->right.B1[i]==p->right.B1[i]);
    }
    for (size_t i=0;i<4*t->k;++i) {
        assert(t->left.C[i]==p->left.C[i] && t->right.C[i]==p->right.C[i]);
        assert(t->pivot_A[i]==p->pivot_A[i] && t->pivot_C[i]==p->pivot_C[i]);
        assert(t->left_sign[i]==p->pivot_left_sign[i] && t->right_sign[i]==p->pivot_right_sign[i]);
    }
}
static void shared_tables(void) {
    SubseqHash3DeletionTables t={0};
    /* Generator/conversion parity at minimum and maximum row capacities. */
    for (size_t k=1;k<=256;k+=255) {
        SubseqHash3PivotTables p={0};
        assert(!subseqhash3_pivot_tables_init(&p,k+1,33,17));
        assert(!subseqhash3_deletion_tables_init(&t,k,33,17));
        compare_parameters(&p,&t); subseqhash3_deletion_tables_destroy(&t);
        p.left.B2[0]=0;
        assert(subseqhash3_deletion_tables_from_pivots(&t,&p)==-1 && !t.k);
        subseqhash3_pivot_tables_destroy(&p);
    }
}
/* Force wide fallback inside the packed rank range, using an allele absent
   from every input. Compare full witnesses in both orders and reuse compact
   workspaces afterwards, including D=1025 partial vectors. */
static void compact_wide_kernels_parity(void) {
    const size_t residues[]={33,1024,1025};
    SubseqHash3DeletionWorkspace *w=subseqhash3_deletion_workspace_create();
    assert(w);
    for (unsigned inward=0;inward<2;++inward)
    for (size_t scenario=0;scenario<2;++scenario) for (size_t di=0;di<3;++di) {
        const size_t n=scenario?50:30,e=scenario?8:4,k=n-e;
        SubseqHash3DeletionTables t={0};
        SubseqHash3DeletionResult compact[256]={{0}},wide[256]={{0}};
        uint8_t s[256];
        for (size_t j=0;j<n;++j) s[j]=(uint8_t)(base()%3);
        assert(!subseqhash3_deletion_tables_init(&t,k,residues[di],37));
        assert(!(inward?subseqhash3_compute_inward_deletion_seeds_with_workspace:
            subseqhash3_compute_deletion_seeds_with_workspace)(s,n,e,&t,NULL,k,compact,w));
        const size_t unused=4*t.d-1; /* allele 3 at the first side rank */
        const int64_t saved=t.left.A[unused]; t.left.A[unused]=128;
        assert(!(inward?subseqhash3_compute_inward_deletion_seeds_with_workspace:
            subseqhash3_compute_deletion_seeds_with_workspace)(s,n,e,&t,NULL,k,wide,w));
        for (size_t j=0;j<k;++j) same(compact+j,wide+j);
        subseqhash3_deletion_results_destroy(wide,k); t.left.A[unused]=saved;
        assert(!(inward?subseqhash3_compute_inward_deletion_seeds_with_workspace:
            subseqhash3_compute_deletion_seeds_with_workspace)(s,n,e,&t,NULL,k,wide,w));
        for (size_t j=0;j<k;++j) same(compact+j,wide+j);
        subseqhash3_deletion_results_destroy(compact,k); subseqhash3_deletion_results_destroy(wide,k);
        subseqhash3_deletion_tables_destroy(&t);
    }
    subseqhash3_deletion_workspace_destroy(w);
}
static void inward_repeats_parity(void) {
    const size_t lengths[]={2,9,30,50},residues[]={7,31,33,1024};
    for (size_t ni=0;ni<4;++ni) for (size_t di=0;di<4;++di) {
        const size_t n=lengths[ni],k=n-1;
        SubseqHash3PivotTables p={0}; SubseqHash3DeletionTables t={0};
        SubseqHash3DeletionResult out[256]={{0}};
        SubseqHash3EditResult reference[256]={{0}};
        uint8_t s[256];
        for (size_t j=0;j<n;++j) s[j]=(uint8_t)base();
        assert(!subseqhash3_pivot_tables_init(&p,n,residues[di],17));
        assert(!subseqhash3_deletion_tables_from_pivots(&t,&p));
        assert(!subseqhash3_compute_inward_deletion_seeds(s,n,1,&t,NULL,k,out));
        assert(!subseqhash3_compute_pivot_seeds(s,n,0,&p,NULL,k,reference));
        for (size_t j=0;j<k;++j) {
            assert(out[j].psi==reference[j].psi && out[j].omega==reference[j].omega);
            validate(s,n,1,j,&t,out+j,1);
        }
        subseqhash3_deletion_results_destroy(out,k);
        subseqhash3_edit_results_destroy(reference,k);
        assert(!subseqhash3_compute_deletion_seeds(s,n,1,&t,NULL,k,out));
        assert(!subseqhash3_compute_outward_pivot_seeds(s,n,0,&p,NULL,k,reference));
        for (size_t j=0;j<k;++j) {
            assert(out[j].psi==reference[j].psi && out[j].omega==reference[j].omega);
            validate(s,n,1,j,&t,out+j,0);
        }
        subseqhash3_deletion_results_destroy(out,k);
        subseqhash3_edit_results_destroy(reference,k);
        subseqhash3_deletion_tables_destroy(&t); subseqhash3_pivot_tables_destroy(&p);
    }
}
int main(void) {
    SubseqHash3DeletionWorkspace *w=subseqhash3_deletion_workspace_create();
    size_t n,e,d,mode,cases=0;
    const size_t residues[]={1,2,7,8,9,31,32,33,65,128,129};
    assert(w);
    for (n=1;n<=9;++n) for (e=0;e<n;++e)
        for (d=0;d<sizeof(residues)/sizeof(*residues);++d) for (mode=0;mode<6;++mode) {
            run(n,e,residues[d],(unsigned)mode,w); ++cases;
        }
    run(30,5,1024,0,w); run(50,8,1025,0,w);
    run(256,1,33,4,w);
    run(256,0,1,0,w); run(256,255,33,0,w); run(256,1,31,0,w);
    {
        uint8_t s[]={0,1,2,3}; size_t p=3;
        SubseqHash3DeletionTables t={0}; SubseqHash3DeletionResult out[4]={{0}};
        assert(!subseqhash3_deletion_tables_init(&t,4,7,12));
        assert(subseqhash3_compute_deletion_seeds(s,4,4,&t,NULL,1,out)==-1);
        assert(subseqhash3_compute_deletion_seeds(s,4,1,&t,&p,1,out)==-1);
        assert(subseqhash3_compute_deletion_seeds(s,4,1,&t,NULL,4,out)==-1);
        assert(subseqhash3_compute_deletion_seeds(s,0,0,&t,NULL,1,out)==-1);
        s[1]=4; assert(subseqhash3_compute_deletion_seeds(s,4,0,&t,NULL,1,out)==-1);
        subseqhash3_deletion_tables_destroy(&t);
    }
    {
        char path[]="/tmp/ssh3-deletion-import-XXXXXX";
        int fd=mkstemp(path); FILE *f;
        SubseqHash3DeletionTables t={0};
        assert(fd>=0); f=fdopen(fd,"w"); assert(f);
        fputs("invalid\n",f); assert(!fclose(f));
        assert(subseqhash3_deletion_tables_import(&t,path,2,7)==-1 && !t.k);
        assert(!unlink(path));
    }
    compact_wide_kernels_parity();
    inward_repeats_parity();
    import_checks();

    shared_tables();
    subseqhash3_deletion_workspace_destroy(w);
    printf("%zu exhaustive cases per traversal (inward and outward) plus large witnesses, ties, unreachable zero, workspace reuse, shared tables and invalid arguments passed\n",cases);
    return 0;
}
