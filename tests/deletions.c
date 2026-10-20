#define _POSIX_C_SOURCE 200809L
#include "ssh3_deletions.h"
#include "ssh3_commands.h"
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
    const SubseqHash3DeletionTables *t,size_t *psi,int64_t *omega) {
    size_t j,lr=0,rr=0;
    int64_t lv=0,rv=0;
    for (j=0;j<p;++j) {
        const size_t row=4*j+s[p-1-j];
        lr=(lr+t->left.C[row])%t->d;
        lv=t->left.B1[row*t->d+lr]*lv+t->left.A[row*t->d+lr];
    }
    for (j=0;j<k-p-1;++j) {
        const size_t row=4*j+s[p+1+j];
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
            evaluate(b->seed,b->k,p,b->t,&psi,&omega);
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
    const SubseqHash3DeletionTables *t,const SubseqHash3DeletionResult *out) {
    size_t j,psi; int64_t omega;
    assert(out->length==n-e && out->deletions==e);
    for (j=0;j<out->length;++j) {
        assert(out->positions[j]<n && (!j || out->positions[j]>out->positions[j-1]));
        assert(out->seed[j]==s[out->positions[j]]);
    }
    evaluate(out->seed,out->length,p,t,&psi,&omega);
    assert(psi==out->psi && omega==out->omega);
}
static void same(const SubseqHash3DeletionResult *a,const SubseqHash3DeletionResult *b) {
    assert(a->length==b->length && a->psi==b->psi && a->omega==b->omega);
    assert(!memcmp(a->seed,b->seed,a->length));
    assert(!memcmp(a->positions,b->positions,a->length*sizeof(size_t)));
}
static void run(size_t n,size_t e,size_t d,unsigned mode,
    SubseqHash3DeletionWorkspace *workspace) {
    SubseqHash3DeletionTables t={0};
    SubseqHash3DeletionResult out[256]={{0}},subset[4]={{0}};
    uint8_t s[256];
    size_t i,k=n-e,pivots[4]={k-1,0,k/2,k-1};
    Brute b={0};
    for (i=0;i<n;++i) s[i]=(uint8_t)(mode==1?0:base());
    assert(!subseqhash3_deletion_tables_init(&t,k,d,71+mode));
    if (mode>=2) {
        for (i=0;i<4*k*d;++i) {
            t.left.A[i]=t.right.A[i]=1;
            t.left.B1[i]=t.right.B1[i]=1;
        }
        for (i=0;i<4*k;++i) {
            t.left.C[i]=t.right.C[i]=0; t.pivot_C[i]=mode==3?d-1:0;
            t.pivot_A[i]=1; t.left_sign[i]=t.right_sign[i]=1;
        }
    }
    assert(!subseqhash3_compute_deletion_seeds_with_workspace(s,n,e,&t,NULL,k,out,workspace));
    if (n<=10) {
        b.s=s; b.n=n; b.k=k; b.t=&t;
        for (i=0;i<k;++i) { b.psi[i]=d; b.omega[i]=INT64_MIN; }
        enumerate(&b,0,0);
    }
    for (i=0;i<k;++i) {
        validate(s,n,e,i,&t,out+i);
        if (n<=10) assert(out[i].psi==b.psi[i] && out[i].omega==b.omega[i]);
    }
    const size_t count=k<4?k:4;
    assert(!subseqhash3_compute_deletion_seeds(s,n,e,&t,pivots,count,subset));
    for (i=0;i<count;++i) same(out+pivots[i],subset+i);
    subseqhash3_deletion_results_destroy(subset,count);
    subseqhash3_deletion_results_destroy(out,k);
    subseqhash3_deletion_tables_destroy(&t);
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
static void cli(void) {
    char dir[]="/tmp/ssh3-deletion-cli-XXXXXX",input[256],csv[256],again[256],row[512],expected[512];
    const uint8_t inputs[2][6]={{0,1,2,3,0,1},{3,3,2,2,1,1}};
    SubseqHash3DeletionTables tables={0};
    SubseqHash3DeletionResult out[3]={{0}};
    FILE *f;
    assert(mkdtemp(dir));
    snprintf(input,sizeof(input),"%s/input",dir); snprintf(csv,sizeof(csv),"%s/out",dir);
    snprintf(again,sizeof(again),"%s/again",dir);
    f=fopen(input,"w"); assert(f); fputs("acgtac\r\nTTGGCC\n",f); assert(!fclose(f));
    char *args[]={"deletions",input,"4","1","9","3",csv,"--seed","42"};
    assert(!ssh3_deletions_main(9,args));
    assert(!subseqhash3_deletion_tables_init(&tables,3,9,42));
    f=fopen(csv,"r"); assert(f);
    assert(fgets(row,sizeof(row),f) && !strcmp(row,"pivot,sequence,window_start,deletions,seed,psi,omega\n"));
    for (size_t seq=0;seq<2;++seq) for (size_t start=0;start<3;++start) {
        assert(!subseqhash3_compute_deletion_seeds(inputs[seq]+start,4,1,&tables,NULL,3,out));
        for (size_t p=0;p<3;++p) {
            char dna[4];
            for (size_t j=0;j<3;++j) dna[j]="ACGT"[out[p].seed[j]];
            dna[3]=0;
            snprintf(expected,sizeof(expected),"%zu,%zu,%zu,1,%s,%zu,%" PRId64 "\n",p,seq,start,dna,out[p].psi,out[p].omega);
            assert(fgets(row,sizeof(row),f) && !strcmp(row,expected));
        }
        subseqhash3_deletion_results_destroy(out,3);
    }
    assert(fgetc(f)==EOF && !ferror(f)); assert(!fclose(f));
    args[6]=input; assert(ssh3_deletions_main(9,args)!=0);
    assert(!link(input,again)); args[6]=again; assert(ssh3_deletions_main(9,args)!=0);
    assert(!unlink(again)); args[6]=csv;
    args[3]="4"; assert(ssh3_deletions_main(9,args)!=0); args[3]="1";
    args[5]="4"; assert(ssh3_deletions_main(9,args)!=0); args[5]="3";
    args[4]="-1"; assert(ssh3_deletions_main(9,args)!=0); args[4]="9";
    args[8]="-1"; assert(ssh3_deletions_main(9,args)!=0); args[8]="42";
    f=fopen(input,"w"); assert(f); fputs("ACNT\n",f); assert(!fclose(f));
    assert(ssh3_deletions_main(9,args)!=0);
    f=fopen(input,"w"); assert(f); fputs("AC\n",f); assert(!fclose(f));
    assert(ssh3_deletions_main(9,args)!=0);
    f=fopen(input,"w"); assert(f); fputs("\n",f); assert(!fclose(f));
    assert(ssh3_deletions_main(9,args)!=0);
    subseqhash3_deletion_tables_destroy(&tables);
    assert(!unlink(csv)); assert(!unlink(input)); assert(!rmdir(dir));
}
int main(void) {
    SubseqHash3DeletionWorkspace *w=subseqhash3_deletion_workspace_create();
    size_t n,e,d,mode,cases=0;
    const size_t residues[]={1,2,7,8,9,31,32,33,65};
    assert(w);
    for (n=1;n<=9;++n) for (e=0;e<n;++e)
        for (d=0;d<sizeof(residues)/sizeof(*residues);++d) for (mode=0;mode<4;++mode) {
            run(n,e,residues[d],(unsigned)mode,w); ++cases;
        }
    run(30,5,1024,0,w); run(50,8,1025,0,w);
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
    import_checks();
    cli();
    subseqhash3_deletion_workspace_destroy(w);
    printf("%zu exhaustive deletion-only cases plus large witnesses, ties, unreachable zero, workspace reuse, CSV/API equivalence and invalid arguments passed\n",cases);
    return 0;
}
