#define _POSIX_C_SOURCE 200809L
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <climits>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <algorithm>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>
#include "ssh3_deletions.h"
#include "subseqhash2seeding.h"

static void write_table(const SubseqHash3DeletionTables &t,const char *path) {
    FILE *f=fopen(path,"w"); assert(f);
    for (int side=0;side<2;++side) {
        const SubseqHash3DeletionSide &s=side?t.left:t.right;
        for (size_t i=0;i<4*t.k*t.d;++i) fprintf(f,"%lld\n",(long long)llabs(s.A[i]));
        for (size_t i=0;i<4*t.k*t.d;++i) fprintf(f,"%d,%d\n",(int)s.B1[i],s.A[i]>0?1:-1);
        for (size_t i=0;i<4*t.k;++i) fprintf(f,"%zu\n",s.C[i]);
    }
    for (size_t i=0;i<4*t.k;++i) fprintf(f,"%d,%d\n",(int)t.right_sign[i],(int)t.left_sign[i]);
    for (size_t i=0;i<4*t.k;++i) fprintf(f,"%lld\n",(long long)llabs(t.pivot_A[i]));
    for (size_t i=0;i<4*t.k;++i) fprintf(f,"%d\n",t.pivot_A[i]>0?1:-1);
    for (size_t i=0;i<4*t.k;++i) fprintf(f,"%zu\n",t.pivot_C[i]);
    for (size_t i=0;i<t.k;++i) fprintf(f,"%zu\n",i);
    assert(!fclose(f));
}
static uint64_t rng=7123,checksum=0;
static FILE *samples=nullptr;
static size_t checked_windows=0;
static uint8_t random_base() {
    rng^=rng<<13; rng^=rng>>7; rng^=rng<<17;
    return (uint8_t)(rng&3);
}
static void original_validate(const std::vector<std::vector<seed>> &old,
    const SubseqHash3DeletionResult *out,const std::vector<uint8_t> &s,size_t k) {
    for (size_t p=0;p<old.size();++p) {
        assert(old[p].size()==1);
        const seed &a=old[p][0];
        if ((int64_t)a.hashval!=out[p].omega || (size_t)a.psi!=out[p].psi) {
            fprintf(stderr,"mismatch n=%zu k=%zu pivot=%zu old=(%d,%lld) new=(%zu,%lld)\n",
                s.size(),k,p,(int)a.psi,(long long)a.hashval,out[p].psi,(long long)out[p].omega);
            abort();
        }
        /* Independently decode the ORIGINAL witness and check the exact seed.
           Different witnesses are allowed on score ties; random 55-bit tables
           in these cases are also expected to produce the same seed. */
        kmer encoded=0;
        for (size_t j=0;j<k;++j) encoded=(encoded<<2)|out[p].seed[j];
        assert(encoded==a.str);
    }
}
static void run(size_t n,size_t e,size_t d,size_t count,bool check_only,const char *fixture=nullptr) {
    const size_t k=n-e,batch=16;
    SubseqHash3DeletionTables tables={},loaded={};
    assert(!subseqhash3_deletion_tables_init(&tables,k,d,12+n+e+d));
    if (fixture) {
        subseqhash3_deletion_tables_destroy(&tables);
        assert(!subseqhash3_deletion_tables_import(&tables,fixture,k,d));
    }
    char path[]="/tmp/ssh3-delete-bench-XXXXXX";
    int fd=mkstemp(path); assert(fd>=0); close(fd);
    write_table(tables,path);
    assert(!subseqhash3_deletion_tables_import(&loaded,path,k,d));
    for (size_t i=0;i<4*k*d;++i) {
        assert(tables.left.A[i]==loaded.left.A[i] && tables.right.A[i]==loaded.right.A[i]);
        assert(tables.left.B1[i]==loaded.left.B1[i] && tables.right.B1[i]==loaded.right.B1[i]);
    }
    for (size_t i=0;i<4*k;++i) {
        assert(tables.left.C[i]==loaded.left.C[i] && tables.right.C[i]==loaded.right.C[i]);
        assert(tables.pivot_A[i]==loaded.pivot_A[i] && tables.pivot_C[i]==loaded.pivot_C[i]);
        assert(tables.left_sign[i]==loaded.left_sign[i] && tables.right_sign[i]==loaded.right_sign[i]);
    }
    std::unique_ptr<subseqhash2seeding> original(new subseqhash2seeding((int)n,(int)k,(int)d,(int)count,INT64_MIN));
    original->init(path); assert(!unlink(path));
    // Exactly one input window. Only n origins can be addressed; allocate
    // n*(n+1)*(k+1)*D cells, rather than the original 500-origin chunk buffer.
    const size_t dimension=(n+1)*(k+1)*d,cells=n*dimension;
    std::vector<DPCell> dp(cells),revdp(cells);
    std::vector<int> h(dimension),revh(dimension);
    std::vector<std::vector<seed>> old(count);
    std::vector<SubseqHash3DeletionResult> out(count);
    std::unique_ptr<SubseqHash3DeletionWorkspace,decltype(&subseqhash3_deletion_workspace_destroy)> workspace(
        subseqhash3_deletion_workspace_create(),subseqhash3_deletion_workspace_destroy);
    assert(workspace);
    std::vector<std::vector<uint8_t>> inputs(batch,std::vector<uint8_t>(n));
    std::vector<std::string> strings(batch,std::string(n,'A'));
    for (size_t b=0;b<batch;++b) {
        for (size_t i=0;i<n;++i) strings[b][i]="ACGT"[inputs[b][i]=random_base()];
        original->getSubseq2Seeds(strings[b],dp.data(),revdp.data(),h.data(),revh.data(),old);
        assert(!subseqhash3_compute_deletion_seeds_with_workspace(inputs[b].data(),n,e,&loaded,nullptr,count,out.data(),workspace.get()));
        original_validate(old,out.data(),inputs[b],k); ++checked_windows;
        subseqhash3_deletion_results_destroy(out.data(),count);
        for (auto &v:old) v.clear();
    }
    if (!check_only) {
        std::vector<double> oldtimes,newtimes;
        const size_t rounds=5,repetitions=4;
        auto measure=[&](bool use_new) {
            const auto start=std::chrono::steady_clock::now();
            for (size_t rep=0;rep<repetitions;++rep) for (size_t b=0;b<batch;++b) {
                if (use_new) {
                    assert(!subseqhash3_compute_deletion_seeds_with_workspace(inputs[b].data(),n,e,&loaded,nullptr,count,out.data(),workspace.get()));
                    checksum+=(uint64_t)out[0].omega;
                    subseqhash3_deletion_results_destroy(out.data(),count);
                } else {
                    original->getSubseq2Seeds(strings[b],dp.data(),revdp.data(),h.data(),revh.data(),old);
                    checksum+=(uint64_t)old[0][0].hashval;
                    for (auto &v:old) v.clear();
                }
            }
            return std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-start).count()/(batch*repetitions);
        };
        for (size_t r=0;r<rounds;++r) {
            if (r&1) { newtimes.push_back(measure(true)); oldtimes.push_back(measure(false)); }
            else { oldtimes.push_back(measure(false)); newtimes.push_back(measure(true)); }
        }
        if (samples) for (size_t j=0;j<rounds;++j)
            fprintf(samples,"%zu,%zu,%zu,%zu,%zu,%.3f,%.3f\n",n,e,d,count,j,oldtimes[j],newtimes[j]);
        std::sort(oldtimes.begin(),oldtimes.end()); std::sort(newtimes.begin(),newtimes.end());
        printf("%zu,%zu,%zu,%zu,%.3f,%.3f,%.2f\n",n,e,d,count,oldtimes[rounds/2],newtimes[rounds/2],oldtimes[rounds/2]/newtimes[rounds/2]);
        fflush(stdout);
    }
    subseqhash3_deletion_tables_destroy(&tables); subseqhash3_deletion_tables_destroy(&loaded);
}
int main(int argc,char **argv) {
    const bool check_only=argc==2 && !strcmp(argv[1],"--check-only");
    if (argc==3 && !strcmp(argv[1],"--samples")) {
        samples=fopen(argv[2],"w"); assert(samples);
        fputs("n,deletions,D,seeds,round,subseqhash2_us,new_us\n",samples);
    }
    if (!check_only) puts("n,deletions,D,seeds,subseqhash2_us,new_us,speedup");
    // Original combine uses a signed 32-bit convolution bitmap: D=11 is
    // representative of its distributed tables, and avoids D>16 shift UB.
    for (size_t count: {size_t(1),size_t(10),size_t(25)}) run(30,5,11,count,check_only);
    for (size_t count: {size_t(1),size_t(10),size_t(42)}) run(50,8,11,count,check_only);
    run(20,1,11,19,check_only); run(20,18,11,2,check_only);
    if (access("subseqhash2/alignment/table_subseq2/25_11",R_OK)==0)
        run(30,5,11,25,true,"subseqhash2/alignment/table_subseq2/25_11");
    if (samples) assert(!fclose(samples));
    fprintf(stderr,"%zu windows match original SubseqHash2 scores and seeds; checksum=%llu\n",checked_windows,(unsigned long long)checksum);
}
