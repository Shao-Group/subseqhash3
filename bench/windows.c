#define _POSIX_C_SOURCE 200809L
#include "ssh3_deletions.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct { uint64_t hash; size_t windows; } Digest;
static void mix(Digest *d,uint64_t x) { d->hash=(d->hash^x)*UINT64_C(1099511628211); }
static int repeats_digest(size_t start,const SubseqHash3EditResult *out,size_t count,void *context) {
    Digest *d=context; ++d->windows; mix(d,start);
    for (size_t a=0;a<count;++a) {
        mix(d,out[a].psi); mix(d,(uint64_t)(int64_t)out[a].omega);
        for (size_t j=0;j<out[a].length;++j) {
            mix(d,out[a].seed[j]); mix(d,out[a].positions[j]); mix(d,out[a].substituted[j]);
        }
    }
    return 0;
}
static int deletions_digest(size_t start,const SubseqHash3DeletionResult *out,size_t count,void *context) {
    Digest *d=context; ++d->windows; mix(d,start);
    for (size_t a=0;a<count;++a) {
        mix(d,out[a].psi); mix(d,(uint64_t)out[a].omega);
        for (size_t j=0;j<out[a].length;++j) { mix(d,out[a].seed[j]); mix(d,out[a].positions[j]); }
    }
    return 0;
}
static double now(void) {
    struct timespec t; assert(!clock_gettime(CLOCK_MONOTONIC,&t));
    return (double)t.tv_sec+(double)t.tv_nsec*1e-9;
}
static Digest run(const uint8_t *s,size_t length,size_t n,size_t e,size_t r,
    SubseqHash3PivotTables *pt,SubseqHash3DeletionTables *dt,int repeats,int inward,int independent) {
    Digest d={UINT64_C(1469598103934665603),0};
    if (!independent) {
        if (repeats) assert(!subseqhash3_seed_repeats_windows(s,length,n,e,pt,1,r,
            (SubseqHash3Traversal)inward,repeats_digest,&d));
        else assert(!subseqhash3_seed_deletion_windows(s,length,n,e,dt,r,
            (SubseqHash3Traversal)inward,deletions_digest,&d));
    } else if (repeats) {
        SubseqHash3PivotWorkspace *w=subseqhash3_pivot_workspace_create(); assert(w);
        SubseqHash3EditResult out[256]={{0}};
        for (size_t start=0;start<=length-n;++start) {
            assert(!(inward?subseqhash3_compute_pivot_repeats_with_workspace:
                subseqhash3_compute_outward_pivot_repeats_with_workspace)(s+start,n,e,pt,1,r,out,w));
            repeats_digest(start,out,r,&d); subseqhash3_edit_results_destroy(out,r);
        }
        subseqhash3_pivot_workspace_destroy(w);
    } else {
        SubseqHash3DeletionWorkspace *w=subseqhash3_deletion_workspace_create(); assert(w);
        SubseqHash3DeletionResult out[256]={{0}};
        for (size_t start=0;start<=length-n;++start) {
            assert(!(inward?subseqhash3_compute_inward_deletion_seeds_with_workspace:
                subseqhash3_compute_deletion_seeds_with_workspace)(s+start,n,e,dt,NULL,r,out,w));
            deletions_digest(start,out,r,&d); subseqhash3_deletion_results_destroy(out,r);
        }
        subseqhash3_deletion_workspace_destroy(w);
    }
    assert(d.windows==length-n+1); return d;
}
int main(int argc,char **argv) {
    if (argc!=7 && argc!=8) {
        fputs("Usage: bench-windows N E D R LENGTH SAMPLES [--independent]\nR=0 means all pivots (separately per algorithm).\n",stderr); return 1;
    }
    size_t n=strtoull(argv[1],NULL,10),e=strtoull(argv[2],NULL,10),d=strtoull(argv[3],NULL,10);
    size_t requested=strtoull(argv[4],NULL,10),length=strtoull(argv[5],NULL,10),samples=strtoull(argv[6],NULL,10);
    assert(n>=2 && n<=256 && e<n && d && length>=n && samples && requested<=n-e && requested<=n-1);
    if (argc==8) assert(!strcmp(argv[7],"--independent"));
    uint8_t *sequence=malloc(length); assert(sequence);
    uint64_t state=823723;
    for (size_t j=0;j<length;++j) { state^=state<<13; state^=state>>7; state^=state<<17; sequence[j]=(uint8_t)(state&3); }
    SubseqHash3PivotTables pt={0}; SubseqHash3DeletionTables dt={0};
    assert(!subseqhash3_pivot_tables_init(&pt,n+1,d,43));
    assert(!subseqhash3_deletion_tables_from_pivots(&dt,&pt));
    puts("algorithm,method,n,E,D,seeds,length,windows,sample,us_per_window,checksum");
    for (int repeats=0;repeats<2;++repeats) {
        const size_t r=requested?requested:n-(repeats?1:e);
        const size_t prefix=length<n+63?length:n+63;
        for (int inward=0;inward<2;++inward) {
            Digest a=run(sequence,prefix,n,e,r,&pt,&dt,repeats,inward,0);
            Digest b=run(sequence,prefix,n,e,r,&pt,&dt,repeats,inward,1);
            assert(a.hash==b.hash && a.windows==b.windows);
        }
        uint64_t hashes[3]={0};
        for (size_t sample=0;sample<samples;++sample) for (int m=0;m<(argc==8?3:2);++m) {
            const int method=(m+(int)(sample%(argc==8?3:2)))%(argc==8?3:2);
            const int inward=method==0,independent=method==2;
            const double begin=now();
            const Digest digest=run(sequence,length,n,e,r,&pt,&dt,repeats,inward,independent);
            const double elapsed=now()-begin;
            if (sample) assert(hashes[method]==digest.hash);
            hashes[method]=digest.hash;
            if (argc==8 && hashes[1] && hashes[2]) assert(hashes[1]==hashes[2]);
            printf("%s,%s,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%.6f,%" PRIu64 "\n",
                repeats?"repeats":"deletions",inward?"inward":independent?"outward-independent":"outward-shared",
                n,e,d,r,length,digest.windows,sample,elapsed*1e6/digest.windows,digest.hash);
            fflush(stdout);
        }
    }
    subseqhash3_pivot_tables_destroy(&pt); subseqhash3_deletion_tables_destroy(&dt); free(sequence);
    return 0;
}
