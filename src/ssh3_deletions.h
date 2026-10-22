#ifndef SSH3_DELETIONS_H
#define SSH3_DELETIONS_H
#include <stddef.h>
#include <stdint.h>
#include "ssh3_pivots.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Original SubseqHash2 orders: reversed left, normal right; minimize psi,
   then maximize omega. A includes B2's sign (and pivot A includes combine3).
   Tables are immutable after initialization. No mixed-edit API is changed. */
typedef struct {
    int64_t *A;
    int8_t *B1;
    size_t *C;
} SubseqHash3DeletionSide;
typedef struct {
    size_t k,d;
    SubseqHash3DeletionSide left,right;
    int64_t *pivot_A;
    int8_t *left_sign,*right_sign;
    size_t *pivot_C;
} SubseqHash3DeletionTables;
typedef struct {
    size_t length,deletions,psi;
    int64_t omega; /* Exact integer, including scores beyond double precision. */
    uint8_t *seed;
    size_t *positions;
} SubseqHash3DeletionResult;
typedef struct SubseqHash3DeletionWorkspace SubseqHash3DeletionWorkspace;

/* Empty table on entry; 1<=k<=256, D>0. Initialization uses the SAME
   parameters as pivot_tables_init(k+1,D,seed), with lossless widening to int64.
   Import reads the ORIGINAL SubseqHash2 text-table format,
   including the final pivot permutation (validated, but API pivots are explicit).
   Imported magnitudes must be <= INT64_MAX/(4*k), ensuring exact arithmetic. */
int subseqhash3_deletion_tables_init(SubseqHash3DeletionTables *,size_t k,
    size_t residues,uint64_t seed);
/* Copy the exact parameters of a common pivot family into owned arrays.
   B2 is folded into A; source remains independently owned and unchanged.
   Side traversal still follows the deletion API's original SubseqHash2 order. */
int subseqhash3_deletion_tables_from_pivots(SubseqHash3DeletionTables *,
    const SubseqHash3PivotTables *);
int subseqhash3_deletion_tables_import(SubseqHash3DeletionTables *,
    const char *path,size_t k,size_t residues);
void subseqhash3_deletion_tables_destroy(SubseqHash3DeletionTables *);
SubseqHash3DeletionWorkspace *subseqhash3_deletion_workspace_create(void);
void subseqhash3_deletion_workspace_destroy(SubseqHash3DeletionWorkspace *);
void subseqhash3_deletion_results_destroy(SubseqHash3DeletionResult *,size_t count);

/* Exactly E deletions, zero substitutions; output length k=n-E. n<=256,
   0<=E<n, tables->k>=k. Input symbols 0..3. 1<=count<=k, NULL pivots means
   0..count-1; explicit pivots must be <k (duplicates allowed). All orders in
   one call share side DPs. Results must be empty; failure clears them.
   Ties prefer the earliest input pivot, then smallest right residue; side
   DP ties prefer deletion. Any tied minimizer is valid in SubseqHash2.
   Workspace is reusable across dimensions but must not be shared concurrently. */
int subseqhash3_compute_deletion_seeds_with_workspace(const uint8_t *,size_t n,
    size_t deletions,const SubseqHash3DeletionTables *,const size_t *pivots,
    size_t count,SubseqHash3DeletionResult *,SubseqHash3DeletionWorkspace *);
int subseqhash3_compute_deletion_seeds(const uint8_t *,size_t n,size_t deletions,
    const SubseqHash3DeletionTables *,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *);
/* Alternative hash family: normal left prefix, reversed right suffix.
   Same minimization, dimensions, pivot and workspace rules as above, but
   scores and seeds differ from the outward SubseqHash2 order. All pivots
   share two side DPs. Wide imported table scores remain exact. */
int subseqhash3_compute_inward_deletion_seeds_with_workspace(const uint8_t *,size_t n,
    size_t deletions,const SubseqHash3DeletionTables *,const size_t *pivots,
    size_t count,SubseqHash3DeletionResult *,SubseqHash3DeletionWorkspace *);
int subseqhash3_compute_inward_deletion_seeds(const uint8_t *,size_t n,size_t deletions,
    const SubseqHash3DeletionTables *,const size_t *pivots,size_t count,
    SubseqHash3DeletionResult *);
/* Same streaming callback/lifetime/position rules as repeats_windows.
   Request the first count pivot orders (1..window-deletions). Wide imported
   coefficients and direct one-deletion cases retain the single-window solver.
   Tables must remain immutable until this call returns. */
typedef int (*SubseqHash3DeletionWindowCallback)(size_t window_start,
    const SubseqHash3DeletionResult *results,size_t count,void *context);
int subseqhash3_seed_deletion_windows(const uint8_t *sequence,size_t length,
    size_t window,size_t deletions,const SubseqHash3DeletionTables *tables,
    size_t count,SubseqHash3Traversal traversal,
    SubseqHash3DeletionWindowCallback callback,void *context);

#ifdef __cplusplus
}
#endif
#endif
