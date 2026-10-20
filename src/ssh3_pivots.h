#ifndef SSH3_PIVOTS_H
#define SSH3_PIVOTS_H

#include "subseqhash3.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ABC^k order family for k=n-1 output bases. Unlike the original paper,
   sides are read from the outer boundaries toward the pivot: the left
   prefix in normal order, the right suffix in reverse order. This permits
   sharing across all pivots of a single input window. */
typedef struct {
    size_t k,d;
    SubseqHash3Tables left,right;
    int8_t *pivot_A,*pivot_left_sign,*pivot_right_sign;
    size_t *pivot_C;
    int8_t *right_reversed_A;
    uint32_t *right_reversed_B1_masks;
} SubseqHash3PivotTables;

typedef struct SubseqHash3PivotWorkspace SubseqHash3PivotWorkspace;
SubseqHash3PivotWorkspace *subseqhash3_pivot_workspace_create(void);
void subseqhash3_pivot_workspace_destroy(SubseqHash3PivotWorkspace *workspace);

/* Initialize an empty table object; 2<=maximum_input_length<=257, D>0. */
int subseqhash3_pivot_tables_init(SubseqHash3PivotTables *tables,
    size_t maximum_input_length,size_t residues,uint64_t seed);
void subseqhash3_pivot_tables_destroy(SubseqHash3PivotTables *tables);

/* Exactly one deletion and E substitutions of retained bases; 0<=E<=n-1.
   Return one length-(n-1) minimizer per requested output pivot. A pivot's
   order minimizes psi, then maximizes omega (not |omega|). Tied scores
   may have multiple minimizing strings; a deterministic witness is returned.
   NULL pivots means pivots 0..count-1; otherwise each pivot must be <n-1.
   Results must be empty on entry. Release them with the existing
   subseqhash3_edit_results_destroy(results,count). Return 0 or -1.
   Positions are offsets into the encoded input sequence (symbols 0..3). */
int subseqhash3_compute_pivot_seeds(const uint8_t *sequence,size_t n,
    size_t substitutions,const SubseqHash3PivotTables *tables,
    const size_t *pivots,size_t count,SubseqHash3EditResult *results);

/* Any positive repeat count: use successive independent table families,
   up to k pivots from each, in family-major then pivot-major order.
   Requires table_count >= ceil(repeats/(n-1)). Families are supplied by
   the caller, so random-table generation is reusable across input windows. */
int subseqhash3_compute_pivot_repeats(const uint8_t *sequence,size_t n,
    size_t substitutions,const SubseqHash3PivotTables *tables,
    size_t table_count,size_t repeats,SubseqHash3EditResult *results);

/* Reuse DP allocations across input windows and table families. Workspaces
   are caller-owned and must not be used concurrently by multiple threads. */
int subseqhash3_compute_pivot_repeats_with_workspace(const uint8_t *sequence,
    size_t n,size_t substitutions,const SubseqHash3PivotTables *tables,
    size_t table_count,size_t repeats,SubseqHash3EditResult *results,
    SubseqHash3PivotWorkspace *workspace);

#ifdef __cplusplus
}
#endif
#endif
