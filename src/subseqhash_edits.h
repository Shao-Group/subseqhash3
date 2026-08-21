#ifndef SUBSEQHASH_EDITS_H
#define SUBSEQHASH_EDITS_H

#include "subseqhash.h"

/*
 * Select exactly n - deletions input positions and change exactly
 * substitutions of their symbols.  A substitution must change the symbol.
 *
 * tables->k must equal n - deletions.  seed has tables->k entries.
 * positions and substituted may be NULL; otherwise they also have
 * tables->k entries.  positions are zero-based indices into s.
 */
int subseqhash_compute_edits(const uint8_t *s, size_t n,
                             size_t deletions, size_t substitutions,
                             const SubseqHashTables *tables,
                             uint8_t *seed, size_t *positions,
                             uint8_t *substituted,
                             size_t *psi, double *omega);

typedef struct {
    size_t deletions, substitutions, length, psi;
    double omega;
    uint8_t *seed;
    size_t *positions;
    uint8_t *substituted;
} SubseqHashEditResult;

/*
 * Compute all total_edits + 1 splits in one DP run.  results must have
 * total_edits + 1 entries.  tables->k must equal n; a result of length k
 * uses ABC rows 0 through k-1.  Call subseqhash_edit_results_destroy()
 * after a successful call.
 */
int subseqhash_compute_all_edit_splits(const uint8_t *s, size_t n,
                                       size_t total_edits,
                                       const SubseqHashTables *tables,
                                       SubseqHashEditResult *results);
void subseqhash_edit_results_destroy(SubseqHashEditResult *results,
                                     size_t count);

#endif
