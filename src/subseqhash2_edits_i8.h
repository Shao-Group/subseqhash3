#ifndef SUBSEQHASH2_EDITS_I8_H
#define SUBSEQHASH2_EDITS_I8_H

#include "subseqhash_edits.h"

/* ABC_K tables for the forward, reverse, and pivot score functions. */
typedef struct {
    size_t k, d;
    int8_t *AF, *BF1, *BF2;
    int8_t *AR, *BR1, *BR2;
    int8_t *AP, *BP1, *BP2;
    size_t *CF, *CR, *CP;
} SubseqHash2I8Tables;

int subseqhash2_i8_tables_init(SubseqHash2I8Tables *tables,
                               size_t maximum_length, size_t d,
                               uint64_t seed);
void subseqhash2_i8_tables_destroy(SubseqHash2I8Tables *tables);

/*
 * Number of results produced per window.  Split s has s substitutions,
 * total_edits-s deletions, length window_length-total_edits+s, and one
 * result for every pivot order in [0,length).  Returns SIZE_MAX on
 * overflow or invalid total_edits.
 */
size_t subseqhash2_i8_results_per_window(size_t window_length,
                                         size_t total_edits);

/*
 * Compute all pivot-order seeds for every length-window_length window of
 * sequence.  Results are window-major, then substitution count 0..E,
 * then order 0..length-1.  result->psi remains the hash residue; use
 * subseqhash2_i8_result_coordinates() to recover the split and order from
 * its per-window result index.  positions are absolute zero-based sequence
 * positions.
 *
 * results must contain
 *   (sequence_length-window_length+1) *
 *   subseqhash2_i8_results_per_window(window_length,total_edits)
 * entries.  Destroy them with subseqhash_edit_results_destroy().
 */
int subseqhash2_compute_long_edits_i8(
    const uint8_t *sequence, size_t sequence_length,
    size_t window_length, size_t total_edits,
    const SubseqHash2I8Tables *tables, SubseqHashEditResult *results);

/* Decode a per-window result index. */
int subseqhash2_i8_result_coordinates(size_t window_length,
                                      size_t total_edits,
                                      size_t result_index,
                                      size_t *substitutions,
                                      size_t *order);

#endif
