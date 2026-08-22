#ifndef SUBSEQHASH_EDITS_A10_H
#define SUBSEQHASH_EDITS_A10_H

#include "subseqhash_edits.h"

/* A values are generated in [1, a_max], where
 * a_max=floor(32766/k).  Thus every omega prefix fits in int16_t;
 * INT16_MAX remains available as the unreachable-state sentinel. */
typedef struct {
    size_t k, d;
    uint16_t *A;
    int8_t *B1, *B2;
    int16_t *signed_A;
    uint16_t a_max;
    uint32_t *B1_masks;
    size_t blocks;
    uint16_t *wrap_indices;
    size_t *C;
} SubseqHashA10Tables;

int subseqhash_a10_tables_init(SubseqHashA10Tables *tables,
                               size_t maximum_length, size_t d,
                               uint64_t seed);
void subseqhash_a10_tables_destroy(SubseqHashA10Tables *tables);

/* Construct an exactly equivalent double table for reference comparisons. */
int subseqhash_a10_tables_make_double(const SubseqHashA10Tables *source,
                                      SubseqHashTables *destination);

int subseqhash_compute_all_edit_splits_a10(
    const uint8_t *s, size_t n, size_t total_edits,
    const SubseqHashA10Tables *tables, SubseqHashEditResult *results);

/* Same DP using the pre-aligned-kernel split-wrap implementation. */
int subseqhash_compute_all_edit_splits_a10_legacy(
    const uint8_t *s, size_t n, size_t total_edits,
    const SubseqHashA10Tables *tables, SubseqHashEditResult *results);

int subseqhash_compute_all_edit_splits_a10_twopass(
    const uint8_t *s, size_t n, size_t total_edits,
    const SubseqHashA10Tables *tables, SubseqHashEditResult *results);

const char *subseqhash_edits_a10_backend(void);

#endif
