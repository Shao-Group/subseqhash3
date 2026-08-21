#ifndef SUBSEQHASH_ADJACENT_H
#define SUBSEQHASH_ADJACENT_H

#include "subseqhash_edits.h"

enum { SUBSEQHASH_ADJACENT_BOS = 4 };

typedef struct {
    size_t k, d;
    double *A;
    int8_t *B1, *B2;
    size_t *C;
} SubseqHashAdjacentTables;

int subseqhash_adjacent_tables_init(SubseqHashAdjacentTables *tables,
                                    size_t maximum_length, size_t d,
                                    uint64_t seed);
void subseqhash_adjacent_tables_destroy(SubseqHashAdjacentTables *tables);

int subseqhash_adjacent_score(const uint8_t *seed, size_t length,
                              const SubseqHashAdjacentTables *tables,
                              size_t *psi, double *omega);

int subseqhash_adjacent_compute_all_edit_splits(
    const uint8_t *s, size_t n, size_t total_edits,
    const SubseqHashAdjacentTables *tables,
    SubseqHashEditResult *results);

int subseqhash_adjacent_compute_edits(
    const uint8_t *s, size_t n, size_t deletions, size_t substitutions,
    const SubseqHashAdjacentTables *tables,
    uint8_t *seed, size_t *positions, uint8_t *substituted,
    size_t *psi, double *omega);

#endif
