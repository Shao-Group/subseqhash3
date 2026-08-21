#ifndef SUBSEQHASH_EDITS_I16_H
#define SUBSEQHASH_EDITS_I16_H

#include "subseqhash_edits.h"

typedef struct {
    size_t k, d;
    int16_t *A;
    int8_t *B1, *B2;
    int16_t *signed_A;
    uint16_t *B1_masks;
    size_t blocks;
    size_t *C;
} SubseqHashI16Tables;

int subseqhash_i16_tables_init(SubseqHashI16Tables *tables,
                               size_t maximum_length, size_t d,
                               uint64_t seed);
void subseqhash_i16_tables_destroy(SubseqHashI16Tables *tables);

/* Construct an exactly equivalent double table for reference comparisons. */
int subseqhash_i16_tables_make_double(const SubseqHashI16Tables *source,
                                      SubseqHashTables *destination);

int subseqhash_compute_all_edit_splits_i16(
    const uint8_t *s, size_t n, size_t total_edits,
    const SubseqHashI16Tables *tables, SubseqHashEditResult *results);

const char *subseqhash_edits_i16_backend(void);

#endif
