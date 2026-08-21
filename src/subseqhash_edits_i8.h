#ifndef SUBSEQHASH_EDITS_I8_H
#define SUBSEQHASH_EDITS_I8_H

#include "subseqhash_edits.h"

typedef struct {
    size_t k, d;
    int8_t *A;
    int8_t *B1, *B2;
    int8_t *signed_A;
    uint32_t *B1_masks;
    size_t blocks;
    size_t *C;
} SubseqHashI8Tables;

int subseqhash_i8_tables_init(SubseqHashI8Tables *tables,
                               size_t maximum_length, size_t d,
                               uint64_t seed);
void subseqhash_i8_tables_destroy(SubseqHashI8Tables *tables);

/* Construct an exactly equivalent double table for reference comparisons. */
int subseqhash_i8_tables_make_double(const SubseqHashI8Tables *source,
                                      SubseqHashTables *destination);

int subseqhash_compute_all_edit_splits_i8(
    const uint8_t *s, size_t n, size_t total_edits,
    const SubseqHashI8Tables *tables, SubseqHashEditResult *results);

const char *subseqhash_edits_i8_backend(void);

#endif
