#ifndef SUBSEQHASH_EDITS_MAX_I8_H
#define SUBSEQHASH_EDITS_MAX_I8_H

#include "subseqhash_edits.h"

#define SUBSEQHASH_MAX_I8_MAX_LENGTH 258u

typedef struct {
    size_t k, d;
    int8_t *A;
    size_t *C;
} SubseqHashMaxI8Tables;

int subseqhash_max_i8_tables_init(SubseqHashMaxI8Tables *tables,
                                  size_t maximum_length, size_t d,
                                  uint64_t seed);
void subseqhash_max_i8_tables_destroy(SubseqHashMaxI8Tables *tables);

int subseqhash_compute_all_edit_splits_max_i8(
    const uint8_t *s, size_t n, size_t total_edits,
    const SubseqHashMaxI8Tables *tables, SubseqHashEditResult *results);

const char *subseqhash_edits_max_i8_backend(void);

#endif
