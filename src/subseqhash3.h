#ifndef SUBSEQHASH3_H
#define SUBSEQHASH3_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t deletions, substitutions, length, psi;
    double omega;
    uint8_t *seed;
    size_t *positions;
    uint8_t *substituted;
} SubseqHash3EditResult;

typedef struct {
    size_t k, d;
    int8_t *A;
    int8_t *B1, *B2;
    int8_t *signed_A;
    uint32_t *B1_masks;
    size_t blocks;
    uint16_t *wrap_indices;
    size_t *C;
    /* Function d occupies f_stride entries.  F[d*f_stride+x] is f_d(x). */
    float *F;
    size_t f_stride;
} SubseqHash3Tables;

int subseqhash3_tables_init(SubseqHash3Tables *tables,
                             size_t maximum_length,size_t residues,
                             uint64_t seed);
void subseqhash3_tables_destroy(SubseqHash3Tables *tables);

int subseqhash3_compute_all_edit_splits(
    const uint8_t *sequence,size_t length,size_t total_edits,
    const SubseqHash3Tables *tables,SubseqHash3EditResult *results);

/* Results are window-major and then substitution count 0..E. Positions are
   absolute offsets in sequence. */
int subseqhash3_compute_long_edit_splits(
    const uint8_t *sequence,size_t sequence_length,size_t window_length,
    size_t total_edits,const SubseqHash3Tables *tables,
    SubseqHash3EditResult *results);

void subseqhash3_edit_results_destroy(SubseqHash3EditResult *results,
                                      size_t count);
const char *subseqhash3_backend(void);

#endif
