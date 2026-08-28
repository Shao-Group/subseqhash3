#ifndef SUBSEQHASH_EDITS_H
#define SUBSEQHASH_EDITS_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t deletions, substitutions, length, psi;
    double omega;
    uint8_t *seed;
    size_t *positions;
    uint8_t *substituted;
} SubseqHashEditResult;

void subseqhash_edit_results_destroy(SubseqHashEditResult *results,
                                     size_t count);

#endif
