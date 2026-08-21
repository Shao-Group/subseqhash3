#ifndef SUBSEQHASH_H
#define SUBSEQHASH_H
#include <stddef.h>
#include <stdint.h>
typedef struct {
    size_t k, d;
    double *A;
    int8_t *B1, *B2;
    size_t *C;
} SubseqHashTables;
int subseqhash_tables_init(SubseqHashTables *tables, size_t k, size_t d,
                           uint64_t seed);
void subseqhash_tables_destroy(SubseqHashTables *tables);
int subseqhash_compute(const uint8_t *s, size_t n,
                       const SubseqHashTables *tables,
                       uint8_t *seed, size_t *positions,
                       size_t *psi, double *omega);
/*
 * Same result as subseqhash_compute(), using rolling score planes and a
 * compact traceback matrix.  The original implementation is retained above
 * as a straightforward reference.
 */
int subseqhash_compute_optimized(const uint8_t *s, size_t n,
                                 const SubseqHashTables *tables,
                                 uint8_t *seed, size_t *positions,
                                 size_t *psi, double *omega);
/*
 * Leftmost-state variant: DP[j,residue] stores the earliest ending position.
 * It minimizes psi exactly, chooses a witness ending as early as possible,
 * and then evaluates that witness using the complete ABC recurrence.
 * Unlike the two functions above, it does not maximize |omega| among all
 * witnesses with the minimum psi.
 */
int subseqhash_compute_leftmost(const uint8_t *s, size_t n,
                                const SubseqHashTables *tables,
                                uint8_t *seed, size_t *positions,
                                size_t *psi, double *omega);
int subseqhash_score(const uint8_t *seed, const SubseqHashTables *tables,
                     size_t *psi, double *omega);
#endif
