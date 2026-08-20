#ifndef SUBSEQUENCE_SCORES_H
#define SUBSEQUENCE_SCORES_H

#include <stddef.h>
#include <stdint.h>

/*
 * Set z[d] to 1 exactly when some length-k subsequence of s has score d.
 *
 * C is a row-major k-by-4 table: C[4*i + symbol].  The caller allocates
 * z with D bytes.  Returns 0 on success and -1 for invalid input or an
 * allocation failure.
 */
int subsequence_scores(const uint8_t *s, size_t n,
                       const size_t *C, size_t k, size_t D,
                       uint8_t *z);

/* Original packed-bit implementation, retained for benchmarking. */
int subsequence_scores_bitset_original(const uint8_t *s, size_t n,
                                       const size_t *C, size_t k, size_t D,
                                       uint8_t *z);

/*
 * Compute z and a witness for its smallest reachable score.  On success,
 * *d_star is the smallest d with z[d] == 1, y contains its k symbols, and
 * positions (when non-NULL) contains their zero-based positions in s.
 */
int subsequence_scores_with_witness(const uint8_t *s, size_t n,
                                    const size_t *C, size_t k, size_t D,
                                    uint8_t *z, size_t *d_star,
                                    uint8_t *y, size_t *positions);

/*
 * Alternative DP using the earliest ending position for each (length, score)
 * state.  It has O(n*4 + k*D*4) time and O(n*4 + D) extra space.
 */
int subsequence_scores_leftmost(const uint8_t *s, size_t n,
                                const size_t *C, size_t k, size_t D,
                                uint8_t *z);

#endif
