#ifndef HASH_COLLISION_ESTIMATOR_H
#define HASH_COLLISION_ESTIMATOR_H
#include <stddef.h>
#include <stdint.h>
typedef struct {
    size_t pairs, valid_pairs, invalid_pairs, collisions;
    double frequency, valid_frequency;
    uint64_t table_seed;
} HashCollisionStats;
int estimate_subseqhash_collisions(const char *filename, size_t k, size_t D,
                                   HashCollisionStats *stats);
int estimate_subseqhash_collisions_original(const char *filename, size_t k,
                                            size_t D,
                                            HashCollisionStats *stats);

/* Backward-compatible name for the SubseqHash estimator. */
int estimate_hash_collisions(const char *filename, size_t k, size_t D,
                             HashCollisionStats *stats);
#endif
