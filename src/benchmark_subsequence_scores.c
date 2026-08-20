#include "subsequence_scores.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

enum {
    BENCHMARK_N = 50,
    BENCHMARK_K = 40,
    BENCHMARK_D = 512000,
    BENCHMARK_CASES = 1000
};

static double now_seconds(void)
{
    struct timeval time;
    assert(gettimeofday(&time, NULL) == 0);
    return (double)time.tv_sec + (double)time.tv_usec / 1000000.0;
}

int main(void)
{
    uint8_t *strings;
    size_t *C;
    uint8_t *bitset_results;
    uint8_t *original_results;
    uint8_t *leftmost_results;
    double start;
    double bitset_seconds;
    double original_seconds;
    double leftmost_seconds;
    double witness_seconds;
    uint8_t witness[BENCHMARK_K];
    size_t witness_positions[BENCHMARK_K];
    size_t d_star;
    uint64_t total_ones = 0;
    size_t minimum_ones = BENCHMARK_D + 1;
    size_t maximum_ones = 0;
    size_t test;
    size_t i;

    strings = malloc((size_t)BENCHMARK_CASES * BENCHMARK_N);
    C = malloc((size_t)BENCHMARK_K * 4 * sizeof(*C));
    bitset_results = malloc(BENCHMARK_D);
    original_results = malloc(BENCHMARK_D);
    leftmost_results = malloc(BENCHMARK_D);
    if (strings == NULL || C == NULL || bitset_results == NULL ||
        original_results == NULL || leftmost_results == NULL) {
        fputs("benchmark allocation failed\n", stderr);
        free(leftmost_results);
        free(bitset_results);
        free(original_results);
        free(C);
        free(strings);
        return EXIT_FAILURE;
    }

    srand(20260820);
    for (i = 0; i < (size_t)BENCHMARK_CASES * BENCHMARK_N; ++i) {
        strings[i] = (uint8_t)(rand() % 4);
    }
    for (i = 0; i < (size_t)BENCHMARK_K * 4; ++i) {
        C[i] = (size_t)rand() % BENCHMARK_D;
    }

    start = now_seconds();
    for (test = 0; test < BENCHMARK_CASES; ++test) {
        assert(subsequence_scores(strings + test * BENCHMARK_N,
                                  BENCHMARK_N, C, BENCHMARK_K, BENCHMARK_D,
                                  bitset_results) == 0);
    }
    bitset_seconds = now_seconds() - start;

    start = now_seconds();
    for (test = 0; test < BENCHMARK_CASES; ++test) {
        assert(subsequence_scores_bitset_original(
                   strings + test * BENCHMARK_N,
                   BENCHMARK_N, C, BENCHMARK_K, BENCHMARK_D,
                   original_results) == 0);
    }
    original_seconds = now_seconds() - start;

    start = now_seconds();
    for (test = 0; test < BENCHMARK_CASES; ++test) {
        assert(subsequence_scores_leftmost(
                   strings + test * BENCHMARK_N,
                   BENCHMARK_N, C, BENCHMARK_K, BENCHMARK_D,
                   leftmost_results) == 0);
    }
    leftmost_seconds = now_seconds() - start;

    start = now_seconds();
    for (test = 0; test < BENCHMARK_CASES; ++test) {
        assert(subsequence_scores_with_witness(
                   strings + test * BENCHMARK_N,
                   BENCHMARK_N, C, BENCHMARK_K, BENCHMARK_D,
                   bitset_results, &d_star, witness,
                   witness_positions) == 0);
    }
    witness_seconds = now_seconds() - start;

    for (test = 0; test < BENCHMARK_CASES; ++test) {
        size_t ones = 0;
        assert(subsequence_scores(strings + test * BENCHMARK_N,
                                  BENCHMARK_N, C, BENCHMARK_K, BENCHMARK_D,
                                  bitset_results) == 0);
        for (i = 0; i < BENCHMARK_D; ++i) {
            ones += bitset_results[i];
        }
        total_ones += ones;
        if (ones < minimum_ones) minimum_ones = ones;
        if (ones > maximum_ones) maximum_ones = ones;
    }

    /* Exact verification outside the measured intervals. */
    for (test = 0; test < 10; ++test) {
        assert(subsequence_scores(strings + test * BENCHMARK_N,
                                  BENCHMARK_N, C, BENCHMARK_K, BENCHMARK_D,
                                  bitset_results) == 0);
        assert(subsequence_scores_bitset_original(
                   strings + test * BENCHMARK_N,
                   BENCHMARK_N, C, BENCHMARK_K, BENCHMARK_D,
                   original_results) == 0);
        assert(subsequence_scores_leftmost(
                   strings + test * BENCHMARK_N,
                   BENCHMARK_N, C, BENCHMARK_K, BENCHMARK_D,
                   leftmost_results) == 0);
        if (memcmp(bitset_results, original_results, BENCHMARK_D) != 0 ||
            memcmp(bitset_results, leftmost_results, BENCHMARK_D) != 0) {
            fputs("error: DP implementations produced different results\n",
                  stderr);
            free(leftmost_results);
            free(bitset_results);
            free(original_results);
            free(C);
            free(strings);
            return EXIT_FAILURE;
        }
    }

    printf("parameters: strings=%d n=%d k=%d D=%d\n",
           BENCHMARK_CASES, BENCHMARK_N, BENCHMARK_K, BENCHMARK_D);
    printf("optimized bitset DP: %.6f s total, %.3f us/string\n",
           bitset_seconds,
           bitset_seconds * 1000000.0 / BENCHMARK_CASES);
    printf("original bitset DP:  %.6f s total, %.3f us/string\n",
           original_seconds,
           original_seconds * 1000000.0 / BENCHMARK_CASES);
    printf("leftmost DP: %.6f s total, %.3f us/string\n",
           leftmost_seconds,
           leftmost_seconds * 1000000.0 / BENCHMARK_CASES);
    printf("bitset DP with witness: %.6f s total, %.3f us/string\n",
           witness_seconds,
           witness_seconds * 1000000.0 / BENCHMARK_CASES);
    printf("bitset speedup: %.2fx\n", original_seconds / bitset_seconds);
    printf("ratio: leftmost DP is %.2fx the optimized bitset DP runtime\n",
           leftmost_seconds / bitset_seconds);
    printf("witness overhead: %.2fx the optimized bitset DP runtime\n",
           witness_seconds / bitset_seconds);
    printf("z ones: %.6f%% aggregate, %.2f average, %zu min, %zu max\n",
           100.0 * (double)total_ones /
               ((double)BENCHMARK_CASES * BENCHMARK_D),
           (double)total_ones / BENCHMARK_CASES,
           minimum_ones, maximum_ones);
    puts("outputs match for 10 exactly verified strings");

    free(leftmost_results);
    free(bitset_results);
    free(original_results);
    free(C);
    free(strings);
    return EXIT_SUCCESS;
}
