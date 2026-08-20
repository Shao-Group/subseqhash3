#include "subsequence_scores.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void brute_visit(const uint8_t *s, size_t n, const size_t *C,
                        size_t k, size_t D, size_t start, size_t chosen,
                        size_t score, uint8_t *z)
{
    size_t position;
    if (chosen == k) {
        z[score] = 1;
        return;
    }
    for (position = start; position + (k - chosen) <= n; ++position) {
        brute_visit(s, n, C, k, D, position + 1, chosen + 1,
                    (score + C[4 * chosen + s[position]]) % D, z);
    }
}

static void brute_smallest_witness(const uint8_t *s, size_t n,
                                   const size_t *C, size_t k, size_t D,
                                   size_t start, size_t chosen, size_t score,
                                   size_t d_star, uint8_t *candidate,
                                   uint8_t *best, int *found)
{
    size_t position;
    if (chosen == k) {
        if (score == d_star &&
            (!*found || memcmp(candidate, best, k) < 0)) {
            memcpy(best, candidate, k);
            *found = 1;
        }
        return;
    }
    for (position = start; position + (k - chosen) <= n; ++position) {
        candidate[chosen] = s[position];
        brute_smallest_witness(
            s, n, C, k, D, position + 1, chosen + 1,
            (score + C[4 * chosen + s[position]]) % D, d_star,
            candidate, best, found);
    }
}

static void one_random_test(size_t n, size_t k, size_t D,
                            int test_brute_force)
{
    uint8_t *s = malloc(n == 0 ? 1 : n * sizeof(*s));
    size_t *C = malloc(k == 0 ? sizeof(*C) : 4 * k * sizeof(*C));
    uint8_t *actual = malloc(D * sizeof(*actual));
    uint8_t *original = malloc(D * sizeof(*original));
    uint8_t *leftmost = malloc(D * sizeof(*leftmost));
    uint8_t *witness_z = malloc(D * sizeof(*witness_z));
    uint8_t *witness = malloc(k == 0 ? 1 : k * sizeof(*witness));
    size_t *witness_positions =
        malloc(k == 0 ? sizeof(*witness_positions)
                      : k * sizeof(*witness_positions));
    size_t d_star;
    uint8_t *expected = test_brute_force
                            ? calloc(D, sizeof(*expected))
                            : NULL;
    size_t i;

    assert(s != NULL && C != NULL && actual != NULL && original != NULL &&
           leftmost != NULL && witness_z != NULL && witness != NULL &&
           witness_positions != NULL);
    assert(!test_brute_force || expected != NULL);
    for (i = 0; i < n; ++i) s[i] = (uint8_t)(rand() % 4);
    for (i = 0; i < 4 * k; ++i) C[i] = (size_t)rand() % D;
    if (test_brute_force) {
        brute_visit(s, n, C, k, D, 0, 0, 0, expected);
    }

    assert(subsequence_scores(s, n, C, k, D, actual) == 0);
    assert(subsequence_scores_bitset_original(s, n, C, k, D, original) == 0);
    assert(subsequence_scores_leftmost(s, n, C, k, D, leftmost) == 0);
    assert(subsequence_scores_with_witness(s, n, C, k, D, witness_z,
                                           &d_star, witness,
                                           witness_positions) == 0);
    assert(memcmp(actual, leftmost, D) == 0);
    assert(memcmp(actual, original, D) == 0);
    assert(memcmp(actual, witness_z, D) == 0);
    if (test_brute_force) {
        assert(memcmp(actual, expected, D) == 0);
        assert(memcmp(leftmost, expected, D) == 0);
    }
    if (k <= n) {
        size_t witness_score = 0;
        for (i = 0; i < D && !actual[i]; ++i) {
        }
        assert(d_star == i);
        for (i = 0; i < k; ++i) {
            assert(witness_positions[i] < n);
            assert(i == 0 || witness_positions[i - 1] < witness_positions[i]);
            assert(witness[i] == s[witness_positions[i]]);
            witness_score =
                (witness_score + C[4 * i + witness[i]]) % D;
        }
        assert(witness_score == d_star);
        if (test_brute_force) {
            uint8_t *candidate = malloc(k == 0 ? 1 : k);
            uint8_t *best = malloc(k == 0 ? 1 : k);
            int found = 0;
            assert(candidate != NULL && best != NULL);
            brute_smallest_witness(s, n, C, k, D, 0, 0, 0, d_star,
                                   candidate, best, &found);
            assert(found);
            assert(memcmp(witness, best, k) == 0);
            free(best);
            free(candidate);
        }
    } else {
        assert(d_star == D);
    }

	if(test_brute_force == 0)
	{
		printf("n=%zu k=%zu D=%zu s=\"", n, k, D);
		for (i = 0; i < n; ++i) {
			printf("%u", (unsigned)s[i]);
		}
		printf("\" z=[");
		for (i = 0; i < D; ++i) {
			printf("%s%u", i == 0 ? "" : "", (unsigned)actual[i]);
		}
		puts("]");
	}

    free(expected);
    free(witness_positions);
    free(witness);
    free(witness_z);
    free(leftmost);
    free(original);
    free(actual);
    free(C);
    free(s);
}

int main(void)
{
    size_t small_n;
    size_t small_k;
    srand(20260820);

    /* Exercise rotations for many non-word-aligned values of D. */
    for (small_n = 0; small_n <= 7; ++small_n) {
        for (small_k = 0; small_k <= small_n + 1; ++small_k) {
            one_random_test(small_n, small_k, 1 + (size_t)rand() % 137, 1);
        }
    }
    one_random_test(7, 4, 64, 1);
    one_random_test(7, 4, 128, 1);
    one_random_test(7, 4, 192, 1);

    size_t n = 50;
    size_t k = 40;
	size_t D = 5120;

    one_random_test(n, k, D, 0);
    return 0;
}
