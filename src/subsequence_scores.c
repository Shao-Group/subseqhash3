#include "subsequence_scores.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { WORD_BITS = 64 };

/* OR source shifted left by shift bits into destination. */
static void or_shift_left(uint64_t *destination, const uint64_t *source,
                          size_t words, size_t shift)
{
    const size_t word_shift = shift / WORD_BITS;
    const unsigned bit_shift = (unsigned)(shift % WORD_BITS);
    size_t i;

    for (i = 0; i + word_shift < words; ++i) {
        const size_t output = i + word_shift;
        destination[output] |= source[i] << bit_shift;
        if (bit_shift != 0 && output + 1 < words) {
            destination[output + 1] |= source[i] >> (WORD_BITS - bit_shift);
        }
    }
}

/* OR source shifted right by shift bits into destination. */
static void or_shift_right(uint64_t *destination, const uint64_t *source,
                           size_t words, size_t shift)
{
    const size_t word_shift = shift / WORD_BITS;
    const unsigned bit_shift = (unsigned)(shift % WORD_BITS);
    size_t i;

    for (i = word_shift; i < words; ++i) {
        const size_t output = i - word_shift;
        destination[output] |= source[i] >> bit_shift;
        if (bit_shift != 0 && output > 0) {
            destination[output - 1] |= source[i] << (WORD_BITS - bit_shift);
        }
    }
}

/* OR a cyclic left rotation in a D-bit bitset into destination. */
static void or_rotate_left(uint64_t *destination, const uint64_t *source,
                           size_t words, size_t D, size_t shift,
                           uint64_t last_word_mask)
{
    size_t i;

    if (shift == 0) {
        for (i = 0; i < words; ++i) {
            destination[i] |= source[i];
        }
    } else {
        /* rotl_D(x, shift) = (x << shift) | (x >> (D - shift)). */
        or_shift_left(destination, source, words, shift);
        or_shift_right(destination, source, words, D - shift);
    }
    destination[words - 1] &= last_word_mask;
}

int subsequence_scores_bitset_original(const uint8_t *s, size_t n,
                                       const size_t *C, size_t k, size_t D,
                                       uint8_t *z)
{
    uint64_t *reachable;
    size_t words;
    uint64_t last_word_mask;
    size_t position;

    if (D == 0 || z == NULL || (n != 0 && s == NULL) ||
        (k != 0 && C == NULL) || k > SIZE_MAX / 4) {
        return -1;
    }

    words = D / WORD_BITS + (D % WORD_BITS != 0);
    if (k == SIZE_MAX || words > SIZE_MAX / (k + 1) ||
        (k + 1) * words > SIZE_MAX / sizeof(*reachable)) {
        return -1;
    }
    last_word_mask = D % WORD_BITS == 0
                         ? UINT64_MAX
                         : (UINT64_C(1) << (D % WORD_BITS)) - 1;

    for (position = 0; position < n; ++position) {
        if (s[position] > 3) {
            return -1;
        }
    }
    for (position = 0; position < 4 * k; ++position) {
        if (C[position] >= D) {
            return -1;
        }
    }

    reachable = calloc((k + 1) * words, sizeof(*reachable));
    if (reachable == NULL) {
        return -1;
    }
    reachable[0] = UINT64_C(1); /* Empty subsequence, score zero. */

    for (position = 0; position < n; ++position) {
        size_t chosen = position + 1 < k ? position + 1 : k;
        const size_t symbol = s[position];

        /* Descending order prevents reusing s[position]. */
        while (chosen > 0) {
            uint64_t *current = reachable + chosen * words;
            const uint64_t *previous = reachable + (chosen - 1) * words;
            const size_t contribution = C[4 * (chosen - 1) + symbol];

            or_rotate_left(current, previous, words, D, contribution,
                           last_word_mask);
            --chosen;
        }
    }

    for (position = 0; position < D; ++position) {
        const uint64_t word = reachable[k * words + position / WORD_BITS];
        z[position] = (uint8_t)((word >> (position % WORD_BITS)) & 1U);
    }
    free(reachable);
    return 0;
}

/*
 * D is a multiple of 64 here.  Compute a complete cyclic rotation in one
 * destination pass and return whether the destination became all ones.
 */
static int or_rotate_left_aligned(uint64_t *destination,
                                  const uint64_t *source,
                                  size_t words, size_t shift)
{
    const size_t word_shift = shift / WORD_BITS;
    const unsigned bit_shift = (unsigned)(shift % WORD_BITS);
    size_t source_index = word_shift == 0 ? 0 : words - word_shift;
    uint64_t all_full = UINT64_MAX;
    size_t output;

    if (bit_shift == 0) {
        const size_t first_count = words - source_index;
        for (output = 0; output < first_count; ++output) {
            const uint64_t merged =
                destination[output] | source[source_index + output];
            destination[output] = merged;
            all_full &= merged;
        }
        for (; output < words; ++output) {
            const uint64_t merged =
                destination[output] | source[output - first_count];
            destination[output] = merged;
            all_full &= merged;
        }
    } else {
        /* Handle source word zero separately; all other ranges are contiguous
         * and friendly to compiler SIMD vectorization. */
        output = 0;
        if (source_index != 0) {
            size_t input;
            for (input = source_index; input < words; ++input, ++output) {
                const uint64_t rotated =
                    (source[input] << bit_shift) |
                    (source[input - 1] >> (WORD_BITS - bit_shift));
                const uint64_t merged = destination[output] | rotated;
                destination[output] = merged;
                all_full &= merged;
            }
        }
        {
            const uint64_t rotated =
                (source[0] << bit_shift) |
                (source[words - 1] >> (WORD_BITS - bit_shift));
            const uint64_t merged = destination[output] | rotated;
            destination[output] = merged;
            all_full &= merged;
            ++output;
        }
        {
            size_t input = 1;
            for (; output < words; ++input, ++output) {
                const uint64_t rotated =
                    (source[input] << bit_shift) |
                    (source[input - 1] >> (WORD_BITS - bit_shift));
                const uint64_t merged = destination[output] | rotated;
                destination[output] = merged;
                all_full &= merged;
            }
        }
    }
    return all_full == UINT64_MAX;
}

int subsequence_scores(const uint8_t *s, size_t n,
                       const size_t *C, size_t k, size_t D,
                       uint8_t *z)
{
    uint64_t *reachable;
    uint8_t *full;
    size_t words;
    size_t position;

    if (D == 0 || z == NULL || (n != 0 && s == NULL) ||
        (k != 0 && C == NULL) || k > SIZE_MAX / 4) {
        return -1;
    }
    for (position = 0; position < n; ++position) {
        if (s[position] > 3) return -1;
    }
    for (position = 0; position < 4 * k; ++position) {
        if (C[position] >= D) return -1;
    }
    if (k > n) {
        memset(z, 0, D);
        return 0;
    }

    /* The original kernel handles the uncommon partial-word modulus case. */
    if (D % WORD_BITS != 0) {
        return subsequence_scores_bitset_original(s, n, C, k, D, z);
    }

    words = D / WORD_BITS;
    if (k == SIZE_MAX || words > SIZE_MAX / (k + 1) ||
        (k + 1) * words > SIZE_MAX / sizeof(*reachable)) {
        return -1;
    }
    reachable = calloc((k + 1) * words, sizeof(*reachable));
    full = calloc(k + 1, sizeof(*full));
    if (reachable == NULL || full == NULL) {
        free(full);
        free(reachable);
        return -1;
    }
    reachable[0] = UINT64_C(1);

    for (position = 1; position <= n; ++position) {
        const size_t symbol = s[position - 1];
        const size_t upper = position < k ? position : k;
        const size_t remaining = n - position;
        const size_t lower = k > remaining ? k - remaining : 1;
        size_t chosen;

        if (lower > upper) continue;
        chosen = upper;
        for (;;) {
            if (!full[chosen]) {
                uint64_t *current = reachable + chosen * words;
                const uint64_t *previous = reachable + (chosen - 1) * words;
                if (full[chosen - 1]) {
                    memset(current, 0xff, words * sizeof(*current));
                    full[chosen] = 1;
                } else {
                    const size_t contribution = C[4 * (chosen - 1) + symbol];
                    full[chosen] = (uint8_t)or_rotate_left_aligned(
                        current, previous, words, contribution);
                }
            }
            if (chosen == lower) break;
            --chosen;
        }
    }

    {
        const uint64_t *answer = reachable + k * words;
        size_t word;
        for (word = 0; word < words; ++word) {
            const uint64_t bits = answer[word];
            uint8_t *output = z + word * WORD_BITS;
            if (bits == 0) {
                memset(output, 0, WORD_BITS);
            } else if (bits == UINT64_MAX) {
                memset(output, 1, WORD_BITS);
            } else {
                unsigned bit;
                for (bit = 0; bit < WORD_BITS; ++bit) {
                    output[bit] = (uint8_t)((bits >> bit) & 1U);
                }
            }
        }
    }
    free(full);
    free(reachable);
    return 0;
}

static unsigned trailing_zeros_u64(uint64_t value)
{
#if defined(__GNUC__) || defined(__clang__)
    return (unsigned)__builtin_ctzll(value);
#else
    unsigned count = 0;
    while ((value & UINT64_C(1)) == 0) {
        value >>= 1;
        ++count;
    }
    return count;
#endif
}

static int subsequence_scores_with_witness_parents(
    const uint8_t *s, size_t n, const size_t *C, size_t k, size_t D,
    uint8_t *z, size_t *d_star, uint8_t *y, size_t *positions)
{
    uint64_t *reachable;
    uint64_t *rotated;
    uint32_t *parent_position;
    size_t words;
    uint64_t last_word_mask;
    size_t position;

    if (D == 0 || z == NULL || d_star == NULL ||
        (n != 0 && s == NULL) || (k != 0 && (C == NULL || y == NULL)) ||
        n >= UINT32_MAX || k > SIZE_MAX / 4) {
        return -1;
    }
    for (position = 0; position < n; ++position) {
        if (s[position] > 3) return -1;
    }
    for (position = 0; position < 4 * k; ++position) {
        if (C[position] >= D) return -1;
    }
    if (k > n) {
        memset(z, 0, D);
        *d_star = D;
        return 0;
    }

    words = D / WORD_BITS + (D % WORD_BITS != 0);
    if (k == SIZE_MAX || words > SIZE_MAX / (k + 1) ||
        (k + 1) * words > SIZE_MAX / sizeof(*reachable) ||
        (k != 0 && D > SIZE_MAX / k) ||
        k * D > SIZE_MAX / sizeof(*parent_position)) {
        return -1;
    }
    last_word_mask = D % WORD_BITS == 0
                         ? UINT64_MAX
                         : (UINT64_C(1) << (D % WORD_BITS)) - 1;
    reachable = calloc((k + 1) * words, sizeof(*reachable));
    rotated = malloc(words * sizeof(*rotated));
    parent_position = calloc(k == 0 ? 1 : k * D,
                             sizeof(*parent_position));
    if (reachable == NULL || rotated == NULL || parent_position == NULL) {
        free(parent_position);
        free(rotated);
        free(reachable);
        return -1;
    }
    reachable[0] = UINT64_C(1);

    for (position = 1; position <= n; ++position) {
        const size_t symbol = s[position - 1];
        const size_t upper = position < k ? position : k;
        const size_t remaining = n - position;
        const size_t lower = k > remaining ? k - remaining : 1;
        size_t chosen;

        if (lower > upper) continue;
        chosen = upper;
        for (;;) {
            uint64_t *current = reachable + chosen * words;
            const uint64_t *previous = reachable + (chosen - 1) * words;
            uint32_t *parents = parent_position + (chosen - 1) * D;
            const size_t contribution = C[4 * (chosen - 1) + symbol];
            size_t word;

            memset(rotated, 0, words * sizeof(*rotated));
            or_rotate_left(rotated, previous, words, D, contribution,
                           last_word_mask);
            for (word = 0; word < words; ++word) {
                uint64_t newly_reached = rotated[word] & ~current[word];
                current[word] |= rotated[word];
                while (newly_reached != 0) {
                    const unsigned bit = trailing_zeros_u64(newly_reached);
                    const size_t score = word * WORD_BITS + bit;
                    if (score < D) parents[score] = (uint32_t)position;
                    newly_reached &= newly_reached - 1;
                }
            }
            if (chosen == lower) break;
            --chosen;
        }
    }

    for (position = 0; position < D; ++position) {
        const uint64_t word = reachable[k * words + position / WORD_BITS];
        z[position] = (uint8_t)((word >> (position % WORD_BITS)) & 1U);
    }
    for (*d_star = 0; *d_star < D && !z[*d_star]; ++*d_star) {
    }

    if (*d_star < D) {
        size_t score = *d_star;
        size_t chosen = k;
        while (chosen > 0) {
            const uint32_t selected_position =
                parent_position[(chosen - 1) * D + score];
            const uint8_t symbol = s[selected_position - 1];
            const size_t contribution = C[4 * (chosen - 1) + symbol];
            y[chosen - 1] = symbol;
            if (positions != NULL) positions[chosen - 1] = selected_position - 1;
            score = score >= contribution
                        ? score - contribution
                        : score + D - contribution;
            --chosen;
        }
    }

    free(parent_position);
    free(rotated);
    free(reachable);
    return 0;
}

static const size_t HISTORY_FULL = SIZE_MAX - 1;

int subsequence_scores_with_witness_history(
    const uint8_t *s, size_t n, const size_t *C, size_t k, size_t D,
    uint8_t *z, size_t *d_star, uint8_t *y, size_t *positions)
{
    uint64_t *reachable;
    uint64_t *history;
    uint8_t *full;
    size_t *history_offset;
    size_t words;
    size_t slots = 0;
    size_t used_slots = 0;
    size_t position;

    /* Keep the general parent implementation as a fallback for partial words. */
    if (D == 0 || D % WORD_BITS != 0) {
        return subsequence_scores_with_witness_parents(
            s, n, C, k, D, z, d_star, y, positions);
    }
    if (z == NULL || d_star == NULL || (n != 0 && s == NULL) ||
        (k != 0 && (C == NULL || y == NULL)) || k > SIZE_MAX / 4 ||
        n == SIZE_MAX || k == SIZE_MAX || n + 1 > SIZE_MAX / (k + 1)) {
        return -1;
    }
    for (position = 0; position < n; ++position) {
        if (s[position] > 3) return -1;
    }
    for (position = 0; position < 4 * k; ++position) {
        if (C[position] >= D) return -1;
    }
    if (k > n) {
        memset(z, 0, D);
        *d_star = D;
        return 0;
    }

    words = D / WORD_BITS;
    for (position = 0; position <= n; ++position) {
        const size_t lower = k > n - position ? k - (n - position) : 0;
        const size_t upper = position < k ? position : k;
        if (lower <= upper) {
            const size_t count = upper - lower + 1;
            if (slots > SIZE_MAX - count) return -1;
            slots += count;
        }
    }
    if (words > SIZE_MAX / (k + 1) ||
        (k + 1) * words > SIZE_MAX / sizeof(*reachable) ||
        slots > SIZE_MAX / words ||
        slots * words > SIZE_MAX / sizeof(*history)) {
        return -1;
    }

    reachable = calloc((k + 1) * words, sizeof(*reachable));
    history = malloc(slots * words * sizeof(*history));
    full = calloc(k + 1, sizeof(*full));
    history_offset = malloc((n + 1) * (k + 1) * sizeof(*history_offset));
    if (reachable == NULL || history == NULL || full == NULL ||
        history_offset == NULL) {
        free(history_offset);
        free(full);
        free(history);
        free(reachable);
        return -1;
    }
    for (position = 0; position < (n + 1) * (k + 1); ++position) {
        history_offset[position] = SIZE_MAX;
    }
    reachable[0] = UINT64_C(1);
    history_offset[0] = used_slots;
    memcpy(history + used_slots++ * words, reachable,
           words * sizeof(*history));

    for (position = 1; position <= n; ++position) {
        const size_t symbol = s[position - 1];
        const size_t upper = position < k ? position : k;
        const size_t remaining = n - position;
        const size_t lower = k > remaining ? k - remaining : 1;
        size_t chosen;

        if (lower <= upper) {
            chosen = upper;
            for (;;) {
                if (!full[chosen]) {
                    uint64_t *current = reachable + chosen * words;
                    const uint64_t *previous =
                        reachable + (chosen - 1) * words;
                    if (full[chosen - 1]) {
                        memset(current, 0xff, words * sizeof(*current));
                        full[chosen] = 1;
                    } else {
                        full[chosen] = (uint8_t)or_rotate_left_aligned(
                            current, previous, words,
                            C[4 * (chosen - 1) + symbol]);
                    }
                }
                if (chosen == lower) break;
                --chosen;
            }
        }

        {
            size_t stored;
            for (stored = lower; stored <= upper; ++stored) {
                const size_t index = position * (k + 1) + stored;
                if (full[stored]) {
                    history_offset[index] = HISTORY_FULL;
                } else {
                    history_offset[index] = used_slots;
                    memcpy(history + used_slots++ * words,
                           reachable + stored * words,
                           words * sizeof(*history));
                }
            }
        }
    }

    {
        const uint64_t *answer = reachable + k * words;
        size_t word;
        for (word = 0; word < words; ++word) {
            const uint64_t bits = answer[word];
            uint8_t *output = z + word * WORD_BITS;
            if (bits == 0) memset(output, 0, WORD_BITS);
            else if (bits == UINT64_MAX) memset(output, 1, WORD_BITS);
            else {
                unsigned bit;
                for (bit = 0; bit < WORD_BITS; ++bit) {
                    output[bit] = (uint8_t)((bits >> bit) & 1U);
                }
            }
        }
    }
    for (*d_star = 0; *d_star < D && !z[*d_star]; ++*d_star) {
    }

    if (*d_star < D) {
        size_t score = *d_star;
        size_t chosen = k;
        position = n;
        while (chosen > 0) {
            const size_t skip_index = (position - 1) * (k + 1) + chosen;
            const size_t offset = history_offset[skip_index];
            int can_skip = offset == HISTORY_FULL;
            if (offset != SIZE_MAX && offset != HISTORY_FULL) {
                can_skip = (int)((history[offset * words + score / WORD_BITS] >>
                                  (score % WORD_BITS)) & 1U);
            }
            if (can_skip) {
                --position;
            } else {
                const uint8_t symbol = s[position - 1];
                const size_t contribution = C[4 * (chosen - 1) + symbol];
                y[chosen - 1] = symbol;
                if (positions != NULL) positions[chosen - 1] = position - 1;
                score = score >= contribution
                            ? score - contribution
                            : score + D - contribution;
                --chosen;
                --position;
            }
        }
    }

    free(history_offset);
    free(full);
    free(history);
    free(reachable);
    return 0;
}

int subsequence_scores_with_witness(const uint8_t *s, size_t n,
                                    const size_t *C, size_t k, size_t D,
                                    uint8_t *z, size_t *d_star,
                                    uint8_t *y, size_t *positions)
{
    uint64_t *suffix;
    uint64_t *history;
    uint64_t *rotated = NULL;
    uint8_t *full;
    size_t *history_offset;
    size_t words;
    size_t slots = 0;
    size_t used_slots = 0;
    uint64_t last_word_mask;
    size_t p;

    if (D == 0 || z == NULL || d_star == NULL ||
        (n != 0 && s == NULL) || (k != 0 && (C == NULL || y == NULL)) ||
        n == SIZE_MAX || k == SIZE_MAX || k > SIZE_MAX / 4 ||
        n + 1 > SIZE_MAX / (k + 1)) {
        return -1;
    }
    for (p = 0; p < n; ++p) {
        if (s[p] > 3) return -1;
    }
    for (p = 0; p < 4 * k; ++p) {
        if (C[p] >= D) return -1;
    }
    if (k > n) {
        memset(z, 0, D);
        *d_star = D;
        return 0;
    }
    if (k == 0) {
        memset(z, 0, D);
        z[0] = 1;
        *d_star = 0;
        return 0;
    }

    words = D / WORD_BITS + (D % WORD_BITS != 0);
    last_word_mask = D % WORD_BITS == 0
                         ? UINT64_MAX
                         : (UINT64_C(1) << (D % WORD_BITS)) - 1;
    for (p = 0; p <= n; ++p) {
        const size_t lower = k > n - p ? k - (n - p) : 0;
        const size_t upper = p < k ? p : k;
        if (lower <= upper) {
            const size_t count = upper - lower + 1;
            if (slots > SIZE_MAX - count) return -1;
            slots += count;
        }
    }
    if (words > SIZE_MAX / (k + 1) ||
        (k + 1) * words > SIZE_MAX / sizeof(*suffix) ||
        slots > SIZE_MAX / words ||
        slots * words > SIZE_MAX / sizeof(*history)) {
        return -1;
    }

    suffix = calloc((k + 1) * words, sizeof(*suffix));
    history = malloc(slots * words * sizeof(*history));
    full = calloc(k + 1, sizeof(*full));
    history_offset = malloc((n + 1) * (k + 1) * sizeof(*history_offset));
    if (D % WORD_BITS != 0) rotated = malloc(words * sizeof(*rotated));
    if (suffix == NULL || history == NULL || full == NULL ||
        history_offset == NULL || (D % WORD_BITS != 0 && rotated == NULL)) {
        free(rotated);
        free(history_offset);
        free(full);
        free(history);
        free(suffix);
        return -1;
    }
    for (p = 0; p < (n + 1) * (k + 1); ++p) {
        history_offset[p] = SIZE_MAX;
    }

    /* At suffix s[n..n), choosing no remaining symbols has score zero. */
    suffix[k * words] = UINT64_C(1);
    history_offset[n * (k + 1) + k] = used_slots;
    memcpy(history + used_slots++ * words, suffix + k * words,
           words * sizeof(*history));

    p = n;
    while (p > 0) {
        const size_t symbol = s[p - 1];
        size_t lower;
        size_t upper;
        size_t chosen;
        --p;
        lower = k > n - p ? k - (n - p) : 0;
        upper = p < k ? p : k;

        /* Ascending order preserves suffix[chosen+1] from before s[p]. */
        for (chosen = lower; chosen <= upper && chosen < k; ++chosen) {
            uint64_t *current = suffix + chosen * words;
            const uint64_t *next = suffix + (chosen + 1) * words;
            const size_t contribution = C[4 * chosen + symbol];
            if (full[chosen]) continue;
            if (full[chosen + 1] && D % WORD_BITS == 0) {
                memset(current, 0xff, words * sizeof(*current));
                full[chosen] = 1;
            } else if (D % WORD_BITS == 0) {
                full[chosen] = (uint8_t)or_rotate_left_aligned(
                    current, next, words, contribution);
            } else {
                size_t word;
                memset(rotated, 0, words * sizeof(*rotated));
                or_rotate_left(rotated, next, words, D, contribution,
                               last_word_mask);
                for (word = 0; word < words; ++word) {
                    current[word] |= rotated[word];
                }
                current[words - 1] &= last_word_mask;
            }
        }

        for (chosen = lower; chosen <= upper; ++chosen) {
            const size_t index = p * (k + 1) + chosen;
            if (full[chosen]) {
                history_offset[index] = HISTORY_FULL;
            } else {
                history_offset[index] = used_slots;
                memcpy(history + used_slots++ * words,
                       suffix + chosen * words,
                       words * sizeof(*history));
            }
        }
    }

    /* suffix[0] is now the complete reachable-score set, so no separate
     * forward Boolean DP is necessary. */
    if (full[0]) {
        memset(z, 1, D);
        *d_star = 0;
    } else {
        const uint64_t *answer = suffix;
        size_t word;
        for (word = 0; word < words; ++word) {
            const uint64_t bits = answer[word];
            const size_t base = word * WORD_BITS;
            const size_t count = D - base < WORD_BITS ? D - base : WORD_BITS;
            uint8_t *output = z + base;
            if (bits == 0) memset(output, 0, count);
            else if (bits == UINT64_MAX && count == WORD_BITS) {
                memset(output, 1, WORD_BITS);
            } else {
                size_t bit;
                for (bit = 0; bit < count; ++bit) {
                    output[bit] = (uint8_t)((bits >> bit) & 1U);
                }
            }
        }
        for (*d_star = 0; *d_star < D && !z[*d_star]; ++*d_star) {
        }
    }

    /* Greedily choose the smallest symbol that admits a feasible completion. */
    {
        size_t chosen;
        size_t start = 0;
        size_t prefix_score = 0;
        for (chosen = 0; chosen < k; ++chosen) {
            size_t symbol;
            int found = 0;
            for (symbol = 0; symbol < 4 && !found; ++symbol) {
                size_t candidate_position;
                for (candidate_position = start;
                     candidate_position < n; ++candidate_position) {
                    size_t new_score;
                    size_t needed;
                    size_t offset;
                    int feasible;
                    if (s[candidate_position] != symbol) continue;
                    new_score = prefix_score + C[4 * chosen + symbol];
                    if (new_score >= D) new_score -= D;
                    needed = *d_star >= new_score
                                 ? *d_star - new_score
                                 : *d_star + D - new_score;
                    offset = history_offset[
                        (candidate_position + 1) * (k + 1) + chosen + 1];
                    feasible = offset == HISTORY_FULL;
                    if (offset != SIZE_MAX && offset != HISTORY_FULL) {
                        feasible = (int)((history[offset * words +
                                                 needed / WORD_BITS] >>
                                          (needed % WORD_BITS)) & 1U);
                    }
                    if (feasible) {
                        y[chosen] = (uint8_t)symbol;
                        if (positions != NULL) positions[chosen] = candidate_position;
                        prefix_score = new_score;
                        start = candidate_position + 1;
                        found = 1;
                        break;
                    }
                }
            }
            if (!found) {
                free(rotated);
                free(history_offset);
                free(full);
                free(history);
                free(suffix);
                return -1;
            }
        }
    }

    free(rotated);
    free(history_offset);
    free(full);
    free(history);
    free(suffix);
    return 0;
}

static int subsequence_scores_leftmost_u16(const uint8_t *s, size_t n,
                                           const size_t *C, size_t k, size_t D,
                                           uint8_t *z)
{
    const uint16_t unreachable = (uint16_t)(n + 1);
    uint16_t *next_position = malloc((n + 2) * 4 * sizeof(*next_position));
    uint16_t *previous = malloc(D * sizeof(*previous));
    uint16_t *current = malloc(D * sizeof(*current));
    size_t *previous_scores = malloc(D * sizeof(*previous_scores));
    size_t *current_scores = malloc(D * sizeof(*current_scores));
    size_t previous_count = 1;
    size_t position;
    size_t chosen;

    if (next_position == NULL || previous == NULL || current == NULL ||
        previous_scores == NULL || current_scores == NULL) {
        free(current_scores);
        free(previous_scores);
        free(current);
        free(previous);
        free(next_position);
        return -1;
    }

    for (position = 0; position < 4; ++position) {
        next_position[4 * n + position] = unreachable;
        next_position[4 * (n + 1) + position] = unreachable;
    }
    position = n;
    while (position > 0) {
        size_t symbol;
        --position;
        for (symbol = 0; symbol < 4; ++symbol) {
            next_position[4 * position + symbol] =
                next_position[4 * (position + 1) + symbol];
        }
        next_position[4 * position + s[position]] = (uint16_t)(position + 1);
    }

    for (position = 0; position < D; ++position) {
        previous[position] = unreachable;
        current[position] = unreachable;
    }
    previous[0] = 0;
    previous_scores[0] = 0;

    /* Sparse startup: avoid D-wide scans while only a few scores exist. */
    for (chosen = 1; chosen <= k && previous_count < D / 4; ++chosen) {
        size_t current_count = 0;
        size_t active;

        for (active = 0; active < previous_count; ++active) {
            const size_t prior_score = previous_scores[active];
            const uint16_t endpoint = previous[prior_score];
            size_t symbol;

            for (symbol = 0; symbol < 4; ++symbol) {
                const uint16_t candidate =
                    next_position[4 * endpoint + symbol];
                size_t score;
                if (candidate > n) {
                    continue;
                }
                score = prior_score + C[4 * (chosen - 1) + symbol];
                if (score >= D) {
                    score -= D;
                }
                if (current[score] == unreachable) {
                    current_scores[current_count++] = score;
                    current[score] = candidate;
                } else if (candidate < current[score]) {
                    current[score] = candidate;
                }
            }
        }

        for (active = 0; active < previous_count; ++active) {
            previous[previous_scores[active]] = unreachable;
        }
        {
            uint16_t *temporary_row = previous;
            size_t *temporary_scores = previous_scores;
            previous = current;
            current = temporary_row;
            previous_scores = current_scores;
            current_scores = temporary_scores;
            previous_count = current_count;
        }
    }

    for (; chosen <= k; ++chosen) {
        const size_t c0 = C[4 * (chosen - 1)];
        const size_t c1 = C[4 * (chosen - 1) + 1];
        const size_t c2 = C[4 * (chosen - 1) + 2];
        const size_t c3 = C[4 * (chosen - 1) + 3];
        size_t score;

        for (score = 0; score < D; ++score) {
            const size_t p0 = score >= c0 ? score - c0 : score + D - c0;
            const size_t p1 = score >= c1 ? score - c1 : score + D - c1;
            const size_t p2 = score >= c2 ? score - c2 : score + D - c2;
            const size_t p3 = score >= c3 ? score - c3 : score + D - c3;
            uint16_t best = next_position[4 * previous[p0]];
            uint16_t candidate = next_position[4 * previous[p1] + 1];

            if (candidate < best) best = candidate;
            candidate = next_position[4 * previous[p2] + 2];
            if (candidate < best) best = candidate;
            candidate = next_position[4 * previous[p3] + 3];
            if (candidate < best) best = candidate;
            current[score] = best;
        }
        {
            uint16_t *temporary = previous;
            previous = current;
            current = temporary;
        }
    }

    for (position = 0; position < D; ++position) {
        z[position] = (uint8_t)(previous[position] <= n);
    }
    free(current_scores);
    free(previous_scores);
    free(current);
    free(previous);
    free(next_position);
    return 0;
}

int subsequence_scores_leftmost(const uint8_t *s, size_t n,
                                const size_t *C, size_t k, size_t D,
                                uint8_t *z)
{
    size_t *next_position;
    size_t *previous;
    size_t *current;
    size_t position;
    size_t chosen;

    if (D == 0 || z == NULL || (n != 0 && s == NULL) ||
        (k != 0 && C == NULL) || n > SIZE_MAX - 2 ||
        n + 2 > SIZE_MAX / (4 * sizeof(*next_position)) ||
        D > SIZE_MAX / sizeof(*previous) || k > SIZE_MAX / 4) {
        return -1;
    }
    for (position = 0; position < n; ++position) {
        if (s[position] > 3) {
            return -1;
        }
    }
    for (position = 0; position < 4 * k; ++position) {
        if (C[position] >= D) {
            return -1;
        }
    }
    if (k > n) {
        memset(z, 0, D * sizeof(*z));
        return 0;
    }
    if (n < UINT16_MAX) {
        return subsequence_scores_leftmost_u16(s, n, C, k, D, z);
    }

    /*
     * next_position[a*(n+2) + p] is the first 1-based occurrence of symbol a
     * strictly after prefix endpoint p.  Index n+1 is a sentinel that maps to
     * itself, allowing the dense DP loop to avoid a reachability branch.
     */
    next_position = malloc((n + 2) * 4 * sizeof(*next_position));
    previous = malloc(D * sizeof(*previous));
    current = malloc(D * sizeof(*current));
    if (next_position == NULL || previous == NULL || current == NULL) {
        free(current);
        free(previous);
        free(next_position);
        return -1;
    }

    {
        const size_t stride = n + 2;
        size_t symbol;
        for (symbol = 0; symbol < 4; ++symbol) {
            size_t *next_for_symbol = next_position + symbol * stride;
            next_for_symbol[n] = n + 1;
            next_for_symbol[n + 1] = n + 1;
            position = n;
            while (position > 0) {
                --position;
                next_for_symbol[position] = next_for_symbol[position + 1];
                if (s[position] == symbol) {
                    next_for_symbol[position] = position + 1;
                }
            }
        }
    }

    for (position = 0; position < D; ++position) {
        previous[position] = n + 1;
    }
    previous[0] = 0; /* The empty subsequence ends before s[1]. */

    for (chosen = 1; chosen <= k; ++chosen) {
        const size_t stride = n + 2;
        size_t symbol;

        /*
         * For each symbol, rotate the previous score row by its contribution.
         * Splitting at D-contribution removes modulo from the hot loop.  The
         * first symbol initializes every current entry; the other three only
         * take minima, so no separate O(D) clearing pass is needed.
         */
        for (symbol = 0; symbol < 4; ++symbol) {
            const size_t contribution = C[4 * (chosen - 1) + symbol];
            const size_t split = D - contribution;
            const size_t *next_for_symbol = next_position + symbol * stride;
            size_t prior_score;

            if (symbol == 0) {
                for (prior_score = 0; prior_score < split; ++prior_score) {
                    current[prior_score + contribution] =
                        next_for_symbol[previous[prior_score]];
                }
                for (; prior_score < D; ++prior_score) {
                    current[prior_score - split] =
                        next_for_symbol[previous[prior_score]];
                }
            } else {
                for (prior_score = 0; prior_score < split; ++prior_score) {
                    const size_t candidate =
                        next_for_symbol[previous[prior_score]];
                    size_t *destination = current + prior_score + contribution;
                    if (candidate < *destination) {
                        *destination = candidate;
                    }
                }
                for (; prior_score < D; ++prior_score) {
                    const size_t candidate =
                        next_for_symbol[previous[prior_score]];
                    size_t *destination = current + prior_score - split;
                    if (candidate < *destination) {
                        *destination = candidate;
                    }
                }
            }
        }
        {
            size_t *temporary = previous;
            previous = current;
            current = temporary;
        }
    }

    for (position = 0; position < D; ++position) {
        z[position] = (uint8_t)(previous[position] <= n);
    }

    free(current);
    free(previous);
    free(next_position);
    return 0;
}
