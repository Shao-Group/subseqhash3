#include "subsequence_scores.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum {
    OUTPUT_N = 50,
    OUTPUT_K = 40,
    OUTPUT_D = 512000,
    OUTPUT_CASES = 1000
};

int main(void)
{
    uint8_t *strings = malloc((size_t)OUTPUT_CASES * OUTPUT_N);
    size_t *C = malloc((size_t)OUTPUT_K * 4 * sizeof(*C));
    uint8_t *z = malloc(OUTPUT_D);
    FILE *output;
    size_t test;
    size_t i;

    if (strings == NULL || C == NULL || z == NULL) {
        fputs("allocation failed\n", stderr);
        free(z);
        free(C);
        free(strings);
        return EXIT_FAILURE;
    }

    srand(20260820);
    for (i = 0; i < (size_t)OUTPUT_CASES * OUTPUT_N; ++i) {
        strings[i] = (uint8_t)(rand() % 4);
    }
    for (i = 0; i < (size_t)OUTPUT_K * 4; ++i) {
        C[i] = (size_t)rand() % OUTPUT_D;
    }

    output = fopen("data/witnesses.txt", "w");
    if (output == NULL) {
        perror("data/witnesses.txt");
        free(z);
        free(C);
        free(strings);
        return EXIT_FAILURE;
    }

    fprintf(output, "# cases=%d n=%d k=%d D=%d seed=20260820\n",
            OUTPUT_CASES, OUTPUT_N, OUTPUT_K, OUTPUT_D);
    for (test = 0; test < OUTPUT_CASES; ++test) {
        const uint8_t *s = strings + test * OUTPUT_N;
        uint8_t witness[OUTPUT_K];
        size_t positions[OUTPUT_K];
        size_t d_star;

        if (subsequence_scores_with_witness(
                s, OUTPUT_N, C, OUTPUT_K, OUTPUT_D, z, &d_star,
                witness, positions) != 0) {
            fprintf(stderr, "witness computation failed for case %zu\n", test);
            fclose(output);
            free(z);
            free(C);
            free(strings);
            return EXIT_FAILURE;
        }

        {
            size_t score = 0;
            for (i = 0; i < OUTPUT_K; ++i) {
                if (positions[i] >= OUTPUT_N ||
                    (i != 0 && positions[i - 1] >= positions[i]) ||
                    witness[i] != s[positions[i]]) {
                    fprintf(stderr, "invalid witness for case %zu\n", test);
                    fclose(output);
                    free(z);
                    free(C);
                    free(strings);
                    return EXIT_FAILURE;
                }
                score = (score + C[4 * i + witness[i]]) % OUTPUT_D;
            }
            if (score != d_star) {
                fprintf(stderr, "invalid witness score for case %zu\n", test);
                fclose(output);
                free(z);
                free(C);
                free(strings);
                return EXIT_FAILURE;
            }
        }

        fprintf(output, "case=%04zu s=", test);
        for (i = 0; i < OUTPUT_N; ++i) fprintf(output, "%u", (unsigned)s[i]);
        fprintf(output, " d*=%zu y=", d_star);
        for (i = 0; i < OUTPUT_K; ++i) {
            fprintf(output, "%u", (unsigned)witness[i]);
        }
        fputs(" positions=[", output);
        for (i = 0; i < OUTPUT_K; ++i) {
            fprintf(output, "%s%zu", i == 0 ? "" : ",", positions[i]);
        }
        fputs("]\n", output);
    }

    if (fclose(output) != 0) {
        perror("data/witnesses.txt");
        free(z);
        free(C);
        free(strings);
        return EXIT_FAILURE;
    }
    free(z);
    free(C);
    free(strings);
    return EXIT_SUCCESS;
}
