# Repeats traversal comparison across residue counts

This repeats the [inward/outward repeats experiment](REPEATS_OUTWARD.md) at
D=32, 128, 1024 and 2048. Both algorithms retain n-1 bases, with exactly
one deletion and E substitutions. Configurations are (n,E,R)=(30,4,10),
(30,4,29), (50,8,10), and (50,8,49). R=29/49 uses all pivots of one family.

## Reproduce

```sh
make BUILD_DIR=build/repeats-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  build/repeats-fast/bench-repeats-traversal
python3 bench/run_repeats_residue_sweep.py \
  --executable build/repeats-fast/bench-repeats-traversal --jobs 8
```

The runner uses only Python's standard library. Eight independent collision
workers run concurrently; each seeding call is single-threaded. Once every
collision worker exits, runtime measurements run sequentially. `--jobs 1`
serializes the whole experiment. Temporary batch files stay under `build/`.
The combined results and machine/parameter metadata are saved under
`bench/repeats-residue-*`.

For individual configurations, the C benchmark accepts:

```sh
build/repeats-fast/bench-repeats-traversal --runtime 128
# FAMILIES, PAIRS_PER_FAMILY, RAW_CSV, D, optional CONFIG and PAIR_OFFSET
build/repeats-fast/bench-repeats-traversal --collisions 32 64 /tmp/d128.csv 128
```

CONFIG 0..3 selects the four configurations in the order listed above.
Offsets support disjoint chunks of the canonical pair stream. Omitting D
preserves the previous benchmark defaults.

## Protocol and interpretation

Runtime is the median of seven samples, each containing 48 window calls on
16 deterministic random inputs. Call order alternates. Both traversals get
identical inputs, table arrays, E, and R within a configuration. Workspaces
are warmed; result allocation, traceback and cleanup are timed, while table
generation and initial workspace allocation are excluded. Scores are
recomputed for each window; overlapping-window sharing is not measured.

Each of the 80 collision conditions has 32 table families and 64 pairs per
family (2,048 trials). The exact same input-pair stream is used across all
four D values and both pivot counts. This stream matches the earlier
D=1024 experiment. Table-generation seed identifiers are shared across D;
the realized tables differ with D. Both traversals use the same realized
tables at a given D.

The five mutation conditions per configuration are E/2 substitutions,
E substitutions, one deletion plus one insertion, two deletions plus two
insertions, and unrelated independent strings. Original DNA is uniform
A/C/G/T. Substitutions affect distinct positions and change their symbols.
Deletion/insertion counts describe generation operations; cancellations
and identical pairs are retained and counted in `nonidentical_pairs`.

The primary probability is at least one identical selected seed at the
same pivot. CSV `any` fields also allow cross-pivot matches, and `per`
fields report individual-pivot match fractions. Seed strings are compared
exactly; numerical hash-score coincidences do not count. Probabilities are
percentages. CSV deltas and intervals are inward minus outward, in
percentage points.

Paired approximate 95% intervals use per-table-family differences:
mean +/- 2.04 * sample_standard_deviation / sqrt(32). Batch counts are merged
within each family before computing intervals. These exploratory intervals
are unadjusted for the multiple conditions. The experiment estimates
sensitivity on these synthetic models, not a universal ranking or a rate
on real reads. Zero observed matches do not establish zero probability.

Consecutive-D intervals use the same paired-family formula, comparing the
next larger D to the preceding D for each traversal. These comparisons
include changes in the realized tables as well as the residue count.
They are recorded separately in `repeats-residue-differences.csv`.

## Results

Measurements use the Xeon Gold 6148 and the AVX-512 flags above. All 20
D=1024 collision conditions reproduce the earlier results exactly, including
rates, counts and paired intervals. Raw batches are checked for gaps and
overlaps before merging; counts are also checked for consistency across
D and for nesting of the first 10 pivots within the all-pivot results.

For n=50, E=8, R=49 and the one-deletion/one-insertion pair model:

| D | Inward us/window | Outward us/window | Inward collision % | Outward collision % |
|---:|---:|---:|---:|---:|
| 32 | 59.344 | 542.557 | 26.42 | 41.41 |
| 128 | 205.113 | 1481.735 | 34.08 | 46.73 |
| 1024 | 892.751 | 12295.166 | 45.46 | 57.03 |
| 2048 | 1704.984 | 26763.758 | 47.75 | 57.91 |

For this model, both traversals show clear increases in collision probability
from D=32 to 128 and from 128 to 1024. Moving from 1024 to 2048 roughly
doubles runtime, while the smaller observed collision gains have intervals
that include zero. This observation is specific to this configuration and
mutation model; the complete results below cover the other conditions.

| Traversal | D change | Collision gain pp | Approx. paired 95% interval |
|:---|:---|---:|:---|
| inward | 32 to 128 | 7.67 | [5.55, 9.78] |
| outward | 32 to 128 | 5.32 | [2.65, 7.99] |
| inward | 128 to 1024 | 11.38 | [8.40, 14.35] |
| outward | 128 to 1024 | 10.30 | [7.76, 12.85] |
| inward | 1024 to 2048 | 2.29 | [-1.40, 5.99] |
| outward | 1024 to 2048 | 0.88 | [-2.61, 4.37] |

## All runtime configurations

These are microseconds per window, median of seven samples. Larger D
increases work and memory. The implementation has a specialized join at
D<=32 and records direct traceback parents at D<=64, so exact scaling need
not be linear across these thresholds.

| n | E | Pivots | D | Inward us | Outward us | Outward / inward |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4 | 10 | 32 | 16.021 | 55.920 | 3.49 |
| 30 | 4 | 10 | 128 | 30.480 | 130.338 | 4.28 |
| 30 | 4 | 10 | 1024 | 136.543 | 947.859 | 6.94 |
| 30 | 4 | 10 | 2048 | 269.240 | 1939.302 | 7.20 |
| 30 | 4 | 29 | 32 | 23.828 | 142.027 | 5.96 |
| 30 | 4 | 29 | 128 | 75.331 | 366.391 | 4.86 |
| 30 | 4 | 29 | 1024 | 292.816 | 2585.171 | 8.83 |
| 30 | 4 | 29 | 2048 | 562.879 | 5391.016 | 9.58 |
| 50 | 8 | 10 | 32 | 19.512 | 118.420 | 6.07 |
| 50 | 8 | 10 | 128 | 56.398 | 323.462 | 5.74 |
| 50 | 8 | 10 | 1024 | 312.446 | 2593.401 | 8.30 |
| 50 | 8 | 10 | 2048 | 647.557 | 5450.703 | 8.42 |
| 50 | 8 | 49 | 32 | 59.344 | 542.557 | 9.14 |
| 50 | 8 | 49 | 128 | 205.113 | 1481.735 | 7.22 |
| 50 | 8 | 49 | 1024 | 892.751 | 12295.166 | 13.77 |
| 50 | 8 | 49 | 2048 | 1704.984 | 26763.758 | 15.70 |

## All mutation-model collision rates

Each cell is **inward % / outward %**, using the at-least-one-match,
same-pivot definition. Per-condition paired intervals, cross-pivot matches
and individual-pivot match rates are in the full CSV.

| n | E | Pivots | Pair generation | D=32 | D=128 | D=1024 | D=2048 |
|---:|---:|---:|:---|:---|:---|:---|:---|
| 30 | 4 | 10 | 2 substitutions | 22.17 / 26.12 | 28.71 / 30.76 | 36.87 / 37.79 | 40.38 / 40.28 |
| 30 | 4 | 10 | 4 substitutions | 1.95 / 2.73 | 3.96 / 3.66 | 6.05 / 5.62 | 5.37 / 5.76 |
| 30 | 4 | 10 | 1 deletion(s) + 1 insertion(s) | 18.07 / 22.36 | 22.51 / 25.98 | 29.05 / 31.25 | 27.59 / 32.08 |
| 30 | 4 | 10 | 2 deletion(s) + 2 insertion(s) | 2.78 / 2.88 | 2.98 / 3.96 | 5.18 / 5.37 | 4.44 / 4.35 |
| 30 | 4 | 29 | 2 substitutions | 48.97 / 58.35 | 61.04 / 67.77 | 71.73 / 75.49 | 74.32 / 77.00 |
| 30 | 4 | 29 | 4 substitutions | 6.10 / 7.37 | 10.11 / 11.28 | 15.33 / 14.70 | 15.14 / 16.80 |
| 30 | 4 | 29 | 1 deletion(s) + 1 insertion(s) | 36.38 / 48.88 | 45.51 / 55.13 | 56.54 / 61.72 | 56.79 / 62.30 |
| 30 | 4 | 29 | 2 deletion(s) + 2 insertion(s) | 6.45 / 7.76 | 7.42 / 10.45 | 10.35 / 11.72 | 11.47 / 10.94 |
| 50 | 8 | 10 | 4 substitutions | 2.29 / 2.69 | 3.22 / 3.56 | 5.13 / 6.01 | 5.81 / 5.47 |
| 50 | 8 | 10 | 8 substitutions | 0.00 / 0.05 | 0.05 / 0.05 | 0.05 / 0.05 | 0.24 / 0.20 |
| 50 | 8 | 10 | 1 deletion(s) + 1 insertion(s) | 8.79 / 12.45 | 11.18 / 15.43 | 15.38 / 19.09 | 17.29 / 19.53 |
| 50 | 8 | 10 | 2 deletion(s) + 2 insertion(s) | 0.59 / 0.98 | 0.98 / 1.12 | 1.95 / 1.42 | 2.25 / 1.27 |
| 50 | 8 | 49 | 4 substitutions | 9.08 / 11.72 | 15.19 / 16.89 | 23.83 / 27.05 | 25.29 / 29.74 |
| 50 | 8 | 49 | 8 substitutions | 0.24 / 0.05 | 0.44 / 0.49 | 0.54 / 0.54 | 0.98 / 0.93 |
| 50 | 8 | 49 | 1 deletion(s) + 1 insertion(s) | 26.42 / 41.41 | 34.08 / 46.73 | 45.46 / 57.03 | 47.75 / 57.91 |
| 50 | 8 | 49 | 2 deletion(s) + 2 insertion(s) | 2.05 / 3.71 | 4.25 / 5.13 | 6.64 / 6.45 | 7.52 / 8.06 |

No collisions were observed in any of the 16 unrelated-string controls.
This is an observation at the tested sample sizes, not a zero-rate bound.

## Artifacts

* [All 16 runtime measurements](repeats-residue-runtime.csv).
* [All 80 collision conditions and traversal intervals](repeats-residue-collisions.csv).
* [All 120 consecutive-D comparisons and intervals](repeats-residue-differences.csv).
* [All 2,560 table-family count rows](repeats-residue-families.csv).
* [Machine and experiment metadata](repeats-residue-metadata.json).

`--aggregate-only` on the Python runner regenerates collision summaries and
D comparisons from saved batch files while retaining the recorded runtimes.
