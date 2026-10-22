# Inward versus outward deletion seeding

Both algorithms choose length-k subsequences with exactly E deletions and
zero substitutions, k=n-E. The default outward order traverses each seed's
left side in reverse and its right side in normal order. `--inward` traverses
the left side normally and the right side in reverse. The same table arrays
are used for both. This is a different hash family, not an equivalent way to
compute the original SubseqHash2 seeds.

```sh
build/simd/subseqhash3 deletions input.seq 30 4 1024 10 outward.csv shared.piv
build/simd/subseqhash3 deletions input.seq 30 4 1024 10 inward.csv shared.piv --inward
```

The inward APIs are `subseqhash3_compute_inward_deletion_seeds()` and
`subseqhash3_compute_inward_deletion_seeds_with_workspace()`, with the same
arguments and result ownership as the outward APIs. Workspaces can alternate
between traversals, table families, dimensions, and compact/wide magnitudes.

## Implementation and correctness

Inward scores are prefix/suffix DPs indexed by consumed input length,
deletion count, and residue. All pivots share these states. Small scores use
16-bit costs `127*j +/- omega`, with reflected right residues for circular
vector joins. Traceback reconstructs deletion ties from stored scores, so
no parent array is needed. Wide table magnitudes use exact 64-bit scores.
The outward algorithm and its tie preferences remain available unchanged.

`make test` checks both orders against independent brute-force subsequence
enumeration: 2,430 configurations per order, including ties, unreachable
zero residues, compact magnitude boundaries, wide scores, duplicate/explicit
pivots and workspace reuse. Additional large witnesses and CLI/API comparisons
cover both modes. For E=1, inward scores also match the independently
implemented repeats algorithm with zero substitutions across 16 configurations.
The compact and wide inward kernels also return identical scores, seeds and
positions in six larger configurations, including D=1025 partial vectors.
Scalar, AVX-512, and address/undefined-behavior sanitizer checks pass.

## Reproduce measurements

```sh
make BUILD_DIR=build/del-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  build/del-fast/bench-deletions-traversal
build/del-fast/bench-deletions-traversal --runtime > bench/deletions-traversal-runtime.csv
build/del-fast/bench-deletions-traversal --collisions 32 64 \
  bench/deletions-traversal-families.csv > bench/deletions-traversal-collisions.csv
```

Runtime uses 16 deterministic random windows, three passes per sample, seven
samples, and rotating order; reported times are medians in microseconds per
window. Both traversals receive identical tables, inputs, E, and requested
pivots 0..R-1. Reusable workspaces are warmed first. Allocation of result
strings, traceback and cleanup are timed; table generation and initial
workspace allocation are excluded. The repeats column is context only:
repeats performs one deletion plus E substitutions and returns length n-1,
so it is not an equivalent task.

## Collision experiment

A collision means **identical selected seed strings**, compared byte for
byte, not matching numerical scores or residues. Each condition uses 32
independent table seeds and 64 independent pairs per table (2,048 pairs).
Both traversals see exactly the same pairs and tables. Random strings have
independent uniform A/C/G/T bases. Mutation models are:

* Substitutions at exactly E/2 or E distinct positions; each replacement
  differs from the original base.
* One or two random deletions followed by the same number of random
  insertions, preserving input length. These are operation counts, not
  guaranteed optimal edit distances; chance cancellations are retained.
  `nonidentical_pairs` records how many pairs actually differ.
* Independently generated unrelated strings as a control.

The main `same` metric counts pairs with at least one matching seed at the
**same output pivot**. `any` allows a match between any two requested pivots.
`per` reports the fraction of matching pivot outputs, rather than the
probability of at least one collision. CSV probability fields are percentages;
`delta` is inward minus outward in percentage points.

Approximate paired 95% intervals use the 32 per-family differences:
mean +/- 2.04 * sample_standard_deviation / sqrt(32). Families are the
independent units, accounting for reuse of each table across 64 pairs.
Per-family counts and table seeds are saved in the raw CSV. These are
exploratory, unadjusted intervals, not a proof of equivalence or a ranking
for all mutation models. Zero observed collisions do not imply zero true
probability. No real reads, GC bias, long repeats or variable-length pair
models were measured in this experiment.

## Measurements (2026-10-21, Xeon Gold 6148, one thread)

At D=1024, inward is 3.4--10.3 times faster than the already optimized
outward kernel in the E=4/8 cases. It is around the repeats runtime when
requesting all pivots; repeats still computes a different edit problem.

| n | E | Pivots | Outward us | Inward us | Speedup | Repeats us |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4 | 10 | 602.598 | 175.974 | 3.42 | 122.871 |
| 30 | 4 | 26 | 1310.432 | 245.659 | 5.33 | 253.896 |
| 50 | 8 | 10 | 1559.619 | 321.073 | 4.86 | 295.495 |
| 50 | 8 | 42 | 5669.658 | 551.840 | 10.27 | 797.337 |

The collision table below uses the same-pivot, at-least-one-match metric.
The main indel differences favor outward, especially at n=50. Substitution
conditions have low or overlapping rates with no consistent inward advantage.

| n | E | Pivots | Pair mutation | Outward % | Inward % | Delta pp | Approx. 95% interval for delta |
|---:|---:|---:|:---|---:|---:|---:|:---|
| 30 | 4 | 10 | 2 substitutions | 14.65 | 14.50 | -0.15 | [-3.08, 2.79] |
| 30 | 4 | 10 | 4 substitutions | 0.24 | 0.68 | 0.44 | [0.00, 0.87] |
| 30 | 4 | 10 | 1 deletions + 1 insertions | 65.23 | 63.33 | -1.90 | [-4.60, 0.79] |
| 30 | 4 | 10 | 2 deletions + 2 insertions | 14.75 | 16.06 | 1.32 | [-0.86, 3.50] |
| 30 | 4 | 26 | 2 substitutions | 30.76 | 29.10 | -1.66 | [-4.78, 1.46] |
| 30 | 4 | 26 | 4 substitutions | 0.93 | 1.22 | 0.29 | [-0.37, 0.96] |
| 30 | 4 | 26 | 1 deletions + 1 insertions | 92.58 | 91.11 | -1.46 | [-3.30, 0.37] |
| 30 | 4 | 26 | 2 deletions + 2 insertions | 32.71 | 32.32 | -0.39 | [-2.72, 1.93] |
| 50 | 8 | 10 | 4 substitutions | 1.17 | 0.78 | -0.39 | [-0.91, 0.13] |
| 50 | 8 | 10 | 8 substitutions | 0.00 | 0.00 | 0.00 | [0.00, 0.00] |
| 50 | 8 | 10 | 1 deletions + 1 insertions | 67.43 | 62.30 | -5.13 | [-8.28, -1.98] |
| 50 | 8 | 10 | 2 deletions + 2 insertions | 17.09 | 14.89 | -2.20 | [-4.05, -0.35] |
| 50 | 8 | 42 | 4 substitutions | 3.42 | 3.32 | -0.10 | [-1.00, 0.81] |
| 50 | 8 | 42 | 8 substitutions | 0.00 | 0.05 | 0.05 | [-0.05, 0.15] |
| 50 | 8 | 42 | 1 deletions + 1 insertions | 99.27 | 96.97 | -2.29 | [-3.10, -1.49] |
| 50 | 8 | 42 | 2 deletions + 2 insertions | 53.71 | 46.14 | -7.57 | [-10.97, -4.17] |

No collisions were observed for unrelated pairs in any of the four
configurations. Allowing cross-pivot matches changes these observed rates
by at most 0.10 percentage points; both metrics are in the full CSV.

These results establish a speed/sensitivity tradeoff on this synthetic
workload. They do not show that inward is always less sensitive. The default
remains outward; `--inward` makes the faster alternative explicit.

Full artifacts: [runtime](deletions-traversal-runtime.csv),
[collision summaries](deletions-traversal-collisions.csv), and
[per-family counts](deletions-traversal-families.csv).
