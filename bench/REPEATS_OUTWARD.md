# Outward versus inward repeats

Both traversals produce length-(n-1) seeds with **exactly one deletion and E
substitutions of retained bases**. They use identical table arrays and
minimize psi, then maximize omega. Inward reads the retained left prefix
normally and the right suffix in reverse; outward reads the left side in
reverse and the right side normally. This changes the hash family and
selected seeds. Outward is the CLI default; `--inward` selects inward.

```sh
build/simd/subseqhash3 repeats input.seq 30 4 1024 10 inward.csv shared.piv --inward
build/simd/subseqhash3 repeats input.seq 30 4 1024 10 outward.csv shared.piv
```

The new APIs are `subseqhash3_compute_outward_pivot_seeds()`,
`subseqhash3_compute_outward_pivot_repeats()` and
`subseqhash3_compute_outward_pivot_repeats_with_workspace()`. Their arguments,
result ownership and table-family batching match the inward APIs. Explicit
pivots may be unsorted or duplicated. Workspaces may alternate traversals,
inputs, dimensions, substitution counts and table families.

## Implementation and verification

Outward traversal reverses each side around an input anchor in a temporary
input, allowing reuse of the same compact SIMD DP kernel, join and score
encoding as inward. Output pivot p needs input anchors p and p+1. A two-slot
cache retains both anchors' side DPs, reusing the overlapping anchor for the
next pivot and retaining scores for traceback. Cached scores are rebuilt for
each new input and table family. Inward builds two shared prefix/suffix DPs;
outward builds side DPs for each needed input anchor. Both paths remain
single-threaded. The cache uses O(n*(E+1)*D) score memory, with a larger
constant than inward, rather than storing all anchors simultaneously.

`make test` now includes independent brute-force checks of both orders: 140
exhaustive edit neighborhoods and 30 larger witnesses per traversal, plus
score limits, unreachable zero residues, explicit/duplicate pivots, workspace
reuse and repeat counts spanning multiple table families. CLI/API tests cover
the outward default, `--inward`, compatibility `--outward`, alternate
traversals in one workspace, partial vectors, n=257 and E=256, and table-file
round trips. With E=0, outward repeats scores also
match the independent outward deletion API at deletion budget 1 across 16
configurations. Scalar, AVX-512 and address/undefined-behavior sanitizer
checks pass.

## Reproduce

```sh
make BUILD_DIR=build/repeats-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test build/repeats-fast/bench-repeats-traversal
build/repeats-fast/bench-repeats-traversal --runtime > bench/repeats-traversal-runtime.csv
build/repeats-fast/bench-repeats-traversal --collisions 32 64 \
  bench/repeats-traversal-families.csv > bench/repeats-traversal-collisions.csv
```

Runtime uses the median of seven samples, 16 deterministic random windows
and three passes per sample, with alternating execution order. Inputs,
tables, E, and requested pivots 0..R-1 are identical. Reusable workspaces are
warmed first. Timing includes result allocations, traceback and cleanup;
it excludes table generation and initial workspace allocation. R=n-1 means
all pivots of one family; larger repeat counts are supported by the API and
tested separately.

## Collision measurements

Each condition uses 32 independent table families and 64 independent input
pairs per family (2,048 pairs). Within each condition, both traversals receive
identical sequences and tables. Input bases are independent uniform A/C/G/T.
Pairs have E/2 or E substitutions at distinct positions, one or two deletions
followed by the same number of insertions, or unrelated independently
generated sequences. Insertions and deletions may cancel by chance: operation
counts are not guaranteed optimal edit distances. `nonidentical_pairs`
records how many pairs actually differ.

A collision means exact equality of selected seed strings, not equal numeric
scores. The main `same` metric is the probability of at least one matching
seed at the same output pivot. `any` allows matches between different output
pivots. `per` is the fraction of matching individual pivot outputs. All
probabilities in the CSV are percentages; `delta` is inward minus outward
in percentage points. All candidate seeds in this experiment have length
n-1, unlike the earlier deletion-only experiment, whose seeds have length
n-E. Its collision rates therefore are not a traversal-only comparison with
these rates.

Approximate paired 95% intervals use the independent table families as
units: mean per-family difference +/- 2.04 * sample_standard_deviation /
sqrt(32). Per-family counts and table seeds are saved for reproduction.
Intervals are exploratory and unadjusted for multiple comparisons. Zero
observed matches do not establish zero probability or equivalence. The
experiment covers synthetic uniform strings, not real reads, GC bias,
long repeats or variable-length pairs.

These measurements hold seed count fixed. An equal-runtime comparison using
additional inward table families would be a separate experiment.

## Results (2026-10-21, Xeon Gold 6148, one thread)

At D=1024 and E=4/8, outward is 6.95--13.69 times slower than inward at
fixed seed count. Both traversals recompute their DPs per input window;
the benchmark does not share scores across overlapping sliding windows.

| n | E substitutions | Pivots | Inward us | Outward us | Outward / inward |
|---:|---:|---:|---:|---:|---:|
| 30 | 4 | 10 | 136.374 | 948.135 | 6.95 |
| 30 | 4 | 29 | 292.988 | 2587.848 | 8.83 |
| 50 | 8 | 10 | 313.531 | 2595.283 | 8.28 |
| 50 | 8 | 49 | 889.511 | 12179.999 | 13.69 |

Outward raises at-least-one-match probability in several conditions,
especially the one-deletion/one-insertion model. It is not a consistent
winner across all mutations. For n=50, E=8, R=49, one deletion plus one
insertion gives 57.03% outward versus 45.46% inward: an outward gain of
11.57 percentage points (approximate paired 95% interval 8.33--14.81).
Two deletions plus two insertions give 6.45% versus 6.64%, an inconclusive
difference. No collisions were observed in the unrelated controls.

| n | E | Pivots | Pair mutation | Inward % | Outward % | Delta (inward - outward) pp | Approx. 95% interval |
|---:|---:|---:|:---|---:|---:|---:|:---|
| 30 | 4 | 10 | 2 substitutions | 36.87 | 37.79 | -0.93 | [-3.76, 1.90] |
| 30 | 4 | 10 | 4 substitutions | 6.05 | 5.62 | 0.44 | [-1.01, 1.89] |
| 30 | 4 | 10 | 1 deletions + 1 insertions | 29.05 | 31.25 | -2.20 | [-4.93, 0.54] |
| 30 | 4 | 10 | 2 deletions + 2 insertions | 5.18 | 5.37 | -0.20 | [-1.27, 0.88] |
| 30 | 4 | 29 | 2 substitutions | 71.73 | 75.49 | -3.76 | [-6.45, -1.07] |
| 30 | 4 | 29 | 4 substitutions | 15.33 | 14.70 | 0.63 | [-1.72, 2.99] |
| 30 | 4 | 29 | 1 deletions + 1 insertions | 56.54 | 61.72 | -5.18 | [-8.03, -2.32] |
| 30 | 4 | 29 | 2 deletions + 2 insertions | 10.35 | 11.72 | -1.37 | [-2.70, -0.03] |
| 50 | 8 | 10 | 4 substitutions | 5.13 | 6.01 | -0.88 | [-2.21, 0.45] |
| 50 | 8 | 10 | 8 substitutions | 0.05 | 0.05 | 0.00 | [-0.14, 0.14] |
| 50 | 8 | 10 | 1 deletions + 1 insertions | 15.38 | 19.09 | -3.71 | [-5.67, -1.75] |
| 50 | 8 | 10 | 2 deletions + 2 insertions | 1.95 | 1.42 | 0.54 | [-0.30, 1.38] |
| 50 | 8 | 49 | 4 substitutions | 23.83 | 27.05 | -3.22 | [-6.17, -0.28] |
| 50 | 8 | 49 | 8 substitutions | 0.54 | 0.54 | 0.00 | [-0.52, 0.52] |
| 50 | 8 | 49 | 1 deletions + 1 insertions | 45.46 | 57.03 | -11.57 | [-14.81, -8.33] |
| 50 | 8 | 49 | 2 deletions + 2 insertions | 6.64 | 6.45 | 0.20 | [-1.13, 1.52] |

The union probability can improve even when individual-pivot rates are
similar. For n=30, E=4, R=29 and two substitutions, per-pivot collision
rates are 4.70% outward and 4.76% inward, while at-least-one-match rates are
75.49% and 71.73%. This is consistent with outward yielding less redundant
matches across the pivot orders; it is not a proof of a general advantage.

Allowing matches between different output pivots changes these observed
rates by at most 0.10 percentage points. The raw CSV includes both definitions.
The approximate intervals and synthetic mutation models do not establish
a universal collision ranking.

Full artifacts: [runtime](repeats-traversal-runtime.csv),
[collision summaries](repeats-traversal-collisions.csv), and
[per-family counts](repeats-traversal-families.csv).

For the expanded comparison at D=32, 128, 1024 and 2048 using matched
input pairs, see the [residue-count sweep](REPEATS_RESIDUES.md).
