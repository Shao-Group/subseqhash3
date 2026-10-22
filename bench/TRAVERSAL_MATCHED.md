# Inward versus outward at matched runtime

Adjusting D changes the conclusion of a fixed-D comparison. On these
synthetic workloads, **inward is the stronger choice for repeats**, whereas
**deletions has no single winner**: outward retains an indel advantage in
some 10-pivot cases, while inward can win when all pivots are requested.
Many deletion differences are too small to distinguish at this sample size.
These conclusions concern the measured configurations and budgets, not a
universal ordering or a proof of optimality over all D.

## Protocol

Measurements were made on 2026-10-21 on the Xeon Gold 6148, using one thread
per seeding call and AVX-512. Deletions retains n-E bases with exactly E
deletions; repeats retains n-1 bases with one deletion and exactly E
substitutions. Each traversal comparison keeps n, E and the requested
pivots fixed. Results compare traversals within an algorithm; deletion and
repeats probabilities describe different seed problems.

For each of the four configurations per algorithm, two budgets start at
inward D=1024 and 2048. A timing-only binary search chooses an outward D,
also checking nearby multiples of 32. Fresh timings validate the chosen
pair. Where SIMD/traceback thresholds leave a gap, inward D is adjusted
using timing alone. The final measured runtime ratios differ from 1 by at
most 3.7%. This is a comparison of matched operating points, not an exhaustive
search over D, and exact timing will vary across runs and machines.

Runtime is the median of seven samples, each containing 48 window calls on
16 deterministic random windows, with alternating traversal order and warmed
workspaces. Allocation, traceback and cleanup are included. Table generation
and initial workspace allocation are excluded. Every window recomputes its
DPs. Calibration runs sequentially before any collision workers start.

Each operating point and mutation model uses 32 independent table families
and 64 input pairs per family (2,048 pairs). Both traversals see identical
pairs. The canonical pair stream is independent of D and preserves the
previous experiments' D=1024 inputs. Family seed identifiers are shared,
but realized tables differ when D differs. Input strings are uniform DNA;
substitutions change distinct positions. The indel models perform one or
two deletions followed by the same number of insertions, preserving length;
chance cancellations and identical pairs are retained. They describe edit
operations, not guaranteed optimal edit distances. Unrelated pairs are
independent uniform strings.

The primary collision probability is at least one byte-identical selected
seed at the same output pivot. Numerical score coincidences do not count.
Raw files also report cross-pivot matches and per-pivot match rates.
Approximate paired 95% intervals use independent table-family differences:
mean +/- 2.04 * sample standard deviation / sqrt(32). These intervals are
exploratory and unadjusted for multiple comparisons; an interval containing
zero does not establish equivalence. No unrelated-pair collisions were
observed, which does not establish a zero underlying probability.

## Matched operating points

D values and times below are inward / outward; times are microseconds per
window. The two rows per configuration are the lower and higher budgets.

| Algorithm | n | E | Pivots | D inward / outward | us inward / outward |
|:---|---:|---:|---:|:---|:---|
| deletions | 30 | 4 | 10 | 1024 / 276 | 177.9 / 180.3 |
| deletions | 30 | 4 | 10 | 2048 / 601 | 378.0 / 369.4 |
| deletions | 30 | 4 | 26 | 1024 / 128 | 241.5 / 235.7 |
| deletions | 30 | 4 | 26 | 2048 / 392 | 530.2 / 533.1 |
| deletions | 50 | 8 | 10 | 1024 / 160 | 321.2 / 309.3 |
| deletions | 50 | 8 | 10 | 2048 / 404 | 643.8 / 642.3 |
| deletions | 50 | 8 | 42 | 1024 / 32 | 534.2 / 525.8 |
| deletions | 50 | 8 | 42 | 2016 / 145 | 1036.2 / 1061.1 |
| repeats | 30 | 4 | 10 | 1024 / 128 | 157.5 / 155.9 |
| repeats | 30 | 4 | 10 | 2048 / 288 | 268.0 / 260.2 |
| repeats | 30 | 4 | 29 | 832 / 64 | 243.8 / 240.0 |
| repeats | 30 | 4 | 29 | 1952 / 192 | 541.8 / 529.9 |
| repeats | 50 | 8 | 10 | 800 / 96 | 247.6 / 247.9 |
| repeats | 50 | 8 | 10 | 2048 / 256 | 643.2 / 620.0 |
| repeats | 50 | 8 | 49 | 640 / 32 | 573.2 / 565.4 |
| repeats | 50 | 8 | 49 | 1888 / 128 | 1543.0 / 1530.2 |

## Deletions

At n=50, E=8, R=10, the lower matched budget favors outward for both indel
models: 66.36% versus 62.30% for one deletion/insertion (outward gain 4.05 pp,
95% interval 1.49 to 6.61), and 18.16% versus 14.89% for two (gain 3.27 pp,
interval 1.14 to 5.41). The higher budget does not resolve a difference in
these cases.

With all 42 pivots at the lower budget, inward wins for two deletions/insertions:
46.14% versus 40.97%, an inward gain of 5.18 pp (interval 2.36 to 7.99).
The one-deletion/insertion rates are nearly equal. At the higher budget,
outward has a small one-deletion/insertion advantage: 98.24% versus 97.22%
(gain 1.03 pp, interval 0.05 to 2.00).

Across 32 non-control deletion comparisons, intervals exclude zero in three
inward directions and three outward directions; 26 remain inconclusive.
Counting these directions is descriptive, not an aggregate significance test.
The practical choice depends on the pivot count, mutation model and budget.

## Repeats

Inward is the better supported choice here after spending its speed advantage
on a larger D. Across 32 non-control comparisons, 17 intervals favor inward,
none favor outward, and 15 include zero. This is strongest for substitutions;
indel cases often have overlapping rates, especially at the higher budget.

For n=50, E=8, R=49 at about 570 us/window, inward D=640 versus outward D=32
has collision probabilities:

* Four substitutions: 23.29% versus 11.72%; gain 11.57 pp [8.88, 14.26].
* One deletion/insertion: 45.41% versus 41.41%; gain 4.00 pp [2.05, 5.96].
* Two deletions/insertions: 5.66% versus 3.71%; gain 1.95 pp [0.97, 2.93].

At about 1,540 us/window, inward D=1888 versus outward D=128 still wins
for four substitutions (23.68% versus 16.89%; gain 6.79 pp [4.33, 9.25]),
while the one-deletion/insertion comparison is inconclusive (48.34% versus
46.73%; gain 1.61 pp [-1.27, 4.50]). Thus an outward advantage at equal D
does not imply an outward advantage at equal runtime.

## All non-control collision probabilities

Each probability cell is inward % / outward %. The full CSV contains the
paired intervals, controls and alternative match definitions. The same
row order and D values as the operating-point table are used below.

| Algorithm | n | E | Pivots | D inward / outward | E/2 substitutions | E substitutions | 1 deletion + insertion | 2 deletions + insertions |
|:---|---:|---:|---:|:---|:---|:---|:---|:---|
| deletions | 30 | 4 | 10 | 1024 / 276 | 14.50 / 13.77 | 0.68 / 0.73 | 63.33 / 64.89 | 16.06 / 14.50 |
| deletions | 30 | 4 | 10 | 2048 / 601 | 14.16 / 13.87 | 0.39 / 0.44 | 62.84 / 64.45 | 15.82 / 14.79 |
| deletions | 30 | 4 | 26 | 1024 / 128 | 29.10 / 27.15 | 1.22 / 1.17 | 91.11 / 91.16 | 32.32 / 29.25 |
| deletions | 30 | 4 | 26 | 2048 / 392 | 31.45 / 30.32 | 1.12 / 1.51 | 91.60 / 91.94 | 32.03 / 31.25 |
| deletions | 50 | 8 | 10 | 1024 / 160 | 0.78 / 0.73 | 0.00 / 0.00 | 62.30 / 66.36 | 14.89 / 18.16 |
| deletions | 50 | 8 | 10 | 2048 / 404 | 0.78 / 0.63 | 0.00 / 0.00 | 65.77 / 66.65 | 16.41 / 16.65 |
| deletions | 50 | 8 | 42 | 1024 / 32 | 3.32 / 1.86 | 0.05 / 0.00 | 96.97 / 97.17 | 46.14 / 40.97 |
| deletions | 50 | 8 | 42 | 2016 / 145 | 3.08 / 2.98 | 0.00 / 0.00 | 97.22 / 98.24 | 48.63 / 48.54 |
| repeats | 30 | 4 | 10 | 1024 / 128 | 36.87 / 30.76 | 6.05 / 3.66 | 29.05 / 25.98 | 5.18 / 3.96 |
| repeats | 30 | 4 | 10 | 2048 / 288 | 40.38 / 34.67 | 5.37 / 5.42 | 27.59 / 29.05 | 4.44 / 4.74 |
| repeats | 30 | 4 | 29 | 832 / 64 | 70.17 / 63.53 | 13.87 / 9.77 | 54.49 / 52.69 | 10.25 / 8.59 |
| repeats | 30 | 4 | 29 | 1952 / 192 | 73.88 / 69.92 | 15.82 / 12.11 | 57.18 / 55.42 | 10.99 / 10.94 |
| repeats | 50 | 8 | 10 | 800 / 96 | 4.64 / 3.17 | 0.15 / 0.05 | 15.19 / 14.21 | 1.90 / 0.98 |
| repeats | 50 | 8 | 10 | 2048 / 256 | 5.81 / 4.59 | 0.24 / 0.10 | 17.29 / 16.80 | 2.25 / 1.22 |
| repeats | 50 | 8 | 49 | 640 / 32 | 23.29 / 11.72 | 0.54 / 0.05 | 45.41 / 41.41 | 5.66 / 3.71 |
| repeats | 50 | 8 | 49 | 1888 / 128 | 23.68 / 16.89 | 0.73 / 0.49 | 48.34 / 46.73 | 6.64 / 5.13 |

## Reproduce and artifacts

```sh
make BUILD_DIR=build/repeats-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  build/repeats-fast/bench-traversal-matched
python3 bench/run_traversal_matched.py \
  --executable build/repeats-fast/bench-traversal-matched --jobs 8
```

Collision workers are independent single-threaded processes; parallelism
only accelerates probability sampling. The runner saves raw files in
`build/traversal-matched/`. `--aggregate-only` regenerates the combined
probability files from these saved counts without rerunning calibration or
sampling. `--skip-calibration` reuses saved operating points, refines any
remaining timing gaps, then reruns sampling. New calibration can select
different D values due to timing variation.

* [Benchmark source](traversal_matched.c) and [runner](run_traversal_matched.py).
* [All calibration timings](traversal-matched-calibration.csv).
* [Matched operating-point runtimes](traversal-matched-runtime.csv).
* [80 collision conditions and paired intervals](traversal-matched-collisions.csv).
* [2,560 per-family count rows](traversal-matched-families.csv).
* [Experiment metadata](traversal-matched-metadata.json).

The matched benchmark reproduces previously saved collision counts at
shared parameters; the validation checks both traversals and multiple D
values. This experiment samples uniform synthetic windows, not real reads,
GC bias, long genomic repeats or variable-length sequence pairs.
