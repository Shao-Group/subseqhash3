# Single-thread seeding runtime improvements

The commands keep the same seeds, scores, tie order, positions, CSV format,
table format and traversal semantics. D need not be a power of two.
The 2x target is reached on several workloads, rather than uniformly across
all parameters. The complete measured range is 1.03x to 7.36x.

## Changes

* For outward sliding windows with D>=128 and nontrivial edit budgets,
  process every use of a global anchor while its two DP tables are hot.
  Retain candidate seeds in a bounded ring until the window is complete.
  This replaces R+1 (sub) or R+E (ind) pairs of cached score tables with one
  reusable pair. Anchors arrive in the original tie order, and callbacks
  still arrive in window/set/pivot order. Keep the previous engine for short
  sequences, small D and the special zero-/one-deletion paths.
* Once a candidate reaches residue zero, skip unreachable-zero candidates'
  searches for a larger residue: they cannot win. Wide imported deletion
  scores fall back before any callback is emitted.
* For CLI input lines containing exactly one window, retain workspaces and
  convert deletion tables once across lines. This removes repeated allocation,
  page faults and table conversion. Mixed short/long inputs are supported.
* Encode each CSV row into a stack buffer and write it once, replacing
  per-base stdio calls.

The anchor ring stores O(S*R*(R+E)*W) seed data (E=1 for sub), independent
of sequence length. Its DP storage no longer grows with the number of
live anchors. Candidate traceback can run more than once when a better
anchor arrives; this tradeoff accounts for some workloads gaining less.

## End-to-end measurements

Intel Xeon Gold 6148; one thread pinned to CPU 2. Both executables use GCC
with identical flags:

```
-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512
```

Baseline: commit `e6f3831c67682823cae0ee08347e379a67ffa134`, preserved at
`build/runtime-before/subseqhash3`. Updated executable:
`build/runtime-after/subseqhash3` (also built at `build/simd/subseqhash3`).
Three alternating before/after pairs per configuration, after one untimed
warmup per executable. Times include process startup, table loading,
conversion, seeding and CSV writing; table generation is excluded. Inputs
are deterministic random DNA. S=1, outward traversal, identical pre-generated
tables. Every full output CSV matches the baseline SHA256 on every run.
Short inputs are 300 lines of length W; long inputs are one line of length
512. R=all means W-1 for sub and W-E for ind.

| W | E | D | Input | R | seeds-sub speedup | seeds-ind speedup |
|---|---|---|---|---|---|---|
| 30 | 4 | 32 | short | 10 | 1.72x | 1.50x |
| 30 | 4 | 32 | short | all | 1.24x | 1.26x |
| 30 | 4 | 32 | long | 10 | 1.04x | 1.08x |
| 30 | 4 | 32 | long | all | 1.05x | 1.08x |
| 30 | 4 | 1024 | short | 10 | 2.03x | 7.36x |
| 30 | 4 | 1024 | short | all | 1.62x | 4.96x |
| 30 | 4 | 1024 | long | 10 | 1.24x | 1.54x |
| 30 | 4 | 1024 | long | all | 2.16x | 1.41x |
| 50 | 8 | 32 | short | 10 | 2.94x | 1.98x |
| 50 | 8 | 32 | short | all | 1.94x | 1.65x |
| 50 | 8 | 32 | long | 10 | 1.03x | 1.09x |
| 50 | 8 | 32 | long | all | 1.03x | 1.10x |
| 50 | 8 | 1024 | short | 10 | 2.63x | 3.50x |
| 50 | 8 | 1024 | short | all | 1.62x | 1.30x |
| 50 | 8 | 1024 | long | 10 | 1.75x | 2.12x |
| 50 | 8 | 1024 | long | all | 2.53x | 2.45x |

At W=50, E=8, D=1024, R=10, short inputs take 2.015 -> 0.767 seconds
for sub (2.63x) and 1.101 -> 0.314 seconds for ind (3.50x).
Long inputs with all pivots take 1.222 -> 0.483 seconds for sub (2.53x)
and 0.504 -> 0.206 seconds for ind (2.45x). Small-D long-sequence workloads
improve only slightly; do not expect a universal 2x improvement. These
measurements apply to the AVX-512 build on this processor.

A separate `/usr/bin/time -f '%M'` check at W=50, E=8, D=1024,
length=512, S=1 and all pivots measured peak RSS of 336,640 -> 10,880 KiB
for sub (30.9x less) and 185,600 -> 11,520 KiB for ind (16.1x less).
These are whole-process peaks; both output files also matched byte for byte.

## Exact typical workload

For W=30, E=5, R=10, S=1, one sequence of length 1000, and outward
traversal, the previous improvement does **not** reach 2x or the later 5x
target. Five alternating paired measurements, after one untimed warmup per
executable, give:

| D | Command | Before (ms) | After (ms) | Speedup |
|---|---|---|---|---|
| 64 | seeds-sub | 28.634 | 27.814 | 1.03x |
| 64 | seeds-ind | 18.781 | 17.205 | 1.09x |
| 1024 | seeds-sub | 194.257 | 162.632 | 1.19x |
| 1024 | seeds-ind | 88.647 | 55.991 | 1.58x |

Each of the four complete CSV outputs is identical to the original on
every run. These results supersede extrapolating the W=30/E=4 or W=50/E=8
measurements to this workload. The implementation in this checkpoint
preserves the original scheme; no experimental replacement scheme is included.

```sh
python3 bench/run_typical_seed_runtime.py --samples 5 --cpu 2
```

* [Typical-workload timing summary](seed-runtime-typical-summary.csv)
* [Individual typical-workload measurements](seed-runtime-typical-samples.csv)
* [Input, table, source and executable hashes](seed-runtime-typical-metadata.json)

## Validation and reproduction

The full AVX-512 and scalar suites pass. Additional ASan/UBSan runs cover
window APIs, multiple sets, callbacks/cancellation, wide-score fallback and
the CLI. Window checks compare scores, seeds, positions and substitution
flags against the independent per-window APIs; 1,032 configurations and
20,964 windows. The R*S suite covers another 336 configurations. The CLI
suite has 147 checks, including repeated single-window inputs interleaved
with long lines, both traversals, multiple table sets and non-power-of-two D.
Leak detection is disabled for the sanitizer run in this environment.

```sh
make BUILD_DIR=build/runtime-after \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' test
python3 bench/run_seed_runtime.py --samples 3 --cpu 2 \
  --before build/runtime-before/subseqhash3 --after build/runtime-after/subseqhash3
```

The runner needs a separately built baseline executable; do not rebuild the
baseline from the modified sources. It pins itself and its sequential child
processes to one CPU and checks every output before accepting timings.

* [Full timing summary](seed-runtime-summary.csv)
* [Individual measurements and CSV hashes](seed-runtime-samples.csv)
* [Machine, binary and source hashes](seed-runtime-metadata.json)
* [Benchmark runner](run_seed_runtime.py)
