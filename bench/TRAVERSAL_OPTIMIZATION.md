# Single-thread traversal optimization

For the current pivot-only CLI and explicit R/S arguments, see
[SEEDING.md](SEEDING.md). Command examples below use the historical interface.

All four paths have been optimized while preserving the existing hash
families, scores, selected seeds and deterministic tie witnesses. The broad
5x target is **not achieved on the earlier D=1024, E=4/8 workloads**. It is
exceeded in the one-deletion/zero-substitution cases at D=1024, and reached
approximately or exceeded in the D=1, E=4/8 cases.

Every ratio below holds n, E, D, requested pivots, input strings and realized
hash tables fixed between the saved implementation and production. Each
seeding call uses one thread; SIMD operates within that thread. Changing D
is not part of an individual before/after speedup.

## Measured speedups

Each range covers n=30/50 with 10 and all pivots. For the first two columns,
E is 4 at n=30 and 8 at n=50. In the last column, deletions uses E=1 and
repeats uses E=0, so both solve exactly one deletion and zero substitutions.

| Algorithm | Traversal | D=1024, E=4/8 | D=1, E=4/8 | D=1024, one deletion, zero substitutions |
|:---|:---|:---|:---|:---|
| deletions | inward | 1.31–1.69x | 4.97–8.88x | 20.61–43.55x |
| deletions | outward | 1.53–1.77x | 10.67–12.12x | 9.89–24.37x |
| repeats | inward | 1.03–1.08x | 5.99–14.49x | 11.75–16.47x |
| repeats | outward | 1.05–1.09x | 21.47–37.20x | 22.49–34.59x |

The 4.97x end of the D=1 inward-deletion range is approximately 5x, rather
than a strict 5x lower bound. These medians vary with CPU frequency and
system activity. The general sweep also covers D=32, 128 and 2048. D=32
repeats is close to parity; the full CSV records every measured case.

For the earlier n=50, E=8, all-pivot workload, times are microseconds per
window and include output allocation, traceback and cleanup:

| Algorithm | Traversal | Pivots | D | Before us | After us | Speedup |
|:---|:---|---:|---:|---:|---:|---:|
| deletions | inward | 42 | 1 | 137.134 | 26.400 | 5.19x |
| deletions | inward | 42 | 1024 | 546.455 | 417.819 | 1.31x |
| deletions | outward | 42 | 1 | 3309.607 | 272.986 | 12.12x |
| deletions | outward | 42 | 1024 | 5865.456 | 3547.480 | 1.65x |
| repeats | inward | 49 | 1 | 280.993 | 38.040 | 7.39x |
| repeats | inward | 49 | 1024 | 884.946 | 859.754 | 1.03x |
| repeats | outward | 49 | 1 | 6660.604 | 179.057 | 37.20x |
| repeats | outward | 49 | 1024 | 12451.813 | 11868.578 | 1.05x |

## Implementation

* **Residue-independent one-deletion solver.** Enumerate the n possible
  deleted positions instead of sweeping D residues. Inward shares unchanged
  prefixes/suffixes and only rescores the part after each deletion. Outward
  updates independent deletion candidates in SIMD lanes while sharing each
  table row. A crossover rule uses this path when D>=128 and D>=4n inward
  or D>=8n outward. Deletions supports exact int64 table values; repeats
  applies this path at zero substitutions. Explicit repeats pivot lists
  longer than 256 retain the general implementation.
* **Compact outward-deletion boundaries and traceback.** Keep boundary scores
  in 16-bit costs through the join, avoiding the previous 64-bit expansion.
  Two masks record strict min/max keep decisions per block of 32 residues.
  At D divisible by 32, both boundary score storage and parent storage are
  roughly one quarter of their previous sizes. The exact int64 fallback remains
  available for table coefficients outside the compact range.
* **Fused deletion validation and packing.** Check magnitudes while packing
  the rank rows actually consumed. Store 127+A in one byte and obtain the
  complementary cost as 254-(127+A). Prepared sign masks and an XOR-based
  complementary blend reduce work inside the vector recurrence. Tables are
  repacked on each call, so edits to externally owned table values remain
  visible.
* **Repeats join and loop improvements.** Retain lane minima across residue
  blocks and candidate edit splits, delaying horizontal reductions. Candidate
  indices preserve the existing strict tie order; the winning residue is
  recovered afterwards. Unroll the two deletion planes and the edit loop.
* **Small-D padding and D=1 kernels.** Cached lane permutations replace
  repeated tiny circular-padding copies. D=1 inward deletions stores
  contiguous scalar DP states. For D=1 repeats with E<32, SIMD lanes carry
  edit counts in a compact DP instead of residues, with the same traceback
  codes. Deletion-free remaining inward traceback prefixes need only copy
  bases and positions.

API signatures and their traversal semantics are unchanged. Both CLI commands
now default to outward; `--inward` selects inward. Repeats retains `--outward`
as an explicit compatibility option. To reproduce repeats output from its
previous default, reuse the same tables and add `--inward`. Each call
recomputes its seed selection. Existing collision-probability measurements
remain applicable to the same parameters; earlier runtime reports describe
the implementation before these optimizations.

Profiling of larger-D repeats identifies DP construction as the dominant
cost. Its dense four-symbol transitions still account for most of the work;
the retained control and join changes provide modest gains in those cases.
A 5x improvement there remains unachieved.

## Measurement and verification

Measurements were made on 2026-10-21 on the Xeon Gold 6148 with CPU affinity
fixed to logical CPU 0. Production and the frozen references use identical
-O3 and AVX-512 compiler flags. Measurements run sequentially, after other
validation jobs finish. Each time is the median of seven samples of 48 calls
on 16 deterministic random windows, with alternating implementation order
and warmed workspaces. Initial workspace allocation and table generation
are excluded; per-call output allocation, reconstruction and destruction
are included. Overlapping-window reuse is not measured.

Before timing each condition, the benchmark compares exact scores, seed
bytes, positions and substitution flags on 16 random windows plus two
homopolymers and a periodic string. The separate reference check covers
118 configurations and 2,242 windows, including D=1, partial vectors,
dispatch crossovers and E=31/32 at D=1. Frozen source snapshots are linked
under renamed symbols in the same process and are not used by production.

The full existing test suite and reference checks pass in scalar, AVX-512,
and AddressSanitizer/UndefinedBehaviorSanitizer builds. Deletion brute-force
coverage now includes D=128/129: 2,970 configurations per traversal, including
wide magnitudes, ties, unreachable zero residues, unsorted/duplicate pivots
and workspace reuse. Repeats additionally checks 140 exhaustive neighborhoods
and 30 larger witnesses per traversal, plus CLI and table-file compatibility.
Leak detection is disabled for the sanitizer run because this execution
environment uses ptrace; address and undefined-behavior checks are enabled.

## CLI defaults

Both commands use outward traversal unless `--inward` is supplied:

```sh
build/subseqhash3 repeats input.seq 30 4 1024 10 outward.csv shared.piv
build/subseqhash3 repeats input.seq 30 4 1024 10 inward.csv shared.piv --inward
build/subseqhash3 deletions input.seq 30 4 1024 10 outward-deletions.csv shared.piv
build/subseqhash3 deletions input.seq 30 4 1024 10 inward-deletions.csv shared.piv --inward
```

CLI regression checks compare CSV seeds and scores against the corresponding
C APIs for the outward default and explicit inward selection. Repeats also
checks its compatibility `--outward` flag and rejects duplicate or conflicting
traversal flags. Both scalar and AVX-512 test suites pass with these defaults.

## Reproduce

```sh
make BUILD_DIR=build/opt \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test check-traversal-speed
python3 bench/run_traversal_speed.py \
  --executable build/opt/bench-traversal-speed --cpu 0 \
  --flags '-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512'
```

`make bench-traversal-speed` runs the D=1024 general and one-deletion
comparisons. The Python runner additionally measures D=1, 32, 128 and 2048
and saves metadata and source hashes. `--residues D` on the C executable
selects one general sweep, and `--single [D]` selects the one-deletion cases.
Every timing process remains single-threaded.

* [General timings: 80 cases](traversal-speed-runtime.csv).
* [One-deletion timings: 16 cases](traversal-single-runtime.csv).
* [Machine, protocol and source hashes](traversal-speed-metadata.json).
* [Benchmark](traversal_speed.c) and [sequential sweep runner](run_traversal_speed.py).
* [Frozen references and provenance](reference/README.md).
