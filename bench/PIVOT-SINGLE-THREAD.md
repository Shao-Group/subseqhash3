# Strictly single-threaded second-round optimization

The requested further **2.5x speedup was not achieved**. The retained change
improves D=31 by about 1.18--1.30x in the final focused run; D=1024 remains
essentially unchanged. Both methods execute on one CPU. No threading code,
threading API, OpenMP dependency, or input-result cache is introduced.

## Algorithm retained

When D<=32, the whole residue space fits in one AVX-512 register. Instead of
reducing every candidate's join to a scalar, the new path maintains minimum
costs and winning candidate IDs in lanes across all input-pivot, pivot-base,
and substitution-allocation candidates. It performs two horizontal reductions
per pivot: one for the best cost, and one for the earliest candidate attaining
that cost. Only then does it recover the first winning residue.

The ID is `((input_pivot-pivot)*4+symbol)*(E+1)+left_substitutions`.
This preserves the previous input/base/allocation priority, even when different
lanes tie. The existing witness scan preserves residue priority. The candidate
ID fits in 16 bits for every supported input length. All finite final costs
are at most `254*256=65024`, below the unreachable sentinel, so adding the
pivot cost in unsigned saturating 16-bit lanes is exact.

If no candidate reaches psi=0, the implementation falls back to the unchanged
general smallest-residue search. It does not assume psi=0 is reachable.
The DP kernel, scalar fallback, D>32 join path, public API, order definition,
edit budgets, seed bytes, source positions, and substitution flags are unchanged.

## Direct before/after comparison

The reference is the **restored first-round optimization**, not the earlier
commit ca35567 and not simple repeating. Its exact source snapshot is
[`reference/pivots_round1.c`](reference/pivots_round1.c); see the accompanying
[reference description and hash](reference/README.md). Both versions link into
one benchmark executable, with the reference exports renamed and the same
compiler options applied to both.

Xeon Gold 6148, GCC 13.3, CPU 0, measured 2026-10-20. Each sample processes
32 random inputs four times; there are seven samples, alternating call order
and warming each method before its timed calls. Values below are per-method
medians in microseconds per input window, including DP allocation, all ten
seeds, traceback, and output destruction. Table initialization is excluded.

| n | E | D | r | Restored baseline (us) | Current (us) | Further speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 5 | 31 | 10 | 15.251 | 11.717 | 1.302x |
| 50 | 8 | 31 | 10 | 23.984 | 20.342 | 1.179x |
| 30 | 5 | 1024 | 10 | 147.372 | 148.251 | 0.994x |
| 50 | 8 | 1024 | 10 | 305.264 | 306.131 | 0.997x |

The sub-1% differences at D=1024 are not evidence of an improvement; that
workload has no retained algorithmic change. Timing varies with CPU frequency
and host load. [Summary CSV](pivot-single-results.csv) and
[all samples](pivot-single-samples.csv) retain the measured values.

## Experiments rejected

The following were tested single-threaded and removed from production code:

- Explicit E=5/8 unrolling: no gain in the focused large-D cases.
- 256-bit DP operations instead of 512-bit operations: slower.
- Rolling score buffers with packed parent history: slower, despite less
  score-history traffic; recording parent decisions added hot-loop work.
- Compile-time D=1024/E=5/8 layout specialization: no material gain.
- Combining small-D circular shifts and sign selection using register
  permutations: no additional reliable gain over the retained join change.

These experiments do not establish that a 2.5x speedup is impossible. They do
not provide evidence for claiming one either. The large-D shared DP still
dominates, and the retained small-D reduction cannot accelerate that kernel.

## Reproduce

On an AVX-512F/BW/VL CPU:

```sh
make BUILD_DIR=build/pivots-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test check-pivots build/pivots-fast/bench-pivots-single
taskset -c 0 build/pivots-fast/bench-pivots-single
```

Use an allowed CPU other than 0 if necessary. `--check-only` skips timing;
`--sweep` additionally measures D=65,257,512. No OpenMP flags or thread-count
environment variables are needed.

The benchmark performs 277 exact baseline comparisons before timing. In
addition to the prior randomized, non-power-of-two, and length-limit cases,
14 adversarial cases cover D=31/32, all-tied scores, n=257, budgets 0/1/256,
costs 0 and 65024, and cases with unreachable psi=0. Before timing it compares
every timed input's complete results against the reference as well.
The unchanged exhaustive-neighborhood/workspace/repeat checks still pass.

Scalar verification:

```sh
make BUILD_DIR=build/pivots-scalar-r2 \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic' \
  test check-pivots check-pivots-single
```

Standalone address/undefined-behavior sanitizer checks:

```sh
make BUILD_DIR=build/pivots-sanitize-r2 \
  CFLAGS='-O1 -g -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -fsanitize=address,undefined -fno-omit-frame-pointer' \
  LDFLAGS='-fsanitize=address,undefined' \
  build/pivots-sanitize-r2/bench-pivots \
  build/pivots-sanitize-r2/bench-pivots-single
ASAN_OPTIONS=detect_leaks=0 build/pivots-sanitize-r2/bench-pivots --check-only
ASAN_OPTIONS=detect_leaks=0 build/pivots-sanitize-r2/bench-pivots-single --check-only
```

Both standalone checks pass. Leak detection is disabled due to the environment's
tracing limitation; address and undefined-behavior checks remain enabled.
The full ordinary scalar/SIMD tests and 365 single-seed equivalence checks pass.
This does not claim the unrelated table-file serialization code is sanitizer-clean.
