# Follow-up pivot optimization

This report preserves the **first optimization round's historical timings**.
Its optimized implementation is frozen in
[`reference/pivots_round1.c`](reference/pivots_round1.c). The current source
also includes the subsequent small-D single-threaded join optimization; see
[the second-round report](PIVOT-SINGLE-THREAD.md). Consequently, running
`bench-pivots-speed` now compares the latest implementation against ca35567,
not the intermediate snapshot measured in the table below. Raw historical
samples are retained unchanged.

The requested further **5x speedup was not achieved**. This implementation
measures 1.45--1.63x faster than the committed pivot approach on the four
focused r=10 cases. These are additional speedups, not comparisons with simple
repeating, and not ratios formed from historical timings on different runs.

## Same-problem measurements

Measured on 2026-10-20, Xeon Gold 6148, CPU affinity 0, GCC 13.3. Both versions
are linked into the same executable using identical compiler flags, inputs,
tables, exact edit counts, requested pivots, and result cleanup. The reference
is `src/ssh3_pivots.c` from commit `ca35567`, with public symbols renamed.
Table initialization is outside timing; per-call workspace allocation,
seeding, traceback, and output destruction are inside timing. No cached input
results or cross-window memoization are used.

Each sample processes 32 random inputs four times (128 windows). Each method
gets one untimed warmup per sample. Seven samples alternate method order;
the table reports each method's median and their ratio. r=10 selects pivots
0 through 9 in one family. Times are microseconds per window, for all 10 seeds.

| n | E substitutions | D | Old pivot (us) | Optimized (us) | Further speedup |
|---:|---:|---:|---:|---:|---:|
| 30 | 5 | 31 | 22.435 | 15.468 | 1.450x |
| 50 | 8 | 31 | 38.714 | 23.727 | 1.632x |
| 30 | 5 | 1024 | 217.911 | 148.366 | 1.469x |
| 50 | 8 | 1024 | 446.646 | 305.995 | 1.460x |

See [summary CSV](pivot-optimization-results.csv) and
[all timing samples](pivot-optimization-samples.csv). Timings vary with
machine load and CPU frequency; the measured ratios are not universal bounds.

## Changes and limits

- Process deletion planes by retained length, improving table/working-set
  locality, and isolate immutable layout metadata from score-array aliasing.
- Name SIMD operands for the three changed bases and unchanged base, removing
  original-base selection branches from the substitution-count loop.
- For D>64, remove all parent recording from the hot DP and reconstruct only
  selected witnesses from predecessor costs. Keep packed parents for D<=64,
  where recomputing traceback was a measured regression.
- Use one vector copy for circular padding at D>=31, replacing the small
  memcpy loop; a 32-lane circular read needs only 31 wrapped elements.
- Evaluate candidate joins without maintaining per-lane argmin indices.
  Resolve the first residue attaining the known winning score only once.

These preserve the API, pivot order construction, psi/omega scores, exact
edit budgets, and deterministic witness ties. The score-only traceback tests
deletion first, ascending changed bases next, then the unchanged base, matching
the previous strict-update order. It does not approximate scores or assume
residue zero is reachable.

Temporary phase timing on the pre-optimization implementation attributed
about 88% of combined DP-build/join time to DP building. Optimizing joins alone
therefore cannot deliver the requested speedup on that workload. The current
algorithm still evaluates O(n*(E+1)*D) shared DP states; none of these changes
removes that dominant state volume. A further 5x has not been demonstrated;
this is not a proof that it is impossible. Hardware performance counters were
unavailable under the host's perf permissions, so no counter-based bottleneck
claim is made.

## Reproduce

Requires repository history containing `ca35567` and an AVX-512F/BW/VL CPU:

```sh
make BUILD_DIR=build/pivots-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test check-pivots build/pivots-fast/bench-pivots-speed
taskset -c 0 build/pivots-fast/bench-pivots-speed
```

Substitute an allowed CPU for 0 when needed. `--sweep` additionally times
D=65, 257, and 512; `--check-only` skips timing.

The new benchmark first checks 260 random/homopolymer configurations with
D in {1,3,7,30,31,32,33,63,64,65,257,1024,1025}, budgets spanning zero through
n-1, and default or explicitly selected (possibly repeated/unsorted) pivots.
Three more comparisons exercise n=257, D=65, E=0/128/256 and boundary/middle
pivots. It compares scores, seed bytes, source positions, flags, and budgets
against the committed implementation. Before timing, each of the 32 inputs
per case also gets an exact output comparison outside the timed interval.
Existing exhaustive-neighborhood, witness, score-limit, workspace, and repeat
batching checks remain unchanged and pass in scalar and SIMD builds.
The current benchmark also includes 14 subsequent adversarial tie/score-limit
cases, bringing its exact comparison count to 277.

Scalar verification:

```sh
make BUILD_DIR=build/pivots-scalar-opt \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic' \
  test check-pivots check-pivots-speed
```

Standalone ASan/UBSan (avoids the unrelated table-file serialization test):

```sh
make BUILD_DIR=build/pivots-sanitize-opt \
  CFLAGS='-O1 -g -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -fsanitize=address,undefined -fno-omit-frame-pointer' \
  LDFLAGS='-fsanitize=address,undefined' \
  build/pivots-sanitize-opt/bench-pivots \
  build/pivots-sanitize-opt/bench-pivots-speed
ASAN_OPTIONS=detect_leaks=0 build/pivots-sanitize-opt/bench-pivots --check-only
ASAN_OPTIONS=detect_leaks=0 build/pivots-sanitize-opt/bench-pivots-speed --check-only
```

Leak detection is disabled because of this environment's tracing limitations;
address and undefined-behavior checks remain enabled. The local DP metadata
copy has explicit 64-byte alignment, including under GCC's ASan fake stack.
