# Exactly one substitution: correctness and timing

For the current pivot-only command line, table generation and R/S contract,
see [SEEDING.md](SEEDING.md). Legacy algorithm and CLI examples below record
the historical interface; the low-level C functions remain available.

For current single-thread `seeds-sub`/`seeds-ind` runtime improvements, exact
output comparisons and end-to-end measurements, see [SEED-RUNTIME.md](SEED-RUNTIME.md).

`subseqhash3_compute_one_substitution(s, n, E, tables, result)` computes
exactly one substitution and exactly `E` deletions. The corresponding
reference call is `subseqhash3_compute_all_edit_splits(s, n, E+1, tables,
results)`, selecting `results[1]`. Both calls use the same input and tables.

The implementation preserves the reference seed, score, positions,
substitution flags, and tie breaking. Existing all-splits functions are
unchanged.

## Reproduce

On an AVX-512F/BW/VL machine:

```sh
make BUILD_DIR=build/fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test bench-one
```

For scalar correctness:

```sh
make test check-one
```

The standalone benchmark performs 365 exact comparisons before timing, plus
32 comparisons for each timed configuration. It covers random and identical
sequences, non-power-of-two residue counts, invalid inputs, extreme scores,
ties, and unreachable high residues. Each timed configuration uses 32 random
sequences, repeated four times per sample. Results are medians of three
samples with alternating call order. Timing includes DP allocation,
traceback, and result cleanup, and excludes table initialization. Both
implementations use the same compiled object and compiler options.

The standalone check also avoids the existing table-file serialization
code, allowing it to run independently under AddressSanitizer/UBSan:

```sh
make BUILD_DIR=build/sanitize-one \
  CFLAGS='-O1 -g -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -fsanitize=address,undefined -fno-omit-frame-pointer' \
  LDFLAGS='-fsanitize=address,undefined' build/sanitize-one/bench-one
ASAN_OPTIONS=detect_leaks=0 build/sanitize-one/bench-one --check-only
```

Leak detection is disabled in this command because LeakSanitizer does not
work under the tracing used by this environment; address and undefined
behavior checks remain enabled.

## Algorithm

The vector path indexes rows by **retained length** `j`, stopping at `n-E`,
and columns by deletion count `h=0..E`. Each cell has two substitution
planes, zero and one. A retained base comes from `(j-1,h)`; a deletion
comes from `(j,h-1)`. Processing residues outside the deletion loop shares
the four symbol-table vectors across all deletion counts.

Scores are represented as two minimum costs:

```text
u = min(127*j + omega)
v = min(127*j - omega)
```

Let `a=A*B2` for the current symbol, retained position, and residue. For
positive `B1`, add `127+a` to `u` and `127-a` to `v`. For negative `B1`,
swap the source costs first. All increments are in `0..254`, so unsigned
16-bit saturating addition preserves the unreachable sentinel `65535`
without separate reachability masks. The existing table length limit
ensures every finite cost is at most `65024`.

For `D>=512`, the last three retained rows use memoized queries starting at the largest
possible output residue. Only the dependencies needed for each queried
residue are evaluated; if it is unreachable, the next smaller residue is
queried. No assumption that the largest residue is reachable is required.
For smaller residue spaces a full vector sweep is cheaper and is used instead.

The vector path supports every `D>=32` through masked partial vectors and
circular padding. Smaller residue counts and builds without AVX-512 use
the previous fixed-split fallback. If only one base remains, a direct
`O(n)` scan replaces the DP entirely.

## Performance limits

Speedup depends strongly on `E`, `n`, and `D`. In measurements on the Xeon
Gold 6148, the `n=50, E=8, D=1024` case remains roughly 3.3--4 times faster
than all-splits, rather than 10 times. Larger deletion budgets can exceed
10 times; `n=64, E=32, D=1024` measured over 25 times faster. These are
single-threaded, exact results. Use the benchmark output to assess the
parameters of interest instead of assuming a universal 10-times gain.

Example run on 2026-10-20 (microseconds per call, median of three samples):

| n | E deletions | D | One substitution | All splits | Speedup |
|---:|---:|---:|---:|---:|---:|
| 20 | 5 | 1024 | 50.397 | 138.374 | 2.75x |
| 30 | 7 | 1024 | 107.719 | 352.718 | 3.27x |
| 50 | 8 | 1024 | 237.421 | 793.963 | 3.34x |
| 50 | 16 | 1024 | 395.838 | 2559.437 | 6.47x |
| 50 | 24 | 1024 | 478.971 | 4311.760 | 9.00x |
| 64 | 32 | 1024 | 850.836 | 24115.928 | 28.34x |
| 50 | 8 | 32 | 7.421 | 36.110 | 4.87x |
| 50 | 8 | 256 | 51.884 | 191.337 | 3.69x |
| 50 | 8 | 1025 | 248.854 | 2041.571 | 8.20x |

# Deletion kernel versus its previous version and repeats

```sh
make BUILD_DIR=build/del-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test bench-deletions-speed
```

`check-deletions-speed` runs only the exact comparisons. No external
SubseqHash2 checkout is required. The benchmark checks 160 random windows
against the immutable `fab085f` deletion kernel, including exact seeds and
positions. Timings use 16 inputs, three passes per sample, five samples,
rotating execution order, and the median. Both deletion implementations
use identical flags, tables, inputs, pivot counts, and reusable workspaces.
Timing includes result allocation, traceback, and cleanup; table creation
is excluded.

Shared pivot tables with side magnitudes at most 127 use 16-bit minimum
costs `127*j + omega` and `127*j - omega`. The AVX-512 kernel updates 32
residues per vector with unsigned saturating addition. It packs table scores
once per side per call, handles the zero-deletion and single-retained-base
states sparsely, and limits each anchor to deletion budgets needed by the
requested pivots. Tables with larger magnitudes retain exact 64-bit scores.
Table magnitudes are checked on every call, so workspace reuse also works
when switching tables or changing their values.

The repeats columns below are context, not an equivalence check: repeats
uses the inward C API for one deletion and E substitutions; deletions uses
the outward C API for E deletions. Both CLI commands now default to outward
and accept `--inward` to select the inward family. The inward
order shares side DPs across pivots. The outward order needs separate side
DPs at each input pivot, so these optimizations do not make deletion timing
match repeats. Changing traversal would change the hash order and outputs.

Measured on 2026-10-21, Xeon Gold 6148, single thread, flags above
(microseconds per call; full results in `deletions-speed-results.csv`):

| n | E deletions | D | Seeds | Previous | Optimized | Speedup | Repeats |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4 | 1024 | 10 | 1761.755 | 651.749 | 2.70x | 143.742 |
| 30 | 4 | 1024 | 26 | 3849.156 | 1417.835 | 2.71x | 283.378 |
| 50 | 8 | 1024 | 10 | 6276.738 | 1688.124 | 3.72x | 342.306 |
| 50 | 8 | 1024 | 42 | 18363.545 | 6126.527 | 3.00x | 848.207 |

For the alternative inward deletion family, see the [implementation, runtime
and paired collision experiment](DELETIONS_INWARD.md).

For outward repeats with one deletion and E substitutions, see the
[implementation, runtime and paired collision experiment](REPEATS_OUTWARD.md).

The [repeats residue sweep](REPEATS_RESIDUES.md) compares runtime and
seed-collision probabilities at D=32, 128, 1024 and 2048.

The [matched-runtime traversal comparison](TRAVERSAL_MATCHED.md) adjusts D
separately for inward and outward at fixed seed count, for both deletions
and repeats.

The [single-thread traversal optimization](TRAVERSAL_OPTIMIZATION.md) records
before/after timings for all four paths, including the exact single-deletion
fast path and the remaining gap to 5x on E=4/8 workloads. Reproduce the
comparison with `make bench-traversal-speed`; verify exact baseline parity
with `make check-traversal-speed`.

The [long-sequence sliding-window comparison](WINDOWS.md) measures inward
against outward with anchor DPs shared across overlapping windows, plus the
previous independent-window loop. Run `make bench-windows` for the sweep;
`make check-windows` checks exact scores and witnesses for the window APIs.

The [outward speed improvements](OUTWARD_SPEED.md) compare the subsequent
optimization against the existing shared-window implementation, covering
fresh and reused single-window workspaces and length-4096 sequences for both
deletions and repeats. Run `python3 bench/run_outward_speed.py --help` for
the paired before/after sweep, including complete output checksum checks.
