# Exactly one substitution: correctness and timing

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
