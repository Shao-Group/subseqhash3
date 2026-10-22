# Sliding windows over long sequences

`seeds-sub` and `seeds-ind` now use public long-sequence APIs. Both commands
emit every length-W window in increasing start order. Outward is the default;
`--inward` selects the alternative hash family. Outward reuses side DPs
across overlapping windows. Inward retains its two shared side DPs within
each window, rebuilding their scores for the next window.

```sh
build/subseqhash3 tables 50 1024 1 shared.piv --seed 42
build/subseqhash3 seeds-sub long.seq 50 8 1024 49 1 outward.csv shared.piv
build/subseqhash3 seeds-sub long.seq 50 8 1024 49 1 inward.csv shared.piv --inward
build/subseqhash3 seeds-ind long.seq 50 8 1024 42 1 outward-deletions.csv shared.piv
build/subseqhash3 seeds-ind long.seq 50 8 1024 42 1 inward-deletions.csv shared.piv --inward
```

The current R-per-set, S-set interface is documented in [SEEDING.md](SEEDING.md).
The C interfaces below remain available with their original semantics.

Input is one A/C/G/T sequence per line, with a separate window stream
for each line. The current CSV includes set and pivot indices. Tables,
edit budgets, scores, seed strings and deterministic witnesses match the
corresponding single-window APIs exactly. Choosing a different traversal
still changes the hash family and may change scores and seeds.

## Public C interfaces

Declared in `src/ssh3_pivots.h` and `src/ssh3_deletions.h`:

```c
int subseqhash3_seed_repeats_windows(
    const uint8_t *sequence, size_t length, size_t window,
    size_t substitutions, const SubseqHash3PivotTables *tables,
    size_t table_count, size_t repeats, SubseqHash3Traversal traversal,
    SubseqHash3RepeatsWindowCallback callback, void *context);

int subseqhash3_seed_deletion_windows(
    const uint8_t *sequence, size_t length, size_t window,
    size_t deletions, const SubseqHash3DeletionTables *tables,
    size_t count, SubseqHash3Traversal traversal,
    SubseqHash3DeletionWindowCallback callback, void *context);
```

Traversal is `SUBSEQHASH3_OUTWARD` or `SUBSEQHASH3_INWARD`. DNA input uses
symbols 0..3. Window and edit limits match the single-window APIs: repeats
uses 2..257 bases and exactly one deletion plus E substitutions; deletions
uses 1..256 bases and E deletions without substitutions. Input length must
be at least the window length. Repeats supports multiple independent table
families and any positive repetition count. Deletions requests the first
`count` pivot orders, with `1<=count<=window-deletions`.

Callbacks receive `(window_start, results, count, context)`, once per window.
Result arrays and their owned buffers are borrowed until the callback
returns; copy anything needed later. Result positions are relative to the
window, matching the existing APIs; add `window_start` for absolute positions.
Returning nonzero stops processing. The APIs return 0 on success and -1 on
invalid input, allocation/seeding failure, or callback cancellation/failure.
Earlier callbacks may already have run when a later failure occurs. All
input symbols are validated before the first callback. Tables must remain
immutable throughout the call. Each invocation owns its workspace and cache.

For example, a repeats callback can forward seeds without retaining the
whole output:

```c
static int consume(size_t start, const SubseqHash3EditResult *seeds,
                   size_t count, void *context) {
    /* Read seeds[0..count-1] here; positions are relative to start. */
    return 0;
}
/* initialized tables and encoded sequence supplied by the caller */
int status = subseqhash3_seed_repeats_windows(
    sequence, length, 50, 8, tables, table_count, 49,
    SUBSEQHASH3_OUTWARD, consume, context);
```

## Sharing across windows

For a fixed global input anchor, outward left states read bases immediately
to its left in reverse order; outward right states read bases immediately to
its right in normal order. Changing the window boundary only changes the
consumed prefix length. The recurrence and rank-indexed table rows for that
prefix are unchanged. We therefore retain every needed consumed-length/edit
state for each active anchor and join the appropriate states for each window.

Repeats needs anchors p and p+1 for output pivot p. For a family requesting
R leading pivots, its ring holds R+1 anchors. Deletions needs anchors p+h,
where h=0..E is the number of left deletions; its ring holds R+E anchors.
Advancing one window discards one old anchor and builds one new anchor per
family. Traceback reads the cached scores or parents and preserves the
single-window tie rules. Joins read cached vectors directly, avoiding a copy
of every boundary vector for every window.

After startup, both traversals build approximately two side DPs per window.
This removes outward's additional factor from rebuilding all requested
anchors separately in every window. Startup fills the anchor ring; its cost
is included in the measurements. Output allocation, traceback and callback
work remain necessary per window for both orders.

The cache trades memory for reuse. For one family, worst-case score storage
is O((R+1)*W*(E+1)*D) for repeats and O((R+E)*W*(E+1)*D) for deletions, compared
with O(W*(E+1)*D) for inward. Repeats caches each required independent family.
Memory does not grow with sequence length; results are freed after each
callback. The current CLI still reads and encodes one complete input line.

Deletion tables with coefficients outside [-127,127] retain the exact int64
single-window solver. Direct cases (zero deletions or one retained base),
and the existing residue-independent one-deletion crossover, also retain
their specialized solvers. Repeats with zero substitutions retains its
existing one-deletion crossover. These paths remain supported but do not use
the dense anchor cache. A single-window repeats input uses its existing
workspace path.

## Benchmark protocol

Measurements use the Xeon Gold 6148, a single thread pinned to CPU 0, with
identical -O3/AVX-512 compiler flags. Each condition fixes input, tables, W,
E, D and requested pivot count between inward and outward. Uniform synthetic
DNA strings have lengths 256 and 4096. Conditions use W=30/E=4 and W=50/E=8,
D=32/1024, and ten or all pivots. "All" means W-E seeds for deletions and W-1
for repeats. No residue-count adjustment is used to equalize runtime.

Each reported time is the median of three complete sequence calls, rotating
method order between samples. Timing includes initial workspace and anchor
cache allocations, startup windows, output allocation, exact traceback,
checksum callback, result destruction and cache cleanup. It excludes table
generation and CSV disk I/O. Reported units are microseconds per emitted
window, so total sequence time is this value times `length-W+1`.

For length 256, the benchmark also times the former outward loop, which
reuses workspace allocations but rebuilds all scores for every window.
Before timing each condition, exact output checksums over the first 64
windows are compared to the independent-window APIs for each traversal.
When the independent loop is timed, its complete-sequence outward checksum
must equal the shared implementation's checksum. Inward and outward are
different families and are not expected to have equal checksums.

## Measured results (2026-10-21)

These are the original shared-window measurements. For the subsequent
outward optimization and its before/after timings, see
[Outward speed improvements](OUTWARD_SPEED.md).

Sharing reduces outward runtime, but the measured difference from inward
does not disappear. On the same length-256 input, sharing accelerates outward
by **2.15–3.26x for deletions** and **2.44–4.50x for repeats**. On length-4096
inputs, shared outward still costs **1.13–2.22x inward for deletions** and
**1.17–2.43x inward for repeats**.

For W=30, D=1024 and ten seeds, the long-sequence gap is only 13% for
deletions and 17% for repeats. With all pivots at W=50, D=1024, the remaining
factors are 2.22x and 2.43x respectively. Long input amortizes startup, while
requesting many pivots retains a much larger active cache. These results
support the predicted removal of outward's extra DP-building factor, but
not equal measured runtimes. Larger working sets, extra cache lookups and
different traceback access patterns remain; identifying their individual
contributions would require profiling.

At W=50/E=8/D=1024 with all pivots, the outward anchor score arrays alone
occupy about 181 MiB for deletions and 363 MiB for repeats, before other
workspace allocations. The corresponding two inward score arrays occupy
about 3.6 MiB and 7.3 MiB. This is a substantial memory tradeoff.

### Length 4096: inward versus shared outward

Times are median microseconds per emitted window; both columns use identical
inputs, tables and parameters within each row.

| Algorithm | W | E | D | Seeds | Inward us | Shared outward us | Outward / inward |
|:---|---:|---:|---:|---:|---:|---:|---:|
| deletions | 30 | 4 | 32 | 10 | 8.37 | 9.98 | 1.19x |
| repeats | 30 | 4 | 32 | 10 | 11.38 | 13.37 | 1.18x |
| deletions | 30 | 4 | 32 | 26 | 17.91 | 24.31 | 1.36x |
| repeats | 30 | 4 | 32 | 29 | 27.36 | 39.58 | 1.45x |
| deletions | 30 | 4 | 1024 | 10 | 103.08 | 116.64 | 1.13x |
| repeats | 30 | 4 | 1024 | 10 | 132.63 | 154.85 | 1.17x |
| deletions | 30 | 4 | 1024 | 26 | 172.29 | 271.04 | 1.57x |
| repeats | 30 | 4 | 1024 | 29 | 284.25 | 589.29 | 2.07x |
| deletions | 50 | 8 | 32 | 10 | 14.92 | 20.56 | 1.38x |
| repeats | 50 | 8 | 32 | 10 | 19.98 | 26.97 | 1.35x |
| deletions | 50 | 8 | 32 | 42 | 43.89 | 67.34 | 1.53x |
| repeats | 50 | 8 | 32 | 49 | 68.10 | 142.84 | 2.10x |
| deletions | 50 | 8 | 1024 | 10 | 205.37 | 344.37 | 1.68x |
| repeats | 50 | 8 | 1024 | 10 | 311.36 | 527.65 | 1.69x |
| deletions | 50 | 8 | 1024 | 42 | 416.06 | 923.86 | 2.22x |
| repeats | 50 | 8 | 1024 | 49 | 876.09 | 2125.38 | 2.43x |

### Length 256: benefit of outward sharing

Both methods seed the same sequence with the same outward family. The old
loop includes reusable allocations; the shared call additionally reuses DP
scores across windows. Cold-cache startup and cleanup are included in both.

| Algorithm | W | E | D | Seeds | Independent outward us | Shared outward us | Sharing speedup |
|:---|---:|---:|---:|---:|---:|---:|---:|
| deletions | 30 | 4 | 32 | 10 | 34.49 | 12.49 | 2.76x |
| repeats | 30 | 4 | 32 | 10 | 51.87 | 16.08 | 3.23x |
| deletions | 30 | 4 | 32 | 26 | 77.20 | 27.95 | 2.76x |
| repeats | 30 | 4 | 32 | 29 | 141.50 | 42.54 | 3.33x |
| deletions | 30 | 4 | 1024 | 10 | 361.70 | 152.21 | 2.38x |
| repeats | 30 | 4 | 1024 | 10 | 882.51 | 196.19 | 4.50x |
| deletions | 30 | 4 | 1024 | 26 | 820.97 | 378.02 | 2.17x |
| repeats | 30 | 4 | 1024 | 29 | 2433.67 | 784.29 | 3.10x |
| deletions | 50 | 8 | 32 | 10 | 91.97 | 28.25 | 3.26x |
| repeats | 50 | 8 | 32 | 10 | 111.80 | 35.87 | 3.12x |
| deletions | 50 | 8 | 32 | 42 | 315.52 | 97.38 | 3.24x |
| repeats | 50 | 8 | 32 | 49 | 536.70 | 220.11 | 2.44x |
| deletions | 50 | 8 | 1024 | 10 | 1061.38 | 494.36 | 2.15x |
| repeats | 50 | 8 | 1024 | 10 | 2552.87 | 673.77 | 3.79x |
| deletions | 50 | 8 | 1024 | 42 | 3547.36 | 1516.56 | 2.34x |
| repeats | 50 | 8 | 1024 | 49 | 12150.54 | 3231.11 | 3.76x |

## Validation

The scalar and AVX-512 full test suites and the address/undefined-behavior
sanitizer suite pass. The window suite compares exact scores, seed bytes,
positions and substitution flags against the independent single-window
implementations for 1,000 configurations and 18,756 windows. Cases include
D=1/2/7/31/32/33/128/129/1024, partial vectors, ring wraparound, random,
homopolymer and periodic inputs, partial/all pivots, multiple repeats
families, E=0 and maximal edit budgets, wide imported magnitudes, the
one-deletion crossover, W=1/256 deletions and W=257 repeats, callback
cancellation and invalid input. Existing CLI/API CSV-equivalence tests now
exercise the long-sequence path. Leak detection is disabled in the sanitizer
run because of this environment's ptrace restriction.

## Reproduce

```sh
make BUILD_DIR=build/opt \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test build/opt/bench-windows
python3 bench/run_windows.py --executable build/opt/bench-windows --cpu 0 \
  --flags '-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512'
```

* [Median runtime results](windows-runtime.csv).
* [Individual samples and output checksums](windows-samples.csv).
* [Protocol, machine and source hashes](windows-metadata.json).
* [Benchmark](windows.c), [runner](run_windows.py) and
  [exact window checks](../tests/windows.c).
