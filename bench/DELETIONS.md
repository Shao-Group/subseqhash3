# Exactly E deletions, up to k SubseqHash2 orders

The separate `src/ssh3_deletions.h` / `src/ssh3_deletions.c` API produces
length-k seeds, where k=n-E, with **exactly E deletions and no substitutions**.
One call returns minimizers for up to k requested pivot orders, sharing its
side DPs across the orders. Existing one-substitution and mixed-edit pivot
algorithms are unchanged. Both commands now use the same binary pivot tables
and default generator; the deletion kernel widens their parameters losslessly
to its internal 64-bit representation.

## CLI

```sh
make simd
# Generate once, reuse the exact same file in both commands.
build/simd/subseqhash3 tables 30 11 10 shared.piv --pivots --seed 42
build/simd/subseqhash3 repeats input.seq 30 5 11 10 repeats.csv shared.piv
build/simd/subseqhash3 deletions input.seq 30 5 11 10 deletions.csv shared.piv
# Or request all k=25 deletion-only orders from that same first family.
build/simd/subseqhash3 deletions input.seq 30 5 11 25 all.csv shared.piv
# Original SubseqHash2 text-table import remains an explicit option.
build/simd/subseqhash3 deletions input.seq 30 5 11 25 original.csv \
  --subseqhash2-tables original-table.txt
```

Use `make` and `build/subseqhash3` for the scalar build. AVX-512F/BW/VL
builds require compatible hardware; there is no runtime CPU dispatch. All
seeding is single-threaded.

```text
subseqhash3 deletions INPUT W E D R OUTPUT [TABLE_FILE] [--seed SEED]
                     [--subseqhash2-tables TEXT_FILE]
```

W is 1..256, 0<=E<W, D>0, and 1<=R<=W-E. R selects numeric pivots
0..R-1. This is a subset of one family, not independent repetitions.
Different orders can produce identical strings. The text table's final
permutation is validated but does not reorder the numeric pivots; to
reproduce the original program's subsampling, pass its selected pivots to
the C API explicitly.

TABLE_FILE is the **same version-1 binary `.piv` format** produced by
`tables --pivots` and accepted by `repeats`. Existing files remain valid.
A file generated with `tables n D R ... --pivots` has n-1 rows per family:

- `repeats` requires n>=W, matching D, and enough families for its repetitions.
- `deletions` requires n-1>=W-E and matching D. It uses only the first family;
  additional families are harmless. Its R can be any value up to W-E,
  independently of the R used when creating the file.
- Zero-deletion windows need W retained rows, so generate n>=W+1. For example,
  `tables 31 11 10 shared.piv --pivots` supports W=30,E=0 deletion seeding.

The file dimension is row capacity, not an edit budget. A table with n=30
supports deletion-only length k<=29, including W=30,E=5 (k=25). The shared
parameters do not change the existing algorithms' different side traversals;
using one file does not make their minimizers identical.

INPUT contains one A/C/G/T sequence per line, either case, with LF or CRLF
endings. Empty lines, non-DNA symbols (including N), and lines shorter than
W are rejected. Windows advance one base. Tables are initialized once and
reused for all windows and input sequences. Without TABLE_FILE, shared tables
are created in `OUTPUT.piv`, just as with `repeats`. The saved file can be
supplied explicitly to reproduce the run. `--seed` selects reproducible
creation; otherwise the randomized generation seed is printed to stderr.
Automatic creation stores one family with n=W for E>=1, and n=W+1 for E=0.
For E>=1, the resulting file also supports `repeats` at the same W.

For example:

```sh
build/simd/subseqhash3 deletions input.seq 30 5 11 10 first.csv --seed 42
build/simd/subseqhash3 deletions input.seq 30 5 11 10 again.csv first.csv.piv
build/simd/subseqhash3 repeats input.seq 30 5 11 10 mixed.csv first.csv.piv
```

To use an original SubseqHash2 text table, supply `--subseqhash2-tables FILE`
with k=W-E and matching D. This replaces the earlier positional-text-file
syntax. TEXT_FILE, a positional shared TABLE_FILE, and `--seed` are mutually
exclusive. The legacy `tables` binary format for seeds/collisions remains
separate and cannot supply pivot orders.

Original text magnitudes must be positive and at most `INT64_MAX/(4*k)`;
signs, residues, final permutation, integer overflow, truncation, and trailing
data are checked. Import preserves the original large parameters without
rescaling. Shared tables use magnitudes 1..127, exactly as `repeats` does.
Consequently the default deletion-only parameters and `--seed` outputs have
changed from the earlier large-magnitude generator; import original text
parameters when reproducing those orders. The recorded large-magnitude
benchmark fixtures remain separate from production table generation.

CSV columns are:

```text
pivot,sequence,window_start,deletions,seed,psi,omega
```

Indices are zero-based. Rows are sequence-major, window-major, then
pivot-major. Seed strings are uppercase and have W-E bases. Each length-L
input line emits `(L-W+1)*R` data rows. Psi and omega are exact integers;
omega is never converted through floating point. The C API additionally
returns the source positions. OUTPUT must not alias INPUT or TABLE_FILE (or TEXT_FILE),
including hard links. Invalid table files are rejected before OUTPUT is
opened; later input or I/O failures can leave a partial CSV. Wall time
printed to stderr includes table initialization/loading and CSV I/O.

## C API

Include `ssh3_deletions.h`. For a 30-base encoded input `sequence` (symbols
0..3), initialize tables with maximum seed length 25, then request all orders:

```c
SubseqHash3DeletionTables tables = {0};
SubseqHash3DeletionResult results[25] = {{0}};
SubseqHash3DeletionWorkspace *workspace =
    subseqhash3_deletion_workspace_create();
int status = -1;
if (workspace && !subseqhash3_deletion_tables_init(&tables, 25, 11, 42))
    status = subseqhash3_compute_deletion_seeds_with_workspace(
        sequence, 30, 5, &tables, NULL, 25, results, workspace);
/* On success, consume results here; keep tables/workspace for more windows. */
subseqhash3_deletion_results_destroy(results, 25);
subseqhash3_deletion_workspace_destroy(workspace);
subseqhash3_deletion_tables_destroy(&tables);
```

Table initialization/import and seed-computation functions return 0 on success
and -1 on invalid input/allocation failure.
Pass empty results and table objects. `deletion_tables_init(k,D,seed)` uses
exactly the parameters of `pivot_tables_init(k+1,D,seed)`. To use an already
loaded shared family, call `subseqhash3_pivot_table_file_open()`, then
`subseqhash3_deletion_tables_from_pivots(&tables, file.tables)`. This creates
independently owned arrays; the shared file object may then be closed.
No score rescaling or change of side signs/residues occurs.

Tables can cover a larger k than the current output; D is fixed by the tables. An explicit pivot array can select
any subset in any order, including duplicate pivots; NULL selects 0..count-1.
The convenience `subseqhash3_compute_deletion_seeds()` owns a temporary
workspace. Reusing the workspace avoids repeated DP allocation; workspaces
must not be shared by concurrent threads.

## Original order definition

This path preserves the original SubseqHash2 side directions, unlike the
existing one-deletion/E-substitution API's documented direction variant.
For a retained string y, pivot p, and tables t:

```text
L = y[p-1], y[p-2], ..., y[0]       (reversed left)
R = y[p+1], y[p+2], ..., y[k-1]     (normal right)
psi_side starts at 0; omega_side starts at 0
psi_side   = (psi_side + C[j,base]) mod D
omega_side = B1[j,base,psi_side] * omega_side + signed_A[j,base,psi_side]
psi_p      = (psi_left + pivot_C[p,y[p]] + psi_right) mod D
omega_p    = left_sign[p,y[p]] * omega_left + pivot_A[p,y[p]]
             + right_sign[p,y[p]] * omega_right
```

The minimizer has the smallest psi, then the largest omega. `signed_A`
includes the original B2 sign; pivot_A includes combine3. Right-side tables
are imported from A/B/C1, left-side tables from revA/revB/C2; combine1 is the
right sign, combine2 the left sign. Scores use int64_t throughout, preserving
the original large integer magnitudes beyond double precision. Bounds ensure
that all finite intermediate and combined values stay below INT64_MAX/4.

On score ties, any minimizer is valid. This implementation prefers the earliest
input pivot, then the smallest normal right residue; side DP ties prefer
skipping the current input base. Scalar and SIMD witnesses agree. No general
claim of identical positions on every tie is made for the original program.
The command emits every window; it does not apply SubseqHash2's cross-window
seed deduplication, reverse-complement handling, or score thresholds.

## Shared algorithm and acceleration

A seed pivot p can originate only at input i=p+h, 0<=h<=E. Its left side
retains p of i bases (h deletions), and its right side retains k-p-1 of
n-i-1 bases (E-h deletions). Each needed input anchor is evaluated once per
side, and its boundary scores serve every requested output pivot.

To preserve the original directions, left DP scans backwards from its
anchor and right DP scans forwards. A state has retained length j, deletion
count h, and residue q. It stores both minimum and maximum omega because B1
can reverse the score ordering. Only 0<=h<=E states are evaluated, and side
lengths are capped by the requested pivot range. The recurrence either skips
the current input base or retains that exact base; it never enumerates four
replacement bases or tracks substitution counts.

Only two retained-length score planes remain live. Completed boundary
scores are saved for joins; intermediate states retain two traceback choice
bits in a byte, rather than full score/parent structs. A reusable workspace
amortizes allocations across windows. For E=0, the fixed sequence is evaluated
directly. For k=1, a direct input scan selects the best retained base.

The AVX-512 path processes eight exact 64-bit residues per vector. Circular
padding converts modular predecessor accesses to contiguous loads. Boundary
right residues are reflected, so joins use two contiguous circular loads
instead of gathers. Joins try psi=0 first. If unreachable for all feasible
input pivots, a circular predecessor-reachability pass finds the exact smallest
psi in O(D); there is no D-squared residue search and no D<=31 bitmap limit.
Vector maxima defer witness-residue selection until the winning input pivot
is known. Arbitrary positive D, including partial vectors and D=1025, works
in both builds.

Worst-case work is O(n*k*(E+1)*D + R*(E+1)*D + R*n), with smaller side work
for subsets of pivots. Boundary scores use O(n*(E+1)*D) 64-bit entries;
traceback uses O(n*k*(E+1)*D) bytes. These retain the original orders rather
than exchanging side directions to reduce the asymptotic number of states.

## Correctness and regression checks

```sh
make test
make BUILD_DIR=build/deletions-simd \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test check-one check-pivots check-pivots-single check-deletions-reference
```

`make test` includes 1,620 exhaustive small cases (n=1..9, every deletion
budget, D=1/2/7/8/9/31/32/33/65, and four random/adversarial configurations).
Every retained subsequence is scored independently of the DP. Tests cover
homopolymers, ties, unreachable psi=0, subsets including duplicate/out-of-order
pivots, workspace reuse across dimensions, n=256, E=0/n-1, D=1024/1025,
64-bit witnesses, malformed/overflowing text tables, and CSV/API equivalence
across multiple sequences and windows. CLI argument and overwrite guards are
also checked. Rejection diagnostics are captured in `BUILD_DIR/deletions-test.log`.

The optional original-source comparison requires a SubseqHash2 checkout at
`subseqhash2/`, or `SUBSEQHASH2_DIR=/path/to/checkout`. It verifies scores and
seed strings on 128 windows over eight configurations. Another 16 windows
use the distributed `subseqhash2/alignment/table_subseq2/25_11` table when
present. The reference sources are compiled unchanged; forced standard-library
includes compensate for their missing standalone header includes. Tests and
benchmarks use temporary files under /tmp and remove them on success.

ASan/UBSan (repeat with AVX-512 flags in CFLAGS to check SIMD):

```sh
ASAN_OPTIONS=detect_leaks=0 make BUILD_DIR=build/deletions-sanitize \
  CFLAGS='-O1 -g -std=c99 -Wall -Wextra -Wpedantic -fsanitize=address,undefined -fno-omit-frame-pointer' \
  LDFLAGS='-fsanitize=address,undefined' check-deletions
```

Shared-table regression checks additionally compare every widened parameter,
use one file with both CLI commands, and compare default/reused-file CSVs
against direct API results. They check n-1>=W-E (including tables smaller than
W), E=0 generation at n=W+1, mismatched D, wrong file formats, automatic-file
hard-link aliases, source lifetime, and minimum/maximum row capacities. At D=1024, an
executable-level test reused one file across both commands: 140 CSV rows per
command matched between scalar and SIMD, and automatic generation reproduced
the explicit file and deletion CSV byte-for-byte.

On 2026-10-20 the scalar and SIMD full test suites, existing single-seed/pivot
checks, and scalar/SIMD ASan/UBSan deletion tests passed. Leak detection was
disabled because LeakSanitizer cannot run under this environment's tracing.
Thirty executable-level CSV comparisons also matched byte-for-byte between
scalar and SIMD at W/E=1/0,12/0,12/1,12/3,12/11,30/5 and
D=1,11,31,33,1024, including original-table CLI loading.

## Measured speed against SubseqHash2

Measurements are single-threaded on an Intel Xeon Gold 6148, CPU 0, with
-O3 and the same AVX-512F/BW/VL/Skylake flags for both implementations.
The original checkout is commit `8d485c66dec44a4f32c4e6f0dd81402f46f82bbb`;
`src/subseqhash2seeding.cpp` SHA-256 is
`da2f31cf69ef458af421abe1e8e6fed900857d50ee5494d04a5de76183640196`.

```sh
make BUILD_DIR=build/deletions-simd \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  build/deletions-simd/bench-deletions
taskset -c 0 build/deletions-simd/bench-deletions \
  --samples bench/deletions-samples.csv > bench/deletions-results.csv
```

These recorded results use the original large-magnitude fixture generator
from commit `3820eaf`. `bench/deletions.cpp` retains that generator privately
to keep these comparisons reproducible after production generation was unified;
these are not timing claims for the new default small-magnitude tables.

Each case uses 16 deterministic input sequences, with four passes per sample
and five alternating old/new samples; times are medians in microseconds per
sequence. Both paths receive identical table parameters and include seed
materialization and destruction. Table loading, first allocation/warmup,
and correctness checks are outside timing; the new DP workspace and original
DP/vector buffers are reused. The original calls getSubseq2Seeds on one
length-n sequence, so deduplication across windows is irrelevant. Its buffers
are reduced to the number of origins actually addressed for that input,
without changing its algorithm or source. The original still builds its DPs
for the entire window when only one order is requested.

D=11 matches the distributed reference tables. Original combine uses a signed
32-bit convolution bitmap and has invalid/truncated shifts for larger D;
we do not extrapolate these speedups to D=1024 or compare against that broken
case. The original k=1 boundary handling is also unsuitable as a correctness
reference; the new k=1 path is checked against exhaustive enumeration instead.

Recorded medians (microseconds per sequence):

| n | E | D | Seeds | SubseqHash2 | New | Speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 5 | 11 | 1 | 1174.192 | 21.136 | 55.55x |
| 30 | 5 | 11 | 10 | 1187.127 | 60.158 | 19.73x |
| 30 | 5 | 11 | 25 | 1185.597 | 120.036 | 9.88x |
| 50 | 8 | 11 | 1 | 4871.593 | 74.252 | 65.61x |
| 50 | 8 | 11 | 10 | 4907.484 | 166.342 | 29.50x |
| 50 | 8 | 11 | 42 | 4987.228 | 464.065 | 10.75x |
| 20 | 1 | 11 | 19 | 116.511 | 25.875 | 4.50x |
| 20 | 18 | 11 | 2 | 35.814 | 14.898 | 2.40x |

Full median results are in [deletions-results.csv](deletions-results.csv),
with all timing samples in [deletions-samples.csv](deletions-samples.csv).
Speedups vary with the deletion budget and requested count; small-window
edge cases have smaller gains than the two main configurations. These are
single-sequence seeding measurements, not whole-genome sliding-window speedups.
