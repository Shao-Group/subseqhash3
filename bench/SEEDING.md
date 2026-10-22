# Pivot-based seeding interface

## R pivots from S independent sets

```text
subseqhash3 tables W D S OUTPUT [--seed SEED]
subseqhash3 seeds-sub INPUT W E D R S OUTPUT [TABLE_FILE] [--inward] [--seed SEED]
subseqhash3 seeds-ind INPUT W E D R S OUTPUT [TABLE_FILE] [--inward] [--seed SEED]
subseqhash3 collisions sub|ind PAIRS_FILE E D R S [TABLE_FILE] [--inward] [--seed SEED] [--jaccard]
```

W sizes the tables for the target window. D is positive. S is the positive
number of independent table sets, generated from distinct deterministic
random streams when a generation seed is specified. R is positive and
selects pivots 0..R-1 from EACH of the first S sets. The total per window is
R*S; the R orders within a set share its two side tables and pivot parameters.
Duplicate output seed strings are possible.

`seeds-sub` uses one deletion and E substitutions; W is 2..257, E is 0..W-1,
and R<=W-1. `seeds-ind` uses E deletions and zero substitutions; W is 1..256,
E is 0..W-1, and R<=W-E. Budgets are exact. Outward is the default;
`--inward` changes the hash traversal and can change scores and seeds.

Input is one nonempty DNA sequence per line; A/C/G/T are case insensitive.
Each line must be at least W bases. Windows advance one base and are emitted
in increasing start order. All lines reuse the same sets of tables.

## Tables

```sh
build/subseqhash3 tables 30 1024 3 shared.piv --seed 42
build/subseqhash3 seeds-sub dna.seq 30 4 1024 5 3 sub.csv shared.piv
build/subseqhash3 seeds-ind dna.seq 30 4 1024 5 3 ind.csv shared.piv
```

Both commands produce 15 seeds per window, from pivots 0..4 of sets 0..2.
Changing R while reusing shared.piv selects a different prefix of pivot
orders in every set, without changing the realized table parameters.

The CLI stores min(W,256) rows per set, enough for either mode at the target
W, including E=0 in seeds-ind. For W=257 only seeds-sub is supported, and its
output needs 256 rows. Larger table capacities can be reused:

* seeds-sub requires at least W-1 rows per set.
* seeds-ind requires at least W-E rows per set.
* Both require matching D and at least S sets. Extra sets are unused.

The version-1 pivot binary layout is retained: a 64-byte little-endian
header followed by S independently owned table sets. The header's internal
`n` remains row_capacity+1, so CLI `tables W` normally stores n=W+1, capped
at 257. This is an internal row-capacity field, not an additional CLI window
parameter. Existing pivot files are accepted when their dimensions suffice;
legacy .abc and SubseqHash2 text files are not CLI table formats.

Omitting TABLE_FILE creates OUTPUT.piv, using the same dimensions as
`tables W D S`. `--seed` is accepted only with automatic creation. The
chosen seed is printed to stderr. Reusing a table file reproduces the same
orders, independently of the time or random state of later invocations.

## CSV contract

```text
set,pivot,sequence,window_start,substitutions,deletions,seed,psi,omega
```

All index columns are zero-based. Rows sort by sequence, window_start, set,
then pivot. Each window contains S contiguous groups of R rows. For W=30,
R=5, S=3, the last group is set=2 and pivot=0..4.

`pivot` is a position in the resulting seed, not a deletion location or an
absolute input coordinate. `sequence` is the input line index. window_start
is an offset into that line. The seed has length W-1 for seeds-sub and W-E
for seeds-ind. psi and omega are the selected hash scores: minimize psi,
then maximize omega. The CSV emits the integer-valued scores exactly within
the supported table magnitudes; the indel C result type is int64_t.

A sequence of length L emits `(L-W+1)*R*S` rows. Output is streamed one
window at a time; the whole sequence line is currently read into memory.

## Collision metrics

`collisions sub` applies the seeds-sub edit model to both strings in every
pair; `collisions ind` applies seeds-ind. There is no sliding-window length
argument: each complete string is one window. R must fit both output lengths
for every pair. Tables must cover the maximum required row count. When
automatically generated, the target W is the maximum input length, and the
file is PAIRS_FILE.piv. Empty files and invalid pairs are rejected.

Each input pair contributes once to each reported probability:

* `any_seed_collision_probability`: any equal seed string across the two
  R*S-element seed sets, including different set/pivot indices.
* `any_same_set_pivot_collision_probability`: at least one aligned
  `(set,pivot)` matches.
* `all_same_set_pivot_collision_probability`: every aligned seed matches.
* `set_F_pivot_P_collision_probability`: that specific aligned order matches.

Strings producing unequal seed lengths do not match. `--jaccard` reports
the mean and population standard deviation of exact edit-neighborhood
Jaccard over input pairs, using the selected edit model. It is independent
of R, S, tables and traversal, and may be expensive.

## C API and lifetime

Include `src/ssh3_seeding.h`. Both entry points have the common argument order:

```c
subseqhash3_seeds_sub(sequence, length, W, E, tables, R, S,
                     SUBSEQHASH3_OUTWARD, sub_callback, context);
subseqhash3_seeds_ind(sequence, length, W, E, tables, R, S,
                     SUBSEQHASH3_OUTWARD, ind_callback, context);
```

`tables` points to S valid SubseqHash3PivotTables objects. DNA is encoded as
0..3. Each callback receives `(window_start, results, R*S, context)`; index
`f*R+p` is set f, pivot p. Results and their buffers are borrowed until the
callback returns. Copy anything retained afterwards. Positions are relative
to the window. Return nonzero from a callback to cancel; APIs return 0 on
success and -1 on invalid input, allocation/seeding failure or cancellation.
Input DNA and all table sets are validated before emitting the first window.
If a later window fails, earlier callbacks may already have run.

Calls own their caches/workspaces, with memory bounded independently of
sequence length. Outward keeps a cache per set; increasing S increases DP
memory accordingly. Indel conversion from shared pivot tables occurs once
per set per call. Passing length==W produces a single window.

`subseqhash3_seed_tables_create(path,W,D,S,seed)` is the matching W-sized
file-generation API. Older APIs remain available, including the low-level
`subseqhash3_pivot_table_file_create(path,n,D,count,seed)` with n-1 row
capacity, and original repeat-batching functions with their old semantics.
The new R-per-set contract applies to the new entry points and commands.

## Migration

* Remove `--pivots` from tables commands; use W D S OUTPUT, where S is an
  actual independent-set count rather than the old requested-repeat count.
* Replace `repeats` with `seeds-sub`, explicitly providing R then S.
  Old repetitions that filled a whole family before starting another are
  replaced by a fixed prefix of R pivots from every selected set.
* Replace `deletions` with `seeds-ind`, adding S (often 1).
* Replace the `repeat,family,pivot` CSV prefix with `set,pivot`.
* Replace legacy `collisions` calls with `collisions sub` or `collisions ind`
  and explicit R and S. The reported metrics now refer to pivot seeds.
* `seeds` and `--subseqhash2-tables` are removed from the CLI. Legacy C
  kernels and the text import function remain available.

## Validation

`make test` includes `tests/seeding.c`, which compares R*S callback batches
against independent per-set C API calls for every tested window, including
exact scores and witnesses, both traversals, partial/all pivots, zero and
maximal budgets, ring wraparound, partial SIMD vectors and callback failures.
`tests/cli.py` checks executable dispatch, CSV ordering and counts, independent
set parameters, auto/reused files, aliases, bounds, removed options, and
collision statistics calculated directly from seeding CSV outputs.

Validated on 2026-10-22 with the full scalar, AVX-512, and
address/undefined-behavior sanitizer suites. The new batch checks cover 336
configurations, and the CLI suite performs 107 executable checks per build.
Existing window checks still match exact outputs for 1,000 configurations
and 18,756 windows; the AVX-512 frozen-reference traversal check also passes
118 configurations and 2,242 windows.

```sh
make test
make BUILD_DIR=build/interface-simd \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test check-traversal-speed
ASAN_OPTIONS=detect_leaks=0 make BUILD_DIR=build/interface-san \
  CFLAGS='-O1 -g -std=c99 -Wall -Wextra -Wpedantic -fsanitize=address,undefined -fno-omit-frame-pointer' \
  LDFLAGS='-fsanitize=address,undefined' test
```
