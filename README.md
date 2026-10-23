# SubseqHash3

SubseqHash3 generates multiple pivot-based DNA seeds with exact edit budgets.
Both seeding commands use outward traversal by default; `--inward` selects
an alternative hash order with different scores and seeds.

## Build and check

```sh
make
make test
```

The normal build uses scalar kernels. For an AVX-512F/BW/VL CPU:

```sh
make simd
make BUILD_DIR=build/simd \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' test
```

The SIMD binary requires compatible hardware; there is no runtime CPU dispatch.
All seeding is single-threaded.

## Commands

```text
subseqhash3 tables W D S OUTPUT [--seed SEED]
subseqhash3 seeds-sub INPUT W E D R S OUTPUT [TABLE_FILE] [--inward] [--seed SEED]
subseqhash3 seeds-ind INPUT W E D R S OUTPUT [TABLE_FILE] [--inward] [--seed SEED]
subseqhash3 collisions sub|ind PAIRS_FILE E D R S [TABLE_FILE] [--inward] [--seed SEED] [--jaccard]
subseqhash3 simulate ...
subseqhash3 test
```

`W` is the window length and `D` is the number of hash residues. `S` is the
number of independent table sets. **`R` is the number of leading pivots used
from each set**, so each window produces **R*S seeds**. Set selection and
pivot selection are separate; increasing R does not create independent sets.
Different orders may select the same seed string.

| Command | Exact edits per window | Window limits | Pivots per set |
|:---|:---|:---|:---|
| `seeds-sub` | One deletion and E substitutions | 2<=W<=257; 0<=E<W | 1<=R<=W-1 |
| `seeds-ind` | E deletions, no substitutions | 1<=W<=256; 0<=E<W | 1<=R<=W-E |

Both commands select pivots 0 through R-1 from every set. Input contains one
ACGT sequence per line, in either case; it is not FASTA/FASTQ. Every length-W
window is emitted, advancing one base at a time. Invalid DNA, empty lines and
sequences shorter than W are rejected.

For W=30, R=5, S=3, every window produces 15 seeds: pivots 0–4 from set 0,
then pivots 0–4 from set 1, then pivots 0–4 from set 2.

## Shared tables and seeding

Generate three independent table sets and reuse the same file for both modes:

```sh
build/subseqhash3 tables 30 1024 3 shared.piv --seed 42
build/subseqhash3 seeds-sub input.seq 30 4 1024 5 3 sub.csv shared.piv
build/subseqhash3 seeds-ind input.seq 30 4 1024 5 3 ind.csv shared.piv
build/subseqhash3 seeds-sub input.seq 30 4 1024 5 3 inward.csv shared.piv --inward
```

`tables` uses W consistently with the seeding commands. Each set reserves
min(W,256) table rows, including capacity for `seeds-ind` with E=0. The
on-disk row-capacity header remains compatible with existing pivot files;
older files can be reused when they contain enough rows, matching D, and at
least S sets. If a file has extra sets, the first S are used.

Omit TABLE_FILE to create `OUTPUT.piv`; add `--seed` for reproducible creation:

```sh
build/subseqhash3 seeds-sub input.seq 30 4 1024 5 3 sub.csv --seed 42
```

`--seed` applies only to table creation and cannot accompany TABLE_FILE.
The former `seeds`, `repeats`, and `deletions` commands, `tables --pivots`,
and the `--subseqhash2-tables` option are removed. Retained low-level C
functions support legacy algorithms and text-table import.

The current single-thread speed changes preserve exact seeds and scores.
See [the measurements against 36de034](bench/EXACT-SEED-SPEED.md) for runtime,
validation and the limits of the improvement.

## Seeding CSV

Both commands write the same columns:

```text
set,pivot,sequence,window_start,substitutions,deletions,seed,psi,omega
```

The first four columns are zero-based. `set` identifies the independent
tables, `pivot` is the distinguished position in the output seed, `sequence`
identifies the input line, and `window_start` is the window offset in that
sequence. The pivot is not the deleted position. Rows are ordered by sequence,
window, set, then pivot. A length-L input sequence produces `(L-W+1)*R*S` rows.
`psi` is minimized, then `omega` is maximized, with deterministic witnesses.

The redundant `repeat` and `family` columns are replaced by explicit `set`
and `pivot` columns. See [the interface guide](bench/SEEDING.md) for migration,
C APIs, exact table requirements and examples.

## Collision evaluation

```sh
build/subseqhash3 collisions sub pairs.csv 4 1024 5 3 shared.piv
build/subseqhash3 collisions ind pairs.csv 4 1024 5 3 shared.piv --inward
```

Each input line contains `STRING1,STRING2`, without a header. Each string is
seeded as one whole window using the same mode, table sets and pivot indices
as the seeding commands. R must be valid for both strings in every pair.
Strings may have unequal lengths; unequal output lengths cannot collide.

The report includes the probability of any seed-string intersection across
all R*S seeds, any/all matching aligned `(set,pivot)` seeds, and each aligned
set/pivot's collision probability. Independent-table trials are not pooled
with pivot orders. `--jaccard` adds exact edit-neighborhood Jaccard statistics;
this can be expensive. Without TABLE_FILE, tables are saved to `PAIRS_FILE.piv`.

## C interface

Include [src/ssh3_seeding.h](src/ssh3_seeding.h) for
`subseqhash3_seeds_sub()` and `subseqhash3_seeds_ind()`. Both accept shared
pivot tables, R, S, and a traversal, and emit one borrowed batch of R*S
results per window. Result index `set*R+pivot` identifies the table set and
pivot. Passing length==W handles a single window. Result positions are
relative to the window; add window_start for absolute positions.
`SUBSEQHASH3_OUTWARD` is zero and `SUBSEQHASH3_INWARD` selects inward traversal.
Both modes share outward anchor DPs across windows separately for every set.

Use `subseqhash3_seed_tables_create()` in
[src/pivot_table_file.h](src/pivot_table_file.h) to generate W-sized table sets.
The older row-capacity, all-splits, single-seed and repeat-batching C APIs
remain available with their existing semantics. They are not CLI commands.

## Validation and benchmarks

`make test` runs scalar/SIMD kernel oracles, table round trips and malformed
file checks, exact window equivalence, R*S callback equivalence to independent
per-set APIs, and executable CLI tests. CLI checks cover both modes and
traversals, table reuse, output ordering, limits, aliases and collision
probabilities against emitted seeds.

The existing benchmark reports retain their recorded measurements and
historical command syntax where noted. Current examples are in
[SEEDING.md](bench/SEEDING.md).

* [Outward single-window and long-sequence speed measurements](bench/OUTWARD_SPEED.md).
* [Sliding-window algorithms and original measurements](bench/WINDOWS.md).
* [Traversal optimization](bench/TRAVERSAL_OPTIMIZATION.md).
* [Deletion algorithms and exact-score implementation](bench/DELETIONS.md).
* [Pivot construction](bench/PIVOTS.md).

The pivot construction is inspired by Xiang Li, Ke Chen, and Mingfu Shao,
*Efficient seeding for error-prone sequences with SubseqHash2*,
Bioinformatics 41(8), 2025, [doi:10.1093/bioinformatics/btaf418](https://doi.org/10.1093/bioinformatics/btaf418).
