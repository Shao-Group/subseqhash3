# SubseqHash3

SubseqHash3 generates and evaluates DNA seeds with exact edit budgets.
The CLI supports `simulate`, `tables`, `seeds`, `repeats`, `collisions`, and
`test`; run `build/subseqhash3 --help` for usage. `repeats` exposes the pivot
approach; the specialized one-substitution function remains a C API.

## Build and check

```sh
make
make test check-one check-pivots
```

For AVX-512F/BW/VL hardware:

```sh
make BUILD_DIR=build/pivots-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test check-one check-pivots
```

The default build uses the scalar fallback. The AVX-512 build requires a
compatible CPU and does not perform runtime CPU dispatch.

## Pivot seeds from the command line

Generate reusable pivot tables, then seed every length-30 window with exactly
one deletion and five substitutions, producing ten repetitions per window:

```sh
make simd
build/simd/subseqhash3 tables 30 1024 10 tables.piv --pivots --seed 42
build/simd/subseqhash3 repeats input.seq 30 5 1024 10 seeds.csv tables.piv
```

`input.seq` contains one DNA sequence per line. Alternatively, omit the last
argument to generate tables automatically in `seeds.csv.piv`:

```sh
build/simd/subseqhash3 repeats input.seq 30 5 1024 10 seeds.csv
```

The default `tables` command still produces legacy files for `seeds` and
`collisions`. These lack pivot parameters and cannot be used with `repeats`;
use `--pivots` for the new format. In pivot mode, R means seed repetitions,
requiring `ceil(R/(W-1))` independent families rather than R full table sets.
All seeding remains single-threaded. See [CLI and table-format details](bench/REPEATS.md).

## Added APIs

| API | Exact substitutions | Exact deletions | Output |
|---|---:|---:|---|
| `subseqhash3_compute_one_substitution()` | 1 | E | One seed of length n-E |
| `subseqhash3_compute_pivot_seeds()` | E | 1 | Up to k=n-1 pivot orders per family |
| `subseqhash3_compute_pivot_repeats()` | E | 1 | Any positive r, batching independent families |

The single-seed API in `src/subseqhash3.h` uses a reduced DP with two
substitution planes. Its SIMD path uses nonnegative cost scores, retained
length rows, shared table vectors, and demand-driven suffix evaluation.
It matches the corresponding result of the existing all-splits function,
including witness tie breaking. See [single-seed details](bench/README.md).

The pivot APIs in `src/ssh3_pivots.h` share prefix/suffix DPs across orders,
use reflected right-side residues for SIMD joins, and support optional
workspace reuse across windows. Their orders follow the ABC^k construction
but use boundary-to-pivot side directions, minimize psi, and maximize omega.
They are a documented variant of the SubseqHash2 paper, not its original
order or table-file format. See [pivot algorithm and usage](bench/PIVOTS.md).

Both APIs return encoded seeds, source positions, substitution flags, and
scores through `SubseqHash3EditResult`. Release result contents with
`subseqhash3_edit_results_destroy()`. Existing all-splits functions are
unchanged. The required `src/ssh3_collisions.h` declaration is also included
in version control so a fresh checkout can build.

## Benchmarks and validation

The dedicated checks include 365 single-seed equivalence cases and 140
exhaustive small edit neighborhoods for the pivot algorithm, plus larger
witnesses, score limits, workspace reuse, and repeat batching. The optimization
benchmark adds 277 exact comparisons against its linked reference, including
all-tied scores, unreachable residue zero, and maximum-length/cost cases.
Standalone ASan/UBSan commands are recorded in the algorithm documents.
The `repeats` addition also fixes the previously recorded legacy table-file
mask-copy over-read caused by misplaced integer division. Historical optimization reports
describe the earlier state. `make test` now includes pivot table-file and CLI
tests; `make check-repeats` runs these independently.

The requested r=10 comparison is reproducible with:

```sh
make BUILD_DIR=build/pivots-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  build/pivots-fast/bench-pivots
taskset -c 0 build/pivots-fast/bench-pivots --r10
```

It tests `(n,E)=(30,5)` and `(50,8)` at both D=31 and D=1024, with 32
inputs, four loops, seven timing samples, and a warmup for every method.
The machine must permit CPU 0 for the affinity command; substitute an
allowed CPU otherwise. See [r=10 results and methodology](bench/R10.md).

The follow-up pivot optimization preserves those exact orders and witnesses.
For its before/after comparison against the original pivot implementation,
see [optimization results](bench/PIVOT-OPTIMIZATION.md) and the
`bench-pivots-speed` / `check-pivots-speed` Makefile targets. This comparison
links the reference source from commit `ca35567` with renamed exports, so it
requires repository history, not just a source archive.

The subsequent strictly single-threaded round accelerates joins for D<=32.
Its [before/after report](bench/PIVOT-SINGLE-THREAD.md) compares against the
restored first optimization, using `bench-pivots-single` and
`check-pivots-single`. It does not achieve the requested further 2.5x gain;
D=1024 performance remains essentially unchanged.

The previous API has the opposite substitution/deletion counts to the
pivot API. Therefore their runtime comparison is not a seed-quality or
same-edit-neighborhood comparison. The benchmark additionally reports
separate calls to the same pivot orders as a same-problem reference.

The pivot construction is inspired by Xiang Li, Ke Chen, and Mingfu Shao,
*Efficient seeding for error-prone sequences with SubseqHash2*,
Bioinformatics 41(8), 2025, [doi:10.1093/bioinformatics/btaf418](https://doi.org/10.1093/bioinformatics/btaf418).
