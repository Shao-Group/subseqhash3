# One deletion, E substitutions, multiple pivot orders

This is a new algorithm and API in `src/ssh3_pivots.h` and
`src/ssh3_pivots.c`. It generates one minimizing length-`k=n-1` string
per requested pivot. Every result uses **exactly one deletion and exactly
E substitutions of retained bases**. A substitution must change its base.
E can be zero through k. Different pivots can produce the same seed.

The design follows the shared two-sided scores and pivot tables in Section
2.1 of `subseqhash2-paper.pdf` and the `combine` routine in
`subseqhash2/src/subseqhash2seeding.cpp`. The linked SubseqHash2 repository
and existing SubseqHash3 algorithms are not modified.

## Defined orders and difference from the paper

The paper hashes the reversed left side and the normal right side of the
pivot. Here the directions are exchanged: hash the left side normally,
starting at the left boundary; hash the right side in reverse, starting at
the right boundary. This deliberate variant permits all pivots in one
window to reuse the same prefix DP and suffix DP.

For an output y of length k and zero-based pivot p, write
`L=y[0..p-1]`, `a=y[p]`, and `R=y[k-1]..y[p+1]`. The tables define

```text
psi_p(y)   = (psi_left(L) + C_pivot[p,a] + psi_right(R)) mod D
omega_p(y) = sign_left[p,a] * omega_left(L)
             + A_pivot[p,a]
             + sign_right[p,a] * omega_right(R)
```

The side psi/omega functions have exactly the ABC recurrences in the
paper. A side's table row is its local retained position, starting at zero.
Each order minimizes psi and then maximizes omega. This differs from the
older SubseqHash3 API's largest-psi, largest-absolute-omega convention.
Scores can tie with finite integer tables. As in the source code's
score comparisons, any tied minimizing string is valid; a deterministic
traceback witness is returned, with no claim of a lexicographic tie break.
Scalar and SIMD implementations return the same witnesses.

Tables use signed pivot magnitudes in 1..127 and side magnitudes in
1..127, permitting exact 16-bit DP arithmetic for k<=256. This is not a
reader for the original program's table files. The k orders in a family
are related, not independent random permutations. Correctness tests verify
their minimizers; collision-probability equivalence to independent hashes
is not established by this implementation or its timing benchmark.

## Shared computation

For every prefix, the forward DP stores both score extrema for each
deletion count h in {0,1}, substitution count e in {0,...,E}, and residue.
The reverse DP applies the same recurrence to the reversed input. Keeping
a base a comes from the previous input prefix with e reduced by whether
a differs from the input base; deleting it comes from deletion count h-1.
Only states with `e<=retained_length` and
`e>=E-(k-retained_length)` are evaluated: the remaining retained positions
must be able to supply all remaining substitutions. Unneeded last
deletion-free layers are skipped as well.

For output pivot p, the input pivot can only be p or p+1:

```text
input pivot i = p + h_left
h_left + h_right = 1
e_left + (pivot_base != input[i]) + e_right = E
```

The optimizer considers both input pivots, all four possible pivot bases,
and all feasible substitution allocations. Each prefix and suffix table is
computed once, with computation trimmed to the requested pivot range.
Traceback reconstructs all seed bases, their original input positions,
and their substitution flags.

Both side extrema are represented as minimum nonnegative costs
`127*retained_length + omega` and `127*retained_length - omega`. All
transition costs are nonnegative, so AVX-512 unsigned saturating addition
preserves the unreachable sentinel without a separate reachability mask.

Right-side residues are stored at `-psi mod D`. For a target output residue
t, a left residue q pairs with the right array at
`(q + C_pivot - t) mod D`. This makes each join a circular, contiguous SIMD
addition followed by a vector minimum reduction. The optimizer first
tries residue zero. If it is unreachable, it finds the smallest feasible
residue in O(D) with a previous-reachable-residue array, then evaluates that
join. It never assumes residue zero is reachable and never enumerates D^2
residue pairs. Masked tails and circular padding support arbitrary positive
D, including 31 and 1025. A scalar fallback implements the same recurrence.

For r<=k requested pivots, the work is O(n*(E+1)*D) for the shared DPs and
O(r*(E+1)*D) for joins, plus O(r*k) traceback. SIMD evaluates 32 residues
at once. Space is O(n*(E+1)*D), plus output storage. More than k repetitions
use ceil(r/k) separately initialized table families, with a partial final
batch. No speedup is guaranteed for every r: at r=1 there is no computation
to share, and at large D the new algorithm can be slower than the earlier
single-seed function.

## C API

For all k pivots in one family (input symbols 0..3):

```c
#include "ssh3_pivots.h"
#include <stdlib.h>

SubseqHash3PivotTables tables = {0};
SubseqHash3EditResult *results = calloc(n - 1, sizeof(*results));
if (results && !subseqhash3_pivot_tables_init(&tables, n, D, random_seed)) {
    if (!subseqhash3_compute_pivot_seeds(sequence, n, E, &tables,
                                      NULL, n - 1, results)) {
        /* results[p] is the minimizer for pivot p. */
    }
}
subseqhash3_edit_results_destroy(results, n - 1);
free(results);
subseqhash3_pivot_tables_destroy(&tables);
```

The caller owns all allocations. Functions return 0 on success and -1 on
invalid input or allocation failure. Pass empty results and destroy their
contents before reusing them. Tables support input lengths 2..257, with
E<=n-1. Tables created for a longer input also support shorter inputs.

To select r<=k orders, pass `NULL,r` for the first r pivots, or pass an
explicit array of r zero-based pivot positions. Requested pivots need not
be contiguous or sorted.

For any r>0, initialize at least `1+(r-1)/(n-1)` families using different
random seeds and call

```c
subseqhash3_compute_pivot_repeats(sequence, n, E, families,
                                family_count, r, results);
```

Results are family-major and then pivot-major. For seeding many windows,
create a workspace once with `subseqhash3_pivot_workspace_create()` and use
`subseqhash3_compute_pivot_repeats_with_workspace(...,workspace)`.
Destroy it with `subseqhash3_pivot_workspace_destroy()` when finished.
The ordinary repeat API already reuses allocations between table families;
the explicit workspace additionally reuses them across input calls. Each
thread must have its own workspace. Tables are read-only during seeding.

## Correctness and benchmark

For the requested r=10 comparison at `(n,E)=(30,5)` and `(50,8)`, see
[the focused benchmark report](R10.md), including raw samples and CPU-pinned
reproduction commands for D=31 and D=1024.

```sh
make test check-pivots
make BUILD_DIR=build/pivots-fast \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' \
  test bench-pivots
```

Example measurements on a Xeon Gold 6148 (2026-10-20), including DP
allocation per call and comparing with r calls to the previous API:

| n | E substitutions | D | r seeds | Shared batch (us) | Previous API x r (us) | Speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 5 | 31 | 1 | 10.889 | 30.535 | 2.80x |
| 30 | 5 | 31 | 8 | 18.998 | 250.025 | 13.16x |
| 30 | 5 | 31 | 29 | 48.649 | 887.658 | 18.25x |
| 30 | 5 | 31 | 61 | 114.931 | 1911.835 | 16.63x |
| 30 | 5 | 1024 | 29 | 410.833 | 2583.269 | 6.29x |
| 50 | 8 | 1024 | 1 | 357.867 | 241.084 | 0.67x |
| 50 | 8 | 1024 | 8 | 404.529 | 1881.613 | 4.65x |
| 50 | 8 | 1024 | 49 | 1073.566 | 11603.661 | 10.81x |
| 50 | 8 | 1024 | 101 | 2475.552 | 24109.869 | 9.74x |

Full measurements, including the same-orders comparison and reusable
workspace API, are in `bench/pivots-results.csv`. For example, the
`n=50,E=8,D=1024,r=49` batch is 15.58x faster than separate single-pivot
calls computing exactly the same neighborhood and orders. The previous
API column uses swapped substitution/deletion counts, as explained below.

The dedicated check independently enumerates all deleted positions and
all strings with exactly E changes for 140 small neighborhoods, evaluates
every pivot's score literally, and checks that each reported score is
optimal. It also validates 30 larger witnesses, single-pivot versus shared
calls, workspace reuse with changing sizes and D, invalid arguments, and
r>k batching.

The default benchmark uses eight random input strings, two loops per timing
sample, and the median of three samples. It rotates method timing order.
Current runs give every method an untimed warmup; the historical snapshot
above predates that uniform warmup. The `--r10` mode uses 32 strings, four
loops, and seven samples, and emits every sample alongside the medians.
Table generation is outside the timed region; result allocation,
traceback, and result cleanup are inside. `shared_us` includes DP buffer
allocation and cleanup per input call. `reuse_us` amortizes DP allocations
through the public workspace API. The reference columns are:

- `separate_pivots_us`: the same new orders and edit neighborhood, obtained
  by calling the single-pivot API r times with matching families/pivots.
- `previous_api_us`: the previously optimized
  `subseqhash3_compute_one_substitution()` called r times with independent
  table sets. That function uses one substitution and E deletions, so this
  column is the requested **timing comparison only**, not a correctness or
  seed-quality comparison for the new neighborhood.

The standalone check can run under ASan/UBSan without invoking the
existing table-file serializer:

```sh
make BUILD_DIR=build/pivots-sanitize \
  CFLAGS='-O1 -g -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -fsanitize=address,undefined -fno-omit-frame-pointer' \
  LDFLAGS='-fsanitize=address,undefined' build/pivots-sanitize/bench-pivots
ASAN_OPTIONS=detect_leaks=0 build/pivots-sanitize/bench-pivots --check-only
```

Leak detection is disabled because LeakSanitizer cannot run under this
environment's tracing. Address and undefined-behavior checks remain active.
