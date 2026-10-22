# Frozen single-threaded benchmark reference

`pivots_round1.c` is the exact restored `src/ssh3_pivots.c` at the start of
the strictly single-threaded second-round request. It includes the earlier
optimization, not the older implementation at commit ca35567. It is only
linked into `bench-pivots-single`, with all exported symbols renamed, and is
not part of the library or CLI. Do not update it when optimizing the current
implementation; that would move the benchmark baseline.

SHA-256: `cc76e256e7d2f98ca432d9e1a8d3200a49d49abb30539e09d766660b920a3f88`

Both implementations are compiled with the same flags, run on one CPU, and
use identical input sequences, initialized tables, requested orders, and edit
budgets. The API header is unchanged by the second-round optimization.

`deletions_round1.c` is the deletion kernel from `fab085f`, before compact
costs, sparse base cases, and per-anchor deletion-budget pruning. It is
compiled under renamed symbols with exactly the same flags as production
by `bench-deletions-speed`.

`traversal_deletions.c` and `traversal_pivots.c` freeze all four traversal
implementations immediately before the October 21 traversal optimization.
`traversal_rename.h` renames their exported symbols so `bench-traversal-speed`
can compare them against production in the same process with identical flags.
The benchmark checks exact scores, seed strings, positions and substitution
flags, including deterministic tie witnesses. Keep these snapshots fixed.
