# SubseqHash3 — notes for Claude

Paper: `recomb_2027_subseqhash3.pdf` (+ `-supp.pdf`). No poppler/pip on this
host; extract text with a stdlib zlib Tj/TJ script. `paper/methods.tex` has
the methods source.

## Model
Neighborhood N_{Es,Ed}(x): delete exactly Ed bases of a length-W window,
substitute exactly Es retained bases; seed length k=W-Ed. k ABC^k orders
(pivot i): psi = (psi_R(left) + C_P[i][z_i] + psi_F(right)) mod d minimized,
then omega = sgnL*omega_R + A_P + sgnR*omega_F maximized. Outward (default):
left side read from pivot leftward, right side from pivot rightward. Inward:
left prefix in normal order, right suffix reversed.

## Commands / code
- `seeds-sub` (Es=E, Ed=1): `src/ssh3_pivots.c`, int8 SubseqHash3PivotTables,
  heavily optimized (compact costs, AVX-512).
- `seeds-ind` (Es=0, Ed=E): `src/ssh3_deletions.c`, int64
  SubseqHash3DeletionTables converted via `subseqhash3_deletion_tables_from_pivots`.
- `seeds` (general Es, Ed): `src/ssh3_general.c`, `subseqhash3_seeds()` in
  `src/ssh3_seeding.h`. Results are SubseqHash3EditResult (sub's layout),
  CLI CSV writer `write_sub`. All three share one order: identical psi/omega
  for the special cases (tied witnesses may differ).
- CLI: `src/ssh3_seed_cli.c` `seeds_main(mode)` (0 sub, 1 ind, 2 general).
- Table file (`tables W D S`) stores pivot tables with k=n-1 rows; general
  needs rows >= W-Ed.

## General solver (ssh3_general.c)
- Side DP `GSide`: block per (len,h,e) = Lc row then Uc row, uint16 compact
  costs (M=127: Lc=M*l+min omega, Uc=M*l-max omega; sign -1 swaps Lc/Uc),
  saturating adds, 0xFFFF sentinel. Row stride holds d residues + 32-entry
  circular tail (row[d+i]=row[i mod d]) + cached row min at d+32 (valid flag
  d+33). len = side input bases consumed, h deletions, e substitutions.
  Retention tries all 4 output bases (delta = base != input); e limited to
  [gemin(l), gemax(l)]. Left side capped at retained length R-1.
- Tables `GTab`: int8 A, sign bitmask `neg`; right side stored at reflected
  residue v=-j mod d so the join is a circular vector sum (v2 = j1+C_P-psi).
- Join: pivot output p, input pivot b=p+h (h=0..Ed), all pivot bases sigma,
  splits eL+eR+delta=Es; row-min lower bounds prune. AVX-512 path (G_SIMD,
  needs F/BW/VL) plus scalar fallback; `-DGPROF_PHASES` prints build/join/
  trace cycle counts.
- Traceback recomputes predecessors from exact costs (deletion preferred).
- Outward: anchor-oriented; for anchor a, windows a-(R-1+Ed)..a; circular
  buffer of R+Ed window slots; window w emitted after anchor w+R-1+Ed.
  Inward: per-window prefix/suffix DPs.
- Optimized 2026-11-11: 9-40x over the int64 scalar version (AVX-512 build,
  identical psi/omega). Remaining ideas: prune infeasible (h,e) states further,
  avoid redundant copies for deletion transitions.

## Tests
`make test` (tests/ and bench/ are gitignored, local only). `tests/general.c`
(`make check-general`) brute-forces the neighborhood for n<=8, all (Es,Ed),
both traversals. `tests/cli.py` checks `seeds` vs `seeds-sub`/`seeds-ind` scores.
Default build is scalar -O2; `make simd` for AVX-512.
