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
- Side DP `GSide`: index ((len*(Ed+1)+h)*(Es+1)+e)*d + j, len = side input
  bases consumed, h deletions, e substitutions, j residue; exact int64 min/max,
  GINF sentinels. Retention tries all 4 output bases (delta = base != input).
  Left side capped at retained length R-1.
- Join (`gjoin`): pivot output p, input pivot b=p+h (h=0..Ed), all pivot
  bases sigma and splits eL+eR+delta=Es. psi search scans targets 0..best
  (worst O(d^2)), omega scan O(d).
- Traceback recomputes predecessors from exact scores (deletion preferred).
- Outward: anchor-oriented; for anchor a, windows a-(R-1+Ed)..a; circular
  buffer of R+Ed window slots; window w emitted after anchor w+R-1+Ed.
  Inward: per-window prefix/suffix DPs.
- Scalar, not yet optimized. Speed ideas: 16-bit compact costs + SIMD as in
  supp. Note 3, reflected-residue O(d) joins, prune infeasible (h,e) states,
  avoid redundant memcpy for deletion transitions.

## Tests
`make test` (tests/ and bench/ are gitignored, local only). `tests/general.c`
(`make check-general`) brute-forces the neighborhood for n<=8, all (Es,Ed),
both traversals. `tests/cli.py` checks `seeds` vs `seeds-sub`/`seeds-ind` scores.
Default build is scalar -O2; `make simd` for AVX-512.
