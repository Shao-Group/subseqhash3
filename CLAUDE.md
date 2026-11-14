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
  circular tail (row[d+i]=row[i mod d]). len = side input bases consumed,
  h deletions, e substitutions. Retention tries all 4 output bases (delta =
  base != input); e limited to [gemin(l), gemax(l)]. Left side capped at
  retained length R-1. Start blocks (h,h,0) are set once per GSide (`ready`).
- Tables `GTab`: int8 A, sign bitmask `neg`; right side stored at reflected
  residue v=-j mod d so the join is a circular vector sum (v2 = j1+C_P-psi).
  SIMD-only `prep` (16-bit M+a per level/chunk/row; kernel derives M-a) and
  `negs`, built once per table set.
- Build kernel `gtarget` (SIMD): one target block per sweep over chunks,
  per-level wrapped offsets `qt` (no wrap branches), two blends for the
  swap. Row minima for join-final levels (l>=minl) go to dense `mins` and
  per-residue-mod-32 minima `mm` (lane-wise min over chunks, stored twice).
- Join: pivot output p, input pivot b=p+h (h=0..Ed), all pivot bases sigma,
  splits eL+eR+delta=Es. Once psi=0 is known, candidates are sorted by an
  omega upper bound (mod-32 class minima paired at shift C_P when 32|d, else
  row minima) and pruned; ties keep (sigma,el) order so output is unchanged.
  psi>0 search `gpsi` uses reachability bitsets (nearest partner, clz/ctz).
  AVX-512 path (G_SIMD, needs F/BW/VL) plus scalar fallback; `-DGPROF_PHASES`
  prints lbuild/build/join/trace cycle counts.
- Small d (SIMD): d<=16 (`GPACK`) packs a block into one vector, Lc lanes
  0..15 and Uc 16..31 (stride 16); `gtarget_packed` does shift+sign swap in
  one `permutexvar` per row via `GTab.pidx`. 16<d<32: `gtarget_small`, rows
  d-periodic over 32 lanes, two-source permute per row; tail 32..63 filled
  only for join-final levels. 10-13x over d76adbf for d<32.
- psi>0: reachability bitsets cached per (side,e) within one join call.
- Traceback recomputes predecessors from exact costs (deletion preferred),
  using `GTab.tt` ((M+a)|sign<<8) and incremental block offsets.
- Outward: anchor-oriented; for anchor a, windows a-(R-1+Ed)..a; circular
  buffer of R+Ed window slots; window w emitted after anchor w+R-1+Ed. If the
  Ed+1 anchors' DPs exceed ~768KB, winners are traced eagerly on improvement
  (one DP per side, stays in L2); else a ring of Ed+1 DPs and one trace per
  winner. Inward: per-window prefix/suffix DPs.
- Optimized 2026-11-11 (7b0206d): 9-40x over the int64 scalar version. Then
  1.2-3.5x more (wall, D=1024 typical 1.5-2.1x; identical CSV output). Build
  is ~2/3 of time and near its limit on this Xeon 6148: zmm int ops issue
  ~1/cycle (no gain over ymm), L1 is 32KB so d=1024 rows stream from L2;
  measured floor ~2.2 TSC per 32-lane transform, in-program ~6. Tried and
  rejected (no gain or slower): uint8 increments, level/plane/column-wide
  chunk-outer kernels (GCC spills), paired targets, software prefetch.

## Tests
`make test` (tests/ and bench/ are gitignored, local only). `tests/general.c`
(`make check-general`) brute-forces the neighborhood for n<=8, all (Es,Ed),
both traversals. `tests/cli.py` checks `seeds` vs `seeds-sub`/`seeds-ind` scores.
Default build is scalar -O2; `make simd` for AVX-512.
