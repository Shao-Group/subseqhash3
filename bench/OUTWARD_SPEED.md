# Outward speed improvements (2026-10-22)

This compares outward seeding before and after the optimization, at identical
input sequences, tables, edit budgets, residue counts and requested pivots.
The baseline is commit `e714ebfab4f47e7e541c78a6e31027de4f0431e5`, which
already includes sharing anchor DPs across overlapping windows.

## Changes

* **Deletions:** pack and reflect the immutable left/right tables once per
  long-sequence call instead of once per new anchor. Preparation stays local
  to the invocation, so public workspace reuse does not cache stale tables.
* **Deletions:** vectorize right boundary reflection and recover the winning
  right residue with a SIMD search that preserves the original witness order.
* **Repeats:** join scores directly from the two anchor DPs. Remove boundary
  copying and the extra pair of DP allocations used only for those copies.
  Reserve the full requested pivot range once, avoiding repeated allocation
  and page faults as the first call advances through anchors.
* **Repeats:** traceback reads the original sequence directly, avoiding a
  transformed-sequence copy for every seed. Parent-based traceback has a
  separate loop that avoids score loads and candidate searches.
* **Both:** exploit empty sides at endpoint pivots. Traceback can copy a fixed
  remaining witness once its edit budget is exhausted.

The hash families, scores, seed bytes, positions, substitution flags and tie
rules are preserved. CLI flags and public C interfaces are unchanged.

## Measurements

See [median results](outward-speed-results.csv),
[individual samples](outward-speed-samples.csv) and
[machine, build flags and source hashes](outward-speed-metadata.json).

Measurements use one thread pinned to CPU 0 on an Intel Xeon Gold 6148 with
the same `-O3`/AVX-512 flags for both binaries. Each entry is the median of
three complete calls. Method order rotates between samples; executable order
rotates between conditions. Timings include workspace allocation, startup,
traceback, output allocation, checksum callbacks and cleanup, and exclude
table generation and disk I/O. Units are microseconds per emitted window.

The sweep covers W=30/E=4 and W=50/E=8, D=32/1024, and ten/all pivots.
All means W-E seeds for deletions and W-1 for repeats. Three scenarios are
reported separately:

* `single-cold`: input length W, with a fresh workspace on every timed call.
  The allocator and CPU caches can remain warm between samples.
* `single-reused`: length-256 input processed with independent single-window
  calls using one reusable workspace; no DP sharing across windows.
* `windows-shared`: length-256/4096 inputs processed with the window API and
  anchor reuse, including cache startup and cleanup.

### Observed speedups

Ranges cover W=30/50 and ten/all pivots. A factor above 1 is faster.
The D=32 long-sequence column uses a separate seven-sample confirmation
sweep; the other columns use the primary three-sample sweep.

| Algorithm | Single, fresh workspace, D=1024 | Single, reused workspace, D=1024 | Length 4096, D=32 | Length 4096, D=1024 |
|:---|---:|---:|---:|---:|
| deletions | 0.99–1.12x | 1.11–1.18x | 1.04–1.14x | 1.19–1.94x |
| repeats | 1.51–2.67x | 1.06–1.07x | 1.08–1.16x | 0.99–1.02x |

**Long-sequence repeats at D=1024 remains near parity.** The measurements
show a consistent improvement for long-sequence repeats at D=32, but do
not establish a consistent improvement at D=1024. Fresh-workspace deletion
calls are also close to parity in some cases. Individual conditions can be
slower; the CSVs include every measured result rather than only wins.

The seven-sample D=32 confirmation has its own
[results](outward-d32-confirm-results.csv),
[samples](outward-d32-confirm-samples.csv) and
[metadata](outward-d32-confirm-metadata.json). It repeats the complete
D=32 sweep, including fresh/reused single-window calls and length-256
inputs. To reproduce this confirmation, append
`--residues 32 --samples 7 --output-prefix bench/outward-d32-confirm`
to the runner command below.

### W=50, E=8, D=1024, all pivots

Times are microseconds per window. All-pivot counts are 42 for deletions
and 49 for repeats.

| Algorithm | Scenario | Before us | After us | Speedup |
|:---|:---|---:|---:|---:|
| deletions | Single, fresh workspace | 3703.86 | 3298.90 | 1.12x |
| deletions | Single, reused workspace | 3609.05 | 3209.21 | 1.12x |
| deletions | Length 4096 | 935.68 | 784.44 | 1.19x |
| repeats | Single, fresh workspace | 18986.97 | 11398.55 | 1.67x |
| repeats | Single, reused workspace | 12134.62 | 11402.04 | 1.06x |
| repeats | Length 4096 | 2162.50 | 2182.41 | 0.99x |

## Correctness

The scalar and AVX-512 full test suites pass, along with scalar and AVX-512
address/undefined-behavior sanitizer suites (`ASAN_OPTIONS=detect_leaks=0`).
The window suite checks exact outputs for 1,000 configurations and 18,756
windows, including ring wraparound, partial vectors, multiple table families,
maximal edits, wide imported tables, callback cancellation and invalid inputs.
The frozen-reference traversal check compares scores, seeds, positions and
substitution flags on another 118 configurations and 2,242 windows per build.
Every benchmark condition also requires equal complete-sequence checksums
between the before/after binaries, for both outward and the inward control.

## Reproduce

Build the baseline in a separate directory so production sources stay intact:

```sh
mkdir -p build/outward-baseline-src
git archive e714ebfab4f47e7e541c78a6e31027de4f0431e5 | tar -x -C build/outward-baseline-src
flags='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512'
make -C build/outward-baseline-src -j4 BUILD_DIR=build CFLAGS="$flags" build/bench-windows
make -j4 BUILD_DIR=build/outward-opt CFLAGS="$flags" test check-traversal-speed build/outward-opt/bench-windows
python3 bench/run_outward_speed.py \
  --before build/outward-baseline-src/build/bench-windows \
  --after build/outward-opt/bench-windows \
  --baseline-ref e714ebfab4f47e7e541c78a6e31027de4f0431e5 \
  --flags "$flags"
```

These measurements describe synthetic uniform DNA on this machine; small
percentage differences are sensitive to CPU frequency and system activity.
