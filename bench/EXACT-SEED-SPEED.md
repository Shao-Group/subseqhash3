# Exact-output improvements after 36de034

This change starts from **36de034** and retains the original full-D hash
family. Seeds, scores, residues, tie order, positions, substitution flags,
seed lengths and collision reports are unchanged. Work is single-threaded.
The rejected experimental schemes have been removed from the active source,
build and CLI. No power-of-two restriction on D is needed.

**The 5x target is not reached.** For the specified workload, the measured
improvement over 36de034 is 1.53x/1.24x for sub/ind at D=64 and
1.03x/1.07x at D=1024. These results preserve sensitivity exactly rather
than trading it for runtime.

## Changes

* Use the anchor-major outward engine at D=64, in addition to D>=128.
  Keeping an anchor's DP tables hot avoids the larger rotating score-table
  cache. Other small and intermediate D settings retain their previous path.
* Retain candidate seed/position/flag buffers across ring reuse instead of
  allocating and freeing them after each callback. Length zero marks an
  inactive candidate; all allocations are freed on success, error or cancellation.
* Round candidate-ring storage to a power of two and use a mask for indexing.
  The logical completion span remains unchanged. This optimizes storage
  addressing and does not change the hash modulus D.
* Separate parent recording and traversal direction into four DP kernels.
  The recurrence stays in one source helper, with these choices resolved
  before the inner loops. Disable the large edit-loop unroll: the generated
  x86 DP kernel text shrinks from about 61 KiB to 17 KiB in this build.
* Reduce 16-bit lane minima using two half-vector minima and the SSE
  minimum-position instruction, avoiding the longer 32-bit reduction sequence.
* Copy fixed witness intervals once all edits are consumed, while generating
  their consecutive positions in a simple loop.
* Format integral CSV fields directly, encode DNA with byte shuffles in the
  SIMD build, and batch rows into an 8 KiB buffer. Scalar formatting is also
  supported. Masked seed-tail loads and reserved row space keep accesses bounded.

The most useful improvement here is D=64 cache locality plus CSV formatting.
The large-D score recurrence still dominates runtime; smaller instruction
changes do not establish a substantial large-D speedup by themselves.

## End-to-end timings

Intel Xeon Gold 6148, GCC 13.3, CPU 2, one thread. Both builds use identical flags:

```
-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512
```

W=30, E=5, R=10, S=1, outward. One deterministic random 1000-base sequence,
971 windows and 9,710 CSV rows. Fifteen alternating before/after pairs after
warmup. Median whole-process wall time includes startup, table loading,
conversion, DP construction, traceback and CSV output; table generation is
excluded. Every output SHA256 matches 36de034 on every run.

|D|Command|36de034 ms|Updated ms|Speedup|
|---|---|---|---|---|
|64|seeds-sub|27.479|17.950|1.53x|
|64|seeds-ind|17.208|13.848|1.24x|
|1024|seeds-sub|159.871|154.671|1.03x|
|1024|seeds-ind|54.794|51.359|1.07x|

A separate batch of 16 deterministic random sequences, **each 1000 bases**,
uses five alternating pairs and the same complete-process protocol:

|D|Command|36de034 ms|Updated ms|Speedup|
|---|---|---|---|---|
|64|seeds-sub|399.768|257.069|1.56x|
|64|seeds-ind|241.923|192.355|1.26x|
|1024|seeds-sub|2488.505|2375.423|1.05x|
|1024|seeds-ind|840.265|785.329|1.07x|

These are measurements for this processor, compiler and input distribution;
the smallest changes are close enough to timing variation that they should
not be treated as universal guarantees. No statistical confidence interval
or genome-wide sensitivity estimate is inferred from the timing samples.

A separate D=64 sub run on 32 local 1000-base sequences measured whole-process
peak RSS at 5,048 -> 1,920 KiB with `/usr/bin/time -f '%M'`; output files were
byte-identical. RSS varies with allocator and environment. Rounded rings add
some candidate storage, while the anchor path removes the larger score cache.

## Correctness and sensitivity

The DP transitions, hash tables, full residue range and strict candidate
ordering are unchanged. These optimizations preserve the selected seeds,
so they preserve which string pairs collide, rather than relying on similar
empirical collision rates from a different scheme.

Validation includes:

* Full scalar and AVX-512 suites, including exhaustive tiny instances, both
  traversals, exact ties, sparse residues and imported wide-score fallback.
* **1,064 window configurations and 23,140 windows** checked against per-window
  APIs for scores, seeds, positions and substitution flags. Added W=30/E=5
  tests cover D=64 random, homopolymer and periodic inputs, plus D=1024.
* **336 R*S configurations**, cancellation and invalid inputs.
* **153 executable checks**, including large CSV batches, full SIMD seed blocks,
  masked tails, maximum window sizes and rejection of the removed scheme options.
* ASan/UBSan window and CLI checks. Leak detection is disabled in this environment.
* **88 complete CSV/collision-report comparisons against the preserved 36de034
  executable**, across three additional table seeds, random/homopolymer/periodic
  1000-base inputs, multiple sets, both traversals and non-power-of-two D.
  Every comparison is byte-identical. Timing samples independently check hashes.

The differential comparisons use the old executable as the reference, avoiding
reliance solely on per-window and sliding APIs that share the optimized helper.

## Reproduce

```sh
make BUILD_DIR=build/simd \
  CFLAGS='-O3 -std=c99 -Wall -Wextra -Wpedantic -mavx512f -mavx512bw -mavx512vl -mtune=skylake-avx512' test
python3 bench/run_exact_seed_runtime.py --before build/typical/previous --samples 15
python3 bench/run_exact_seed_runtime.py --before build/typical/previous \
  --lines 16 --samples 5 --output bench/seed-runtime-exact-batch
python3 bench/check_exact_seed_outputs.py --before build/typical/previous
```

The runners expect a separately preserved or rebuilt 36de034 executable.
`--before` and `--after` override executable paths. Build 36de034 from an
isolated checkout/archive with the same flags; do not rebuild the reference
from the modified source. Commands and output formats remain unchanged.

* [Single-sequence medians](seed-runtime-exact-summary.csv), [samples](seed-runtime-exact-samples.csv), [metadata and hashes](seed-runtime-exact-metadata.json)
* [Batch medians](seed-runtime-exact-batch-summary.csv), [samples](seed-runtime-exact-batch-samples.csv), [metadata](seed-runtime-exact-batch-metadata.json)
* [Differential check results and hashes](seed-runtime-exact-checks.json)
* [Timing runner](run_exact_seed_runtime.py), [differential runner](check_exact_seed_outputs.py)
* [Earlier checkpoint measurements](SEED-RUNTIME.md)
