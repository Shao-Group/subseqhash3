#!/usr/bin/env python3
"""Compare exact CLI output and single-thread runtime against a saved baseline."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import tempfile
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, default=Path('build/runtime-before/subseqhash3'))
    parser.add_argument('--after', type=Path, default=Path('build/runtime-after/subseqhash3'))
    parser.add_argument('--cpu', type=int, default=2)
    parser.add_argument('--samples', type=int, default=5)
    parser.add_argument('--short-lines', type=int, default=300)
    parser.add_argument('--long-length', type=int, default=512)
    parser.add_argument('--output', type=Path, default=Path('bench/seed-runtime'))
    args = parser.parse_args()
    assert args.samples > 0 and args.short_lines > 0 and args.long_length >= 100
    os.sched_setaffinity(0, {args.cpu})  # Children inherit a single CPU.
    executables = {'before': args.before.resolve(), 'after': args.after.resolve()}
    samples, summary = [], []
    fields = ['mode', 'W', 'E', 'D', 'R', 'S', 'length', 'lines']

    def save():
        for suffix, rows in [('samples', samples), ('summary', summary)]:
            if rows:
                with Path(f'{args.output}-{suffix}.csv').open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
                    writer.writeheader()
                    writer.writerows(rows)

    with tempfile.TemporaryDirectory(prefix='ssh3-runtime-') as temporary:
        root = Path(temporary)
        for n, e in [(30, 4), (50, 8)]:
            for d in [32, 1024]:
                table = root / 'shared.piv'
                subprocess.run([str(executables['before']), 'tables', str(n), str(d), '1',
                                str(table), '--seed', '42'], check=True, capture_output=True)
                for length, lines in [(n, args.short_lines), (args.long_length, 1)]:
                    rng = random.Random(125 + n + length)
                    dna = root / 'dna.seq'
                    dna.write_text(''.join(''.join(rng.choices('ACGT', k=length)) + '\n'
                                           for _ in range(lines)))
                    for mode in ['sub', 'ind']:
                        for r in [10, n - (1 if mode == 'sub' else e)]:
                            config = dict(zip(fields, [mode, n, e, d, r, 1, length, lines]))
                            times = {'before': [], 'after': []}
                            checksum = None
                            # One untimed warmup per executable, followed by alternating pairs.
                            for sample in range(-1, args.samples):
                                order = ['before', 'after'] if sample % 2 == 0 else ['after', 'before']
                                for version in order:
                                    output = root / f'{version}.csv'
                                    command = [str(executables[version]), f'seeds-{mode}', str(dna),
                                               str(n), str(e), str(d), str(r), '1', str(output), str(table)]
                                    begin = time.perf_counter()
                                    subprocess.run(command, check=True, capture_output=True)
                                    elapsed = time.perf_counter() - begin
                                    current = digest(output)
                                    if checksum is None:
                                        checksum = current
                                    assert current == checksum, (config, version, sample)
                                    if sample >= 0:
                                        times[version].append(elapsed)
                                        samples.append(config | {'version': version, 'sample': sample,
                                                                 'seconds': elapsed, 'sha256': current})
                            a, b = [statistics.median(times[v]) for v in ['before', 'after']]
                            summary.append(config | {'before_seconds': a, 'after_seconds': b,
                                                     'speedup': a / b, 'sha256': checksum})
                            save()
                            print(f'{config}: {a / b:.2f}x ({a:.4f}s -> {b:.4f}s)', flush=True)
    model = next(line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                 if line.startswith('model name'))
    metadata = {
        'cpu': args.cpu, 'cpu_model': model, 'platform': platform.platform(),
        'samples': args.samples, 'executables': {k: str(v) for k, v in executables.items()},
        'executable_sha256': {k: digest(v) for k, v in executables.items()},
        'protocol': 'One CPU; no worker threads; one untimed warmup, alternating before/after pairs; '
                    'median whole-process wall time including table loading, conversion, seeding and CSV I/O; '
                    'shared pre-generated tables; SHA256 equality of complete CSV on every run.',
        'sources_sha256': {p: digest(Path(p)) for p in ['src/ssh3_pivots.c', 'src/ssh3_deletions.c',
                                                       'src/ssh3_seed_cli.c', __file__]},
    }
    Path(f'{args.output}-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
