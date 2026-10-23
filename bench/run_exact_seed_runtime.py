#!/usr/bin/env python3
"""Exact-output comparison against 36de034 at the typical geometry."""
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


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, default=Path('build/typical/previous'))
    parser.add_argument('--after', type=Path, default=Path('build/simd/subseqhash3'))
    parser.add_argument('--cpu', type=int, default=2)
    parser.add_argument('--lines', type=int, default=1)
    parser.add_argument('--samples', type=int, default=5)
    parser.add_argument('--output', type=Path, default=Path('bench/seed-runtime-exact'))
    args = parser.parse_args()
    assert args.samples > 0 and args.lines > 0
    os.sched_setaffinity(0, {args.cpu})
    executables = {name: path.resolve() for name, path in [('before', args.before), ('after', args.after)]}
    samples, summary = [], []
    with tempfile.TemporaryDirectory(prefix='ssh3-typical-') as directory:
        root = Path(directory)
        rng = random.Random(562)
        dna = root / 'input.seq'
        dna.write_text(''.join(''.join(rng.choices('ACGT', k=1000)) + '\n' for _ in range(args.lines)))
        input_sha = sha256(dna)
        table_hashes = {}
        for d in [64, 1024]:
            table = root / 'table.piv'
            subprocess.run([str(executables['before']), 'tables', '30', str(d), '1',
                            str(table), '--seed', '42'], check=True, capture_output=True)
            table_hashes[str(d)] = sha256(table)
            for mode in ['sub', 'ind']:
                times = {name: [] for name in executables}
                checksum = None
                for sample in range(-1, args.samples):
                    order = ['before', 'after'] if sample % 2 == 0 else ['after', 'before']
                    for version in order:
                        output = root / f'{version}.csv'
                        command = [str(executables[version]), f'seeds-{mode}', str(dna), '30',
                                   '5', str(d), '10', '1', str(output), str(table)]
                        begin = time.perf_counter()
                        subprocess.run(command, check=True, capture_output=True)
                        elapsed = time.perf_counter() - begin
                        current = sha256(output)
                        if checksum is None:
                            checksum = current
                        assert current == checksum, (d, mode, version, sample)
                        if sample >= 0:
                            times[version].append(elapsed)
                            samples.append({'D': d, 'mode': mode, 'version': version, 'sample': sample,
                                            'seconds': elapsed, 'sha256': current})
                a, b = [statistics.median(times[version]) for version in ['before', 'after']]
                summary.append({'D': d, 'mode': mode, 'before_seconds': a,
                                'after_seconds': b, 'speedup': a / b, 'sha256': checksum})
                print(f'D={d}, seeds-{mode}: {a:.6f}s -> {b:.6f}s ({a / b:.2f}x)', flush=True)
    write_csv(Path(f'{args.output}-samples.csv'), samples)
    write_csv(Path(f'{args.output}-summary.csv'), summary)
    model = next(line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                 if line.startswith('model name'))
    metadata = {
        'cpu': args.cpu, 'cpu_model': model, 'platform': platform.platform(), 'samples': args.samples,
        'parameters': {'W': 30, 'E': 5, 'D': [64, 1024], 'R': 10, 'S': 1, 'length': 1000,
                       'lines': args.lines, 'traversal': 'outward', 'sequence_seed': 562, 'table_seed': 42},
        'protocol': 'One CPU, single-threaded commands, one untimed warmup per executable; alternating '
                    'before/after pairs; median whole-process wall time including table loading, conversion, '
                    'seeding and CSV I/O; excludes table generation; complete CSV hashes checked every run.',
        'input_sha256': input_sha, 'table_sha256': table_hashes,
        'executables': {name: str(path) for name, path in executables.items()},
        'executable_sha256': {name: sha256(path) for name, path in executables.items()},
        'source_sha256': {path: sha256(path) for path in ['src/ssh3_pivots.c', 'src/ssh3_deletions.c',
                                                        'src/ssh3_seed_cli.c', __file__]},
    }
    Path(f'{args.output}-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
