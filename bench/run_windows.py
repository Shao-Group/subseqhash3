#!/usr/bin/env python3
"""Sequential same-parameter sliding-window traversal benchmarks."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
from datetime import datetime, timezone


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, default=Path('build/opt/bench-windows'))
    parser.add_argument('--cpu', type=int, default=0)
    parser.add_argument('--lengths', nargs='+', type=int, default=[256, 4096])
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--flags', default='')
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    samples, results = [], []
    for length in args.lengths:
        for n, e in [(30, 4), (50, 8)]:
            for d in [32, 1024]:
                for count in [10, 0]:
                    command = [str(args.executable), str(n), str(e), str(d),
                               str(count), str(length), str(args.samples)]
                    # Include the old independent-window loop on the short input.
                    if length == min(args.lengths):
                        command += ['--independent']
                    completed = subprocess.run(command, check=True, text=True, capture_output=True)
                    rows = list(csv.DictReader(io.StringIO(completed.stdout)))
                    samples.extend(rows)
                    groups = {}
                    for row in rows:
                        key = (row['algorithm'], row['method'])
                        groups.setdefault(key, []).append(row)
                    for (_, _), group in groups.items():
                        row = group[0]
                        results.append({key: row[key] for key in
                                        ['algorithm', 'method', 'n', 'E', 'D', 'seeds', 'length', 'windows']} |
                                       {'median_us_per_window': statistics.median(float(r['us_per_window']) for r in group),
                                        'samples': len(group), 'checksum': row['checksum']})
                    write_csv(Path('bench/windows-samples.csv'), samples)
                    write_csv(Path('bench/windows-runtime.csv'), results)
                    print(f'finished length={length}, W={n}, E={e}, D={d}, seeds={count or "all"}', flush=True)
    sources = ['src/ssh3_pivots.c', 'src/ssh3_deletions.c', 'src/ssh3_pivots.h',
               'src/ssh3_deletions.h', 'bench/windows.c', 'bench/run_windows.py']
    model = next((line.split(':', 1)[1].strip() for line in
                  Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')), '')
    metadata = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'machine': platform.uname()._asdict(), 'cpu_model': model, 'cpu': args.cpu,
        'compiler_flags': args.flags, 'lengths': args.lengths, 'samples': args.samples,
        'protocol': 'single thread, sequential, median of alternating method order samples; entire sequence call including cold anchor-cache allocation, warm-up windows, output allocation, traceback, checksum callback and cleanup; excludes table generation and disk I/O',
        'equivalence': 'exact output checksum of first 64 windows versus independent-window APIs before timing; full-sequence outward checksums also compared when independent mode is timed',
        'source_sha256': {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources},
    }
    Path('bench/windows-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
