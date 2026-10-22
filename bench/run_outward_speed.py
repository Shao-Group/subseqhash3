#!/usr/bin/env python3
"""Compare frozen and optimized window benchmarks, including exact checksums."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
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
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--baseline-ref', required=True)
    parser.add_argument('--cpu', type=int, default=0)
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--residues', type=int, nargs='+', default=[32, 1024])
    parser.add_argument('--flags', required=True)
    parser.add_argument('--output-prefix', type=Path, default=Path('bench/outward-speed'))
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    samples, results = [], []
    condition = 0
    for n, e in [(30, 4), (50, 8)]:
        for d in args.residues:
            for count in [10, 0]:
                for length in [n, 256, 4096]:
                    versions = {}
                    # Rotate executable order as well as the benchmark's method order.
                    order = ['before', 'after'] if condition % 2 == 0 else ['after', 'before']
                    for version in order:
                        command = [str(getattr(args, version)), str(n), str(e), str(d),
                                   str(count), str(length), str(args.samples)]
                        if length == 256:
                            command.append('--independent')
                        completed = subprocess.run(command, check=True, text=True, capture_output=True)
                        rows = list(csv.DictReader(io.StringIO(completed.stdout)))
                        groups = {}
                        for row in rows:
                            groups.setdefault((row['algorithm'], row['method']), []).append(row)
                            samples.append({'version': version} | row)
                        versions[version] = groups
                    for key, before in versions['before'].items():
                        after = versions['after'][key]
                        checksums = {r['checksum'] for r in before + after}
                        if len(checksums) != 1:
                            raise RuntimeError(f'Output mismatch: {n,e,d,count,length,key}: {checksums}')
                        # Keep inward samples as a control, report outward results.
                        if key[1] == 'inward':
                            continue
                        old = statistics.median(float(r['us_per_window']) for r in before)
                        new = statistics.median(float(r['us_per_window']) for r in after)
                        scenario = ('single-cold' if length == n else 'single-reused'
                                    if key[1] == 'outward-independent' else 'windows-shared')
                        results.append({k: before[0][k] for k in
                                        ['algorithm', 'n', 'E', 'D', 'seeds', 'length', 'windows']} |
                                       {'scenario': scenario, 'before_us': old, 'after_us': new,
                                        'speedup': old / new, 'samples': args.samples,
                                        'checksum': before[0]['checksum']})
                    write_csv(Path(str(args.output_prefix) + '-samples.csv'), samples)
                    write_csv(Path(str(args.output_prefix) + '-results.csv'), results)
                    print(f'checked W={n} E={e} D={d} seeds={count or "all"} length={length}', flush=True)
                    condition += 1
    sources = ['src/ssh3_pivots.c', 'src/ssh3_deletions.c', 'bench/windows.c',
               'bench/run_outward_speed.py']
    metadata = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'baseline_ref': args.baseline_ref, 'cpu': args.cpu,
        'cpu_model': next(line.split(':', 1)[1].strip() for line in
                          Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),
        'compiler_flags': args.flags, 'samples': args.samples, 'residues': args.residues,
        'protocol': 'One thread; median complete calls, alternating method and executable order; '
                    'cold single-window APIs, reused-workspace independent window loop at length 256, '
                    'shared windows at lengths 256/4096. Includes allocation, traceback, checksum callback '
                    'and cleanup; excludes table generation and disk I/O. Exact before/after checksums '
                    'required for every algorithm, method and condition, including inward controls.',
        'source_sha256': {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources},
        'executable_sha256': {v: hashlib.sha256(getattr(args, v).read_bytes()).hexdigest()
                              for v in ['before', 'after']},
    }
    Path(str(args.output_prefix) + '-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
