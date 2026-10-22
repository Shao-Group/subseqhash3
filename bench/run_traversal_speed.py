#!/usr/bin/env python3
"""Measure all four traversal paths sequentially against the frozen baseline."""
import argparse
import csv
import datetime
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, default=Path('build/opt/bench-traversal-speed'))
    parser.add_argument('--cpu', type=int, default=min(os.sched_getaffinity(0)))
    parser.add_argument('--residues', type=int, nargs='+', default=[1, 32, 128, 1024, 2048])
    parser.add_argument('--flags', help='compiler flags used to build both implementations')
    args = parser.parse_args()
    assert args.executable.is_file() and args.cpu in os.sched_getaffinity(0)
    assert all(d > 0 for d in args.residues) and len(set(args.residues)) == len(args.residues)
    command = ['taskset', '-c', str(args.cpu), str(args.executable.resolve())]
    runtime = []
    for d in args.residues:
        print(f'Measuring D={d} alone on CPU {args.cpu}', flush=True)
        result = subprocess.run(command+['--residues', str(d)], check=True, capture_output=True, text=True)
        rows = list(csv.DictReader(io.StringIO(result.stdout)))
        assert len(rows) == 16 and all(int(row['D']) == d for row in rows)
        runtime.extend(rows)
    write_csv(Path('bench/traversal-speed-runtime.csv'), runtime)
    print('Measuring one-deletion, zero-substitution cases at D=1024', flush=True)
    result = subprocess.run(command+['--single'], check=True, capture_output=True, text=True)
    single = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(single) == 16
    write_csv(Path('bench/traversal-single-runtime.csv'), single)
    sources = ['src/ssh3_deletions.c', 'src/ssh3_pivots.c', 'src/single_deletion.h',
               'bench/reference/traversal_deletions.c', 'bench/reference/traversal_pivots.c']
    metadata = dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        machine=platform.uname()._asdict(), cpu=args.cpu, compiler_flags=args.flags,
        cpu_model=next(line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                       if line.startswith('model name')),
        executable=str(args.executable.resolve()), residues=args.residues,
        runtime='median of 7 samples, 48 calls per sample, 16 deterministic random windows; alternating implementation order; warmed workspaces; allocations, traceback and cleanup included; table generation excluded',
        equivalence='scores, seeds, positions, substitution flags on 16 random and 3 structured inputs before timing each condition',
        threading='one CPU, one thread per process; all measurements sequential',
        source_sha256={path: hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in sources})
    Path('bench/traversal-speed-metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(f'Saved {len(runtime)} general cases and {len(single)} one-deletion cases.', flush=True)


if __name__ == '__main__':
    main()
