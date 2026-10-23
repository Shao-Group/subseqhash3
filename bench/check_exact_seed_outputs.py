#!/usr/bin/env python3
"""Compare complete seeds and collision reports against checkpoint 36de034."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import subprocess
import tempfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, default=Path('build/typical/previous'))
    parser.add_argument('--after', type=Path, default=Path('build/simd/subseqhash3'))
    parser.add_argument('--output', type=Path, default=Path('bench/seed-runtime-exact-checks.json'))
    args = parser.parse_args()
    programs = [args.before.resolve(), args.after.resolve()]
    checks = []
    with tempfile.TemporaryDirectory(prefix='ssh3-exact-check-') as directory:
        root = Path(directory)
        for seed in [14731, 73191, 287551]:
            for d in [64, 1024]:
                table = root / 'tables.piv'
                subprocess.run([str(programs[0]), 'tables', '30', str(d), '2', str(table), '--seed', str(seed)], check=True, capture_output=True)
                for pattern in ['random', 'homopolymer', 'periodic']:
                    rng = random.Random(seed)
                    sequence = ''.join(rng.choices('ACGT', k=1000)) if pattern == 'random' else 'A'*1000 if pattern == 'homopolymer' else 'ACGT'*250
                    dna = root / 'dna.seq'
                    dna.write_text(sequence + '\n')
                    for mode in ['sub', 'ind']:
                        hashes = []
                        for index, program in enumerate(programs):
                            output = root / f'out{index}.csv'
                            subprocess.run([str(program), f'seeds-{mode}', str(dna), '30', '5', str(d), '10', '1', str(output), str(table)], check=True, capture_output=True)
                            hashes.append(sha(output))
                        assert hashes[0] == hashes[1], (seed, d, pattern, mode)
                        checks.append(dict(kind='seeds', W=30, E=5, D=d, R=10, S=1, length=1000, seed=seed, pattern=pattern, mode=mode, sha256=hashes[0]))
                # Multiple sets, partial pivots, both traversals, and a mixed input.
                dna.write_text('ACGT'*12 + '\n' + 'A'*30 + '\n' + 'TGCA'*10 + '\n')
                for mode in ['sub', 'ind']:
                    for inward in [False, True]:
                        hashes = []
                        for index, program in enumerate(programs):
                            output = root / f'out{index}.csv'
                            flags = ['--inward'] if inward else []
                            subprocess.run([str(program), f'seeds-{mode}', str(dna), '30', '5', str(d), '10', '2', str(output), str(table), *flags], check=True, capture_output=True)
                            hashes.append(sha(output))
                        assert hashes[0] == hashes[1]
                        checks.append(dict(kind='multi_set_seeds', D=d, seed=seed, mode=mode, inward=inward, sha256=hashes[0]))
                rng = random.Random(seed+1)
                pairs = root / 'pairs.csv'
                lines = []
                for i in range(40):
                    a = ''.join(rng.choices('ACGT', k=30))
                    b = list(a)
                    for position in rng.sample(range(30), i%6):
                        b[position] = rng.choice([x for x in 'ACGT' if x != b[position]])
                    lines.append(a + ',' + ''.join(b) + '\n')
                pairs.write_text(''.join(lines))
                for mode in ['sub', 'ind']:
                    reports = [subprocess.run([str(p), 'collisions', mode, str(pairs), '5', str(d), '10', '2', str(table)], check=True, capture_output=True).stdout for p in programs]
                    assert reports[0] == reports[1]
                    checks.append(dict(kind='collisions', D=d, seed=seed, mode=mode, pairs=40, sha256=hashlib.sha256(reports[0]).hexdigest()))
        # Non-power-of-two source dimensions keep their original behavior.
        for d in [1, 7, 31, 32, 33, 65, 128, 129]:
            dna = root / 'dna.seq'
            dna.write_text('ACGT'*14 + '\n')
            table = root / 'tables.piv'
            subprocess.run([str(programs[0]), 'tables', '12', str(d), '1', str(table), '--seed', '9182'], check=True, capture_output=True)
            for mode in ['sub', 'ind']:
                hashes = []
                for index, program in enumerate(programs):
                    output = root / f'out{index}.csv'
                    subprocess.run([str(program), f'seeds-{mode}', str(dna), '12', '3', str(d), '8', '1', str(output), str(table)], check=True, capture_output=True)
                    hashes.append(sha(output))
                assert hashes[0] == hashes[1]
                checks.append(dict(kind='other_dimensions', D=d, mode=mode, sha256=hashes[0]))
    args.output.write_text(json.dumps(dict(checkpoint='36de034', executables={str(p): sha(p) for p in programs}, checks=checks), indent=2) + '\n')
    print(f'{len(checks)} complete CSV/collision comparisons match checkpoint 36de034')


if __name__ == '__main__':
    main()
