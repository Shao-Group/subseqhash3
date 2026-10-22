#!/usr/bin/env python3
"""Executable-level regression checks for the pivot-only R-by-S interface."""
import argparse
import csv
import math
from pathlib import Path
import struct
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    exe = parser.parse_args().executable.resolve()
    calls = 0

    def run(*args, ok=True):
        nonlocal calls
        calls += 1
        result = subprocess.run([str(exe), *map(str, args)], capture_output=True, text=True)
        assert (result.returncode == 0) == ok, (args, result.returncode, result.stdout, result.stderr)
        return result.stdout

    def rows(path):
        with path.open() as stream:
            reader = csv.DictReader(stream)
            assert reader.fieldnames == ['set', 'pivot', 'sequence', 'window_start',
                                         'substitutions', 'deletions', 'seed', 'psi', 'omega']
            return list(reader)

    with tempfile.TemporaryDirectory(prefix='ssh3-cli-') as directory:
        root = Path(directory)
        dna, table, out = (root / name for name in ['dna', 'table.piv', 'out.csv'])
        inputs = ['ACGTACGTTACGA', 'TTTTACGTACGTA']
        dna.write_text(inputs[0].lower() + '\r\n' + inputs[1] + '\n')
        run('tables', 7, 33, 3, table, '--seed', 42)
        header = struct.unpack('<8s7Q', table.read_bytes()[:64])
        assert header[2:5] == (8, 33, 3)  # n=row capacity+1, D, S
        help_text = run('--help')
        assert 'seeds-sub' in help_text and 'seeds-ind' in help_text
        for old in ['seeds', 'repeats', 'deletions']:
            run(old, '--help', ok=False)
        run('tables', 7, 33, 3, table, '--pivots', ok=False)
        for mode, e in [('sub', 2), ('ind', 2), ('ind', 0)]:
            for inward in [False, True]:
                flag = ['--inward'] if inward else []
                run('seeds-' + mode, dna, 7, e, 33, 2, 3, out, table, *flag)
                actual = rows(out)
                assert len(actual) == 2 * 7 * 2 * 3
                keys = [(int(a['sequence']), int(a['window_start']), int(a['set']), int(a['pivot']))
                        for a in actual]
                assert keys == [(seq, start, f, p) for seq in range(2) for start in range(7)
                                for f in range(3) for p in range(2)]
                assert all(len(a['seed']) == 7 - (1 if mode == 'sub' else e) for a in actual)
                for f in range(3):
                    single = root / f'single-{f}.piv'
                    seed = (42 + f * 0x9e3779b97f4a7c15) % 2**64
                    run('tables', 7, 33, 1, single, '--seed', seed)
                    one_out = root / 'one.csv'
                    run('seeds-' + mode, dna, 7, e, 33, 2, 1, one_out, single, *flag)
                    selected = [dict(a, set='0') for a in actual if int(a['set']) == f]
                    assert selected == rows(one_out)
                # Automatic and explicit table generation produce identical files and rows.
                auto = root / 'auto.csv'
                run('seeds-' + mode, dna, 7, e, 33, 2, 3, auto, '--seed', 42, *flag)
                assert rows(auto) == actual
                assert Path(str(auto) + '.piv').read_bytes() == table.read_bytes()
                # R changes only the leading pivots selected in EACH set.
                run('seeds-' + mode, dna, 7, e, 33, 1, 3, auto, table, *flag)
                assert rows(auto) == [a for a in actual if a['pivot'] == '0']
        # Validate dimensions, unsupported options and aliases without truncating OUTPUT.
        out.write_text('KEEP\n')
        for mode in ['sub', 'ind']:
            for r, s in [(0, 1), (7, 1), (2, 0), (2, 4), (2, 2**64-1)]:
                run('seeds-' + mode, dna, 7, 1, 33, r, s, out, table, ok=False)
                assert out.read_text() == 'KEEP\n'
            for flags in [['--inward', '--inward'], ['--subseqhash2-tables', table], ['--unknown']]:
                run('seeds-' + mode, dna, 7, 1, 33, 2, 3, out, table, *flags, ok=False)
            run('seeds-' + mode, dna, 7, 1, 33, 2, 3, table, table, ok=False)
            run('seeds-' + mode, dna, 7, 1, 33, 2, 3, dna, table, ok=False)
            bad = root / 'bad.piv'
            bad.write_text('bad\n')
            run('seeds-' + mode, dna, 7, 1, 33, 2, 3, out, bad, ok=False)
            assert out.read_text() == 'KEEP\n'
        # Extreme budgets/window sizes and R limits.
        for mode, w, e, r in [('sub', 2, 1, 1), ('sub', 257, 256, 1),
                               ('ind', 1, 0, 1), ('ind', 256, 0, 1), ('ind', 7, 6, 1)]:
            dna.write_text('A' * (w+1) + '\n')
            run('seeds-' + mode, dna, w, e, 1, r, 2, out, '--seed', 7)
            assert len(rows(out)) == 2*r*2
        dna.write_text('ACNTACGT\n')
        run('seeds-sub', dna, 7, 1, 33, 2, 3, out, table, ok=False)
        # Output hard links and generated-table aliases must preserve input.
        dna.write_text('ACGTACGT\n')
        alias = root / 'alias'
        alias.hardlink_to(dna)
        run('seeds-sub', dna, 7, 1, 33, 2, 3, alias, table, ok=False)
        auto = root / 'alias-auto.csv'
        auto_table = Path(str(auto) + '.piv')
        auto_table.hardlink_to(dna)
        run('seeds-sub', dna, 7, 1, 33, 2, 3, auto, '--seed', 42, ok=False)
        assert dna.read_text() == 'ACGTACGT\n'
        auto_table.unlink()
        auto.symlink_to(auto_table)
        run('seeds-sub', dna, 7, 1, 33, 2, 3, auto, '--seed', 42, ok=False)
        assert auto_table.read_bytes()[:8] == b'SQH3PIV\0'
        # Collision summaries must match seeds emitted by the actual seeding commands.
        pairs = [('ACGT', 'ACGT'), ('ACGT', 'TGCA'), ('AAAA', 'CCCC')]
        pair_file = root / 'pairs.csv'
        pair_file.write_text(''.join(a + ',' + b + '\n' for a, b in pairs))
        dna.write_text(''.join(a + '\n' + b + '\n' for a, b in pairs))
        run('tables', 4, 7, 2, table, '--seed', 99)
        for mode in ['sub', 'ind']:
            for inward in [False, True]:
                flag = ['--inward'] if inward else []
                run('seeds-' + mode, dna, 4, 1, 7, 2, 2, out, table, *flag)
                data = rows(out)
                report = dict(line.split(' ', 1) for line in
                              run('collisions', mode, pair_file, 1, 7, 2, 2, table, *flag).splitlines())
                seeds = {seq: [a['seed'] for a in data if int(a['sequence']) == seq] for seq in range(6)}
                aligned = [[a == b for a, b in zip(seeds[2*i], seeds[2*i+1])] for i in range(3)]
                expected = {
                    'any_seed_collision_probability': sum(bool(set(seeds[2*i]) & set(seeds[2*i+1])) for i in range(3))/3,
                    'any_same_set_pivot_collision_probability': sum(any(a) for a in aligned)/3,
                    'all_same_set_pivot_collision_probability': sum(all(a) for a in aligned)/3,
                }
                expected.update({f'set_{i//2}_pivot_{i%2}_collision_probability':
                                 sum(a[i] for a in aligned)/3 for i in range(4)})
                assert report['seeds_per_string'] == '4'
                for key, value in expected.items():
                    assert math.isclose(float(report[key]), value, abs_tol=1e-9), (key, report, value)
        # Unequal resulting lengths cannot collide; exact neighborhoods also have Jaccard zero.
        pair_file.write_text('ACG,ACGT\n')
        for mode in ['sub', 'ind']:
            report = dict(line.split(' ', 1) for line in
                          run('collisions', mode, pair_file, 0, 7, 1, 2, '--seed', 99, '--jaccard').splitlines())
            assert float(report['any_seed_collision_probability']) == 0
            assert float(report['average_jaccard']) == 0
        pair_file.write_text('ACGT,ACGT\n')
        report = run('collisions', 'sub', pair_file, 1, 7, 2, 2, table, '--jaccard')
        assert 'average_jaccard 1.000000000' in report
        for text in ['', 'ACNT,ACGT\n', 'ACGT;ACGT\n', 'A,ACGT\n']:
            pair_file.write_text(text)
            run('collisions', 'sub', pair_file, 1, 7, 2, 2, table, ok=False)
    print(f'pivot CLI: {calls} executable checks passed (R*S, both modes/traversals, tables, CSV, collisions and aliases)')


if __name__ == '__main__':
    main()
