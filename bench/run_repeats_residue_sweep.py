#!/usr/bin/env python3
"""Run paired repeats experiments at several residue counts; aggregate raw counts."""
import argparse
import concurrent.futures
import csv
import datetime
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess

KEYS = ['n', 'substitutions', 'D', 'seeds', 'mutation', 'changes']
COUNTS = ['pairs', 'outward_same', 'inward_same', 'outward_any', 'inward_any',
          'outward_per', 'inward_per', 'nonidentical_pairs', 'inward_only', 'outward_only']
MUTATIONS = {'substitutions': 0, 'delete_insert': 1, 'unrelated': 2}


def read_csv(path):
    with open(path, newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, fields, rows):
    with open(path, 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def order(row):
    return tuple(int(row[k]) for k in KEYS[:4]) + (MUTATIONS[row['mutation']], int(row['changes']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, default=Path('build/repeats-fast/bench-repeats-traversal'))
    parser.add_argument('--output-prefix', type=Path, default=Path('bench/repeats-residue'))
    parser.add_argument('--parts-dir', type=Path, default=Path('build/repeats-residue-parts'))
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--aggregate-only', action='store_true', help='reuse saved collision batches; retain existing runtime measurements')
    parser.add_argument('--families', type=int, default=32)
    parser.add_argument('--pairs-per-family', type=int, default=64)
    parser.add_argument('--chunk', type=int, default=16)
    parser.add_argument('--residues', type=int, nargs='+', default=[32, 128, 1024, 2048])
    args = parser.parse_args()
    assert args.executable.is_file() and args.jobs > 0 and args.families >= 32
    assert args.pairs_per_family > 0 and args.chunk > 0
    assert all(d > 0 for d in args.residues) and len(set(args.residues)) == len(args.residues)
    args.parts_dir.mkdir(parents=True, exist_ok=True)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    executable = str(args.executable.resolve())
    jobs = [(d, c, offset, min(args.chunk, args.pairs_per_family-offset))
            for d in args.residues for c in range(4)
            for offset in range(0, args.pairs_per_family, args.chunk)]
    # Schedule longer jobs first. Each child is single-threaded; only collision
    # sampling runs concurrently. Runtime measurements run after all children exit.
    jobs.sort(key=lambda job: job[0] * [0.95, 2.6, 2.6, 12.2][job[1]], reverse=True)

    def run(job):
        d, config, offset, size = job
        stem = args.parts_dir / f'd{d}-c{config}-offset{offset}'
        command = [executable, '--collisions', str(args.families), str(size), str(stem)+'.raw.csv',
                   str(d), str(config), str(offset)]
        with open(str(stem)+'.summary.csv', 'w') as out, open(str(stem)+'.log', 'w') as err:
            subprocess.run(command, stdout=out, stderr=err, check=True)
        return stem

    stems = []
    if args.aggregate_only:
        stems = [args.parts_dir / f'd{d}-c{config}-offset{offset}' for d, config, offset, _ in jobs]
        assert all(Path(str(stem)+'.raw.csv').is_file() for stem in stems)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(run, job): job for job in jobs}
            for future in concurrent.futures.as_completed(futures):
                stems.append(future.result())
                d, config, offset, size = futures[future]
                print(f'Collision batches {len(stems)}/{len(jobs)}: D={d}, config={config}, pairs={offset}:{offset+size}', flush=True)

    merged = {}
    offsets = {}
    for stem in stems:
        for row in read_csv(str(stem)+'.raw.csv'):
            key = tuple(row[k] for k in KEYS+['family'])
            offset, size = int(row['pair_offset']), int(row['pairs'])
            if key not in merged:
                merged[key] = {**{k: row[k] for k in KEYS+['family', 'table_seed']},
                               **{field: 0 for field in COUNTS}, 'pair_offset': 0}
                offsets[key] = []
            result = merged[key]
            assert result['table_seed'] == row['table_seed']
            offsets[key].append((offset, size))
            for field in COUNTS:
                result[field] += int(row[field])
    for key, row in merged.items():
        expected = 0
        for offset, size in sorted(offsets[key]):
            assert offset == expected
            expected += size
        assert expected == args.pairs_per_family == row['pairs']
        for mode in ['inward', 'outward']:
            assert 0 <= row[mode+'_same'] <= row[mode+'_any'] <= row['pairs']
            assert row[mode+'_same'] <= row[mode+'_per'] <= row['pairs']*int(row['seeds'])
        assert row['inward_same']-row['outward_same'] == row['inward_only']-row['outward_only']

    groups = {}
    for row in merged.values():
        groups.setdefault(tuple(row[k] for k in KEYS), []).append(row)
    summaries = []
    for key, group in groups.items():
        assert len(group) == args.families
        result = dict(zip(KEYS, key))
        total = sum(row['pairs'] for row in group)
        result.update(families=args.families, pairs=total,
                      nonidentical_pairs=sum(row['nonidentical_pairs'] for row in group))
        for field in ['same', 'any']:
            for mode in ['outward', 'inward']:
                result[mode+'_'+field+'_pct'] = 100*sum(row[mode+'_'+field] for row in group)/total
            differences = [100*(row['inward_'+field]-row['outward_'+field])/row['pairs'] for row in group]
            delta = statistics.mean(differences)
            half = 2.04*statistics.stdev(differences)/math.sqrt(len(group))
            result['delta_'+field+'_pp'] = delta
            result['delta_'+field+'_low95'] = delta-half
            result['delta_'+field+'_high95'] = delta+half
        for mode in ['outward', 'inward']:
            result[mode+'_per_pct'] = 100*sum(row[mode+'_per'] for row in group)/(total*int(result['seeds']))
            result[mode+'_only'] = sum(row[mode+'_only'] for row in group)
        summaries.append(result)
    summaries.sort(key=order)
    assert len(summaries) == 20*len(args.residues)
    family_rows = sorted(merged.values(), key=lambda row: order(row)+(int(row['family']),))
    family_fields = KEYS+['family', 'table_seed', 'pairs', 'outward_same', 'inward_same',
                         'outward_any', 'inward_any', 'outward_per', 'inward_per',
                         'nonidentical_pairs', 'inward_only', 'outward_only', 'pair_offset']
    write_csv(str(args.output_prefix)+'-families.csv', family_fields, family_rows)
    fields = KEYS+['families', 'pairs', 'nonidentical_pairs', 'outward_same_pct', 'inward_same_pct',
                  'delta_same_pp', 'delta_same_low95', 'delta_same_high95', 'outward_any_pct', 'inward_any_pct',
                  'delta_any_pp', 'delta_any_low95', 'delta_any_high95', 'outward_per_pct', 'inward_per_pct',
                  'inward_only', 'outward_only']
    write_csv(str(args.output_prefix)+'-collisions.csv', fields,
              [{k: f'{v:.4f}' if isinstance(v, float) else v for k, v in row.items()} for row in summaries])

    # Verify the same source-pair stream is used across D and pivot counts.
    pair_counts = {}
    lookup = {}
    for row in family_rows:
        key = tuple(row[k] for k in ['n', 'substitutions', 'mutation', 'changes', 'family'])
        counts = (row['pairs'], row['nonidentical_pairs'])
        assert pair_counts.setdefault(key, counts) == counts
        lookup[tuple(row[k] for k in KEYS+['family'])] = row
    for row in family_rows:
        if row['seeds'] != '10':
            continue
        key = tuple(str(int(row['n'])-1) if k == 'seeds' else row[k] for k in KEYS+['family'])
        larger = lookup[key]
        assert all(larger[mode+'_'+field] >= row[mode+'_'+field]
                   for mode in ['inward', 'outward'] for field in ['same', 'any', 'per'])

    baseline = Path('bench/repeats-traversal-collisions.csv')
    if 1024 in args.residues and args.families == 32 and args.pairs_per_family == 64 and baseline.is_file():
        old = {tuple(row[k] for k in KEYS): row for row in read_csv(baseline)}
        for row in summaries:
            if row['D'] != '1024':
                continue
            previous = old[tuple(row[k] for k in KEYS)]
            assert all(abs(float(row[k])-float(previous[k])) < .000051 for k in fields if k not in KEYS)
        print('D=1024 reproduces all 20 previous collision conditions exactly.', flush=True)

    # Consecutive-D comparisons use the same input stream and table-family
    # identifiers. Realized tables vary with D; family pairing captures that
    # variation as part of the observed D effect.
    contrasts = []
    residues = sorted(args.residues)
    for source_d, target_d in zip(residues, residues[1:]):
        for summary in summaries:
            if int(summary['D']) != source_d:
                continue
            source_key = tuple(summary[k] for k in KEYS)
            for mode in ['inward', 'outward']:
                before, after = [], []
                for family in range(args.families):
                    origin = lookup[source_key+(str(family),)]
                    target_key = tuple(str(target_d) if k == 'D' else summary[k] for k in KEYS)+(str(family),)
                    target = lookup[target_key]
                    before.append(100*origin[mode+'_same']/origin['pairs'])
                    after.append(100*target[mode+'_same']/target['pairs'])
                differences = [b-a for a, b in zip(before, after)]
                mean = statistics.mean(differences)
                half = 2.04*statistics.stdev(differences)/math.sqrt(args.families)
                contrasts.append({**{k: summary[k] for k in KEYS if k != 'D'},
                    'from_D': source_d, 'to_D': target_d, 'traversal': mode,
                    'from_pct': statistics.mean(before), 'to_pct': statistics.mean(after),
                    'difference_pp': mean, 'low95': mean-half, 'high95': mean+half})
    if contrasts:
        write_csv(str(args.output_prefix)+'-differences.csv', list(contrasts[0]),
                  [{k: f'{v:.4f}' if isinstance(v, float) else v for k, v in row.items()} for row in contrasts])
    if args.aggregate_only:
        print(f'Aggregated {len(summaries)} conditions and {len(contrasts)} consecutive-D comparisons; existing runtimes retained.')
        return

    runtimes = []
    for d in sorted(args.residues):
        print(f'Measuring runtime alone: D={d}', flush=True)
        completed = subprocess.run([executable, '--runtime', str(d)], capture_output=True, text=True, check=True)
        import io
        runtimes.extend(csv.DictReader(io.StringIO(completed.stdout)))
    write_csv(str(args.output_prefix)+'-runtime.csv', list(runtimes[0]), runtimes)
    metadata = dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    machine=platform.uname()._asdict(), executable=executable,
                    residues=args.residues, families=args.families, pairs_per_family=args.pairs_per_family,
                    jobs=args.jobs, chunk=args.chunk, canonical_pair_stream='D-independent; preserves previous D=1024 inputs',
                    runtime='sequential, after all collision workers finish; median of 7 samples, 48 window calls per sample',
                    confidence='paired table-family differences, mean +/- 2.04*sample_stdev/sqrt(families)',
                    consecutive_D_intervals='paired table-family differences, target D minus source D')
    Path(str(args.output_prefix)+'-metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(f'Finished: {len(summaries)} collision conditions, {len(family_rows)} family rows, {len(runtimes)} runtime rows.', flush=True)


if __name__ == '__main__':
    main()
