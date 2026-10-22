#!/usr/bin/env python3
"""Calibrate D at fixed seed count, then sample paired collision probabilities."""
import argparse
import concurrent.futures
import csv
import datetime
import io
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess


def rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def write(path, data):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, list(data[0]))
        w.writeheader()
        w.writerows(data)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable', default='build/bench-traversal-matched')
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--families', type=int, default=32)
    p.add_argument('--pairs', type=int, default=64)
    p.add_argument('--aggregate-only', action='store_true')
    p.add_argument('--calibrate-only', action='store_true')
    p.add_argument('--skip-calibration', action='store_true')
    a = p.parse_args()
    assert a.jobs > 0 and a.families >= 32 and a.pairs > 0
    directory = Path('build/traversal-matched')
    directory.mkdir(parents=True, exist_ok=True)
    prefix = 'bench/traversal-matched'
    configs = [(alg, n, e, r, di) for alg in ['deletions', 'repeats']
               for n, e in [(30, 4), (50, 8)] for r in [10, n-(e if alg == 'deletions' else 1)]
               for di in [1024, 2048]]
    trials, calibrated = [], []
    if not a.aggregate_only and not a.skip_calibration:
        for alg, n, e, r, di in configs:
            cache = {}
            def measure(d):
                result = subprocess.run([a.executable, 'runtime', alg, str(n), str(e), str(r), str(di), str(d)],
                                        check=True, text=True, capture_output=True)
                t = {k: float(v) for k, v in next(csv.DictReader(io.StringIO(result.stdout))).items()}
                cache[d] = t
                trials.append(dict(algorithm=alg, n=n, E=e, seeds=r, inward_D=di, outward_D=d, **t))
                return t['outward_us']/t['inward_us']
            lo, hi = 1, di
            # Find a D close to the budget. Performance can jump at SIMD and
            # traceback thresholds, so also inspect nearby vector boundaries.
            while hi-lo > 1:
                mid = (lo+hi)//2
                if measure(mid) <= 1: lo = mid
                else: hi = mid
            for d in sorted({lo, hi, max(1, lo//32*32), (lo//32+1)*32}):
                if d not in cache: measure(d)
            feasible = [d for d, t in cache.items() if t['outward_us'] <= t['inward_us']]
            assert feasible, (alg, n, e, r, di)
            selected = max(feasible)
            # Fresh validation; choose using calibration only, never collisions.
            validation = measure(selected)
            t = cache[selected]
            calibrated.append(dict(algorithm=alg, n=n, E=e, seeds=r, inward_D=di,
                                   outward_D=selected, **t, outward_over_inward=validation))
            write(prefix+'-calibration.csv', trials)
            write(prefix+'-runtime.csv', calibrated)
            print('Calibrated', calibrated[-1], flush=True)
    else:
        calibrated = rows(prefix+'-runtime.csv')
    if not a.aggregate_only:
        # Close residual gaps at outward SIMD/parent-storage thresholds by
        # tuning inward D. This selection still uses runtime alone.
        for c in calibrated:
            c['budget_inward_D'] = c.get('budget_inward_D', c['inward_D'])
            best = dict(c)
            tried = {int(c['inward_D'])}
            for _ in range(4):
                ratio = float(c['outward_over_inward'])
                if .96 <= ratio <= 1.04: break
                di = max(32, round(int(c['inward_D'])*ratio/32)*32)
                if di in tried: break
                tried.add(di)
                result = subprocess.run([a.executable, 'runtime', c['algorithm'],
                    *[str(c[k]) for k in ['n', 'E', 'seeds']], str(di), str(c['outward_D'])],
                    check=True, text=True, capture_output=True)
                t = {k: float(v) for k, v in next(csv.DictReader(io.StringIO(result.stdout))).items()}
                c.update(inward_D=di, **t, outward_over_inward=t['outward_us']/t['inward_us'])
                trials.append(dict(algorithm=c['algorithm'], n=c['n'], E=c['E'], seeds=c['seeds'],
                                   inward_D=di, outward_D=c['outward_D'], **t))
                if abs(math.log(float(c['outward_over_inward']))) < abs(math.log(float(best['outward_over_inward']))):
                    best = dict(c)
            c.update(best)
            print('Refined', c, flush=True)
        write(prefix+'-runtime.csv', calibrated)
        if trials:
            if a.skip_calibration:
                trials = rows(prefix+'-calibration.csv')+trials
            write(prefix+'-calibration.csv', trials)
    if a.calibrate_only: return
    def work(job):
        i, conf, kind, changes = job
        path = directory/f'c{i}-k{kind}-e{changes}.csv'
        if not a.aggregate_only:
            with open(path, 'w') as out:
                subprocess.run([a.executable, 'collisions', conf['algorithm'], *[str(conf[k]) for k in
                                ['n', 'E', 'seeds', 'inward_D', 'outward_D']], str(a.families), str(a.pairs),
                                str(kind), str(changes)], stdout=out, check=True)
        return conf, kind, changes, rows(path)
    jobs = [(i, c, kind, changes) for i, c in enumerate(calibrated)
            for kind, changes in [(0, int(c['E'])//2), (0, int(c['E'])), (1, 1), (1, 2), (2, 0)]]
    raw, summaries = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for c, kind, changes, data in pool.map(work, jobs):
            conf = {k: c[k] for k in ['algorithm', 'n', 'E', 'seeds', 'inward_D', 'outward_D']}
            conf.update(mutation=['substitutions', 'delete_insert', 'unrelated'][kind], changes=changes)
            assert len(data) == a.families
            assert [int(row['family']) for row in data] == list(range(a.families))
            for row in data:
                assert int(row['pairs']) == a.pairs
                assert int(row['inward_same'])-int(row['outward_same']) == int(row['inward_only'])-int(row['outward_only'])
                for mode in ['inward', 'outward']:
                    assert 0 <= int(row[mode+'_same']) <= int(row[mode+'_any']) <= a.pairs
                raw.append(dict(**conf, **row))
            summary = dict(**conf, families=a.families, pairs=a.families*a.pairs)
            for metric in ['same', 'any', 'per']:
                scale = 100/a.pairs/(int(c['seeds']) if metric == 'per' else 1)
                for mode in ['inward', 'outward']:
                    summary[mode+'_'+metric+'_pct'] = statistics.mean(int(row[mode+'_'+metric])*scale for row in data)
                delta = [(int(row['inward_'+metric])-int(row['outward_'+metric]))*scale for row in data]
                mean = statistics.mean(delta)
                half = 2.04*statistics.stdev(delta)/math.sqrt(a.families)
                summary.update({f'delta_{metric}_pp': mean, f'delta_{metric}_low95': mean-half,
                                f'delta_{metric}_high95': mean+half})
            summaries.append(summary)
            print('Completed', conf, flush=True)
    write(prefix+'-families.csv', raw)
    write(prefix+'-collisions.csv', summaries)
    metadata_path = Path(prefix+'-metadata.json')
    provenance = json.loads(metadata_path.read_text()) if a.aggregate_only and metadata_path.exists() else dict(
        timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        machine=platform.uname()._asdict(), executable=str(Path(a.executable).resolve()))
    provenance.update(families=a.families, pairs_per_family=a.pairs,
        budget_inward_D=[1024,2048], runtime='median 7 samples of 48 windows, warmed workspaces, allocations/traceback included; calibration sequential before collisions',
        calibration='binary search plus nearby multiples of 32; largest measured feasible outward D; fresh timing validation; refine inward D to close threshold gaps',
        confidence='paired family mean +/- 2.04 stdev / sqrt(families); exploratory',
        pair_stream='canonical D-independent stream from repeats residue sweep')
    metadata_path.write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
