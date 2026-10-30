#!/usr/bin/env python3
"""Seal the held-out set (report §4.1 O5, §5.2 item 9).

Usage:
  heldout-seal.py seal <build-dir> <src-root> <ctest-baseline.log> <secret-dir> <public-dir>
  heldout-seal.py lane-tree <src-root> <secret-dir> <out-dir>

`seal` draws the sealed set and writes it only to <secret-dir> (owned by the
oracle, never committed); <public-dir> gets the method, the per-stratum
counts, the exclusions and a sha256 commitment to the secret file.

Units are test directories, not ctest names: a directory any quick test
uses is never sealed, because the lane tree removes whole directories.
Pool: every fem/tests directory with no quick test, plus every elmerice
test directory. Excluded from the pool:
  - D12: directories whose tests compile a test-local Fortran plugin;
  - pristine failures: any of its tests failed or timed out in the
    pristine ctest baseline (it cannot be graded against a pass).
fem directories are stratified by (serial|parallel, family label) and 30%
of each stratum is sealed, ranked by sha256(salt + name) with a secret
random salt; every eligible elmerice directory is sealed (§5.2).

`lane-tree` copies fem/tests and elmerice/Tests without the sealed
directories, for candidate lanes.
"""
import collections
import hashlib
import json
import os
import re
import secrets
import shutil
import subprocess
import sys

FRACTION = 0.30
MIN_SEALED = 21  # a 10-point gap must mean more than ~2 tests (§4.1 O5)
GENERIC = {'serial', 'parallel', 'quick', 'slow', 'benchmark'}


def registry(oracle, build, suite):
    out = subprocess.run([oracle, 'registry', build, '--suite', suite],
                         check=True, capture_output=True, text=True).stdout
    return [l.split('\t') for l in out.splitlines()]


def baseline_failures(log):
    bad = set()
    for line in open(log, encoding='latin-1'):
        m = re.search(r'Test +#\d+: (\S+) \.+\s*\*+(\w+)', line)
        if m:
            bad.add(m.group(1))
    return bad


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()


def seal(build, src, log, secret_dir, public_dir):
    oracle = os.environ.get('ELMER_ORACLE', 'elmer-oracle')
    failing = baseline_failures(log)
    dirs = collections.defaultdict(lambda: {'tests': [], 'labels': set(), 'np': set()})
    for suite, rel in (('fem/tests', 'fem/tests'), ('elmerice', 'elmerice/Tests')):
        for name, d, np_, labels in registry(oracle, build, suite):
            e = dirs[(rel, d)]
            e['tests'].append(name)
            e['labels'].update(labels.split(';'))
            e['np'].add(int(np_))
    pool, excluded = [], collections.defaultdict(list)
    for (rel, d), e in sorted(dirs.items()):
        if rel == 'fem/tests' and 'quick' in e['labels']:
            continue
        key = '%s/%s' % (rel, d)
        bdir = os.path.join(build, rel, d)
        if os.path.isdir(bdir) and any(f.endswith('.so') for f in os.listdir(bdir)):
            excluded['d12-test-local-plugin'].append(key)
            continue
        if any(t in failing for t in e['tests']):
            excluded['pristine-failing'].append(key)
            continue
        pool.append((rel, d, e))

    freq = collections.Counter(l for rel, d, e in pool if rel == 'fem/tests'
                               for l in e['labels'] - GENERIC)
    strata = collections.defaultdict(list)
    for rel, d, e in pool:
        if rel == 'elmerice/Tests':
            strata['elmerice'].append((rel, d))
            continue
        kind = 'parallel' if max(e['np']) > 1 else 'serial'
        fams = sorted(e['labels'] - GENERIC, key=lambda l: (-freq[l], l))
        strata['%s/%s' % (kind, fams[0] if fams else 'unlabelled')].append((rel, d))

    salt = secrets.token_hex(16)
    sealed, counts = [], {}
    for s in sorted(strata):
        members = sorted(strata[s], key=lambda m: hashlib.sha256((salt + m[1]).encode()).hexdigest())
        k = len(members) if s == 'elmerice' else int(FRACTION * len(members) + 0.5)
        counts[s] = {'pool': len(members), 'sealed': k}
        sealed += [{'suite': rel, 'dir': d, 'stratum': s, 'tests': sorted(dirs[(rel, d)]['tests'])}
                   for rel, d in members[:k]]
    if len(sealed) < MIN_SEALED:
        sys.exit('sealed set too small: %d < %d' % (len(sealed), MIN_SEALED))

    secret = {'version': 1, 'salt': salt, 'sealed': sealed}
    if os.path.exists(os.path.join(secret_dir, 'sealed.json')):
        sys.exit('a seal already exists in %s; resealing needs an ADR' % secret_dir)
    os.makedirs(secret_dir, exist_ok=True)
    blob = canonical(secret)
    with open(os.path.join(secret_dir, 'sealed.json'), 'wb') as f:
        f.write(blob)
    public = {
        'version': 1,
        'method': 'directories with no quick test; fem strata (serial|parallel, most common family label) sealed at %.0f%% (rounded half up) ranked by sha256(salt+dir) with a secret salt; all eligible elmerice dirs sealed' % (100 * FRACTION),
        'commitment_sha256': hashlib.sha256(blob).hexdigest(),
        'secret_location': os.path.join(secret_dir, 'sealed.json'),
        'pool_dirs': len(pool),
        'sealed_dirs': len(sealed),
        'sealed_tests': sum(len(s['tests']) for s in sealed),
        'sealed_elmerice_dirs': sum(1 for s in sealed if s['suite'] == 'elmerice/Tests'),
        'min_sealed': MIN_SEALED,
        'strata': counts,
        'excluded': {k: sorted(v) for k, v in excluded.items()},
    }
    os.makedirs(public_dir, exist_ok=True)
    with open(os.path.join(public_dir, 'seal.json'), 'w') as f:
        json.dump(public, f, indent=1, sort_keys=True)
        f.write('\n')
    print('pool %d sealed %d dirs (%d tests), elmerice %d; excluded %s' % (
        len(pool), len(sealed), public['sealed_tests'], public['sealed_elmerice_dirs'],
        {k: len(v) for k, v in excluded.items()}))


def lane_tree(src, secret_dir, out):
    secret = json.load(open(os.path.join(secret_dir, 'sealed.json')))
    sealed = {(s['suite'], s['dir']) for s in secret['sealed']}
    if os.path.exists(out):
        sys.exit('%s exists; build the lane tree into a fresh directory' % out)
    for rel in ('fem/tests', 'elmerice/Tests'):
        root = os.path.join(src, rel)
        drop = {d for r, d in sealed if r == rel}
        shutil.copytree(root, os.path.join(out, rel), symlinks=True,
                        ignore=lambda path, names, root=root, drop=drop:
                        [n for n in names if path == root and n in drop])
    print('lane tree at %s without %d sealed dirs' % (out, len(sealed)))


if __name__ == '__main__':
    if len(sys.argv) == 7 and sys.argv[1] == 'seal':
        seal(*sys.argv[2:7])
    elif len(sys.argv) == 5 and sys.argv[1] == 'lane-tree':
        lane_tree(*sys.argv[2:5])
    else:
        sys.exit(__doc__)
