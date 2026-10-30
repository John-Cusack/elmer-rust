#!/usr/bin/env python3
"""Build the committed golden index and grading-class table from a capture.

Usage: golden-index.py <tests-dir> <golden-root> <arm-build-dir> <tolerances.json> <out-dir>

Writes <out-dir>/quick-release.tsv (one row per registered test: verdict,
runner mode, result files, sha256 of its golden manifest) and
<out-dir>/classes.tsv (one row per test and class with its grading mode).
"""
import glob
import hashlib
import json
import os
import re
import sys

tests_dir, golden, build, tol_path, out = sys.argv[1:6]
reg = [l.split('\t') for l in open(os.path.join(golden, 'registry.tsv')).read().splitlines()]
tol = json.load(open(tol_path))

worst = {}
structure = set()
for e in tol['entries']:
    s = e['sigma_f']
    if s is None:
        structure.add(e['test'])
        s = float('inf')
    worst[e['test']] = max(worst.get(e['test'], 0.0), s)

def sif_texts(d):
    for f in sorted(glob.glob(os.path.join(tests_dir, d, '*.sif'))):
        yield open(f, 'rb').read().decode('latin-1')

def active(pattern, text):
    return re.search(r'(?im)^[^!#\n]*' + pattern, text) is not None

rows, classes = [], []
for name, d, np_, labels in reg:
    run = os.path.join(golden, name, 'run')
    verdict = open(os.path.join(golden, name, 'verdict.tsv')).read().split('\t')
    results = []
    for root, _, files in os.walk(run):
        for f in files:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, run)
            if rel.endswith('.pos') or os.path.exists(os.path.join(tests_dir, d, rel)):
                continue
            head = open(p, 'rb').read(16).lstrip()
            if head.startswith(b'BINARY 3.') or head.startswith(b'ASCII 3'):
                results.append(rel)
    results.sort()
    man = open(os.path.join(golden, name, 'golden.manifest'), 'rb').read()
    rows.append([name, d, np_, verdict[1], verdict[2], str(len(results)), ';'.join(results),
                 hashlib.sha256(man).hexdigest()])

    texts = list(sif_texts(d))
    tols = []
    for t in texts:
        for m in re.finditer(r'(?im)^[^!#\n]*reference\s+norm\s+tolerance\s*=\s*(?:real\s+)?([0-9.eEdD+-]+)', t):
            try:
                tols.append(float(m.group(1).replace('d', 'e').replace('D', 'e')))
            except ValueError:
                pass
    reftol = min(tols) if tols else 1e-5
    eigen = any(active(r'eigen\s+analysis\s*=\s*(logical\s+)?true', t) for t in texts)
    adaptive = any(active(r'adaptive\s+mesh\s+refinement\s*=\s*(logical\s+)?true', t) for t in texts)
    def cls(c, grading, why):
        classes.append([name, c, grading, why])
    if eigen:
        cls('eigen', 'eigen-mode', 'Eigen Analysis: sort eigenvalues, normalise signs, compare clusters by principal angles')
    if adaptive:
        cls('adaptive', 'O1+O4', 'Adaptive Mesh Refinement: path-dependent remeshing')
    if name in ('HeatSolveVecCondMatrix', 'StatCurrentVecResMatrix'):
        cls('lumped-matrix', 'matrix-file', 'golden is the lumped matrix .dat the solver writes; no field file')
    elif name not in worst:
        cls('no-field-output', 'O1+logs', 'no field variables in any result file; golden is its logs and SaveScalars output')
    if name in structure:
        cls('structure-differs', 'O1+O4+O6', 'a field differs in structure between sensitivity arms')
    elif not eigen and name in worst and worst[name] > 0.1 * reftol:
        cls('tier-a-unfit', 'O1+O4+O6', 'sigma_F %.3g > 0.1 x RefNorm tolerance %.1g' % (worst[name], reftol))
    if glob.glob(os.path.join(build, 'fem', 'tests', d, '*.so')):
        cls('test-local-plugin', 'D12-fixture', 'compiles a test-local Fortran plugin; Rust fixture port needed')
    if int(np_) > 1:
        cls('mpi', 'per-rank', 'np=%s; result files per partition' % np_)
    if verdict[2] == 'replay':
        cls('replay', 'cmake-replay', 'runtest.cmake has custom logic; captured through cmake -P')

os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'quick-release.tsv'), 'w') as f:
    f.write('test\tdir\tnp\tverdict\tmode\tresults\tresult_files\tgolden_manifest_sha256\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
with open(os.path.join(out, 'classes.tsv'), 'w') as f:
    f.write('test\tclass\tgrading\treason\n')
    for c in classes:
        f.write('\t'.join(c) + '\n')
