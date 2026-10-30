#!/usr/bin/env python3
"""Pristine-vs-pristine determinism (report §4.3 item 9, §5.2 item 9).

Usage: determinism-report.py <out.tsv> <name>=<run-a>:<run-b>:<sigma-per-test.tsv> ...

Each pair is two captures of the same tests on the same arm
(sweep-capture.sh) and the per-test output of
`elmer-oracle sigma-table <tests> <json> <run-a> b=<run-b>`. A test is:
  - deterministic: both runs pass and every field agrees bitwise
    (sigma 0, after the Perm mapping);
  - NONDETERMINISTIC: both pass but some field differs (the worst
    relative distance is reported, inf for structure);
  - no-fields: both pass, no field variables to compare;
  - pristine-fail: either run fails or times out (both verdicts listed).
"""
import collections
import sys


def matrix(run):
    rows = {}
    for line in open(run + '/matrix.tsv').read().splitlines()[1:]:
        c = line.split('\t')
        rows[c[0]] = c[1]
    return rows


out_path = sys.argv[1]
rows = []
totals = collections.Counter()
for spec in sys.argv[2:]:
    name, rest = spec.split('=', 1)
    a, b, sig = rest.split(':')
    va, vb = matrix(a), matrix(b)
    worst = {}
    for line in open(sig).read().splitlines():
        t, w, _n = line.split('\t')
        worst[t] = float(w)
    for t in sorted(va):
        x, y = va[t], vb.get(t, 'MISSING')
        if x != 'pass' or y != 'pass':
            verdict, detail = 'pristine-fail', '%s/%s' % (x, y)
        elif t not in worst:
            verdict, detail = 'no-fields', ''
        elif worst[t] == 0.0:
            verdict, detail = 'deterministic', ''
        else:
            verdict, detail = 'NONDETERMINISTIC', '%.3g' % worst[t]
        totals[(name, verdict)] += 1
        rows.append((name, t, verdict, detail))
with open(out_path, 'w') as f:
    f.write('set\ttest\tverdict\tdetail\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
for (name, verdict), n in sorted(totals.items()):
    print('%s\t%s\t%d' % (name, verdict, n))
