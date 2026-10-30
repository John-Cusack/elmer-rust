#!/usr/bin/env python3
"""Build the capability matrix (report §5.2 item 6) from a capability sweep.

Usage: capability-matrix.py <tests-dir> <sweep-out> <classes.tsv> <out-dir>

Per registered test, from the LISTUSAGE arm's logs:
  - procedures: every `AddEquationBasics: Using procedure: <lib> <fn>`,
    marked `injected` when no sif names that library in a `Procedure`
    (legacy `Equation = Heat Equation` solvers, the VTU writer `Post File`
    adds, ...);
  - keywords present in the sifs, per section kind, split into used and
    unused by Elmer's own DEVEL_LISTUSAGE report;
  - queried_absent: keywords the code asked for that no sif set (the
    ListFind-miss trace): the defaults a capability check over present
    keywords would miss;
  - classes from classes.tsv.
Writes quick.json, keywords.tsv and modules.tsv into <out-dir>.
"""
import collections
import glob
import json
import os
import re
import sys

tests_dir, sweep, classes_path, out = sys.argv[1:5]
reg = [l.split('\t') for l in open(os.path.join(sweep, 'registry.tsv')).read().splitlines()]

classes = collections.defaultdict(list)
for line in open(classes_path).read().splitlines()[1:]:
    t, c = line.split('\t')[:2]
    classes[t].append(c)

SECTION_KINDS = ['boundary condition', 'initial condition', 'body force', 'simulation',
                 'constants', 'equation', 'material', 'solver', 'body', 'component',
                 'boundary', 'header', 'run control']


def strip_comment(line):
    q = False
    for i, c in enumerate(line):
        if c == '"':
            q = not q
        elif c in '!#' and not q:
            return line[:i]
    return line


def section_kind(name):
    n = name.strip().lower()
    for k in SECTION_KINDS:
        if n.startswith(k):
            return k.title()
    return None


def opens_table(code):
    v = code.rsplit('=', 1)[-1].lower().split()
    while v and v[-1] in ('cubic', 'monotone', 'harmonic', 'linear'):
        v.pop()
    return bool(v) and v[-1] in ('real', 'integer', 'logical', 'string', 'file')


def sif_keywords(text):
    """(section kind, keyword) pairs, following LoadInputFile's top level."""
    out = []
    section = None
    in_table = in_comment = False
    for raw in text.splitlines():
        code = strip_comment(raw).strip()
        low = code.lower()
        if in_comment:
            in_comment = low != '*/'
            continue
        if in_table:
            in_table = not low.startswith('end')
            continue
        if section is not None:
            if low == 'end':
                section = None
                continue
            if '=' in code:
                key = ' '.join(code.split('=', 1)[0].split()).lower()
                key = re.sub(r'\(.*\)$', '', key).strip()
                if section:
                    out.append((section, key))
                in_table = opens_table(code)
            continue
        if not low or low.startswith('$'):
            continue
        if '::' in low:
            sec, rest = low.split('::', 1)
            kind = section_kind(sec)
            if kind and '=' in rest:
                key = ' '.join(rest.split('=', 1)[0].split())
                out.append((kind, re.sub(r'\(.*\)$', '', key).strip()))
            continue
        if low.startswith(('include', 'echo ', 'check ')) or low in ('run', 'stop'):
            continue
        if low == '/*':
            in_comment = True
            continue
        section = section_kind(low) or ''
    return out


def sif_procedures(text):
    libs = set()
    for raw in text.splitlines():
        code = strip_comment(raw)
        m = re.match(r'(?i)\s*procedure\s*=\s*(?:file\s+)?"([^"]+)"', code)
        if m:
            libs.add(os.path.basename(m.group(1)).lower())
    return libs


tests = {}
kw_present = collections.Counter()
kw_unused = collections.Counter()
kw_absent = collections.Counter()
mod_tests = collections.defaultdict(set)
mod_injected = collections.defaultdict(set)

for name, d, np_, labels in reg:
    base = os.path.join(sweep, name)
    rc = open(os.path.join(base, 'rc')).read().strip() if os.path.exists(os.path.join(base, 'rc')) else 'missing'
    run = os.path.join(base, 'run')
    src_sifs = sorted(glob.glob(os.path.join(tests_dir, d, '*.sif')))
    present = set()
    declared = set()
    original_text = ''
    for f in src_sifs:
        t = open(f, 'rb').read().decode('latin-1')
        original_text += t
        present.update(sif_keywords(t))
        declared |= sif_procedures(t)
    logs = sorted(glob.glob(os.path.join(run, '*stdout*.log')))
    procs, unused, absent = set(), set(), set()
    for lg in logs:
        in_unused = False
        for line in open(lg, 'rb').read().decode('latin-1').splitlines():
            if line.startswith('AddEquationBasics: Using procedure:'):
                parts = line.split(':', 2)[2].split()
                if len(parts) >= 2:
                    procs.add((os.path.basename(parts[0]), parts[1]))
            elif line.startswith('ListFindMiss: '):
                absent.add(line[len('ListFindMiss: '):].strip().lower())
            elif line.strip() == 'Unused keywords:':
                in_unused = True
            elif in_unused:
                if line.startswith('       ') and len(line) > 29:
                    sec = line[7:29].strip()
                    kind = section_kind(sec)
                    unused.add((kind or sec, line[29:].strip().lower()))
                else:
                    in_unused = False
    # The sweep itself raised Max Output Level; never report it.
    if 'max output level' not in original_text.lower():
        unused.discard(('Simulation', 'max output level'))
    used = sorted(p for p in present if p not in unused)
    injected = sorted({lib for lib, _ in procs if lib.lower() not in declared})
    tests[name] = {
        'dir': d,
        'np': int(np_),
        'rc': rc,
        'classes': classes.get(name, []),
        'procedures': sorted([lib, fn] for lib, fn in procs),
        'injected_libraries': injected,
        'keywords_used': [list(p) for p in used],
        'keywords_unused': sorted(list(p) for p in unused),
        'queried_absent': sorted(absent),
    }
    for p in present:
        kw_present[p] += 1
    for p in unused:
        kw_unused[p] += 1
    for k in absent:
        kw_absent[k] += 1
    for lib, _ in procs:
        mod_tests[lib].add(name)
    for lib in injected:
        mod_injected[lib].add(name)

os.makedirs(out, exist_ok=True)
# Queried-absent names are shared across tests: store each test's list as
# ids into queried-absent.tsv (id = row order, most common first).
absent_order = sorted(kw_absent.items(), key=lambda x: (-x[1], x[0]))
absent_id = {k: i for i, (k, _) in enumerate(absent_order)}
with open(os.path.join(out, 'quick.json'), 'w') as f:
    f.write('{"version": 1, "arm": "listusage", "queried_absent_ids": "row ids of queried-absent.tsv", "tests": {\n')
    rows = []
    for name in sorted(tests):
        t = dict(tests[name])
        t['queried_absent'] = sorted(absent_id[k] for k in t['queried_absent'])
        rows.append('%s: %s' % (json.dumps(name), json.dumps(t, sort_keys=True, separators=(',', ':'))))
    f.write(',\n'.join(rows))
    f.write('\n}}\n')
with open(os.path.join(out, 'keywords.tsv'), 'w') as f:
    f.write('section\tkeyword\ttests_present\ttests_unused\n')
    for (sec, k), n in sorted(kw_present.items(), key=lambda x: (-x[1], x[0])):
        f.write('%s\t%s\t%d\t%d\n' % (sec, k, n, kw_unused.get((sec, k), 0)))
with open(os.path.join(out, 'queried-absent.tsv'), 'w') as f:
    f.write('id\tkeyword\ttests_querying_without_setting\n')
    for i, (k, n) in enumerate(absent_order):
        f.write('%d\t%s\t%d\n' % (i, k, n))
with open(os.path.join(out, 'modules.tsv'), 'w') as f:
    f.write('library\ttests\ttests_injected\n')
    for lib in sorted(mod_tests, key=lambda l: (-len(mod_tests[l]), l)):
        f.write('%s\t%d\t%d\n' % (lib, len(mod_tests[lib]), len(mod_injected.get(lib, ()))))
print('tests', len(tests), 'rc', collections.Counter(t['rc'] for t in tests.values()))
print('keywords present', len(kw_present), 'queried-absent', len(kw_absent), 'libraries', len(mod_tests))
print('tests with injected libraries', sum(1 for t in tests.values() if t['injected_libraries']))
