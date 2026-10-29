import os, random
random.seed(11)
out = 'case3d'

def write_mesh(d, coords, elems, body):
    os.makedirs(d, exist_ok=True)
    types = {}
    for t, _ in elems: types[t] = types.get(t, 0) + 1
    with open(f'{d}/mesh.header', 'w') as f:
        f.write(f'{len(coords)} {len(elems)} 0\n{len(types)}\n')
        for t, c in sorted(types.items()): f.write(f'{t} {c}\n')
    with open(f'{d}/mesh.nodes', 'w') as f:
        for i, (x, y, z) in enumerate(coords, 1): f.write(f'{i} -1 {x!r} {y!r} {z!r}\n')
    with open(f'{d}/mesh.elements', 'w') as f:
        for i, (t, nodes) in enumerate(elems, 1):
            f.write(f'{i} {body} {t} ' + ' '.join(str(n + 1) for n in nodes) + '\n')
    open(f'{d}/mesh.boundary', 'w').close()

def grid(lo, hi, n, jitter):
    c = []
    for k in range(n[2] + 1):
        for j in range(n[1] + 1):
            for i in range(n[0] + 1):
                p = [lo[d] + (hi[d] - lo[d]) * idx / n[d] for d, idx in enumerate((i, j, k))]
                if 0 < i < n[0] and 0 < j < n[1] and 0 < k < n[2]:
                    p = [x + random.uniform(-jitter, jitter) for x in p]
                c.append(tuple(p))
    nid = lambda i, j, k: (k * (n[1] + 1) + j) * (n[0] + 1) + i
    cells = []
    for k in range(n[2]):
        for j in range(n[1]):
            for i in range(n[0]):
                cells.append((i, [nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                                  nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1)]))
    return c, cells

# src: hexes for the left half, wedges on the right
c, cells = grid((0, 0, 0), (1, 1, 1), (4, 3, 3), 0.05)
el = []
for i, h in cells:
    if i < 2: el.append((808, h))
    else:
        a, b, cc, d, e, f, g, hh = h
        el.append((706, [a, b, cc, e, f, g])); el.append((706, [a, cc, d, e, g, hh]))
write_mesh(f'{out}/src', c, el, 1)

# dst: tets from a shifted/smaller grid, partly outside src
c, cells = grid((0.1, 0.1, 0.05), (1.15, 0.95, 1.0), (5, 4, 4), 0.03)
el = []
for _, h in cells:
    for t in ([0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]):
        el.append((504, [h[q] for q in t]))
write_mesh(f'{out}/dst', c, el, 2)
