import os, math, random
random.seed(7)
out = 'case2d'

def write_mesh(d, coords, elems, body):
    os.makedirs(d, exist_ok=True)
    types = {}
    for t, _ in elems: types[t] = types.get(t, 0) + 1
    with open(f'{d}/mesh.header', 'w') as f:
        f.write(f'{len(coords)} {len(elems)} 0\n{len(types)}\n')
        for t, c in sorted(types.items()): f.write(f'{t} {c}\n')
    with open(f'{d}/mesh.nodes', 'w') as f:
        for i, (x, y) in enumerate(coords, 1): f.write(f'{i} -1 {x!r} {y!r} 0.0\n')
    with open(f'{d}/mesh.elements', 'w') as f:
        for i, (t, nodes) in enumerate(elems, 1):
            f.write(f'{i} {body} {t} ' + ' '.join(str(n + 1) for n in nodes) + '\n')
    open(f'{d}/mesh.boundary', 'w').close()

# src: 7 x 6 nodes on the unit square, jittered, quads on the left, triangles right
nx, ny = 6, 5
src = []
for j in range(ny + 1):
    for i in range(nx + 1):
        x, y = i / nx, j / ny
        if 0 < i < nx and 0 < j < ny:
            x += random.uniform(-0.06, 0.06); y += random.uniform(-0.06, 0.06)
        src.append((x, y))
id_ = lambda i, j: j * (nx + 1) + i
se = []
for j in range(ny):
    for i in range(nx):
        a, b, c, d = id_(i, j), id_(i + 1, j), id_(i + 1, j + 1), id_(i, j + 1)
        if i < 3: se.append((404, [a, b, c, d]))
        else:
            se.append((303, [a, b, c])); se.append((303, [a, c, d]))
write_mesh(f'{out}/src', src, se, 1)

# dst: 10 x 9 nodes on [0.15,1.25] x [-0.05,0.95], triangles, nodes poke outside src
mx, my = 9, 8
dst = []
for j in range(my + 1):
    for i in range(mx + 1):
        x = 0.15 + 1.10 * i / mx; y = -0.05 + 1.0 * j / my
        if 0 < i < mx and 0 < j < my:
            x += random.uniform(-0.02, 0.02); y += random.uniform(-0.02, 0.02)
        dst.append((x, y))
idd = lambda i, j: j * (mx + 1) + i
de = []
for j in range(my):
    for i in range(mx):
        a, b, c, d = idd(i, j), idd(i + 1, j), idd(i + 1, j + 1), idd(i, j + 1)
        de.append((303, [a, b, c])); de.append((303, [a, c, d]))
write_mesh(f'{out}/dst', dst, de, 2)
