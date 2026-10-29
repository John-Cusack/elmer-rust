import struct, sys, re
vtu, out = sys.argv[1], sys.argv[2]
raw = open(vtu, 'rb').read()
head = raw[:raw.index(b'<AppendedData')].decode()
npts = int(re.search(r'NumberOfPoints="(\d+)"', head).group(1))
i = raw.index(b'<AppendedData')
i = raw.index(b'_', i) + 1
n = struct.unpack('<I', raw[i:i+4])[0]
vals = struct.unpack(f'<{n//8}d', raw[i+4:i+4+n])
assert len(vals) == npts
open(out, 'w').write('\n'.join(repr(v) for v in vals) + '\n')
print(npts, min(vals), max(vals), sum(1 for v in vals if v == -1000.0))
