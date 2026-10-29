import struct, sys, re
vtu = sys.argv[1]
raw = open(vtu, 'rb').read()
head = raw[:raw.index(b'<AppendedData')].decode()
body = raw[raw.index(b'_', raw.index(b'<AppendedData')) + 1:]
arrs = {}
for m in re.finditer(r'<DataArray type="(\w+)"(?: Name="(\w+)")?(?: NumberOfComponents="(\d+)")? format="appended" offset="(\d+)"', head):
    typ, name, nc, off = m.group(1), m.group(2) or 'points', m.group(3), int(m.group(4))
    n = struct.unpack('<I', body[off:off+4])[0]
    fmt = {'Float64': 'd', 'Int32': 'i'}[typ]
    sz = struct.calcsize(fmt)
    arrs[name] = struct.unpack(f'<{n//sz}{fmt}', body[off+4:off+4+n])
for k, v in arrs.items(): print(k, len(v))
if len(sys.argv) > 2:
    open(sys.argv[2] + '_values.txt', 'w').write('\n'.join(repr(x) for x in arrs['tempa']) + '\n')
    open(sys.argv[2] + '_points.txt', 'w').write('\n'.join(repr(x) for x in arrs['points']) + '\n')
