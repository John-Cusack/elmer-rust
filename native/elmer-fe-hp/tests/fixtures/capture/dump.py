import gdb, json, re, os, struct
OUT = open(os.environ.get("DUMP_OUT", "/tmp/dump.jsonl"), "w")
INF = None

def flat(v):
    s = str(v).replace("<optimized out>", "nan")
    return [[float(x) for x in g.split(",") if x.strip()] for g in re.findall(r"\(([^()]*)\)", s)]

def rd(addr, n):
    return bytes(INF.read_memory(addr, n))

def read_desc(dptr):
    """gfortran array descriptor -> list of columns (one list per last-index value)."""
    if dptr == 0:
        return None
    base = struct.unpack("<Q", rd(dptr, 8))[0]
    rank = rd(dptr + 28, 1)[0]
    dims = []
    for r in range(rank):
        st, lb, ub = struct.unpack("<qqq", rd(dptr + 40 + 24 * r, 24))
        dims.append((st, lb, ub))
    if rank == 1:
        st, lb, ub = dims[0]
        n = ub - lb + 1
        return [[struct.unpack("<d", rd(base + 8 * st * i, 8))[0] for i in range(n)]]
    st0, lb0, ub0 = dims[0]
    st1, lb1, ub1 = dims[1]
    cols = []
    for j in range(ub1 - lb1 + 1):
        cols.append([struct.unpack("<d", rd(base + 8 * (st0 * i + st1 * j), 8))[0] for i in range(ub0 - lb0 + 1)])
    return cols

gdb.execute("set pagination off"); gdb.execute("set print elements 0"); gdb.execute("set print repeats 0")
gdb.execute("set print frame-arguments none")
gdb.execute("set breakpoint pending on")
BP = gdb.Breakpoint("__eleminfo_MOD_elementinfo")
gdb.execute("run case.sif")
INF = gdb.selected_inferior()
while True:
    try:
        f = gdb.selected_frame()
    except gdb.error:
        break
    if f.name() is None or "elementinfo" not in f.name().lower():
        break
    try:
        caller = f.older().name() or ""
        if "localmatrix" in caller.lower():
            cfa = int(f.read_register("rbp")) + 16
            if os.environ.get("NOEDGE"):
                INF.write_memory(cfa + 48, (0).to_bytes(8, "little"))
                INF.write_memory(cfa + 56, (0).to_bytes(8, "little"))
            ptrs = struct.unpack("<10Q", rd(cfa, 80))
            el = f.read_var("element")
            try:
                pd = el["pdefs"].dereference()
                pinfo = {"tetratype": int(pd["tetratype"]), "p": int(pd["p"]), "sd": bool(pd["serendipity"])}
            except Exception as e:
                pinfo = {"err": str(e)}
            rec = {"pinfo": pinfo, "code": int(el["type"].dereference()["elementcode"]),
                   "nodeindexes": [int(x) for x in re.findall(r"-?\d+", str(el["nodeindexes"]).replace("(","").replace(")",""))],
                   "xyz": flat(f.read_var("nodes")["xyz"]),
                   "u": float(f.read_var("u")), "v": float(f.read_var("v")), "w": float(f.read_var("w"))}
            # detJ is the 6th (register) argument: read after return from the saved address
            detj_addr = int(f.read_register("r9")) if False else None
            detj = f.read_var("detj")
            detj_a = int(detj.address)
            BP.enabled = False   # recursive ElementInfo calls must not interrupt the finish
            gdb.execute("finish", to_string=True)
            BP.enabled = True
            rec["detj"] = struct.unpack("<d", rd(detj_a, 8))[0]
            for name, k in (("basis", 0), ("dbasisdx", 1), ("edgebasis", 6), ("rotbasis", 7)):
                rec[name] = read_desc(ptrs[k])
            OUT.write(json.dumps(rec) + "\n"); OUT.flush()
    except Exception as e:
        import traceback; OUT.write(json.dumps({"error": traceback.format_exc()}) + "\n"); OUT.flush()
    try:
        gdb.execute("continue", to_string=True)
    except gdb.error:
        break
    try:
        gdb.selected_frame()
    except gdb.error:
        break
gdb.execute("quit")
