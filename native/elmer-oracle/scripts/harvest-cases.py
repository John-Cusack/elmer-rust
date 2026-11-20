#!/usr/bin/env python3
"""Harvest the ABI campaign's staged .sif differentials (report §3.6, §5.2 item 11).

Each campaign proof (`crates/elmer-fem-b/ports/<solver>/proof-<solver>.sh`
and the ElmerIce USF proofs) materialised its cases in a lane scratch dir:
`<lane>/scratch/proof-<solver>/case*.sif` next to the mesh source (an
ElmerGrid `.grd`, or a mesh directory copied from an upstream test). This
script reads those dirs (never writes them) and turns each case into an
oracle test directory that `elmer-oracle capture` plans like a registered
test:

    <out>/<id>/CMakeLists.txt  ADD_ELMER_TEST(<id> LABELS harvest)
    <out>/<id>/runtest.cmake   ElmerGrid 1 2 <mesh> (when a .grd), RUN_ELMER_TEST()
    <out>/<id>/ELMERSOLVER_STARTINFO
    <out>/<id>/<case>.sif      verbatim
    <out>/<id>/<mesh>.grd      verbatim, or
    <out>/<id>/<mesh>  ->  <out>/_meshes/<mesh>-<sha12>/  (mesh dirs, stored once)
    <out>/<id>/<other>         files the sif names that sit beside it

It also writes <out>/cases.tsv: id, lane, proof dir, sif, sha256 of the
sif, mesh source, procedures, and why a case was skipped.

usage: harvest-cases.py <lane-data-root> <out-dir>
"""

import hashlib
import os
import re
import shutil
import sys

LANES = ["rs-elmer-t020", "rs-elmer-t021a", "rs-elmer-t021b", "rs-elmer-t022a", "rs-elmer-t023", "rs-elmer-t025"]
MESH_FILES = ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary")


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def code_lines(text):
    """The sif without `!` comments (quoted strings kept)."""
    out = []
    for line in text.splitlines():
        q = False
        for i, c in enumerate(line):
            if c == '"':
                q = not q
            elif c == "!" and not q:
                line = line[:i]
                break
        out.append(line)
    return out


def referenced_files(text, here):
    """Quoted names in the sif that are files or directories beside it
    (data tables, restart inputs, OpenFOAM case dirs), excluding the mesh
    and procedure libraries. A file's `.pos`/`.names` companions (binary
    positions, column names) come with it."""
    names = set()
    for line in code_lines(text):
        low = line.lower()
        if "procedure" in low or "mesh db" in low:
            continue
        cands = re.findall(r'"([^"]+)"', line)
        m = re.search(r"=\s*(?:file\s+)?([A-Za-z0-9_./-]+\.(?:dat|txt|csv|xyz|ep|result))\s*$", line, re.I)
        if m:
            cands.append(m.group(1))
        for q in cands:
            if q.startswith("/") or ".." in q:
                continue
            p = os.path.join(here, q)
            if os.path.isfile(p) or os.path.isdir(p):
                names.add(q)
                for ext in (".pos", ".names"):
                    if os.path.isfile(p + ext):
                        names.add(q + ext)
    return sorted(names)


def main():
    root, out = sys.argv[1], sys.argv[2]
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)
    rows = []
    seen = {}
    for lane in LANES:
        scratch = os.path.join(root, lane, "scratch")
        if not os.path.isdir(scratch):
            continue
        for d in sorted(os.listdir(scratch)):
            here = os.path.join(scratch, d)
            if not d.startswith("proof-") or not os.path.isdir(here):
                continue
            solver = d[len("proof-"):]
            for sif in sorted(f for f in os.listdir(here) if re.fullmatch(r"case[A-Za-z0-9_]*\.sif", f)):
                case = sif[:-4]
                path = os.path.join(here, sif)
                text = open(path, errors="replace").read()
                digest = sha256(path)
                base = f"{solver}-{case}"
                # One copy per distinct sif: a later lane re-staging the
                # same bytes adds nothing.
                if (base, digest) in seen:
                    continue
                cid = base if base not in {k[0] for k in seen} else f"{lane[len('rs-elmer-'):]}-{base}"
                seen[(base, digest)] = cid
                procs = sorted({p.lower() for p in re.findall(r'Procedure\s*=\s*(?:File\s+)?"([^"]+)"', text, re.I)})
                m = re.search(r'Mesh\s+DB\s+"([^"]*)"\s+"([^"]*)"', text, re.I)
                row = {"id": cid, "lane": lane, "proof": d, "sif": sif, "sha256": digest,
                       "mesh": "", "procedures": ",".join(procs), "skip": ""}
                rows.append(row)
                if not m:
                    row["skip"] = "no Mesh DB"
                    continue
                mdir, mesh = m.group(1), m.group(2)
                if mdir not in (".", "./"):
                    row["skip"] = f"Mesh DB directory {mdir!r}"
                    continue
                grd = os.path.join(here, mesh + ".grd")
                mpath = os.path.join(here, mesh)
                dst = os.path.join(out, cid)
                if os.path.exists(dst):
                    shutil.rmtree(dst)
                os.makedirs(dst)
                if os.path.isfile(grd):
                    shutil.copy2(grd, dst)
                    row["mesh"] = f"{mesh}.grd"
                    grid = f"execute_process(COMMAND ${{ELMERGRID_BIN}} 1 2 {mesh})\n"
                elif all(os.path.isfile(os.path.join(mpath, f)) for f in MESH_FILES):
                    # Mesh directories (mostly upstream test meshes) are
                    # stored once, by content, and linked; the runner's
                    # staging follows the link.
                    h = hashlib.sha256()
                    for top, dirs, files in sorted(os.walk(mpath)):
                        dirs.sort()
                        for f in sorted(files):
                            p = os.path.join(top, f)
                            with open(p, "rb") as fh:
                                h.update(os.path.relpath(p, mpath).encode() + b"\0" + fh.read())
                    shared = os.path.join(out, "_meshes", f"{mesh}-{h.hexdigest()[:12]}")
                    if not os.path.isdir(shared):
                        shutil.copytree(mpath, shared)
                    os.symlink(os.path.relpath(shared, dst), os.path.join(dst, mesh))
                    row["mesh"] = f"{mesh}/"
                    grid = ""
                else:
                    shutil.rmtree(dst)
                    row["skip"] = f"no mesh source for {mesh!r}"
                    continue
                shutil.copy2(path, dst)
                extra = referenced_files(text, here)
                for f in extra:
                    os.makedirs(os.path.dirname(os.path.join(dst, f)) or dst, exist_ok=True)
                    if os.path.isdir(os.path.join(here, f)):
                        shutil.copytree(os.path.join(here, f), os.path.join(dst, f), dirs_exist_ok=True)
                    else:
                        shutil.copy2(os.path.join(here, f), os.path.join(dst, f))
                copies = " ".join(["ELMERSOLVER_STARTINFO", row["mesh"].rstrip("/")] + extra)
                with open(os.path.join(dst, "CMakeLists.txt"), "w") as f:
                    f.write("INCLUDE(test_macros)\nINCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)\n\n"
                            f"CONFIGURE_FILE({sif} {sif} COPYONLY)\n\n"
                            f'file(COPY {copies} DESTINATION "${{CMAKE_CURRENT_BINARY_DIR}}/")\n\n'
                            f"ADD_ELMER_TEST({cid} LABELS harvest)\n")
                with open(os.path.join(dst, "runtest.cmake"), "w") as f:
                    f.write("include(test_macros)\n" + grid + "RUN_ELMER_TEST()\n")
                with open(os.path.join(dst, "ELMERSOLVER_STARTINFO"), "w") as f:
                    f.write(sif + "\n")
    cols = ["id", "lane", "proof", "sif", "sha256", "mesh", "procedures", "skip"]
    with open(os.path.join(out, "cases.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in cols) + "\n")
    kept = sum(1 for r in rows if not r["skip"])
    print(f"{len(rows)} distinct staged cases, {kept} harvested, {len(rows) - kept} skipped", file=sys.stderr)


if __name__ == "__main__":
    main()
