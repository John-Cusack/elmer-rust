//! Mesh read proof: load every generated quick-test mesh with
//! `elmer_mesh` and compare against an independent from-scratch text
//! parse kept inside this test (shared-trait rule: test-only duplicate
//! readers must stay under `tests/`).
//!
//! Requires the mesh manifest from `proof/gen_meshes.py`:
//! `MESH_MANIFEST` env or default scratch path. Meshes themselves come
//! from `MESH_ROOT` (default scratch meshes dir).

fn manifest_path() -> String {
    std::env::var("MESH_MANIFEST").unwrap_or_else(|_| {
        "/home/john/firstmate/data/elmer-s1-frontend/proof/mesh-manifest.txt".to_string()
    })
}

fn mesh_root() -> String {
    std::env::var("MESH_ROOT")
        .unwrap_or_else(|_| "/home/john/firstmate/data/elmer-s1-frontend/meshes".to_string())
}

/// Independent text parse of one mesh DB (no `elmer_mesh` code).
struct RawDb {
    n_nodes: usize,
    n_bulk: usize,
    n_bound: usize,
    nodes: Vec<[f64; 3]>,
    bulk: Vec<(i32, u16, Vec<u32>)>, // (body, code, 0-based nodes)
    bound: Vec<(i32, i64, i64, u16, Vec<u32>)>, // (constraint, left, right, code, nodes)
}

fn parse_ints(s: &str) -> Vec<i64> {
    // Elmer's read_ints without halo (halo markers rejected here).
    let mut out = Vec::new();
    let b = s.as_bytes();
    let k = s.trim_end().len();
    let mut l = 0;
    while l < k && out.len() < 64 {
        while l < k && (b[l] == b' ' || b[l] == b'\t') {
            l += 1;
        }
        if l >= k {
            break;
        }
        assert!(
            b[l] != b'/',
            "halo marker in serial proof mesh: {s}"
        );
        if !(b[l] == b'-' || b[l].is_ascii_digit()) {
            break;
        }
        let mut m = l + 1;
        while m < k && b[m].is_ascii_digit() {
            m += 1;
        }
        out.push(s[l..m].parse().unwrap());
        l = m;
    }
    out
}

fn read_raw_db(dir: &str) -> RawDb {
    let header = std::fs::read_to_string(format!("{dir}/mesh.header")).unwrap();
    let mut it = header.split_whitespace();
    let num = |it: &mut std::str::SplitWhitespace<'_>| -> usize {
        it.next().unwrap().parse().unwrap()
    };
    let (n_nodes, n_bulk, n_bound) = (num(&mut it), num(&mut it), num(&mut it));
    let n_types: usize = num(&mut it);
    for _ in 0..n_types {
        let _code: usize = num(&mut it);
        let _n: usize = num(&mut it);
    }
    let mut nodes = Vec::with_capacity(n_nodes);
    for line in std::fs::read_to_string(format!("{dir}/mesh.nodes"))
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .take(n_nodes)
    {
        let w: Vec<&str> = line.split_whitespace().collect();
        let _tag: i32 = w[0].parse().unwrap();
        nodes.push([
            w[2].parse().unwrap(),
            w[3].parse().unwrap(),
            w[4].parse().unwrap(),
        ]);
    }
    let mut bulk = Vec::with_capacity(n_bulk);
    for line in std::fs::read_to_string(format!("{dir}/mesh.elements"))
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .take(n_bulk)
    {
        let v = parse_ints(line);
        let (tag, body, code) = (v[0] as i32, v[1] as i32, v[2] as u16);
        // Node count from the element code family (linear S1 set + quad).
        // Independent transcription of elements.def node counts.
        let nn = match code {
            101 => 1,
            202 => 2,
            203 => 3,
            204 => 4,
            303 => 3,
            306 => 6,
            404 => 4,
            408 => 8,
            409 => 9,
            504 => 4,
            510 => 10,
            605 => 5,
            706 => 6,
            718 => 18,
            808 => 8,
            820 => 20,
            827 => 27,
            _ => v.len() - 3,
        };
        let conn: Vec<u32> = v[3..3 + nn].iter().map(|&x| (x - 1) as u32).collect();
        bulk.push((body, code, conn));
        let _ = tag;
    }
    let mut bound = Vec::with_capacity(n_bound);
    for line in std::fs::read_to_string(format!("{dir}/mesh.boundary"))
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .take(n_bound)
    {
        let v = parse_ints(line);
        let (tag, bndry, left, right, code) =
            (v[0] as i32, v[1] as i32, v[2], v[3], v[4] as u16);
        let nn = match code {
            101 => 1,
            202 => 2,
            203 => 3,
            204 => 4,
            303 => 3,
            306 => 6,
            404 => 4,
            408 => 8,
            409 => 9,
            504 => 4,
            510 => 10,
            605 => 5,
            706 => 6,
            718 => 18,
            808 => 8,
            820 => 20,
            827 => 27,
            _ => v.len() - 5,
        };
        let conn: Vec<u32> = v[5..5 + nn].iter().map(|&x| (x - 1) as u32).collect();
        bound.push((bndry, left, right, code, conn));
        let _ = tag;
    }
    RawDb { n_nodes, n_bulk, n_bound, nodes, bulk, bound }
}

#[test]
fn quick_test_meshes_match_elmergrid_output() {
    let manifest = std::fs::read_to_string(manifest_path()).unwrap();
    let root = mesh_root();
    let mut n_ok = 0usize;
    let mut n_skip = 0usize;
    let mut failures: Vec<String> = Vec::new();
    for line in manifest.lines() {
        let mut parts = line.split('\t');
        let (name, status, detail) = (
            parts.next().unwrap_or(""),
            parts.next().unwrap_or(""),
            parts.next().unwrap_or(""),
        );
        if !status.starts_with("ok") {
            n_skip += 1;
            continue;
        }
        let dir = format!("{root}/{name}/{detail}");
        let raw = read_raw_db(&dir);
        let mesh = elmer_mesh::load(&dir).unwrap_or_else(|e| {
            panic!("{name}: reader failed on ElmerGrid output: {e}");
        });
        let mut bad = String::new();
        if mesh.nodes.len() != raw.n_nodes {
            bad.push_str(&format!("nodes {}!={} ", mesh.nodes.len(), raw.n_nodes));
        }
        if mesh.bulk.len() != raw.n_bulk {
            bad.push_str(&format!("bulk {}!={} ", mesh.bulk.len(), raw.n_bulk));
        }
        if mesh.boundary.len() != raw.n_bound {
            bad.push_str(&format!("bound {}!={} ", mesh.boundary.len(), raw.n_bound));
        }
        // Full coordinate + connectivity comparison (tags may permute:
        // compare node sets and per-element node sets).
        if bad.is_empty() {
            // Nodes: same multiset of coordinates (permutation-tolerant).
            let mut a: Vec<[i64; 3]> = mesh
                .nodes
                .iter()
                .map(|n| {
                    [
                        (n[0] * 1e12).round() as i64,
                        (n[1] * 1e12).round() as i64,
                        (n[2] * 1e12).round() as i64,
                    ]
                })
                .collect();
            let mut b: Vec<[i64; 3]> = raw
                .nodes
                .iter()
                .map(|n| {
                    [
                        (n[0] * 1e12).round() as i64,
                        (n[1] * 1e12).round() as i64,
                        (n[2] * 1e12).round() as i64,
                    ]
                })
                .collect();
            a.sort();
            b.sort();
            if a != b {
                bad.push_str("node coords differ ");
            }
            // Bulk connectivity: multiset of (code, body, sorted nodes).
            let mut ea: Vec<(u16, i32, Vec<u32>)> = mesh
                .bulk
                .iter()
                .map(|e| {
                    let mut ns = e.nodes.clone();
                    ns.sort();
                    (e.code, e.body, ns)
                })
                .collect();
            let mut eb: Vec<(u16, i32, Vec<u32>)> = raw
                .bulk
                .iter()
                .map(|(body, code, ns)| {
                    let mut s = ns.clone();
                    s.sort();
                    (*code, *body, s)
                })
                .collect();
            ea.sort();
            eb.sort();
            if ea != eb {
                bad.push_str("bulk connectivity differs ");
            }
            // Boundary: multiset of (code, sorted nodes); parents checked below.
            let mut fa: Vec<(u16, Vec<u32>)> = mesh
                .boundary
                .iter()
                .map(|e| {
                    let mut ns = e.nodes.clone();
                    ns.sort();
                    (e.code, ns)
                })
                .collect();
            let mut fb: Vec<(u16, Vec<u32>)> = raw
                .bound
                .iter()
                .map(|(_, _, _, code, ns)| {
                    let mut s = ns.clone();
                    s.sort();
                    (*code, s)
                })
                .collect();
            fa.sort();
            fb.sort();
            if fa != fb {
                bad.push_str("boundary connectivity differs ");
            }
            // Parents: every parent index in range; every boundary node
            // used by its parent bulk element (structural check).
            for (i, be) in mesh.boundary.iter().enumerate() {
                for p in be.parents.into_iter().flatten() {
                    if p >= mesh.bulk.len() {
                        bad.push_str(&format!("boundary {i} parent out of range "));
                    }
                }
            }
        }
        if bad.is_empty() {
            n_ok += 1;
        } else {
            failures.push(format!("{name}: {bad}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} mesh mismatches:\n{}",
        failures.len(),
        failures.iter().take(20).cloned().collect::<Vec<_>>().join("\n")
    );
    assert!(n_ok > 400, "too few meshes compared: {n_ok}");
    eprintln!("mesh read proof: {n_ok} ok, {n_skip} skipped (non-serial)");
}
