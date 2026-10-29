//! `MeshParser`: the Mesh2D input (`.mif`) reader, with the tokenizer
//! behaviour of the C++ stream (`ignoreComment`, `readString` quoting).

pub struct NodeToken {
    pub tag: i32,
    pub boundary_tag: i32,
    pub x: f64,
    pub y: f64,
    pub delta: f64,
}

pub struct EdgeToken {
    pub tag: i32,
    pub boundary_tag: i32,
    pub nodes: Vec<i32>,
    pub seg_count: i32,
    pub delta: f64,
}

pub struct LoopToken {
    pub direction: i32,
    pub edges: Vec<i32>,
}

#[derive(Default)]
pub struct BgMeshToken {
    pub kind: String,
    /// `(x, y, delta)`.
    pub nodes: Vec<(f64, f64, f64)>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum LayerKind {
    VoronoiVertex,
    /// `MovingFront`.
    VoronoiSegment,
}

pub struct LayerToken {
    pub tag: i32,
    pub kind: LayerKind,
    pub loops: Vec<LoopToken>,
    pub fixed_nodes: Vec<NodeToken>,
    pub bg: BgMeshToken,
    pub delta: f64,
}

pub struct BodyToken {
    pub tag: i32,
    pub layers: Vec<LayerToken>,
    pub delta: f64,
    pub parabolic: bool,
}

pub struct Parsed {
    pub delta: f64,
    pub scale: f64,
    pub nodes: Vec<NodeToken>,
    pub edges: Vec<EdgeToken>,
    pub bodies: Vec<BodyToken>,
}

struct Reader<'a> {
    b: &'a [u8],
    p: usize,
}

type R<T> = Result<T, String>;

impl<'a> Reader<'a> {
    fn eof(&self) -> bool {
        self.p >= self.b.len()
    }

    fn peek(&self) -> Option<u8> {
        self.b.get(self.p).copied()
    }

    /// `while (ignoreComment());`
    fn skip(&mut self) {
        loop {
            while self.peek().is_some_and(|c| c.is_ascii_whitespace()) {
                self.p += 1;
            }
            if matches!(self.peek(), Some(b'!' | b'#')) {
                while let Some(c) = self.peek() {
                    self.p += 1;
                    if c == b'\n' {
                        break;
                    }
                }
            } else {
                return;
            }
        }
    }

    fn token(&mut self) -> &'a str {
        let s = self.p;
        while self.peek().is_some_and(|c| !c.is_ascii_whitespace()) {
            self.p += 1;
        }
        std::str::from_utf8(&self.b[s..self.p]).unwrap_or("")
    }

    fn word(&mut self, want: &str) -> R<()> {
        self.skip();
        let w = self.token();
        if w == want {
            Ok(())
        } else {
            Err(format!("expected {want}, found {w:?}"))
        }
    }

    fn int(&mut self) -> R<i32> {
        self.skip();
        let s = self.p;
        if matches!(self.peek(), Some(b'+' | b'-')) {
            self.p += 1;
        }
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            self.p += 1;
        }
        std::str::from_utf8(&self.b[s..self.p])
            .unwrap_or("")
            .parse()
            .map_err(|_| "expected an integer".to_string())
    }

    fn double(&mut self) -> R<f64> {
        self.skip();
        let s = self.p;
        if matches!(self.peek(), Some(b'+' | b'-')) {
            self.p += 1;
        }
        while self.peek().is_some_and(|c| c.is_ascii_digit() || c == b'.') {
            self.p += 1;
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            let save = self.p;
            self.p += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.p += 1;
            }
            if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    self.p += 1;
                }
            } else {
                self.p = save;
            }
        }
        std::str::from_utf8(&self.b[s..self.p])
            .unwrap_or("")
            .parse()
            .map_err(|_| "expected a number".to_string())
    }

    /// `readString`: a word with `"` quoting and `!`/`#` comments.
    fn string(&mut self) -> R<String> {
        let mut out = String::new();
        let (mut head, mut inquote, mut incomment) = (true, false, false);
        while !self.eof() {
            let ch = self.b[self.p];
            self.p += 1;
            if incomment {
                if ch != b'\n' {
                    continue;
                }
                incomment = false;
            }
            if ch == b'"' {
                inquote = !inquote;
            } else if !inquote && ch.is_ascii_whitespace() {
                if head {
                    continue;
                }
                break;
            } else if !inquote && (ch == b'!' || ch == b'#') {
                incomment = true;
            } else {
                head = false;
                out.push(ch as char);
            }
        }
        if out.is_empty() && self.eof() {
            return Err("unexpected end of file".into());
        }
        Ok(out)
    }

    /// The optional `H:` / `R:` size, `-1` if absent.
    fn size(&mut self, delta: f64) -> R<f64> {
        self.skip();
        match self.peek() {
            Some(b'H') => {
                self.word("H:")?;
                self.double()
            }
            Some(b'R') => {
                self.word("R:")?;
                Ok(self.double()? * delta)
            }
            _ => Ok(-1.0),
        }
    }
}

fn node(r: &mut Reader, word: &str, delta: f64) -> R<NodeToken> {
    r.word(word)?;
    let tag = r.int()?;
    let boundary_tag = r.int()?;
    let d = r.size(delta)?;
    let x = r.double()?;
    let y = r.double()?;
    Ok(NodeToken {
        tag,
        boundary_tag,
        x,
        y,
        delta: d,
    })
}

fn explicit_bg(r: &mut Reader, bg: &mut BgMeshToken) -> R<()> {
    let len = r.int()?;
    for _ in 0..len {
        let x = r.double()?;
        let y = r.double()?;
        let d = r.double()?;
        bg.nodes.push((x, y, d));
    }
    Ok(())
}

/// `MeshParser::readInputFile`. `external_bg` is the `--bgmesh=` mesh, used
/// by every layer that does not name its own.
pub fn parse(text: &str, external_bg: Option<&BgMeshToken>) -> R<Parsed> {
    let mut r = Reader {
        b: text.as_bytes(),
        p: 0,
    };
    r.word("Geometry2D:")
        .map_err(|e| format!("Error reading header: {e}"))?;
    r.word("H:")?;
    let delta = r.double()?;
    r.word("MeshScalingFactor:")?;
    let scale = r.double()?;
    r.word("Nodes:")?;
    let node_count = r.int()?;
    r.word("Edges:")?;
    let edge_count = r.int()?;
    r.word("Bodies:")?;
    let body_count = r.int()?;

    let mut nodes = Vec::new();
    for _ in 0..node_count {
        nodes
            .push(node(&mut r, "NodeId:", delta).map_err(|e| format!("Error reading nodes: {e}"))?);
    }

    let mut edges = Vec::new();
    for _ in 0..edge_count {
        r.word("EdgeId:")
            .map_err(|e| format!("Error reading edges: {e}"))?;
        let tag = r.int()?;
        let boundary_tag = r.int()?;
        r.skip();
        let (d, seg_count) = match r.peek() {
            Some(b'H') => {
                r.word("H:")?;
                (r.double()?, 0)
            }
            Some(b'R') => {
                r.word("R:")?;
                (r.double()? * delta, 0)
            }
            Some(b'N') => {
                r.word("N:")?;
                (-1.0, r.int()?)
            }
            _ => (-1.0, 0),
        };
        let len = r.int()?;
        let mut ns = Vec::new();
        for _ in 0..len {
            ns.push(r.int()?);
        }
        edges.push(EdgeToken {
            tag,
            boundary_tag,
            nodes: ns,
            seg_count,
            delta: d,
        });
    }

    let mut bodies = Vec::new();
    for _ in 0..body_count {
        bodies.push(
            body(&mut r, delta, external_bg).map_err(|e| format!("Error reading bodies: {e}"))?,
        );
    }
    Ok(Parsed {
        delta,
        scale,
        nodes,
        edges,
        bodies,
    })
}

fn body(r: &mut Reader, delta: f64, external_bg: Option<&BgMeshToken>) -> R<BodyToken> {
    r.word("BodyId:")?;
    let tag = r.int()?;
    let bdelta = r.size(delta)?;
    r.word("ElementOrder:")?;
    let order = r.string()?;
    let parabolic = match order.as_str() {
        "Parabolic" => true,
        "Linear" => false,
        o => return Err(format!("unknown element order {o}")),
    };
    r.word("Layers:")?;
    let layer_count = r.int()?;
    let mut layers = Vec::new();
    for _ in 0..layer_count {
        layers.push(layer(r, delta, external_bg)?);
    }
    Ok(BodyToken {
        tag,
        layers,
        delta: bdelta,
        parabolic,
    })
}

fn layer(r: &mut Reader, delta: f64, external_bg: Option<&BgMeshToken>) -> R<LayerToken> {
    r.word("LayerId:")?;
    let tag = r.int()?;
    let ldelta = r.size(delta)?;
    r.word("LayerType:")?;
    let kind = match r.string()?.as_str() {
        "VoronoiVertex" => LayerKind::VoronoiVertex,
        "MovingFront" => LayerKind::VoronoiSegment,
        other => return Err(format!("unsupported layer type {other}")),
    };
    r.word("FixedNodes:")?;
    let fixed_count = r.int()?;
    let mut fixed_nodes = Vec::new();
    for _ in 0..fixed_count {
        fixed_nodes.push(node(r, "NodeId:", delta)?);
    }
    r.word("BGMesh:")?;
    let mut bg = BgMeshToken {
        kind: r.string()?,
        nodes: Vec::new(),
    };
    match bg.kind.as_str() {
        "Explicit" => explicit_bg(r, &mut bg)?,
        "External" => return Err("External background meshes are not supported".into()),
        _ => {}
    }
    let bg = match external_bg {
        Some(e) => BgMeshToken {
            kind: e.kind.clone(),
            nodes: e.nodes.clone(),
        },
        None => bg,
    };
    r.word("Loops:")?;
    let loop_count = r.int()?;
    let mut loops = Vec::new();
    for _ in 0..loop_count {
        r.word("LoopId:")?;
        let _tag = r.int()?;
        r.word("Direction:")?;
        let direction = r.int()?;
        if direction != 1 && direction != -1 {
            return Err("loop direction must be 1 or -1".into());
        }
        r.word("Edges:")?;
        let n = r.int()?;
        let mut edges = Vec::new();
        for _ in 0..n {
            edges.push(r.int()?);
        }
        loops.push(LoopToken { direction, edges });
    }
    Ok(LayerToken {
        tag,
        kind,
        loops,
        fixed_nodes,
        bg,
        delta: ldelta,
    })
}
