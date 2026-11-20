//! `Particle_t` and the location routines of ParticleUtils.F90:
//! `LocateParticles` and `LocateParticleInMeshMarch` (the robust march
//! through the faces, `Particle Locate Robust` on by default).

use elmer_core::MeshView;

use crate::faces::{
    center, element_nodes, segment_element_inside, segment_element_intersection2, Faces, V3,
};

pub(crate) const PARTICLE_ALLOCATED: i32 = 1;
pub(crate) const PARTICLE_INITIATED: i32 = 3;
pub(crate) const PARTICLE_LOCATED: i32 = 4;
pub(crate) const PARTICLE_FACEBOUNDARY: i32 = 6;
pub(crate) const PARTICLE_HIT: i32 = 8;
pub(crate) const PARTICLE_READY: i32 = 9;
pub(crate) const PARTICLE_FIXEDCOORD: i32 = 10;
pub(crate) const PARTICLE_FIXEDVELO: i32 = 11;
pub(crate) const PARTICLE_WALLBOUNDARY: i32 = 12;
pub(crate) const PARTICLE_LOST: i32 = 13;

/// A particle variable (`Particles % Variables`): one value per particle.
#[derive(Clone, Debug)]
pub(crate) struct ParticleVar {
    pub name: String,
    pub values: Vec<f64>,
}

/// `Particle_t` as far as ParticleDynamics uses it. Element indexes are
/// 1-based with 0 for none, as in the Fortran.
#[derive(Default)]
pub(crate) struct Particles {
    pub dim: usize,
    pub time_order: i32,
    pub number: usize,
    pub coord: Vec<V3>,
    pub velo: Vec<V3>,
    pub force: Vec<V3>,
    pub prev_coord: Vec<V3>,
    pub status: Vec<i32>,
    pub element: Vec<usize>,
    pub face: Vec<usize>,
    pub local_min: V3,
    pub local_max: V3,
    pub global_min: V3,
    pub global_max: V3,
    pub dtime: f64,
    pub time: f64,
    pub vars: Vec<ParticleVar>,
}

impl Particles {
    pub fn var(&self, name: &str) -> Option<usize> {
        self.vars.iter().position(|v| v.name == name)
    }

    /// `ParticleVariableCreate` (one dof).
    pub fn var_create(&mut self, name: &str) -> usize {
        if let Some(i) = self.var(name) {
            return i;
        }
        self.vars.push(ParticleVar {
            name: name.to_string(),
            values: vec![0.0; self.number],
        });
        self.vars.len() - 1
    }

    /// `AllocateParticles` for a fresh set of `n` particles.
    pub fn allocate(&mut self, n: usize) {
        self.coord = vec![[0.0; 3]; n];
        self.velo = vec![[0.0; 3]; n];
        self.force = vec![[0.0; 3]; n];
        self.prev_coord = vec![[0.0; 3]; n];
        self.status = vec![PARTICLE_ALLOCATED; n];
        self.element = vec![0; n];
        self.face = vec![0; n];
    }

    /// `DeleteLostParticles`: packs the particles that are not lost.
    pub fn delete_lost(&mut self) {
        let prev = self.number;
        if prev == 0 {
            return;
        }
        let keep: Vec<usize> = (0..prev)
            .filter(|&i| self.status[i] != PARTICLE_LOST)
            .collect();
        if keep.len() == prev {
            return;
        }
        let n2 = keep.len();
        macro_rules! pack {
            ($v:expr, $zero:expr) => {{
                let packed: Vec<_> = keep.iter().map(|&i| $v[i]).collect();
                $v[..n2].copy_from_slice(&packed);
                for x in &mut $v[n2..prev] {
                    *x = $zero;
                }
            }};
        }
        pack!(self.coord, [0.0; 3]);
        pack!(self.velo, [0.0; 3]);
        pack!(self.force, [0.0; 3]);
        pack!(self.prev_coord, [0.0; 3]);
        pack!(self.status, PARTICLE_ALLOCATED);
        pack!(self.face, 0);
        pack!(self.element, 0);
        for v in &mut self.vars {
            pack!(v.values, 0.0);
        }
        self.number = n2;
    }
}

/// The solver keywords the march reads once (`SAVE`d in the Fortran).
pub(crate) struct MarchParams {
    pub stop_at_face: bool,
    pub max_trials: usize,
    pub eps: f64,
}

/// `LocateParticleInMeshMarch` (robust): follows `rinit -> rfin` from
/// element `elem` (1-based, 0 to start from element 1 at its centre) face
/// by face. Returns the element (0 when lost), the status, the stop face
/// (1-based, 0 for none) and updates `rfin`, `velo` at a wall.
#[allow(clippy::too_many_arguments)]
pub(crate) fn march(
    mesh: &dyn MeshView,
    faces: &Faces,
    mp: &MarchParams,
    elem: &mut usize,
    rinit: V3,
    rfin: &mut V3,
    init: bool,
    status: &mut i32,
    accurate_at_face: bool,
    stop_face: &mut usize,
    velo: &mut V3,
) {
    let mut do_init = init;
    let use_center = if *elem == 0 {
        do_init = true;
        *elem = 1;
        true
    } else {
        !accurate_at_face
    };
    if !use_center {
        let mut ds2 = 0.0;
        for d in 0..3 {
            ds2 += (rinit[d] - rfin[d]).powi(2);
        }
        if ds2 < f64::EPSILON {
            return;
        }
    }
    let mut prev_element = 0usize;
    *status = PARTICLE_LOST;
    *stop_face = 0;
    let mut element = *elem;
    let mut rtmp = [0.0f64; 3];
    let mut min_lambda = 0.0f64;
    for i in 1..=mp.max_trials {
        if use_center {
            rtmp = center(&element_nodes(mesh, element));
        } else if i == 1 {
            rtmp = rinit;
        } else {
            for d in 0..3 {
                rtmp[d] += min_lambda * (rfin[d] - rtmp[d]);
            }
        }
        let (face, lambda) = segment_element_intersection2(mesh, faces, element, rtmp, *rfin);
        min_lambda = lambda;
        let Some(f) = face else {
            let mut ds2 = 0.0;
            for d in 0..3 {
                ds2 += (rtmp[d] - rfin[d]).powi(2);
            }
            if ds2 < f64::EPSILON {
                *status = PARTICLE_HIT;
            }
            break;
        };
        if min_lambda > 1.0 - mp.eps {
            *status = PARTICLE_HIT;
            break;
        }
        // The face has no boundary constraint (FindMeshEdges2D sets none),
        // so the particle passes to the other parent or stops at a wall.
        let fc = &faces.faces[f];
        let next;
        if fc.left != 0 && fc.right != 0 {
            next = if element == fc.left {
                fc.right
            } else {
                fc.left
            };
            if mp.stop_at_face && !do_init {
                *status = PARTICLE_FACEBOUNDARY;
            }
        } else {
            next = 0;
            *status = PARTICLE_WALLBOUNDARY;
        }
        if *status == PARTICLE_WALLBOUNDARY || *status == PARTICLE_FACEBOUNDARY {
            *stop_face = f + 1;
            for d in 0..3 {
                rfin[d] = rtmp[d] + min_lambda * (rfin[d] - rtmp[d]);
            }
            *velo = [0.0; 3];
            break;
        }
        if next == prev_element {
            if segment_element_inside(mesh, faces, next, *rfin) {
                element = next;
                *status = PARTICLE_HIT;
                break;
            }
            if segment_element_inside(mesh, faces, element, *rfin) {
                *status = PARTICLE_HIT;
                break;
            }
            *status = PARTICLE_LOST;
            break;
        }
        prev_element = element;
        element = next;
    }
    *elem = if *status == PARTICLE_LOST { 0 } else { element };
}

/// `LocateParticles` in one partition with `Particle Accurate At Face` off.
pub(crate) fn locate_particles(
    mesh: &dyn MeshView,
    faces: &Faces,
    mp: &MarchParams,
    particles: &mut Particles,
    accurate_always: bool,
) {
    // `FaceIndex` carries over from one particle to the next where the march
    // returns before setting it, as the Fortran's local does.
    let mut face_index = 0usize;
    for no in 0..particles.number {
        let mut status = particles.status[no];
        if status >= PARTICLE_LOST || status < PARTICLE_INITIATED {
            continue;
        }
        if status == PARTICLE_WALLBOUNDARY || status == PARTICLE_FIXEDCOORD {
            continue;
        }
        let init_location = status < PARTICLE_LOCATED;
        let accurate_now = if status == PARTICLE_INITIATED {
            false
        } else {
            accurate_always
        };
        let mut element = particles.element[no];
        let mut rfin = particles.coord[no];
        let mut velo = particles.velo[no];
        let rinit = particles.prev_coord[no];
        let mut speed = 0.0;
        for v in velo {
            speed += v.powi(2);
        }
        let speed = speed.sqrt();
        if speed < f64::EPSILON && status > PARTICLE_INITIATED {
            continue;
        }
        march(
            mesh,
            faces,
            mp,
            &mut element,
            rinit,
            &mut rfin,
            init_location,
            &mut status,
            accurate_now,
            &mut face_index,
            &mut velo,
        );
        if element == 0 {
            // The Fortran restores `Rfin0`, `Velo0` that this path never set;
            // the particle is lost and nothing reads them again.
            velo = [0.0; 3];
        }
        particles.element[no] = element;
        particles.status[no] = status;
        particles.face[no] = face_index;
        particles.coord[no] = rfin;
        particles.velo[no] = velo;
    }
}
