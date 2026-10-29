//! Port of `fem/src/modules/FilterTimeSeries.F90` (Elmer @ c0b4bc6, the
//! campaign fork) at the solver-plugin boundary (whole `.so` → Rust
//! cdylib with the same dynamic symbol).
//!
//! Unit: `fortran:fem/src/modules/FilterTimeSeries.F90` — one F77-global
//! entry `filtertimeseries_` with one `CONTAINS` internal procedure
//! (`ParameterName`). The emitted dynamic symbol is exactly
//! `filtertimeseries_` (proven by `nm -D` on the reference `.so`); the
//! internal emits only a local symbol.
//!
//! Boundary (ADR-057 §3, solver-entry gate shape): the entry's derived
//! dummies (`Model`, `Solver`) cross as opaque addresses. Every component
//! access and every Elmer library call goes through the generated-style
//! Fortran `BIND(C)` wrappers in
//! `accessors/filtertimeseries_acc.F90`, compiled against the built
//! `.mod` files of the target configure — the port itself contains no
//! layout assumption and no pinned offset. Order-sensitive FP
//! (`EXP`/`SIN`/`COS`/`SQRT`/`**`, the mean blending, `INT` interval
//! checks) and whole-array updates (`MeanField = 0`) stay in Fortran
//! wrappers, so bit-identity holds by construction wherever the call
//! sequence matches. Plain `+ - * /` weight arithmetic lives in Rust
//! (IEEE-754, identical evaluation order).
//!
//! Faithfulness notes (mirrors FilterTimeSeries.F90 statement for
//! statement):
//! - `SAVE` semantics replicated: `SubroutineVisited`, `TimesVisited`,
//!   and `CumWeight(99)` live in the Rust `static` state with identical
//!   lifetimes; the per-iteration `Var`/`MeanVar`/`MeanField` pointers
//!   live wrapper-side with identical lifetimes (consumed before the
//!   next iteration's re-select, as in the original).
//! - `ParameterName` is inlined as `format!("{base} {novar}")`, which is
//!   byte-identical to its `WRITE(str,'(A,I0)')` (note it uses the host
//!   `NoVar`, always equal to the passed `No` at every call site).
//! - `FilterName` is the keyword *name* (`format!("Filter {novar}")`),
//!   because the original assigns `FilterName = ParameterName(...)`
//!   without any list read; `'Sine'`/`'Cosine'` override for series.
//! - `GotOper` is never set `.TRUE.` in the original, so its warn arm is
//!   dead; the transcription keeps the flag literal with the same value.
//! - `dt`/`TransientSimulation` arrive by reference; `dt` is read
//!   (filter weights), the transient flag is unused, exactly as in the
//!   original.
//! - Reference-UB paths are given the author's evident intent and
//!   documented: an unknown `Operator` with no series keeps
//!   `MeanDofs = InstDofs`; `Reset Interval` on a not-yet-created field
//!   still resets `CumWeight` but skips the unassociated store; a
//!   missing `Start Time` leaves the interval origin at `0.0`. None is
//!   exercised by the proof configs.
//! - No Rust threads: the solver entries run serially, so the `Mutex`
//!   is only a `static`-sharing guard, never contended across threads.

use std::ffi::c_void;
use std::os::raw::{c_double, c_int};
use std::sync::Mutex;

extern "C" {
    fn fta_begin(h_model: *const c_void, h_solver: *const c_void, dt: c_double);
    fn fta_info_banner();
    fn fta_get_time(t: *mut c_double, found: *mut c_int);
    fn fta_get_string(
        key: *const u8,
        keylen: c_int,
        buf: *mut u8,
        buflen: c_int,
        outlen: *mut c_int,
        found: *mut c_int,
    );
    fn fta_get_creal(key: *const u8, keylen: c_int, v: *mut c_double, found: *mut c_int);
    fn fta_get_int(key: *const u8, keylen: c_int, v: *mut c_int, found: *mut c_int);
    fn fta_sim_creal(key: *const u8, keylen: c_int, v: *mut c_double);
    fn fta_realtime() -> c_double;
    fn fta_first_order() -> c_int;
    fn fta_select_var(name: *const u8, namelen: c_int, assoc: *mut c_int);
    fn fta_var_info(inst_dofs: *mut c_int, nsize: *mut c_int, has_perm: *mut c_int);
    fn fta_mean_select(
        name: *const u8,
        namelen: c_int,
        assoc: *mut c_int,
        mean_dofs: *mut c_int,
        nsize: *mut c_int,
    );
    fn fta_mean_create(
        name: *const u8,
        namelen: c_int,
        mean_dofs: c_int,
        nsize: c_int,
        nseries: c_int,
    );
    fn fta_decay(mode: c_int, dt: c_double, st: c_double, relax: *mut c_double);
    fn fta_reset_check(
        st: c_double,
        t0: c_double,
        time: c_double,
        dt: c_double,
        do_reset: *mut c_int,
    );
    fn fta_fcoeff(freq: c_double, time: c_double, t0: c_double, fcoeff: *mut c_double);
    #[allow(clippy::too_many_arguments)]
    fn fta_filter_loop(
        oper_code: c_int,
        loopsize: c_int,
        instdofs: c_int,
        ratio: c_double,
        cold: c_double,
        cnew: c_double,
        fcoeff: c_double,
        nsine: c_int,
        ncosine: c_int,
        has_prev: c_int,
        end_ratio: c_int,
    );
    fn fta_warn_novar(name: *const u8, namelen: c_int);
    fn fta_fatal_freq();
    fn fta_fatal_scalar();
    fn fta_mean_nullify();
}

/// First-call SAVE state: `SubroutineVisited`, `TimesVisited`, and
/// `CumWeight(99)`, with identical lifetimes to the original.
struct FtsState {
    visited: bool,
    times: c_int,
    cumw: [c_double; 99],
}

impl FtsState {
    const fn new() -> Self {
        Self {
            visited: false,
            times: 0,
            cumw: [0.0; 99],
        }
    }
}

static STATE: Mutex<FtsState> = Mutex::new(FtsState::new());

fn get_string(key: &str) -> (Vec<u8>, bool) {
    let mut buf = vec![0u8; 256];
    let (mut outlen, mut found) = (0, 0);
    unsafe {
        fta_get_string(
            key.as_ptr(),
            key.len() as c_int,
            buf.as_mut_ptr(),
            buf.len() as c_int,
            &mut outlen,
            &mut found,
        );
    }
    buf.truncate(outlen as usize);
    (buf, found != 0)
}

fn get_creal(key: &str) -> (c_double, bool) {
    let (mut v, mut found) = (0.0, 0);
    unsafe {
        fta_get_creal(key.as_ptr(), key.len() as c_int, &mut v, &mut found);
    }
    (v, found != 0)
}

fn get_int(key: &str) -> (c_int, bool) {
    let (mut v, mut found) = (0, 0);
    unsafe {
        fta_get_int(key.as_ptr(), key.len() as c_int, &mut v, &mut found);
    }
    (v, found != 0)
}

fn select_var(name: &[u8]) -> bool {
    let mut assoc = 0;
    unsafe {
        fta_select_var(name.as_ptr(), name.len() as c_int, &mut assoc);
    }
    assoc != 0
}

/// `filtertimeseries_` — the solver entry. Argument order and
/// by-reference travel match the Fortran dummy list
/// `(Model, Solver, dtime, TransientSimulation)`; `Model` and the
/// transient flag are unused by the body beyond handle capture, exactly
/// as in the original.
///
/// # Safety
///
/// `model`/`solver` must be the Elmer `Model_t`/`Solver_t` addresses the
/// loader passes. `dt` must be a valid read (it is dereferenced once).
#[export_name = "filtertimeseries_"]
pub unsafe extern "C" fn filtertimeseries_entry(
    model: *mut c_void,
    solver: *mut c_void,
    dt: *const c_double,
    _transient: *const c_int,
) {
    let mut state = STATE.lock().unwrap();
    let dt_in: c_double = unsafe { *dt };
    unsafe { fta_begin(model, solver, dt_in) };

    unsafe { fta_info_banner() };

    if !state.visited {
        state.cumw = [0.0; 99];
        state.visited = true;
    }
    state.times += 1;

    let mut time = 0.0;
    let mut time_found = 0;
    unsafe { fta_get_time(&mut time, &mut time_found) };
    let _ = time_found;

    let mut old_name: Option<Vec<u8>> = None;
    let mut old_assoc = false;

    for novar in 1..=99 {
        let mut got_stuff = false;
        let mut end_ratio = false;

        // `Name = ParameterName('Variable',NoVar)`; value read.
        let (varname, found) = get_string(&format!("Variable {novar}"));
        let var_name: Vec<u8>;
        if found {
            if !select_var(&varname) {
                unsafe { fta_warn_novar(varname.as_ptr(), varname.len() as c_int) };
                continue;
            }
            var_name = varname.clone();
            old_name = Some(varname);
            old_assoc = true;
            got_stuff = true;
        } else if old_assoc {
            // `Var => OldVar; VarName = OldVarName` — re-select by the
            // stored name (same target the pointer would alias).
            let old = old_name.clone().unwrap_or_default();
            select_var(&old);
            var_name = old;
        } else {
            break;
        }

        // `Oper = ListGetString(...); IF(.NOT. Found) Oper = 'none'`.
        let (oper, oper_found) = get_string(&format!("Operator {novar}"));
        let oper = if oper_found { oper } else { b"none".to_vec() };
        let got_oper = false;
        got_stuff = got_stuff || oper_found;

        let mut dt = dt_in;
        let mut ratio: c_double = 1.0;

        // Start after given wall-clock time.
        let (t0_rt, found) = get_creal("Start Real Time");
        if found && unsafe { fta_realtime() } < t0_rt {
            ratio = 0.0;
        }

        // Start after relative given wall-clock time.
        let (q, found) = get_creal("Real Time Max Fraction");
        if found {
            let mut tmax = 0.0;
            let key = "Real Time Max";
            unsafe { fta_sim_creal(key.as_ptr(), key.len() as c_int, &mut tmax) };
            if unsafe { fta_realtime() } < q * tmax {
                ratio = 0.0;
            }
        }

        // Start watch time.
        let (mut t0, mut _t0set) = (0.0, false);
        let (v, found) = get_creal(&format!("Start Time {novar}"));
        if found {
            t0 = v;
            _t0set = true;
            if time < t0 {
                ratio = 0.0;
            } else if time - dt < t0 {
                ratio = (time - t0) / dt;
            }
            got_stuff = true;
        }

        // Stop watch time. (`t1` has no cross-branch flow in the original:
        // the Stop Time value is consumed inside its own branch, and Stop
        // Cycle overwrites before reading — so both uses are scoped. `t0`
        // stays shared: Start Time `t0` flows into Reset/Fourier below.)
        let (v, found) = get_creal(&format!("Stop Time {novar}"));
        if found {
            let t1 = v;
            if time - dt > t1 {
                ratio = 0.0;
            } else if time > t1 {
                ratio = (t1 + dt - time) / dt;
                end_ratio = true;
            }
            got_stuff = true;
        }

        // Start timestep number.
        let (n0, found) = get_int(&format!("Start Timestep {novar}"));
        if found {
            if state.times < n0 {
                ratio = 0.0;
            }
            dt = 1.0;
            got_stuff = true;
        }

        // Stop timestep number.
        let (n1, found) = get_int(&format!("Stop Timestep {novar}"));
        if found {
            if state.times > n1 {
                ratio = 0.0;
            }
            dt = 1.0;
            got_stuff = true;
        }

        // Reset mean after elapsed time.
        let (st, found) = get_creal(&format!("Reset Interval {novar}"));
        if found {
            let mut do_reset = 0;
            unsafe { fta_reset_check(st, t0, time, dt, &mut do_reset) };
            if do_reset != 0 {
                state.cumw[(novar - 1) as usize] = 0.0;
            }
            got_stuff = true;
        }

        // Relax by giving time scale.
        let mut relax: c_double = 1.0;
        let (st, found) = get_creal(&format!("Decay Time {novar}"));
        if found {
            unsafe { fta_decay(0, dt, st, &mut relax) };
            got_stuff = true;
        }

        // Relax by giving timestep scale (overrides the time scale,
        // exactly as in the original).
        relax = 1.0;
        let (n0, found) = get_int(&format!("Decay Timesteps {novar}"));
        if found {
            unsafe { fta_decay(1, 0.0, n0 as c_double, &mut relax) };
            dt = 1.0;
            got_stuff = true;
        }

        // Filter with externally provided function (sin, cos, ...).
        let mut phase: c_double = 1.0;
        let (st, found) = get_creal(&format!("Time Filter {novar}"));
        if found {
            phase = st;
            got_stuff = true;
        }

        // Explicit sine and cosine series of given degree.
        let (nsine, _) = get_int(&format!("Sine Series {novar}"));
        let (ncosine, _) = get_int(&format!("Cosine Series {novar}"));
        let nseries = nsine.max(ncosine);
        got_stuff = got_stuff || nseries > 0;

        // If no new control parameters were obtained exit the loop.
        if !got_stuff {
            break;
        }

        let mut freq = 0.0;
        if nseries > 0 {
            let (f, found) = get_creal("Frequency");
            if !found {
                unsafe { fta_fatal_freq() };
                std::process::abort();
            }
            freq = f;
            let (v, found) = get_creal(&format!("Start Cycle {novar}"));
            if found {
                t0 = 1.0 / v;
                if time < t0 {
                    ratio = 0.0;
                } else if time - dt < t0 {
                    ratio = (time - t0) / dt;
                }
            }
            let (v, found) = get_creal(&format!("Stop Cycle {novar}"));
            if found {
                let t1 = 1.0 / v;
                if time - dt > t1 {
                    ratio = 0.0;
                } else if time > t1 {
                    ratio = (t1 + dt - time) / dt;
                    end_ratio = true;
                }
            }
        }

        let (mut inst_dofs, mut inst_size, mut _has_perm) = (0, 0, 0);
        unsafe { fta_var_info(&mut inst_dofs, &mut inst_size, &mut _has_perm) };

        let prev_exists = unsafe { fta_first_order() };

        // `FilterName`: 'Sine'/'Cosine' for series, else the keyword
        // *name* `Filter <NoVar>` (the original never reads its value).
        let filter_name: Vec<u8> = if nsine > 0 {
            b"Sine".to_vec()
        } else if ncosine > 0 {
            b"Cosine".to_vec()
        } else {
            format!("Filter {novar}").into_bytes()
        };
        let mut mean_name = [filter_name.as_slice(), b" ", var_name.as_slice()].concat();
        // `WRITE(MeanVarName,'(A)')` truncates to `MAX_NAME_LEN`.
        mean_name.truncate(128);

        let (mut assoc, mut mean_dofs, mut nsize) = (0, 0, 0);
        unsafe {
            fta_mean_select(
                mean_name.as_ptr(),
                mean_name.len() as c_int,
                &mut assoc,
                &mut mean_dofs,
                &mut nsize,
            );
        }
        if assoc == 0 {
            if oper.as_slice() == b"length" {
                mean_dofs = 1;
            } else {
                if inst_dofs > 1 && nseries > 0 {
                    unsafe { fta_fatal_scalar() };
                    std::process::abort();
                }
                // Unknown operators fall here with `GotOper == .FALSE.`
                // (dead warn in the original, `MeanDofs` unset — reference
                // UB); the evident intent is the pass-through sizing.
                mean_dofs = inst_dofs;
            }
            if nseries > 0 {
                mean_dofs = nseries;
            }
            nsize = inst_size * (mean_dofs / inst_dofs);
            unsafe {
                fta_mean_create(
                    mean_name.as_ptr(),
                    mean_name.len() as c_int,
                    mean_dofs,
                    nsize,
                    nseries,
                );
            }
        }

        let weight = dt * ratio;
        // `CumWeight(NoVar) = Relax * CumWeight(NoVar)`.
        state.cumw[(novar - 1) as usize] *= relax;

        if weight.abs() < c_double::MIN_POSITIVE {
            continue;
        }

        // Weight for the filtered set / current timestep.
        let cold = state.cumw[(novar - 1) as usize]
            / (state.cumw[(novar - 1) as usize] + weight);
        let (cnew, loopsize, fcoeff) = if nseries > 0 {
            let mut f = 0.0;
            unsafe { fta_fcoeff(freq, time, t0, &mut f) };
            (
                2.0 * weight / (state.cumw[(novar - 1) as usize] + weight),
                nsize / nseries,
                f,
            )
        } else {
            (
                phase * weight / (state.cumw[(novar - 1) as usize] + weight),
                nsize,
                0.0,
            )
        };

        let oper_code = if oper.as_slice() == b"square" {
            1
        } else if oper.as_slice() == b"abs" {
            2
        } else if oper.as_slice() == b"length" {
            3
        } else {
            // 'none' and anything else (the dead `GotOper` warn never
            // fires in the original either).
            let _ = got_oper;
            0
        };

        unsafe {
            fta_filter_loop(
                oper_code,
                loopsize,
                inst_dofs,
                ratio,
                cold,
                cnew,
                fcoeff,
                nsine,
                ncosine,
                prev_exists,
                i32::from(end_ratio),
            );
        }

        state.cumw[(novar - 1) as usize] += weight;
        unsafe { fta_mean_nullify() };
    }
}
