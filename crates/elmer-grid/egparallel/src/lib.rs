// Faithful Rust transcription of elmergrid/src/egparallel.c @ 9f6af2f85.
// Workspace layout: shared declarations live in eg-prelude; this crate
// appends only the unit body.
// Config scope: exe-only unit (not in elmergrid_core); USE_METIS=1 branches
// live (exe config), HAVE_UNISTD_H-undefined include is a no-op (unconditional
// unistd.h include above it).
// SAFETY: single-threaded utility; raw-pointer arithmetic mirrors C 1:1.
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(static_mut_refs)]
#![allow(dead_code)]
#![allow(non_snake_case)]

use eg_prelude::*;

// ---- METIS ABI surface (pinned to bundled metis-5.1.0) -------------------
// IDXTYPEWIDTH=32 -> idx_t is int32 (== c_int); REALTYPEWIDTH=32 ->
// real_t is float, but every call below passes NULL weights, so the width
// is unobservable. Enum values counted from metis.h (zero-based).
pub type idx_t = c_int;
pub type real_t = c_float;

pub const METIS_NOPTIONS: usize = 40;
pub const METIS_OPTION_PTYPE: usize = 0;
pub const METIS_OPTION_OBJTYPE: usize = 1;
pub const METIS_OPTION_CTYPE: usize = 2;
pub const METIS_OPTION_IPTYPE: usize = 3;
pub const METIS_OPTION_RTYPE: usize = 4;
pub const METIS_OPTION_DBGLVL: usize = 5;
pub const METIS_OPTION_NCUTS: usize = 7;
pub const METIS_OPTION_SEED: usize = 8;
pub const METIS_OPTION_MINCONN: usize = 10;
pub const METIS_OPTION_CONTIG: usize = 11;
pub const METIS_OPTION_NUMBERING: usize = 17;

pub const METIS_CTYPE_SHEM: idx_t = 1;
pub const METIS_IPTYPE_GROW: idx_t = 0;
pub const METIS_RTYPE_GREEDY: idx_t = 1;
pub const METIS_OBJTYPE_VOL: idx_t = 1;

// ---- sibling-unit imports (declared in headers, resolved at final link) --
extern "C" {
    fn AddExtension(fname1: *const c_char, fname2: *mut c_char, newext: *const c_char);
    fn bigerror(error_text: *const c_char);
    fn smallerror(error_text: *const c_char);
    fn dvector(nl: c_int, nh: c_int) -> *mut c_double;
    fn ivector(nl: c_int, nh: c_int) -> *mut c_int;
    fn imatrix(nrl: c_int, nrh: c_int, ncl: c_int, nch: c_int) -> *mut *mut c_int;
    fn free_dvector(v: *mut c_double, nl: c_int, nh: c_int);
    fn free_ivector(v: *mut c_int, nl: c_int, nh: c_int);
    fn free_imatrix(m: *mut *mut c_int, nrl: c_int, nrh: c_int, ncl: c_int, nch: c_int);
    fn next_int(start: *mut *mut c_char) -> c_int;
    fn next_real(start: *mut *mut c_char) -> Real;
    fn SortIndex(n: c_int, key: *mut c_double, ord: *mut c_int);
    fn CalculateIndexwidth(data: *mut FemType, indxis: c_int, indx: *mut c_int) -> c_int;
    fn CreateDualGraph(data: *mut FemType, unconnected: c_int, info: c_int) -> c_int;
    fn CreateInverseTopology(data: *mut FemType, info: c_int) -> c_int;
    fn CreateNodalGraph(data: *mut FemType, full: c_int, info: c_int) -> c_int;
    fn GetBoundaryElement(
        sideind: c_int,
        bound: *mut BoundaryType,
        data: *mut FemType,
        ind: *mut c_int,
        sideelemtype: *mut c_int,
    );
    fn GetElementSide(
        element: c_int,
        side: c_int,
        normal: c_int,
        data: *mut FemType,
        ind: *mut c_int,
        sideelemtype: *mut c_int,
    );
    fn GetMaxElementType(data: *mut FemType) -> c_int;
    fn SetConnectedElements(data: *mut FemType, info: c_int) -> c_int;
    fn METIS_SetDefaultOptions(options: *mut idx_t) -> c_int;
    fn METIS_PartGraphRecursive(
        nvtxs: *mut idx_t,
        ncon: *mut idx_t,
        xadj: *mut idx_t,
        adjncy: *mut idx_t,
        vwgt: *mut idx_t,
        vsize: *mut idx_t,
        adjwgt: *mut idx_t,
        nparts: *mut idx_t,
        tpwgts: *mut real_t,
        ubvec: *mut real_t,
        options: *mut idx_t,
        edgecut: *mut idx_t,
        part: *mut idx_t,
    ) -> c_int;
    fn METIS_PartGraphKway(
        nvtxs: *mut idx_t,
        ncon: *mut idx_t,
        xadj: *mut idx_t,
        adjncy: *mut idx_t,
        vwgt: *mut idx_t,
        vsize: *mut idx_t,
        adjwgt: *mut idx_t,
        nparts: *mut idx_t,
        tpwgts: *mut real_t,
        ubvec: *mut real_t,
        options: *mut idx_t,
        edgecut: *mut idx_t,
        part: *mut idx_t,
    ) -> c_int;
    fn METIS_PartMeshDual(
        ne: *mut idx_t,
        nn: *mut idx_t,
        eptr: *mut idx_t,
        eind: *mut idx_t,
        vwgt: *mut idx_t,
        vsize: *mut idx_t,
        ncommon: *mut idx_t,
        nparts: *mut idx_t,
        tpwgts: *mut real_t,
        options: *mut idx_t,
        objval: *mut idx_t,
        epart: *mut idx_t,
        npart: *mut idx_t,
    ) -> c_int;
    fn METIS_PartMeshNodal(
        ne: *mut idx_t,
        nn: *mut idx_t,
        eptr: *mut idx_t,
        eind: *mut idx_t,
        vwgt: *mut idx_t,
        vsize: *mut idx_t,
        nparts: *mut idx_t,
        tpwgts: *mut real_t,
        options: *mut idx_t,
        objval: *mut idx_t,
        epart: *mut idx_t,
        npart: *mut idx_t,
    ) -> c_int;
    fn METIS_NodeND(
        nvtxs: *mut idx_t,
        xadj: *mut idx_t,
        adjncy: *mut idx_t,
        vwgt: *mut idx_t,
        options: *mut idx_t,
        perm: *mut idx_t,
        iperm: *mut idx_t,
    ) -> c_int;
}
// NOTE: C file-local LONGLINE (2048) as a Rust const.
pub const FUSE_LONGLINE: usize = 2048;

/// # Safety
///
/// partitioned-solution fusion; `goto skip` becomes a guarded block;
/// res is block-local (assigned iff novctrs != 0, read iff the same).
#[no_mangle]
pub unsafe extern "C" fn FuseSolutionElmerPartitioned(
    prefix: *mut c_char,
    outfile: *mut c_char,
    decimals: c_int,
    parts: c_int,
    minstep: c_int,
    maxstep: c_int,
    dstep: c_int,
    info: c_int,
) -> c_int {
    
    
    #[allow(unused_assignments)]
    let (mut novctrs, mut elemcode): (c_int, c_int) = (0, 0);
    let (totknots, totelements, mut sumknots, mut sumelements): (c_int, c_int, c_int, c_int);
    // NOTE: timesteps is an sscanf-style &mut target below (next_int writes
    // through &mut cp, but timesteps itself is assigned by value); the 0
    // init is dead by analysis on the nofiles>=2 path yet required for the
    // possibly-uninitialized proof — allowed, not deleted. Same for cp/xyz.
    #[allow(unused_assignments)]
    let (mut timesteps, mut i, mut j, mut k, mut l, mut step): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    ) = (0, 0, 0, 0, 0, 0);
    let mut ind: [c_int; MAXNODESD3] = [0; MAXNODESD3];
    let (mut nofiles, mut activestep): (c_int, c_int);
    // NOTE: assigned iff novctrs != 0 below, read iff the same (the
    // `goto skip` shape); block-local so the compiler proves it.
    let mut inn: [*mut FILE; MAXPARTITIONS + 1] = [core::ptr::null_mut(); MAXPARTITIONS + 1];
    let (mut intest, out): (*mut FILE, *mut FILE);
    let mut line: [c_char; FUSE_LONGLINE] = [0; FUSE_LONGLINE];
    let mut filename: [c_char; MAXFILESIZE] = [0; MAXFILESIZE];
    let mut text: [c_char; MAXNAMESIZE] = [0; MAXNAMESIZE];
    let mut outstyle: [c_char; MAXFILESIZE] = [0; MAXFILESIZE];
    // NOTE: cp/charend/xyz are &mut- and return-slot handoffs to fgets and
    // the next_* parsers; inits are analysis-dead but proof-required.
    // charend is write-only in C (never read back).
    #[allow(unused_assignments)]
    let (mut cp, mut _charend): (*mut c_char, *mut c_char) =
        (core::ptr::null_mut(), core::ptr::null_mut());
    #[allow(unused_assignments)]
    let (mut x, mut y, mut z): (Real, Real, Real) = (0.0, 0.0, 0.0);
    if (minstep != 0 || maxstep != 0 || dstep != 0)
        && info != 0 {
            printf(
                c"Saving results in the interval from %d to %d with step %d\n".as_ptr()
                    as *const c_char,
                minstep,
                maxstep,
                dstep,
            );
        }
    i = 0;
    loop {
        sprintf(
            filename.as_mut_ptr(),
            c"%s.ep.%d".as_ptr() as *const c_char,
            prefix,
            i,
        );
        intest = fopen(filename.as_ptr(), c"r".as_ptr() as *const c_char);
        if intest.is_null() {
            break;
        }
        if i <= MAXPARTITIONS as c_int {
            inn[i as usize] = intest;
        }
        i += 1;
    }
    if i > MAXPARTITIONS as c_int {
        printf(
            c"**********************************************************\n".as_ptr()
                as *const c_char,
        );
        printf(
            c"Only data for %d partitions is fused (%d)\n".as_ptr() as *const c_char,
            MAXPARTITIONS as c_int,
            i,
        );
        printf(
            c"**********************************************************\n".as_ptr()
                as *const c_char,
        );
        i = MAXPARTITIONS as c_int;
    }
    nofiles = i;
    if nofiles < 2 {
        printf(
            c"Opening of partitioned data from file %s wasn't successful!\n".as_ptr()
                as *const c_char,
            filename.as_ptr(),
        );
        return 2;
    } else {
        if parts > 0 {
            nofiles = c_min_int(parts, nofiles);
        }
        if info != 0 {
            printf(
                c"Loading Elmer results from %d partitions.\n".as_ptr() as *const c_char,
                nofiles,
            );
        }
    }
    if (minstep != 0 || maxstep != 0 || dstep != 0)
        && info != 0 {
            printf(
                c"Saving results in the interval from %d to %d with step %d\n".as_ptr()
                    as *const c_char,
                minstep,
                maxstep,
                dstep,
            );
        }
    // NOTE: Ivector macro family (0-based here) = ivector.
    let noknots: *mut c_int = ivector(0, nofiles - 1);
    let noelements: *mut c_int = ivector(0, nofiles - 1);
    sumknots = 0;
    sumelements = 0;
    i = 0;
    while i < nofiles {
        _charend = fgets(line.as_mut_ptr(), FUSE_LONGLINE as c_int, inn[i as usize]);
        if i == 0 {
            cp = line.as_mut_ptr();
            *noknots.wrapping_offset(i as isize) = next_int(&mut cp);
            *noelements.wrapping_offset(i as isize) = next_int(&mut cp);
            novctrs = next_int(&mut cp);
            timesteps = next_int(&mut cp);
        } else {
            sscanf(
                line.as_ptr(),
                c"%d %d".as_ptr() as *const c_char,
                &mut *noknots.wrapping_offset(i as isize),
                &mut *noelements.wrapping_offset(i as isize),
            );
        }
        sumknots += *noknots.wrapping_offset(i as isize);
        sumelements += *noelements.wrapping_offset(i as isize);
        i += 1;
    }
    totknots = sumknots;
    totelements = sumelements;
    // NOTE: `res` lives in the novctrs-guarded block below (see SAFETY).
    if info != 0 {
        printf(
            c"There are %d nodes, %d elements and %d vectors.\n".as_ptr() as *const c_char,
            totknots,
            sumelements,
            novctrs,
        );
    }
    AddExtension(
        outfile,
        filename.as_mut_ptr(),
        c"ep".as_ptr() as *const c_char,
    );
    if info != 0 {
        printf(
            c"Saving ElmerPost data to %s.\n".as_ptr() as *const c_char,
            filename.as_ptr(),
        );
    }
    out = fopen(filename.as_ptr(), c"w".as_ptr() as *const c_char);
    if out.is_null() {
        printf(c"opening of file was not successful\n".as_ptr() as *const c_char);
        return 3;
    }
    i = timesteps;
    if minstep != 0 || maxstep != 0 || dstep != 0 {
        if dstep > 1 {
            i = 0;
            step = minstep;
            while step <= maxstep {
                if (step - minstep) % dstep == 0 {
                    i += 1;
                }
                step += 1;
            }
        } else {
            i = maxstep - minstep + 1;
        }
    }
    fprintf(
        out,
        c"%d %d %d %d %s %s".as_ptr() as *const c_char,
        totknots,
        totelements,
        novctrs + 1,
        i,
        c"scalar: Partition".as_ptr() as *const c_char,
        cp,
    );
    if info != 0 {
        printf(
            c"Reading and writing %d coordinates.\n".as_ptr() as *const c_char,
            totknots,
        );
    }
    sprintf(
        outstyle.as_mut_ptr(),
        c"%%.%dg %%.%dg %%.%dg\n".as_ptr() as *const c_char,
        decimals,
        decimals,
        decimals,
    );
    j = 0;
    while j < nofiles {
        i = 1;
        while i <= *noknots.wrapping_offset(j as isize) {
            loop {
                _charend = fgets(line.as_mut_ptr(), FUSE_LONGLINE as c_int, inn[j as usize]);
                if line[0] != '#' as c_char {
                    break;
                }
            }
            sscanf(
                line.as_ptr(),
                c"%le %le %le".as_ptr() as *const c_char,
                &mut x,
                &mut y,
                &mut z,
            );
            fprintf(out, outstyle.as_ptr(), x, y, z);
            i += 1;
        }
        j += 1;
    }
    if info != 0 {
        printf(
            c"Reading and writing %d element topologies.\n".as_ptr() as *const c_char,
            totelements,
        );
    }
    sumknots = 0;
    j = 0;
    while j < nofiles {
        i = 1;
        while i <= *noelements.wrapping_offset(j as isize) {
            loop {
                _charend = fgets(line.as_mut_ptr(), FUSE_LONGLINE as c_int, inn[j as usize]);
                if line[0] != '#' as c_char {
                    break;
                }
            }
            sscanf(
                line.as_ptr(),
                c"%s".as_ptr() as *const c_char,
                text.as_mut_ptr(),
            );
            cp = strstr(line.as_ptr(), c" ".as_ptr() as *const c_char);
            elemcode = next_int(&mut cp);
            k = 0;
            while k < elemcode % 100 {
                // NOTE: strspn returns size_t; l is int in C (narrowing).
                l = strspn(cp, c" ".as_ptr() as *const c_char) as c_int;
                if l == 0 {
                    _charend = fgets(line.as_mut_ptr(), FUSE_LONGLINE as c_int, inn[j as usize]);
                    cp = line.as_mut_ptr();
                }
                ind[k as usize] = next_int(&mut cp);
                k += 1;
            }
            if elemcode == 102 {
                elemcode = 101;
            }
            fprintf(
                out,
                c"%s %d".as_ptr() as *const c_char,
                text.as_ptr(),
                elemcode,
            );
            k = 0;
            while k < elemcode % 100 {
                fprintf(
                    out,
                    c" %d".as_ptr() as *const c_char,
                    ind[k as usize] + sumknots,
                );
                k += 1;
            }
            fprintf(out, c"\n".as_ptr() as *const c_char);
            i += 1;
        }
        sumknots += *noknots.wrapping_offset(j as isize);
        j += 1;
    }
    if info != 0 {
        printf(
            c"Reading and writing %d degrees of freedom.\n".as_ptr() as *const c_char,
            novctrs,
        );
    }
    // NOTE: `goto skip` — the dof section runs iff novctrs != 0.
    if novctrs != 0 {
        let res: *mut Real = dvector(1, novctrs);
        sprintf(
            outstyle.as_mut_ptr(),
            c"%%.%dg ".as_ptr() as *const c_char,
            decimals,
        );
        activestep = FALSE;
        if maxstep != 0 {
            timesteps = c_min_int(timesteps, maxstep);
        }
        step = 1;
        while step <= timesteps {
            if step >= minstep {
                if dstep > 0 {
                    activestep = ((step - minstep) % dstep == 0) as c_int;
                } else {
                    activestep = TRUE;
                }
            }
            k = 0;
            while k < nofiles {
                i = 1;
                while i <= *noknots.wrapping_offset(k as isize) {
                    loop {
                        _charend =
                            fgets(line.as_mut_ptr(), FUSE_LONGLINE as c_int, inn[k as usize]);
                        if activestep != 0
                            && k == 0
                                && !strstr(line.as_ptr(), c"#time".as_ptr() as *const c_char)
                                    .is_null()
                            {
                                fprintf(out, c"%s".as_ptr() as *const c_char, line.as_ptr());
                                fprintf(stderr, c"%s".as_ptr() as *const c_char, line.as_ptr());
                            }
                        if line[0] != '#' as c_char {
                            break;
                        }
                    }
                    if activestep != 0 {
                        cp = line.as_mut_ptr();
                        j = 1;
                        while j <= novctrs {
                            *res.wrapping_offset(j as isize) = next_real(&mut cp);
                            j += 1;
                        }
                        fprintf(out, c"%d ".as_ptr() as *const c_char, k + 1);
                        j = 1;
                        while j <= novctrs {
                            fprintf(out, outstyle.as_ptr(), *res.wrapping_offset(j as isize));
                            j += 1;
                        }
                        fprintf(out, c"\n".as_ptr() as *const c_char);
                    }
                    i += 1;
                }
                k += 1;
            }
            step += 1;
        }
    }
    // NOTE: `skip:` label target.
    i = 0;
    while i < nofiles {
        fclose(inn[i as usize]);
        i += 1;
    }
    fclose(out);
    if info != 0 {
        printf(c"Successfully fused partitioned Elmer results\n".as_ptr() as *const c_char);
    }
    0
}
// SAFETY: read-only table build over valid elements (C contract).
unsafe fn CreatePartitionTable(data: *mut FemType, info: c_int) -> c_int {
    let (mut i, mut j, mut k, mut m): (c_int, c_int, c_int, c_int);
    let (noelements, noknots, partitions, mut nonodes, periodic): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (mut maxneededtimes, mut sharings, mut part, mut ind, mut hit, mut notinany): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    // NOTE: assigned iff periodic below; reads are periodic-guarded (the
    // flag never changes mid-function), so the null init is unobservable.
    let mut indxper: *mut c_int = core::ptr::null_mut();
    let mut debug: c_int;
    printf(
        c"Creating a table showing all parenting partitions of nodes.\n".as_ptr()
            as *const c_char,
    );
    if (*data).maxpartitiontable != 0 {
        printf(c"The partition table already exists!\n".as_ptr() as *const c_char);
        smallerror(c"Partition table not done".as_ptr() as *const c_char);
        return 0;
    }
    partitions = (*data).nopartitions;
    if partitions <= 1 {
        bigerror(
            c"Cannot do partition table for less than two partitions".as_ptr() as *const c_char,
        );
    }
    noelements = (*data).noelements;
    periodic = (*data).periodicexist;
    noknots = (*data).noknots;
    if periodic != 0 {
        indxper = (*data).periodic;
    }
    maxneededtimes = 0;
    sharings = 0;
    i = 1;
    while i <= noelements {
        part = *(*data).elempart.wrapping_offset(i as isize);
        nonodes = *(*data).elementtypes.wrapping_offset(i as isize) % 100;
        j = 0;
        while j < nonodes {
            ind = *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize);
            if periodic != 0 {
                ind = *indxper.wrapping_offset(ind as isize);
            }
            debug = FALSE;
            if debug != 0 {
                printf(
                    c"ind=%d i=%d j=%d part=%d\n".as_ptr() as *const c_char,
                    ind,
                    i,
                    j,
                    part,
                );
            }
            hit = 0;
            k = 1;
            while k <= maxneededtimes {
                if *(*data).partitiontable[k as usize].wrapping_offset(ind as isize) == 0 {
                    hit = k;
                }
                if *(*data).partitiontable[k as usize].wrapping_offset(ind as isize) == part {
                    hit = -k;
                }
                if hit != 0 {
                    break;
                }
                k += 1;
            }
            if hit > 0 {
                *(*data).partitiontable[hit as usize].wrapping_offset(ind as isize) = part;
                if hit == 2 {
                    sharings += 1;
                }
            } else if hit == 0 {
                maxneededtimes += 1;
                (*data).partitiontable[maxneededtimes as usize] = ivector(1, noknots);
                m = 1;
                while m <= noknots {
                    *(*data).partitiontable[maxneededtimes as usize].wrapping_offset(m as isize) =
                        0;
                    m += 1;
                }
                *(*data).partitiontable[maxneededtimes as usize].wrapping_offset(ind as isize) =
                    part;
                if maxneededtimes == 2 {
                    sharings += 1;
                }
            }
            if debug != 0 {
                printf(c"hit = %d\n".as_ptr() as *const c_char, hit);
            }
            j += 1;
        }
        i += 1;
    }
    notinany = 0;
    i = 1;
    while i <= noknots {
        if periodic != 0
            && i != *indxper.wrapping_offset(i as isize) {
                i += 1;
                continue;
            }
        hit = FALSE;
        k = 1;
        while k <= maxneededtimes {
            if *(*data).partitiontable[k as usize].wrapping_offset(i as isize) == 0 {
                break;
            }
            if *(*data).partitiontable[k as usize].wrapping_offset(i as isize)
                == *(*data).nodepart.wrapping_offset(i as isize)
            {
                hit = k;
                break;
            }
            k += 1;
        }
        if hit > 1 {
            *(*data).partitiontable[hit as usize].wrapping_offset(i as isize) =
                *(*data).partitiontable[1].wrapping_offset(i as isize);
            *(*data).partitiontable[1].wrapping_offset(i as isize) =
                *(*data).nodepart.wrapping_offset(i as isize);
        } else if hit == 0 {
            // NOTE: `if(0)` diagnostic omitted (runtime-dead).
            notinany += 1;
            *(*data).nodepart.wrapping_offset(i as isize) =
                *(*data).partitiontable[1].wrapping_offset(i as isize);
        }
        i += 1;
    }
    if periodic != 0 {
        i = 1;
        while i <= noknots {
            ind = *indxper.wrapping_offset(i as isize);
            if ind == i {
                i += 1;
                continue;
            }
            k = 1;
            while k <= maxneededtimes {
                if *(*data).partitiontable[k as usize].wrapping_offset(ind as isize) == 0 {
                    break;
                }
                *(*data).partitiontable[k as usize].wrapping_offset(i as isize) =
                    *(*data).partitiontable[k as usize].wrapping_offset(ind as isize);
                k += 1;
            }
            *(*data).nodepart.wrapping_offset(i as isize) =
                *(*data).nodepart.wrapping_offset(ind as isize);
            i += 1;
        }
    }
    if info != 0 {
        printf(
            c"Nodes belong to %d partitions in maximum\n".as_ptr() as *const c_char,
            maxneededtimes,
        );
        printf(
            c"There are %d shared nodes which is %.2f %% of all nodes.\n".as_ptr()
                as *const c_char,
            sharings,
            100.0 * sharings as Real / noknots as Real,
        );
        printf(
            c"The initial owner was not any of the elements for %d nodes\n".as_ptr()
                as *const c_char,
            notinany,
        );
    }
    // NOTE: `if(0)` dump loop omitted (runtime-dead).
    (*data).maxpartitiontable = maxneededtimes;
    (*data).partitiontableexists = TRUE;
    0
}

// SAFETY: two-pass (omstart) majority vote over valid elements (C contract).
unsafe fn PartitionElementsByNodes(data: *mut FemType, info: c_int) -> c_int {
    let (mut i, mut j): (c_int, c_int);
    let (noelements, nopartitions, mut part, mut maxpart, mut maxpart2, mut minpart): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (elempart, nodepart, nodesinpart, cuminpart): (
        *mut c_int,
        *mut c_int,
        *mut c_int,
        *mut c_int,
    );
    let (knows, cumknows): (*mut *mut c_int, *mut *mut c_int);
    let mut set: c_int;
    if (*data).partitionexist == 0 {
        return 1;
    }
    noelements = (*data).noelements;
    nopartitions = (*data).nopartitions;
    elempart = (*data).elempart;
    nodepart = (*data).nodepart;
    nodesinpart = ivector(1, nopartitions);
    cuminpart = ivector(1, nopartitions);
    j = 1;
    while j <= nopartitions {
        *cuminpart.wrapping_offset(j as isize) = 0;
        j += 1;
    }
    knows = imatrix(1, nopartitions, 1, nopartitions);
    cumknows = imatrix(1, nopartitions, 1, nopartitions);
    i = 1;
    while i <= nopartitions {
        j = 1;
        while j <= nopartitions {
            *(*knows.wrapping_offset(i as isize)).wrapping_offset(j as isize) = 0;
            *(*cumknows.wrapping_offset(i as isize)).wrapping_offset(j as isize) = 0;
            j += 1;
        }
        i += 1;
    }
    set = FALSE;
    // NOTE: `omstart:` backward goto becomes a two-pass loop.
    loop {
        i = 1;
        while i <= noelements {
            j = 1;
            while j <= nopartitions {
                *nodesinpart.wrapping_offset(j as isize) = 0;
                j += 1;
            }
            j = 0;
            while j < *(*data).elementtypes.wrapping_offset(i as isize) % 100 {
                part = *nodepart.wrapping_offset(
                    *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize)
                        as isize,
                );
                *nodesinpart.wrapping_offset(part as isize) += 1;
                j += 1;
            }
            maxpart = 1;
            maxpart2 = 1;
            j = 1;
            while j <= nopartitions {
                if *nodesinpart.wrapping_offset(j as isize)
                    > *nodesinpart.wrapping_offset(maxpart as isize)
                {
                    maxpart = j;
                }
                j += 1;
            }
            if maxpart == 1 {
                maxpart2 = 2;
            }
            j = 1;
            while j <= nopartitions {
                if j == maxpart {
                    j += 1;
                    continue;
                }
                if *nodesinpart.wrapping_offset(j as isize)
                    > *nodesinpart.wrapping_offset(maxpart2 as isize)
                {
                    maxpart2 = j;
                }
                j += 1;
            }
            if *nodesinpart.wrapping_offset(maxpart as isize)
                > *nodesinpart.wrapping_offset(maxpart2 as isize)
            {
                if set != 0 {
                    *elempart.wrapping_offset(i as isize) = maxpart;
                } else {
                    *cuminpart.wrapping_offset(maxpart as isize) += 1;
                }
            } else {
                if set != 0 {
                    *(*cumknows.wrapping_offset(maxpart as isize))
                        .wrapping_offset(maxpart2 as isize) += 1;
                    if *(*cumknows.wrapping_offset(maxpart as isize))
                        .wrapping_offset(maxpart2 as isize)
                        > *(*knows.wrapping_offset(maxpart as isize))
                            .wrapping_offset(maxpart2 as isize)
                            / 2
                    {
                        *elempart.wrapping_offset(i as isize) = maxpart2;
                        *cuminpart.wrapping_offset(maxpart2 as isize) += 1;
                    } else {
                        *elempart.wrapping_offset(i as isize) = maxpart;
                        *cuminpart.wrapping_offset(maxpart as isize) += 1;
                    }
                } else {
                    *(*knows.wrapping_offset(maxpart as isize))
                        .wrapping_offset(maxpart2 as isize) += 1;
                }
            }
            i += 1;
        }
        if set != 0 {
            break;
        }
        set = TRUE;
    }
    minpart = *cuminpart.wrapping_offset(1);
    maxpart = *cuminpart.wrapping_offset(1);
    j = 1;
    while j <= nopartitions {
        minpart = c_min_int(minpart, *cuminpart.wrapping_offset(j as isize));
        maxpart = c_max_int(maxpart, *cuminpart.wrapping_offset(j as isize));
        j += 1;
    }
    if info != 0 {
        printf(
            c"Set the element partitions by the dominating nodal partition\n".as_ptr()
                as *const c_char,
        );
        printf(
            c"There are from %d to %d elements in the %d partitions.\n".as_ptr() as *const c_char,
            minpart,
            maxpart,
            nopartitions,
        );
    }
    free_ivector(nodesinpart, 1, nopartitions);
    free_ivector(cuminpart, 1, nopartitions);
    free_imatrix(knows, 1, nopartitions, 1, nopartitions);
    free_imatrix(cumknows, 1, nopartitions, 1, nopartitions);
    0
}
// SAFETY: two-pass (omstart) majority vote over valid nodes (C contract).
unsafe fn PartitionNodesByElements(data: *mut FemType, info: c_int) -> c_int {
    let (mut i, mut j, mut k): (c_int, c_int, c_int);
    let (noknots, nopartitions, mut part, mut minpart, mut maxpart): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let mut maxpart2: c_int;
    let (cuminpart, knows, cumknows): (*mut c_int, *mut *mut c_int, *mut *mut c_int);
    let (elempart, nodepart, nodesinpart): (*mut c_int, *mut c_int, *mut c_int);
    let (invrow, invcol): (*mut c_int, *mut c_int);
    let mut set: c_int;
    if (*data).partitionexist == 0 {
        return 1;
    }
    CreateInverseTopology(data, info);
    invrow = (*data).invtopo.rows;
    invcol = (*data).invtopo.cols;
    noknots = (*data).noknots;
    nopartitions = (*data).nopartitions;
    elempart = (*data).elempart;
    nodepart = (*data).nodepart;
    if info != 0 {
        printf(
            c"Number of nodal partitions: %d\n".as_ptr() as *const c_char,
            nopartitions,
        );
    }
    nodesinpart = ivector(1, nopartitions);
    cuminpart = ivector(1, nopartitions);
    j = 1;
    while j <= nopartitions {
        *cuminpart.wrapping_offset(j as isize) = 0;
        j += 1;
    }
    knows = imatrix(1, nopartitions, 1, nopartitions);
    cumknows = imatrix(1, nopartitions, 1, nopartitions);
    i = 1;
    while i <= nopartitions {
        j = 1;
        while j <= nopartitions {
            *(*knows.wrapping_offset(i as isize)).wrapping_offset(j as isize) = 0;
            *(*cumknows.wrapping_offset(i as isize)).wrapping_offset(j as isize) = 0;
            j += 1;
        }
        i += 1;
    }
    set = FALSE;
    // NOTE: `omstart:` backward goto becomes a two-pass loop.
    loop {
        i = 1;
        while i <= noknots {
            j = 1;
            while j <= nopartitions {
                *nodesinpart.wrapping_offset(j as isize) = 0;
                j += 1;
            }
            j = *invrow.wrapping_offset((i - 1) as isize);
            while j < *invrow.wrapping_offset(i as isize) {
                k = *invcol.wrapping_offset(j as isize) + 1;
                part = *elempart.wrapping_offset(k as isize);
                if part > 0 {
                    *nodesinpart.wrapping_offset(part as isize) += 1;
                }
                j += 1;
            }
            maxpart = 0;
            maxpart2 = 0;
            j = 1;
            while j <= nopartitions {
                if *nodesinpart.wrapping_offset(j as isize) > 0
                    && (maxpart == 0
                        || *nodesinpart.wrapping_offset(j as isize)
                            > *nodesinpart.wrapping_offset(maxpart as isize))
                {
                    maxpart = j;
                }
                j += 1;
            }
            j = 1;
            while j <= nopartitions {
                if *nodesinpart.wrapping_offset(j as isize) > 0 {
                    if j == maxpart {
                        j += 1;
                        continue;
                    }
                    if maxpart2 == 0
                        || *nodesinpart.wrapping_offset(j as isize)
                            > *nodesinpart.wrapping_offset(maxpart2 as isize)
                    {
                        maxpart2 = j;
                    }
                }
                j += 1;
            }
            if maxpart == 0 && maxpart2 == 0 {
                printf(
                    c"Node is not included in any partition: %d\n".as_ptr() as *const c_char,
                    i,
                );
            } else if maxpart2 == 0 {
                *nodepart.wrapping_offset(i as isize) = maxpart;
                *cuminpart.wrapping_offset(maxpart as isize) += 1;
            } else if *nodesinpart.wrapping_offset(maxpart as isize)
                > *nodesinpart.wrapping_offset(maxpart2 as isize)
            {
                if set != 0 {
                    *nodepart.wrapping_offset(i as isize) = maxpart;
                } else {
                    *cuminpart.wrapping_offset(maxpart as isize) += 1;
                }
            } else {
                if set != 0 {
                    *(*cumknows.wrapping_offset(maxpart as isize))
                        .wrapping_offset(maxpart2 as isize) += 1;
                    if *(*cumknows.wrapping_offset(maxpart as isize))
                        .wrapping_offset(maxpart2 as isize)
                        > *(*knows.wrapping_offset(maxpart as isize))
                            .wrapping_offset(maxpart2 as isize)
                            / 2
                    {
                        *nodepart.wrapping_offset(i as isize) = maxpart2;
                        *cuminpart.wrapping_offset(maxpart2 as isize) += 1;
                    } else {
                        *nodepart.wrapping_offset(i as isize) = maxpart;
                        *cuminpart.wrapping_offset(maxpart as isize) += 1;
                    }
                } else {
                    *(*knows.wrapping_offset(maxpart as isize))
                        .wrapping_offset(maxpart2 as isize) += 1;
                }
            }
            i += 1;
        }
        if set != 0 {
            break;
        }
        set = TRUE;
    }
    minpart = *cuminpart.wrapping_offset(1);
    maxpart = *cuminpart.wrapping_offset(1);
    j = 1;
    while j <= nopartitions {
        minpart = c_min_int(minpart, *cuminpart.wrapping_offset(j as isize));
        maxpart = c_max_int(maxpart, *cuminpart.wrapping_offset(j as isize));
        j += 1;
    }
    if info != 0 {
        printf(
            c"Set the node partitions by the dominating element partition.\n".as_ptr()
                as *const c_char,
        );
        printf(
            c"There are from %d to %d nodes in the %d partitions.\n".as_ptr() as *const c_char,
            minpart,
            maxpart,
            nopartitions,
        );
    }
    free_ivector(nodesinpart, 1, nopartitions);
    free_ivector(cuminpart, 1, nopartitions);
    free_imatrix(knows, 1, nopartitions, 1, nopartitions);
    free_imatrix(cumknows, 1, nopartitions, 1, nopartitions);
    0
}
/// # Safety
///
/// recursive geometric partitioning over valid meshes (C contract).
/// `goto skippart` becomes an if/else split (the skipped region frees nothing
/// on that path in C either — same leak, same behavior).
#[no_mangle]
pub unsafe extern "C" fn PartitionSimpleElements(
    data: *mut FemType,
    eg: *mut ElmergridType,
    bound: *mut BoundaryType,
    dimpart: *mut c_int,
    dimper: *mut c_int,
    partorder: c_int,
    corder: *mut Real,
    parttol: Real,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k, mut l): (c_int, c_int, c_int, c_int);
    let (mut ind, mut minpart, mut maxpart): (c_int, c_int, c_int);
    let (noknots, noelements, mut nonodes, mut elemsinpart, periodic): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (partitions1, partitions2, mut partitions3, mut partitions): (c_int, c_int, c_int, c_int);
    let (mut vpartitions1, mut vpartitions2, mut vpartitions3, vpartitions): (
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (mut noelements0, noelements1, mut noparts0): (c_int, c_int, c_int);
    // NOTE: NULL-init in C; assigned iff a connect mesh exists, read under
    // the same flag (never reassigned mid-function) — unobservable init.
    let mut elemconnect: *mut c_int = core::ptr::null_mut();
    let (indx, nopart, inpart): (*mut c_int, *mut c_int, *mut c_int);
    let arrange: *mut Real;
    let (mut x, mut y, mut z, cx, cy, cz): (Real, Real, Real, Real, Real, Real);
    if info != 0 {
        printf(c"PartitionSimpleElements\n".as_ptr() as *const c_char);
    }
    noelements = (*data).noelements;
    noknots = (*data).noknots;
    partitions1 = *dimpart.wrapping_offset(0);
    partitions2 = *dimpart.wrapping_offset(1);
    partitions3 = *dimpart.wrapping_offset(2);
    if (*data).dim < 3 {
        partitions3 = 1;
    }
    partitions = partitions1 * partitions2 * partitions3;
    if partitions1 < 2 && partitions2 < 2 && partitions3 < 2 && (*eg).connect == 0 {
        printf(
            c"Some of the divisions must be larger than one: %d %d %d\n".as_ptr()
                as *const c_char,
            partitions1,
            partitions2,
            partitions3,
        );
        bigerror(c"Partitioning not performed".as_ptr() as *const c_char);
    }
    if partitions >= noelements {
        printf(
            c"There must be fever partitions than elements (%d vs %d)!\n".as_ptr()
                as *const c_char,
            partitions,
            noelements,
        );
        bigerror(c"Partitioning not performed".as_ptr() as *const c_char);
    }
    if (*eg).partbcz > 1 || (*eg).partbcr != 0 {
        PartitionConnectedElements1D(data, bound, eg, info);
    } else if (*eg).partbcmetis > 1 {
        PartitionConnectedElementsMetis(data, bound, eg, (*eg).partbcmetis, 3, info);
    } else if (*eg).connect != 0 {
        PartitionConnectedElementsStraight(data, bound, eg, info);
    }
    if (*data).nodeconnectexist != 0 {
        ExtendBoundaryPartitioning(data, bound, (*eg).partbclayers, info);
    }
    if (*data).partitionexist == 0 {
        (*data).partitionexist = TRUE;
        (*data).elempart = ivector(1, noelements);
        (*data).nodepart = ivector(1, noknots);
        (*data).nopartitions = partitions;
    }
    inpart = (*data).elempart;
    // NOTE: C stores vpartitions{1,2,3} = partitions{1,2,3} here and again
    // below; the first stores are never read and dropped.
    // NOTE: `if(0)` diagnostic omitted (runtime-dead).
    noparts0 = 0;
    noelements0 = 0;
    if (*data).elemconnectexist != 0 || (*data).nodeconnectexist != 0 {
        elemconnect = (*data).elemconnect;
        noparts0 = 0;
        i = 1;
        while i <= (*data).noelements {
            if *elemconnect.wrapping_offset(i as isize) != 0 {
                noelements0 += 1;
            }
            noparts0 = c_max_int(noparts0, *elemconnect.wrapping_offset(i as isize));
            i += 1;
        }
        if info != 0 {
            printf(
                c"There are %d initial partitions in the connected mesh\n".as_ptr()
                    as *const c_char,
                noparts0,
            );
        }
        if info != 0 {
            printf(
                c"There are %d initial elements in the connected mesh\n".as_ptr()
                    as *const c_char,
                noelements0,
            );
        }
    }
    noelements1 = noelements - noelements0;
    vpartitions1 = partitions1;
    vpartitions2 = partitions2;
    vpartitions3 = partitions3;
    periodic = *dimper.wrapping_offset(0) | *dimper.wrapping_offset(1) | *dimper.wrapping_offset(2);
    if periodic != 0 {
        if *dimper.wrapping_offset(0) != 0 && partitions1 > 1 {
            vpartitions1 *= 2;
        }
        if *dimper.wrapping_offset(1) != 0 && partitions2 > 1 {
            vpartitions2 *= 2;
        }
        if *dimper.wrapping_offset(2) != 0 && partitions3 > 1 {
            vpartitions3 *= 2;
        }
    }
    vpartitions = vpartitions1 * vpartitions2 * vpartitions3;
    nopart = ivector(1, vpartitions + noparts0);
    // NOTE: `goto skippart` — the hybrid fast path skips the whole
    // partitioning region below (frees included, as in C).
    if vpartitions == 1 && (*data).elemconnectexist != 0 {
        if info != 0 {
            printf(
                c"Only one regular partitions requested, skipping 2nd part of hybrid partitioning\n"
                    .as_ptr() as *const c_char,
            );
        }
        j = 1;
        while j <= noelements {
            if *elemconnect.wrapping_offset(j as isize) > 0 {
                *inpart.wrapping_offset(j as isize) = *elemconnect.wrapping_offset(j as isize);
            } else {
                *inpart.wrapping_offset(j as isize) = noparts0 + 1;
            }
            j += 1;
        }
        partitions = noparts0 + 1;
        (*data).nopartitions = partitions;
    } else {
        if info != 0 {
            printf(
                c"Making a simple partitioning for %d elements in %d-dimensions.\n".as_ptr()
                    as *const c_char,
                noelements1,
                (*data).dim,
            );
        }
        arrange = dvector(1, noelements);
        indx = ivector(1, noelements);
        if partorder != 0 {
            cx = *corder.wrapping_offset(0);
            cy = *corder.wrapping_offset(1);
            cz = *corder.wrapping_offset(2);
        } else {
            cx = 1.0;
            cy = 0.0001;
            cz = cy * cy;
        }
        // NOTE: C stores z = 0.0 here; every read is dominated by a later
        // store, so it is dropped for warning-freedom.
        i = 1;
        while i <= noelements {
            *inpart.wrapping_offset(i as isize) = 1;
            i += 1;
        }
        if vpartitions1 > 1 {
            if info != 0 {
                printf(
                    c"Ordering 1st direction with (%.3lg*x + %.3lg*y + %.3lg*z)\n".as_ptr()
                        as *const c_char,
                    cx,
                    cy,
                    cz,
                );
            }
            j = 1;
            while j <= noelements {
                if (*data).elemconnectexist != 0
                    && *elemconnect.wrapping_offset(j as isize) > 0 {
                        *arrange.wrapping_offset(j as isize) = 1.0e9;
                        j += 1;
                        continue;
                    }
                nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
                x = 0.0;
                y = 0.0;
                z = 0.0;
                i = 0;
                while i < nonodes {
                    k = *(*(*data).topology.wrapping_offset(j as isize))
                        .wrapping_offset(i as isize);
                    x += *(*data).x.wrapping_offset(k as isize);
                    y += *(*data).y.wrapping_offset(k as isize);
                    z += *(*data).z.wrapping_offset(k as isize);
                    i += 1;
                }
                *arrange.wrapping_offset(j as isize) = (cx * x + cy * y + cz * z) / nonodes as Real;
                j += 1;
            }
            SortIndex(noelements, arrange, indx);
            i = 1;
            while i <= noelements1 {
                ind = *indx.wrapping_offset(i as isize);
                k = (i * vpartitions1 - 1) / noelements1 + 1;
                *inpart.wrapping_offset(ind as isize) = k;
                i += 1;
            }
        }
        if vpartitions2 > 1 {
            if info != 0 {
                printf(c"Ordering in the 2nd direction.\n".as_ptr() as *const c_char);
            }
            j = 1;
            while j <= noelements {
                if (*data).elemconnectexist != 0
                    && *elemconnect.wrapping_offset(j as isize) > 0 {
                        *arrange.wrapping_offset(j as isize) = 1.0e9;
                        j += 1;
                        continue;
                    }
                nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
                x = 0.0;
                y = 0.0;
                z = 0.0;
                i = 0;
                while i < nonodes {
                    k = *(*(*data).topology.wrapping_offset(j as isize))
                        .wrapping_offset(i as isize);
                    x += *(*data).x.wrapping_offset(k as isize);
                    y += *(*data).y.wrapping_offset(k as isize);
                    z += *(*data).z.wrapping_offset(k as isize);
                    i += 1;
                }
                *arrange.wrapping_offset(j as isize) =
                    (-cy * x + cx * y + cz * z) / nonodes as Real;
                j += 1;
            }
            SortIndex(noelements, arrange, indx);
            i = 1;
            while i <= vpartitions {
                *nopart.wrapping_offset(i as isize) = 0;
                i += 1;
            }
            elemsinpart = noelements1 / (vpartitions1 * vpartitions2);
            i = 1;
            while i <= noelements1 {
                j = 0;
                ind = *indx.wrapping_offset(i as isize);
                loop {
                    j += 1;
                    k = (*inpart.wrapping_offset(ind as isize) - 1) * vpartitions2 + j;
                    if !(*nopart.wrapping_offset(k as isize) >= elemsinpart && j < vpartitions2) {
                        break;
                    }
                }
                *nopart.wrapping_offset(k as isize) += 1;
                *inpart.wrapping_offset(ind as isize) =
                    (*inpart.wrapping_offset(ind as isize) - 1) * vpartitions2 + j;
                i += 1;
            }
        }
        if vpartitions3 > 1 {
            if info != 0 {
                printf(c"Ordering in the 3rd direction.\n".as_ptr() as *const c_char);
            }
            j = 1;
            while j <= noelements {
                if (*data).elemconnectexist != 0
                    && *elemconnect.wrapping_offset(j as isize) > 0 {
                        *arrange.wrapping_offset(j as isize) = 1.0e9;
                        j += 1;
                        continue;
                    }
                nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
                x = 0.0;
                y = 0.0;
                z = 0.0;
                i = 0;
                while i < nonodes {
                    k = *(*(*data).topology.wrapping_offset(j as isize))
                        .wrapping_offset(i as isize);
                    x += *(*data).x.wrapping_offset(k as isize);
                    y += *(*data).y.wrapping_offset(k as isize);
                    z += *(*data).z.wrapping_offset(k as isize);
                    i += 1;
                }
                *arrange.wrapping_offset(j as isize) =
                    (-cz * x - cy * y + cx * z) / nonodes as Real;
                j += 1;
            }
            SortIndex(noelements, arrange, indx);
            i = 1;
            while i <= vpartitions {
                *nopart.wrapping_offset(i as isize) = 0;
                i += 1;
            }
            elemsinpart = noelements1 / (vpartitions1 * vpartitions2 * vpartitions3);
            i = 1;
            while i <= noelements1 {
                j = 0;
                ind = *indx.wrapping_offset(i as isize);
                loop {
                    j += 1;
                    k = (*inpart.wrapping_offset(ind as isize) - 1) * vpartitions3 + j;
                    if !(*nopart.wrapping_offset(k as isize) >= elemsinpart && j < vpartitions3) {
                        break;
                    }
                }
                *nopart.wrapping_offset(k as isize) += 1;
                *inpart.wrapping_offset(ind as isize) =
                    (*inpart.wrapping_offset(ind as isize) - 1) * vpartitions3 + j;
                i += 1;
            }
        }
        if parttol > 1.0e-20 {
            // NOTE: `if(0)` count dumps omitted (runtime-dead).
            j = 1;
            while j <= noelements {
                if (*data).elemconnectexist != 0
                    && *elemconnect.wrapping_offset(j as isize) > 0 {
                        *arrange.wrapping_offset(j as isize) = 1.0e9;
                        j += 1;
                        continue;
                    }
                nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
                x = 0.0;
                y = 0.0;
                z = 0.0;
                i = 0;
                while i < nonodes {
                    k = *(*(*data).topology.wrapping_offset(j as isize))
                        .wrapping_offset(i as isize);
                    x += *(*data).x.wrapping_offset(k as isize);
                    y += *(*data).y.wrapping_offset(k as isize);
                    z += *(*data).z.wrapping_offset(k as isize);
                    i += 1;
                }
                *arrange.wrapping_offset(j as isize) = (cx * x + cy * y + cz * z) / nonodes as Real;
                j += 1;
            }
            SortIndex(noelements, arrange, indx);
            // NOTE: kprev/l are stride-block locals (kprev needs its own
            // binding; l reuses the function loop counter as in C).
            let mut kprev: c_int = *inpart.wrapping_offset(*indx.wrapping_offset(1) as isize);
            l = 0;
            i = 2;
            while i <= noelements1 {
                ind = *indx.wrapping_offset(i as isize);
                k = *inpart.wrapping_offset(ind as isize);
                if k != kprev
                    && fabs(
                        *arrange.wrapping_offset(i as isize)
                            - *arrange.wrapping_offset((i - 1) as isize),
                    ) < parttol
                    {
                        *nopart.wrapping_offset(k as isize) -= 1;
                        *nopart.wrapping_offset(kprev as isize) += 1;
                        *inpart.wrapping_offset(ind as isize) = kprev;
                        k = kprev;
                        // NOTE: `if(0)` diagnostic omitted (runtime-dead).
                        l += 1;
                    }
                kprev = k;
                i += 1;
            }
            if info != 0 {
                printf(
                    c"Number of partition changes due to tolerance: %d\n".as_ptr()
                        as *const c_char,
                    l,
                );
            }
        }
        free_dvector(arrange, 1, noelements);
        free_ivector(indx, 1, noelements);
        if (*data).elemconnectexist != 0 {
            j = 1;
            while j <= noelements {
                if *elemconnect.wrapping_offset(j as isize) > 0 {
                    *inpart.wrapping_offset(j as isize) = *elemconnect.wrapping_offset(j as isize);
                } else {
                    *inpart.wrapping_offset(j as isize) =
                        *inpart.wrapping_offset(j as isize) + noparts0;
                }
                j += 1;
            }
            partitions += noparts0;
            (*data).nopartitions = partitions;
        }
        if periodic != 0 {
            // NOTE: block-scoped C locals become a Rust block.
            {
                
                let (mut p1, mut p2, mut p3, mut q1, mut q2, mut q3): (
                    c_int,
                    c_int,
                    c_int,
                    c_int,
                    c_int,
                    c_int,
                );
                // NOTE: P/Q match the C locals (allowed by non_snake_case).
                let (mut P, mut Q): (c_int, c_int);
                if (*data).elemconnectexist != 0 {
                    bigerror(
                        c"Cannot use connect flag with periodic systems\n".as_ptr()
                            as *const c_char,
                    );
                }
                // NOTE: C pre-inits p1=p2=p3=1, but the for-loops below
                // reassign first; the stores are dead and dropped.
                let partmap: *mut c_int = ivector(1, vpartitions);
                i = 1;
                while i <= vpartitions {
                    *partmap.wrapping_offset(i as isize) = 0;
                    i += 1;
                }
                p1 = 1;
                while p1 <= vpartitions1 {
                    q1 = p1;
                    if *dimper.wrapping_offset(0) != 0 && vpartitions1 > 1 {
                        if q1 == vpartitions1 {
                            q1 = 0;
                        }
                        q1 = q1 / 2 + 1;
                    }
                    p2 = 1;
                    while p2 <= vpartitions2 {
                        q2 = p2;
                        if *dimper.wrapping_offset(1) != 0 && vpartitions2 > 1 {
                            if q2 == vpartitions2 {
                                q2 = 0;
                            }
                            q2 = q2 / 2 + 1;
                        }
                        p3 = 1;
                        while p3 <= vpartitions3 {
                            q3 = p3;
                            if *dimper.wrapping_offset(2) != 0 && vpartitions3 > 1 {
                                if q3 == vpartitions3 {
                                    q3 = 0;
                                }
                                q3 = q3 / 2 + 1;
                            }
                            P = vpartitions3 * vpartitions2 * (p1 - 1)
                                + vpartitions3 * (p2 - 1)
                                + p3;
                            Q = partitions3 * partitions2 * (q1 - 1) + partitions3 * (q2 - 1) + q3;
                            *partmap.wrapping_offset(P as isize) = Q;
                            p3 += 1;
                        }
                        p2 += 1;
                    }
                    p1 += 1;
                }
                i = 1;
                while i <= noelements {
                    *inpart.wrapping_offset(i as isize) =
                        *partmap.wrapping_offset(*inpart.wrapping_offset(i as isize) as isize);
                    i += 1;
                }
                free_ivector(partmap, 1, vpartitions);
            }
        }
    }
    // NOTE: `skippart:` label target.
    i = 1;
    while i <= partitions {
        *nopart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    i = 1;
    while i <= noelements {
        *nopart.wrapping_offset(*inpart.wrapping_offset(i as isize) as isize) += 1;
        i += 1;
    }
    minpart = *nopart.wrapping_offset(1);
    maxpart = *nopart.wrapping_offset(1);
    i = 1;
    while i <= partitions {
        minpart = c_min_int(minpart, *nopart.wrapping_offset(i as isize));
        maxpart = c_max_int(maxpart, *nopart.wrapping_offset(i as isize));
        i += 1;
    }
    free_ivector(nopart, 1, vpartitions);
    PartitionNodesByElements(data, info);
    if info != 0 {
        printf(
            c"Successfully made a partitioning with %d to %d elements.\n".as_ptr()
                as *const c_char,
            minpart,
            maxpart,
        );
    }
    0
}
/// # Safety
///
/// non-recursive geometric partitioning over valid meshes (C contract).
#[no_mangle]
pub unsafe extern "C" fn PartitionSimpleElementsNonRecursive(
    data: *mut FemType,
    dimpart: *mut c_int,
    dimper: *mut c_int,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k): (c_int, c_int, c_int);
    let mut minpart: c_int;
    let mut maxpart: c_int;
    let (noknots, noelements, mut nonodes, periodic): (c_int, c_int, c_int, c_int);
    let (partitions1, partitions2, mut partitions3, partitions): (c_int, c_int, c_int, c_int);
    let (mut indx1, mut indy1, mut indz1): (c_int, c_int, c_int);
    let (nopart, inpart): (*mut c_int, *mut c_int);
    let (mut x, mut y, mut z, mut maxx, mut minx, mut maxy, mut miny, mut maxz, mut minz): (
        Real,
        Real,
        Real,
        Real,
        Real,
        Real,
        Real,
        Real,
        Real,
    );
    partitions1 = *dimpart.wrapping_offset(0);
    partitions2 = *dimpart.wrapping_offset(1);
    partitions3 = *dimpart.wrapping_offset(2);
    if (*data).dim < 3 {
        partitions3 = 1;
    }
    partitions = partitions1 * partitions2 * partitions3;
    if partitions1 < 2 && partitions2 < 2 && partitions3 < 2 {
        printf(
            c"Some of the divisions must be larger than one: %d %d %d\n".as_ptr()
                as *const c_char,
            partitions1,
            partitions2,
            partitions3,
        );
        bigerror(c"Partitioning not performed".as_ptr() as *const c_char);
    }
    if (*data).partitionexist == 0 {
        (*data).partitionexist = TRUE;
        (*data).elempart = ivector(1, (*data).noelements);
        (*data).nodepart = ivector(1, (*data).noknots);
        (*data).nopartitions = partitions;
    }
    inpart = (*data).elempart;
    noelements = (*data).noelements;
    noknots = (*data).noknots;
    periodic = *dimper.wrapping_offset(0) | *dimper.wrapping_offset(1) | *dimper.wrapping_offset(2);
    if periodic != 0 {
        printf(
            c"Implement periodicity for this partitioning routine!\n".as_ptr() as *const c_char,
        );
    }
    if info != 0 {
        printf(
            c"Making a simple partitioning for %d elements in %d-dimensions.\n".as_ptr()
                as *const c_char,
            noelements,
            (*data).dim,
        );
        printf(
            c"There can be at maximum %d partitions\n".as_ptr() as *const c_char,
            partitions,
        );
    }
    nopart = ivector(1, partitions);
    i = 1;
    while i <= partitions {
        *nopart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    // NOTE: C stores z = 0.0 here; dropped (see SimpleElements note).
    // NOTE: C initializes IndZ (not IndX/IndY); indz1 is read unconditionally.
    indz1 = 1;
    i = 1;
    while i <= noelements {
        *inpart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    maxx = *(*data).x.wrapping_offset(1);
    minx = *(*data).x.wrapping_offset(1);
    maxy = *(*data).y.wrapping_offset(1);
    miny = *(*data).y.wrapping_offset(1);
    maxz = *(*data).z.wrapping_offset(1);
    minz = *(*data).z.wrapping_offset(1);
    i = 1;
    while i <= noknots {
        x = *(*data).x.wrapping_offset(i as isize);
        y = *(*data).y.wrapping_offset(i as isize);
        z = *(*data).z.wrapping_offset(i as isize);
        maxx = c_max_real(maxx, x);
        minx = c_min_real(minx, x);
        maxy = c_max_real(maxy, y);
        miny = c_min_real(miny, y);
        maxz = c_max_real(maxz, z);
        minz = c_min_real(minz, z);
        i += 1;
    }
    if info != 0 {
        printf(
            c"Range in x-direction: %12.5le %12.5le\n".as_ptr() as *const c_char,
            minx,
            maxx,
        );
        printf(
            c"Range in y-direction: %12.5le %12.5le\n".as_ptr() as *const c_char,
            miny,
            maxy,
        );
        printf(
            c"Range in z-direction: %12.5le %12.5le\n".as_ptr() as *const c_char,
            minz,
            maxz,
        );
    }
    j = 1;
    while j <= noelements {
        nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
        x = 0.0;
        y = 0.0;
        z = 0.0;
        i = 0;
        while i < nonodes {
            k = *(*(*data).topology.wrapping_offset(j as isize)).wrapping_offset(i as isize);
            x += *(*data).x.wrapping_offset(k as isize);
            y += *(*data).y.wrapping_offset(k as isize);
            z += *(*data).z.wrapping_offset(k as isize);
            i += 1;
        }
        x /= nonodes as Real;
        y /= nonodes as Real;
        z /= nonodes as Real;
        indx1 = ceil(partitions1 as Real * (x - minx) / (maxx - minx)) as c_int;
        indy1 = ceil(partitions2 as Real * (y - miny) / (maxy - miny)) as c_int;
        if (*data).dim == 3 {
            indz1 = ceil(partitions3 as Real * (z - minz) / (maxz - minz)) as c_int;
        }
        k = (indz1 - 1) * partitions1 * partitions2 + (indy1 - 1) * partitions1 + indx1;
        *inpart.wrapping_offset(j as isize) = k;
        *nopart.wrapping_offset(k as isize) += 1;
        j += 1;
    }
    maxpart = 0;
    minpart = noelements;
    j = 0;
    if info != 0 {
        printf(c"Renumbering the partitions (showing only 64):\n".as_ptr() as *const c_char);
    }
    i = 1;
    while i <= partitions {
        k = *nopart.wrapping_offset(i as isize);
        if k != 0 {
            maxpart = c_max_int(k, maxpart);
            minpart = c_min_int(k, minpart);
            j += 1;
            *nopart.wrapping_offset(i as isize) = j;
            if info != 0 && j <= 64 {
                printf(c"%d -> %d\n".as_ptr() as *const c_char, i, j);
            }
        }
        i += 1;
    }
    if info != 0 {
        printf(
            c"There are %d active partitions out of %d possible ones\n".as_ptr() as *const c_char,
            j,
            partitions,
        );
    }
    (*data).nopartitions = j;
    i = 1;
    while i <= noelements {
        j = *inpart.wrapping_offset(i as isize);
        *inpart.wrapping_offset(i as isize) = *nopart.wrapping_offset(j as isize);
        i += 1;
    }
    free_ivector(nopart, 1, partitions);
    PartitionNodesByElements(data, info);
    if info != 0 {
        printf(
            c"Successfully made a partitioning with %d to %d elements.\n".as_ptr()
                as *const c_char,
            minpart,
            maxpart,
        );
    }
    0
}
// NOTE: C function-local MAXCATEGORY as a Rust const.
pub const MAXCATEGORY: c_int = 500;

/// # Safety
///
/// rotational geometric partitioning over valid meshes (C contract).
#[no_mangle]
pub unsafe extern "C" fn PartitionSimpleElementsRotational(
    data: *mut FemType,
    dimpart: *mut c_int,
    dimper: *mut c_int,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k): (c_int, c_int, c_int);
    let (mut minpart, mut maxpart, dim, mut hit): (c_int, c_int, c_int, c_int);
    let (noknots, noelements, mut nonodes, periodic): (c_int, c_int, c_int, c_int);
    let (partitions1, partitions2, mut partitions3, mut partitions): (c_int, c_int, c_int, c_int);
    let (mut indr, mut indf, mut indz, mut connect): (c_int, c_int, c_int, c_int);
    let (nopart, inpart, cumf, cumr): (*mut c_int, *mut c_int, *mut c_int, *mut c_int);
    // NOTE: NULL-init in C; assigned iff dim==3, read under the same flag
    // (never reassigned mid-function) — unobservable init.
    let mut cumz: *mut c_int = core::ptr::null_mut();
    // NOTE: NULL-init in C; assigned iff nodeconnectexist, read under the
    // same flag — unobservable init.
    let mut elemconnect: *mut c_int = core::ptr::null_mut();
    let (mut x, mut y, mut z, mut r, mut f): (Real, Real, Real, Real, Real);
    let (mut maxr, mut minr, mut maxf, mut minf, mut maxz, mut minz): (
        Real,
        Real,
        Real,
        Real,
        Real,
        Real,
    );
    dim = (*data).dim;
    partitions1 = *dimpart.wrapping_offset(0);
    partitions2 = *dimpart.wrapping_offset(1);
    partitions3 = *dimpart.wrapping_offset(2);
    if dim < 3 {
        partitions3 = 1;
    }
    partitions = partitions1 * partitions2 * partitions3;
    connect = FALSE;
    if (*data).nodeconnectexist != 0 {
        SetConnectedElements(data, info);
        connect = TRUE;
        elemconnect = (*data).elemconnect;
        partitions += partitions3;
    }
    if partitions1 < 2 && partitions2 < 2 && partitions3 < 2 {
        printf(
            c"Some of the divisions must be larger than one: %d %d %d\n".as_ptr()
                as *const c_char,
            partitions1,
            partitions2,
            partitions3,
        );
        bigerror(c"Partitioning not performed".as_ptr() as *const c_char);
    }
    if (*data).partitionexist == 0 {
        (*data).partitionexist = TRUE;
        (*data).elempart = ivector(1, (*data).noelements);
        (*data).nodepart = ivector(1, (*data).noknots);
        (*data).nopartitions = partitions;
    }
    inpart = (*data).elempart;
    noelements = (*data).noelements;
    noknots = (*data).noknots;
    periodic = *dimper.wrapping_offset(0) | *dimper.wrapping_offset(1) | *dimper.wrapping_offset(2);
    if periodic != 0 {
        printf(
            c"Implement periodicity for this partitioning routine!\n".as_ptr() as *const c_char,
        );
    }
    if info != 0 {
        printf(
            c"Making a simple rotational partitioning for %d elements in %d-dimensions.\n"
                .as_ptr() as *const c_char,
            noelements,
            (*data).dim,
        );
        printf(
            c"There can be at maximum %d partitions\n".as_ptr() as *const c_char,
            partitions,
        );
    }
    nopart = ivector(1, partitions);
    i = 1;
    while i <= partitions {
        *nopart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    // NOTE: C stores z = 0.0 here; dropped (see SimpleElements note).
    indz = 1;
    i = 1;
    while i <= noelements {
        *inpart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    x = *(*data).x.wrapping_offset(1);
    y = *(*data).y.wrapping_offset(1);
    z = *(*data).z.wrapping_offset(1);
    r = sqrt(x * x + y * y);
    f = 180.0 * atan2(y, x) / FM_PI;
    if f < 0.0 {
        f += 360.0;
    }
    maxr = r;
    minr = r;
    maxf = f;
    minf = f;
    maxz = z;
    minz = z;
    i = 1;
    while i <= noknots {
        x = *(*data).x.wrapping_offset(i as isize);
        y = *(*data).y.wrapping_offset(i as isize);
        z = *(*data).z.wrapping_offset(i as isize);
        r = sqrt(x * x + y * y);
        f = 180.0 * atan2(y, x) / FM_PI;
        if f < 0.0 {
            f += 360.0;
        }
        maxr = c_max_real(maxr, r);
        minr = c_min_real(minr, r);
        maxf = c_max_real(maxf, f);
        minf = c_min_real(minf, f);
        maxz = c_max_real(maxz, z);
        minz = c_min_real(minz, z);
        i += 1;
    }
    if info != 0 {
        printf(
            c"Range in r-direction: %12.5le %12.5le\n".as_ptr() as *const c_char,
            minr,
            maxr,
        );
        printf(
            c"Range in f-direction: %12.5le %12.5le\n".as_ptr() as *const c_char,
            minf,
            maxf,
        );
        printf(
            c"Range in z-direction: %12.5le %12.5le\n".as_ptr() as *const c_char,
            minz,
            maxz,
        );
    }
    if maxf - minf > 180.0 {
        maxf = 360.0;
        minf = 0.0;
    }
    cumr = ivector(0, MAXCATEGORY);
    cumf = ivector(0, MAXCATEGORY);
    if dim == 3 {
        cumz = ivector(0, MAXCATEGORY);
    }
    i = 0;
    while i <= MAXCATEGORY {
        *cumf.wrapping_offset(i as isize) = 0;
        *cumr.wrapping_offset(i as isize) = 0;
        if dim == 3 {
            *cumz.wrapping_offset(i as isize) = 0;
        }
        i += 1;
    }
    j = 1;
    while j <= noelements {
        nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
        x = 0.0;
        y = 0.0;
        z = 0.0;
        i = 0;
        while i < nonodes {
            k = *(*(*data).topology.wrapping_offset(j as isize)).wrapping_offset(i as isize);
            x += *(*data).x.wrapping_offset(k as isize);
            y += *(*data).y.wrapping_offset(k as isize);
            z += *(*data).z.wrapping_offset(k as isize);
            i += 1;
        }
        x /= nonodes as Real;
        y /= nonodes as Real;
        z /= nonodes as Real;
        r = sqrt(x * x + y * y);
        f = 180.0 * atan2(y, x) / FM_PI;
        if f < 0.0 {
            f += 360.0;
        }
        indr = ceil(MAXCATEGORY as Real * (r - minr) / (maxr - minr)) as c_int;
        indf = ceil(MAXCATEGORY as Real * (f - minf) / (maxf - minf)) as c_int;
        if dim == 3 {
            indz = ceil(MAXCATEGORY as Real * (z - minz) / (maxz - minz)) as c_int;
        }
        if !(1..=MAXCATEGORY).contains(&indr) {
            printf(
                c"IndR out of bounds : %d\n".as_ptr() as *const c_char,
                indr,
            );
            indr = c_min_int(c_max_int(indr, 1), MAXCATEGORY);
        }
        if !(1..=MAXCATEGORY).contains(&indf) {
            printf(
                c"IndF out of bounds : %d\n".as_ptr() as *const c_char,
                indf,
            );
            indf = c_min_int(c_max_int(indf, 1), MAXCATEGORY);
        }
        if dim == 3 {
            if !(1..=MAXCATEGORY).contains(&indz) {
                printf(
                    c"IndZ out of bounds : %d\n".as_ptr() as *const c_char,
                    indz,
                );
                indz = c_min_int(c_max_int(indz, 1), MAXCATEGORY);
            }
            *cumz.wrapping_offset(indz as isize) += 1;
        }
        if connect != 0 && *elemconnect.wrapping_offset(j as isize) < 0 {
            j += 1;
            continue;
        }
        *cumr.wrapping_offset(indr as isize) += 1;
        *cumf.wrapping_offset(indf as isize) += 1;
        j += 1;
    }
    i = 1;
    while i <= MAXCATEGORY {
        *cumr.wrapping_offset(i as isize) =
            *cumr.wrapping_offset(i as isize) + *cumr.wrapping_offset((i - 1) as isize);
        *cumf.wrapping_offset(i as isize) =
            *cumf.wrapping_offset(i as isize) + *cumf.wrapping_offset((i - 1) as isize);
        if dim == 3 {
            *cumz.wrapping_offset(i as isize) =
                *cumz.wrapping_offset(i as isize) + *cumz.wrapping_offset((i - 1) as isize);
        }
        i += 1;
    }
    j = noelements - (*data).elemconnectexist;
    i = 1;
    while i <= MAXCATEGORY {
        *cumr.wrapping_offset(i as isize) = ceil(
            1.0 * partitions1 as Real * *cumr.wrapping_offset(i as isize) as Real
                / noelements as Real,
        ) as c_int;
        *cumf.wrapping_offset(i as isize) =
            ceil(1.0 * partitions2 as Real * *cumf.wrapping_offset(i as isize) as Real / j as Real)
                as c_int;
        if dim == 3 {
            *cumz.wrapping_offset(i as isize) = ceil(
                1.0 * partitions3 as Real * *cumz.wrapping_offset(i as isize) as Real / j as Real,
            ) as c_int;
        }
        i += 1;
    }
    j = 1;
    while j <= noelements {
        nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
        x = 0.0;
        y = 0.0;
        z = 0.0;
        i = 0;
        while i < nonodes {
            k = *(*(*data).topology.wrapping_offset(j as isize)).wrapping_offset(i as isize);
            x += *(*data).x.wrapping_offset(k as isize);
            y += *(*data).y.wrapping_offset(k as isize);
            z += *(*data).z.wrapping_offset(k as isize);
            i += 1;
        }
        x /= nonodes as Real;
        y /= nonodes as Real;
        z /= nonodes as Real;
        r = sqrt(x * x + y * y);
        f = 180.0 * atan2(y, x) / FM_PI;
        if f < 0.0 {
            f += 360.0;
        }
        indr = ceil(MAXCATEGORY as Real * (r - minr) / (maxr - minr)) as c_int;
        indf = ceil(MAXCATEGORY as Real * (f - minf) / (maxf - minf)) as c_int;
        if dim == 3 {
            indz = ceil(MAXCATEGORY as Real * (z - minz) / (maxz - minz)) as c_int;
        }
        indr = c_min_int(c_max_int(indr, 1), MAXCATEGORY);
        indf = c_min_int(c_max_int(indf, 1), MAXCATEGORY);
        if dim == 3 {
            indz = c_min_int(c_max_int(indz, 1), MAXCATEGORY);
        }
        indr = *cumr.wrapping_offset(indr as isize);
        indf = *cumf.wrapping_offset(indf as isize);
        if dim == 3 {
            indz = *cumz.wrapping_offset(indz as isize);
        } else {
            indz = 1;
        }
        k = (indz - 1) * partitions1 * partitions2 + (indf - 1) * partitions1 + indr;
        if connect != 0 {
            if *elemconnect.wrapping_offset(j as isize) < 0 {
                k = indz;
            } else {
                k += 1;
            }
        }
        *inpart.wrapping_offset(j as isize) = k;
        *nopart.wrapping_offset(k as isize) += 1;
        j += 1;
    }
    maxpart = 0;
    minpart = noelements;
    j = 0;
    hit = FALSE;
    i = 1;
    while i <= partitions {
        if *nopart.wrapping_offset(i as isize) == 0 {
            hit = TRUE;
        }
        i += 1;
    }
    if hit != 0 {
        if info != 0 {
            printf(c"Renumbering the partitions (showing only 64):\n".as_ptr() as *const c_char);
        }
        i = 1;
        while i <= partitions {
            k = *nopart.wrapping_offset(i as isize);
            if k != 0 {
                maxpart = c_max_int(k, maxpart);
                minpart = c_min_int(k, minpart);
                j += 1;
                *nopart.wrapping_offset(i as isize) = j;
                if info != 0 && i != j && j <= 64 {
                    printf(c"%d -> %d %d\n".as_ptr() as *const c_char, i, j, k);
                }
            }
            i += 1;
        }
        if info != 0 {
            printf(
                c"There are %d active partitions out of %d possible ones\n".as_ptr()
                    as *const c_char,
                j,
                partitions,
            );
        }
        (*data).nopartitions = j;
        i = 1;
        while i <= noelements {
            j = *inpart.wrapping_offset(i as isize);
            *inpart.wrapping_offset(i as isize) = *nopart.wrapping_offset(j as isize);
            i += 1;
        }
    } else {
        if info != 0 {
            printf(
                c"All possible %d out of %d partitions are active\n".as_ptr() as *const c_char,
                j,
                partitions,
            );
        }
        (*data).nopartitions = partitions;
    }
    free_ivector(nopart, 1, partitions);
    free_ivector(cumr, 0, MAXCATEGORY);
    free_ivector(cumf, 0, MAXCATEGORY);
    if dim == 3 {
        free_ivector(cumz, 0, MAXCATEGORY);
    }
    PartitionNodesByElements(data, info);
    if info != 0 {
        printf(
            c"Successfully made a partitioning with %d to %d elements.\n".as_ptr()
                as *const c_char,
            minpart,
            maxpart,
        );
    }
    0
}
/// # Safety
///
/// direct BC-to-partition linking over valid meshes (C contract).
/// NOTE: C-dead locals (dim, allocated, partz, IndZ, noconnect, maxbcelem,
/// cumz, partelems, nodeconnect) are dropped — declared/assigned, never read
/// in the C source either.
#[no_mangle]
pub unsafe extern "C" fn PartitionConnectedElementsStraight(
    data: *mut FemType,
    bound: *mut BoundaryType,
    eg: *mut ElmergridType,
    info: c_int,
) -> c_int {
    // NOTE: C-dead locals j/debug/noknots/sidenodes dropped (assigned,
    // never read in the C source either).
    let (mut i, mut k, mut l): (c_int, c_int, c_int);
    let (mut hit, mut bctype): (c_int, c_int);
    let (noelements, mut bcelem, mut bc): (c_int, c_int, c_int);
    let mut totpartelems: c_int;
    let mut sideelemtype: c_int = 0;
    
    let mut sideind: [c_int; MAXNODESD2] = [0; MAXNODESD2];
    // NOTE: C stores debug = FALSE here but never reads it (dropped).
    if info != 0 {
        printf(c"Link the connected set directly to the partition\n".as_ptr() as *const c_char);
    }
    if (*data).nodeconnectexist == 0 {
        printf(c"There are no connected boundary nodes?\n".as_ptr() as *const c_char);
        return 1;
    }
    // NOTE: C stores noknots here but never reads it (dropped).
    noelements = (*data).noelements;
    totpartelems = 0;
    bcelem = 0;
    (*data).elemconnect = ivector(1, noelements);
    let elemconnect: *mut c_int = (*data).elemconnect;
    i = 1;
    while i <= noelements {
        *elemconnect.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    bc = 0;
    while bc < MAXBOUNDARIES as c_int {
        if (*bound.wrapping_offset(bc as isize)).created == FALSE {
            bc += 1;
            continue;
        }
        if (*bound.wrapping_offset(bc as isize)).nosides == 0 {
            bc += 1;
            continue;
        }
        i = 1;
        while i <= (*bound.wrapping_offset(bc as isize)).nosides {
            GetElementSide(
                *(*bound.wrapping_offset(bc as isize))
                    .parent
                    .wrapping_offset(i as isize),
                *(*bound.wrapping_offset(bc as isize))
                    .side
                    .wrapping_offset(i as isize),
                *(*bound.wrapping_offset(bc as isize))
                    .normal
                    .wrapping_offset(i as isize),
                data,
                sideind.as_mut_ptr(),
                &mut sideelemtype,
            );
            // NOTE: C stores sidenodes here but never reads it (dropped).
            hit = FALSE;
            k = 1;
            while k <= (*eg).connect {
                bctype = *(*eg)
                    .connectbounds
                    .as_mut_ptr()
                    .wrapping_offset((k - 1) as isize);
                hit = (*(*bound.wrapping_offset(bc as isize))
                    .types
                    .wrapping_offset(i as isize)
                    == bctype) as c_int;
                if hit != 0 {
                    break;
                }
                k += 1;
            }
            if hit == 0 {
                i += 1;
                continue;
            }
            bcelem += 1;
            k = *(*bound.wrapping_offset(bc as isize))
                .parent
                .wrapping_offset(i as isize);
            l = *(*bound.wrapping_offset(bc as isize))
                .parent2
                .wrapping_offset(i as isize);
            if k != 0
                && *elemconnect.wrapping_offset(k as isize) == 0 {
                    *elemconnect.wrapping_offset(k as isize) = 1;
                    totpartelems += 1;
                }
            if l != 0
                && *elemconnect.wrapping_offset(l as isize) == 0 {
                    *elemconnect.wrapping_offset(l as isize) = 1;
                    totpartelems += 1;
                }
            i += 1;
        }
        bc += 1;
    }
    (*data).elemconnectexist = totpartelems;
    if info != 0 {
        printf(
            c"Number of constrained boundary elements = %d\n".as_ptr() as *const c_char,
            bcelem,
        );
    }
    if info != 0 {
        printf(
            c"Number of constrained bulk elements = %d\n".as_ptr() as *const c_char,
            totpartelems,
        );
    }
    if info != 0 {
        printf(
            c"Successfully made connected elements to a partition\n".as_ptr() as *const c_char,
        );
    }
    0
}

/// # Safety
///
/// 1D BC partitioning with two-pass (omstart) category mapping.
/// NOTE: C-dead locals (partelems, nodeconnect, maxbcelem store) dropped as above.
#[no_mangle]
pub unsafe extern "C" fn PartitionConnectedElements1D(
    data: *mut FemType,
    bound: *mut BoundaryType,
    eg: *mut ElmergridType,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k, mut l): (c_int, c_int, c_int, c_int);
    let (dim, mut allocated, debug): (c_int, c_int, c_int);
    let (partz, partr, parts, mut hit, mut bctype): (c_int, c_int, c_int, c_int, c_int);
    let (noknots, noelements, mut bcelem, mut bc): (c_int, c_int, c_int, c_int);
    let (mut indz, mut noconnect, mut totpartelems): (c_int, c_int, c_int);
    let mut sideelemtype: c_int = 0;
    // NOTE: C declares sidehits/nohits but never uses them.
    let mut sidenodes: c_int;
    // NOTE: elemconnect is loop-carried (assigned end of pass 1, read from
    // pass 2 under the allocated flag); null init is unobservable.
    
    let mut elemconnect: *mut c_int = core::ptr::null_mut();
    let mut sideind: [c_int; MAXNODESD2] = [0; MAXNODESD2];
    let (mut val, mut z, mut maxz, mut minz): (Real, Real, Real, Real);
    debug = FALSE;
    partz = (*eg).partbcz;
    partr = (*eg).partbcr;
    if partz == 0 && partr == 0 {
        return 0;
    }
    parts = c_max_int(partz, partr);
    if info != 0 {
        if partz != 0 {
            printf(
                c"Making a simple 1D partitioning in z for the connected elements only\n".as_ptr()
                    as *const c_char,
            );
        } else {
            printf(
                c"Making a simple 1D partitioning in r for the connected elements only\n".as_ptr()
                    as *const c_char,
            );
        }
    }
    if (*data).nodeconnectexist == 0 {
        printf(c"There are no connected boundary nodes?\n".as_ptr() as *const c_char);
        return 1;
    }
    dim = (*data).dim;
    if dim < 3 {
        printf(
            c"PartitionConnectedElements1D is only applicable in 3D\n".as_ptr() as *const c_char,
        );
        return 2;
    }
    noknots = (*data).noknots;
    noelements = (*data).noelements;
    totpartelems = 0;
    if partr != 0 {
        z = sqrt(
            *(*data).x.wrapping_offset(1) * *(*data).z.wrapping_offset(1)
                + *(*data).y.wrapping_offset(1) * *(*data).y.wrapping_offset(1),
        );
    } else {
        z = *(*data).z.wrapping_offset(1);
    }
    maxz = z;
    minz = z;
    i = 1;
    while i <= noknots {
        if partr != 0 {
            z = sqrt(
                *(*data).x.wrapping_offset(i as isize) * *(*data).x.wrapping_offset(i as isize)
                    + *(*data).y.wrapping_offset(i as isize)
                        * *(*data).y.wrapping_offset(i as isize),
            );
        } else {
            z = *(*data).z.wrapping_offset(i as isize);
        }
        maxz = c_max_real(maxz, z);
        minz = c_min_real(minz, z);
        i += 1;
    }
    if info != 0 {
        printf(
            c"Range in coordinate extent: %12.5le %12.5le\n".as_ptr() as *const c_char,
            minz,
            maxz,
        );
    }
    let cumz: *mut c_int = ivector(0, MAXCATEGORY);
    i = 0;
    while i <= MAXCATEGORY {
        *cumz.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    allocated = FALSE;
    bcelem = 0;
    // NOTE: `omstart:` backward goto becomes a two-pass loop (first pass
    // counts into cumz, second assigns elemconnect).
    loop {
        bc = 0;
        while bc < MAXBOUNDARIES as c_int {
            if (*bound.wrapping_offset(bc as isize)).created == FALSE {
                bc += 1;
                continue;
            }
            if (*bound.wrapping_offset(bc as isize)).nosides == 0 {
                bc += 1;
                continue;
            }
            i = 1;
            while i <= (*bound.wrapping_offset(bc as isize)).nosides {
                GetElementSide(
                    *(*bound.wrapping_offset(bc as isize))
                        .parent
                        .wrapping_offset(i as isize),
                    *(*bound.wrapping_offset(bc as isize))
                        .side
                        .wrapping_offset(i as isize),
                    *(*bound.wrapping_offset(bc as isize))
                        .normal
                        .wrapping_offset(i as isize),
                    data,
                    sideind.as_mut_ptr(),
                    &mut sideelemtype,
                );
                sidenodes = sideelemtype % 100;
                hit = FALSE;
                k = 1;
                while k <= (*eg).connect {
                    bctype = *(*eg)
                        .connectbounds
                        .as_mut_ptr()
                        .wrapping_offset((k - 1) as isize);
                    if bctype > 0 {
                        if *(*bound.wrapping_offset(bc as isize))
                            .types
                            .wrapping_offset(i as isize)
                            == bctype
                        {
                            hit = TRUE;
                        }
                    } else if bctype == -1 {
                        if *(*bound.wrapping_offset(bc as isize))
                            .parent
                            .wrapping_offset(i as isize)
                            != 0
                        {
                            hit = TRUE;
                        }
                    } else if bctype == -2 {
                        if *(*bound.wrapping_offset(bc as isize))
                            .parent
                            .wrapping_offset(i as isize)
                            != 0
                            && *(*bound.wrapping_offset(bc as isize))
                                .parent2
                                .wrapping_offset(i as isize)
                                != 0
                        {
                            hit = TRUE;
                        }
                    } else if bctype == -3
                        && *(*bound.wrapping_offset(bc as isize))
                            .parent
                            .wrapping_offset(i as isize)
                            != 0
                            && *(*bound.wrapping_offset(bc as isize))
                                .parent2
                                .wrapping_offset(i as isize)
                                == 0
                        {
                            hit = TRUE;
                        }
                    if hit != 0 {
                        break;
                    }
                    k += 1;
                }
                if hit == 0 {
                    i += 1;
                    continue;
                }
                z = 0.0;
                j = 0;
                while j < sidenodes {
                    k = sideind[j as usize];
                    if partr != 0 {
                        val = sqrt(
                            *(*data).x.wrapping_offset(k as isize)
                                * *(*data).x.wrapping_offset(k as isize)
                                + *(*data).y.wrapping_offset(k as isize)
                                    * *(*data).y.wrapping_offset(k as isize),
                        );
                    } else {
                        val = *(*data).z.wrapping_offset(k as isize);
                    }
                    z += val;
                    j += 1;
                }
                z /= sidenodes as Real;
                indz = ceil(MAXCATEGORY as Real * (z - minz) / (maxz - minz)) as c_int;
                indz = c_min_int(c_max_int(indz, 1), MAXCATEGORY);
                if allocated != 0 {
                    indz = *cumz.wrapping_offset(indz as isize);
                    // NOTE: `if(0)` counter omitted (runtime-dead).
                    k = *(*bound.wrapping_offset(bc as isize))
                        .parent
                        .wrapping_offset(i as isize);
                    l = *(*bound.wrapping_offset(bc as isize))
                        .parent2
                        .wrapping_offset(i as isize);
                    if k != 0 {
                        *elemconnect.wrapping_offset(k as isize) = indz;
                        totpartelems += 1;
                    }
                    if l != 0 {
                        *elemconnect.wrapping_offset(l as isize) = indz;
                        totpartelems += 1;
                    }
                } else {
                    bcelem += 1;
                    *cumz.wrapping_offset(indz as isize) += 1;
                }
                i += 1;
            }
            bc += 1;
        }
        if allocated == 0 {
            // NOTE: C stores maxbcelem = bcelem here but never reads it.
            if debug != 0 {
                printf(
                    c"Number of constrained boundary elements = %d\n".as_ptr() as *const c_char,
                    bcelem,
                );
                printf(
                    c"Differential categories (showing only 20 active ones from %d)\n".as_ptr()
                        as *const c_char,
                    MAXCATEGORY,
                );
                k = 0;
                j = 0;
                while j <= MAXCATEGORY {
                    if *cumz.wrapping_offset(j as isize) > 0 {
                        k += 1;
                        printf(
                            c"%d : %d\n".as_ptr() as *const c_char,
                            j,
                            *cumz.wrapping_offset(j as isize),
                        );
                        if k == 20 {
                            break;
                        }
                    }
                    j += 1;
                }
            }
            i = 1;
            while i <= MAXCATEGORY {
                *cumz.wrapping_offset(i as isize) =
                    *cumz.wrapping_offset(i as isize) + *cumz.wrapping_offset((i - 1) as isize);
                i += 1;
            }
            if debug != 0 {
                printf(c"Cumulative categories\n".as_ptr() as *const c_char);
                j = 0;
                while j <= MAXCATEGORY {
                    printf(
                        c"%d : %d\n".as_ptr() as *const c_char,
                        j,
                        *cumz.wrapping_offset(j as isize),
                    );
                    j += 1;
                }
            }
            noconnect = bcelem;
            i = 1;
            while i <= MAXCATEGORY {
                *cumz.wrapping_offset(i as isize) = ceil(
                    1.0 * parts as Real * *cumz.wrapping_offset(i as isize) as Real
                        / noconnect as Real,
                ) as c_int;
                i += 1;
            }
            if debug != 0 {
                printf(c"Partition categories\n".as_ptr() as *const c_char);
                j = 0;
                while j <= MAXCATEGORY {
                    printf(
                        c"%d : %d\n".as_ptr() as *const c_char,
                        j,
                        *cumz.wrapping_offset(j as isize),
                    );
                    j += 1;
                }
            }
            (*data).elemconnect = ivector(1, noelements);
            elemconnect = (*data).elemconnect;
            i = 1;
            while i <= noelements {
                *elemconnect.wrapping_offset(i as isize) = 0;
                i += 1;
            }
            allocated = TRUE;
            if info != 0 {
                printf(c"Allocated and now setting the entries\n".as_ptr() as *const c_char);
            }
            // NOTE: `goto omstart` via loop-continue.
            continue;
        }
        break;
    }
    (*data).elemconnectexist = totpartelems;
    free_ivector(cumz, 0, MAXCATEGORY);
    if info != 0 {
        printf(
            c"Successfully made a partitioning 1D with %d elements\n".as_ptr() as *const c_char,
            totpartelems,
        );
    }
    0
}
// SAFETY: METIS option setup; indices/values pinned to bundled metis-5.1.0.
unsafe fn SetMetisOptions(options: *mut idx_t, eg: *mut ElmergridType, info: c_int) -> c_int {
    if info != 0 {
        printf(c"Setting default Metis options!\n".as_ptr() as *const c_char);
    }
    METIS_SetDefaultOptions(options);
    if info != 0 {
        printf(c"Setting user defined Metis options!\n".as_ptr() as *const c_char);
    }
    *options.wrapping_add(METIS_OPTION_NUMBERING) = 0;
    *options.wrapping_add(METIS_OPTION_DBGLVL) = 3;
    *options.wrapping_add(METIS_OPTION_CTYPE) = METIS_CTYPE_SHEM;
    *options.wrapping_add(METIS_OPTION_IPTYPE) = METIS_IPTYPE_GROW;
    *options.wrapping_add(METIS_OPTION_RTYPE) = METIS_RTYPE_GREEDY;
    if (*eg).metis_contig != 0 {
        *options.wrapping_add(METIS_OPTION_CONTIG) = 1;
    }
    if (*eg).metis_seed != 0 {
        *options.wrapping_add(METIS_OPTION_SEED) = (*eg).metis_seed;
    }
    if (*eg).metis_ncuts > 1 {
        *options.wrapping_add(METIS_OPTION_NCUTS) = (*eg).metis_ncuts;
    }
    if (*eg).metis_volcut != 0 {
        *options.wrapping_add(METIS_OPTION_OBJTYPE) = METIS_OBJTYPE_VOL;
    }
    if (*eg).metis_minconn != 0 {
        *options.wrapping_add(METIS_OPTION_MINCONN) = 1;
    }
    0
}

/// # Safety
///
/// dual-graph build + METIS call over valid connected BCs.
/// NOTE: C reads maxbcelem uninitialized on the first step-2 growth (stack
/// garbage sizing the Ivector); transcribed as a deterministic 0-init
/// (fail-safe exit(1) via the ivector guard instead of silent UB).
#[no_mangle]
pub unsafe extern "C" fn PartitionConnectedElementsMetis(
    data: *mut FemType,
    bound: *mut BoundaryType,
    eg: *mut ElmergridType,
    mut nparts: c_int,
    metisopt: c_int,
    info: c_int,
) -> c_int {
    // NOTE: C-dead locals n/partelems/nodepart/dim/set dropped (never read).
    let (mut i, mut j, mut k, mut l, mut m): (c_int, c_int, c_int, c_int, c_int);
    let (noknots, noelements, mut sidenodes, mut ind, mut nohits, mut totpartelems): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (mut dualmaxcon, mut invmaxcon, mut totcon, mut step, mut bc, mut bcelem, mut bcelem2): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (mut hit, mut noconnect): (c_int, c_int);
    // NOTE: minpartelems/maxpartelems declared but never used in C.
    let mut sideelemtype: c_int = 0;
    let mut con: c_int;
    // NOTE: C reads maxbcelem uninitialized on the first step-2 growth;
    // deterministic 0-init (fail-safe ivector guard, see SAFETY above).
    let mut maxbcelem: c_int = 0;
    let (elemconnect, nodeperm, neededby): (*mut c_int, *mut c_int, *mut c_int);
    let mut bcdualgraph: [*mut c_int; MAXCONNECTIONS] = [core::ptr::null_mut(); MAXCONNECTIONS];
    let mut bcinvtopo: [*mut c_int; MAXCONNECTIONS] = [core::ptr::null_mut(); MAXCONNECTIONS];
    let mut sideind: [c_int; MAXNODESD1] = [0; MAXNODESD1];
    let mut nn: c_int;
    let mut options: [idx_t; METIS_NOPTIONS] = [0; METIS_NOPTIONS];
    let (xadj, adjncy, vwgt, adjwgt): (*mut idx_t, *mut idx_t, *mut idx_t, *mut idx_t);
    // NOTE: ncon is indeterminate in C (never assigned, read by METIS);
    // transcribed as 1 (the only sane constraint count, as in MetisGraph).
    // edgecut is METIS-written, never read by C — 0-init for the handoff.
    // NOTE: C assigns numflag = 0 but never reads it (the METIS calls take
    // &ncon); dropped.
    let (mut edgecut, mut ncon): (idx_t, idx_t) = (0, 1);
    let mut wgtflag: idx_t;
    
    if info != 0 {
        printf(
            c"Making a Metis partitioning for the connected BC elements only\n".as_ptr()
                as *const c_char,
        );
    }
    if (*data).nodeconnectexist == 0 {
        printf(c"There are no connected elements\n".as_ptr() as *const c_char);
        return 1;
    }
    SetMetisOptions(options.as_mut_ptr(), eg, info);
    // NOTE: C stores dim = data->dim here but never reads it.
    noknots = (*data).noknots;
    noelements = (*data).noelements;
    noconnect = 0;
    i = 1;
    while i <= noknots {
        if *(*data).nodeconnect.wrapping_offset(i as isize) != 0 {
            noconnect += 1;
        }
        i += 1;
    }
    if noconnect == 0 {
        printf(c"There are really no connected nodes?\n".as_ptr() as *const c_char);
        return 2;
    }
    if info != 0 {
        printf(
            c"Number of connected nodes is %d\n".as_ptr() as *const c_char,
            noconnect,
        );
    }
    nodeperm = ivector(1, noknots);
    i = 1;
    while i <= noknots {
        *nodeperm.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    noconnect = 0;
    i = 1;
    while i <= noknots {
        if *(*data).nodeconnect.wrapping_offset(i as isize) != 0 {
            noconnect += 1;
            *nodeperm.wrapping_offset(i as isize) = noconnect;
        }
        i += 1;
    }
    if info != 0 {
        printf(c"Permuted the connected nodes\n".as_ptr() as *const c_char);
    }
    invmaxcon = 0;
    dualmaxcon = 0;
    totcon = 0;
    neededby = ivector(1, noconnect);
    i = 1;
    while i <= noconnect {
        *neededby.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    step = 1;
    while step <= 2 {
        bcelem = 0;
        bc = 0;
        while bc < MAXBOUNDARIES as c_int {
            if (*bound.wrapping_offset(bc as isize)).created == FALSE {
                bc += 1;
                continue;
            }
            if (*bound.wrapping_offset(bc as isize)).nosides == 0 {
                bc += 1;
                continue;
            }
            i = 1;
            while i <= (*bound.wrapping_offset(bc as isize)).nosides {
                GetBoundaryElement(
                    i,
                    bound.wrapping_offset(bc as isize),
                    data,
                    sideind.as_mut_ptr(),
                    &mut sideelemtype,
                );
                // NOTE: commented-out GetElementSide alternative omitted.
                sidenodes = sideelemtype % 100;
                nohits = 0;
                j = 0;
                while j < sidenodes {
                    if *nodeperm.wrapping_offset(sideind[j as usize] as isize) != 0 {
                        nohits += 1;
                    }
                    j += 1;
                }
                if nohits < sidenodes {
                    i += 1;
                    continue;
                }
                bcelem += 1;
                j = 0;
                while j < sidenodes {
                    ind = *nodeperm.wrapping_offset(sideind[j as usize] as isize);
                    if step == 1 {
                        *neededby.wrapping_offset(ind as isize) += 1;
                        k = *neededby.wrapping_offset(ind as isize);
                        if k > invmaxcon {
                            invmaxcon += 1;
                            bcinvtopo[invmaxcon as usize] = ivector(1, noconnect);
                            // NOTE: `if(0)` diagnostic omitted.
                            m = 1;
                            while m <= noconnect {
                                *bcinvtopo[invmaxcon as usize].wrapping_offset(m as isize) = 0;
                                m += 1;
                            }
                        }
                        *bcinvtopo[k as usize].wrapping_offset(ind as isize) = bcelem;
                    } else if step == 2 {
                        k = 1;
                        while k <= invmaxcon {
                            bcelem2 = *bcinvtopo[k as usize].wrapping_offset(ind as isize);
                            if bcelem2 == 0 {
                                break;
                            }
                            if bcelem2 == bcelem {
                                k += 1;
                                continue;
                            }
                            hit = FALSE;
                            l = 1;
                            while l <= dualmaxcon {
                                if *bcdualgraph[l as usize].wrapping_offset(bcelem as isize)
                                    == bcelem2
                                {
                                    hit = TRUE;
                                }
                                if *bcdualgraph[l as usize].wrapping_offset(bcelem as isize) == 0 {
                                    break;
                                }
                                l += 1;
                            }
                            if hit == 0 {
                                if l > dualmaxcon {
                                    dualmaxcon += 1;
                                    if l >= MAXCONNECTIONS as c_int {
                                        printf(
                                            c"Number of connections %d vs. static limit %d\n"
                                                .as_ptr()
                                                as *const c_char,
                                            l,
                                            MAXCONNECTIONS as c_int - 1,
                                        );
                                        bigerror(
                                            c"Maximum of connections in dual graph larger than the static limit!"
                                                .as_ptr() as *const c_char,
                                        );
                                    }
                                    // NOTE: `if(0)` diagnostic omitted.
                                    bcdualgraph[dualmaxcon as usize] = ivector(1, maxbcelem);
                                    m = 1;
                                    while m <= maxbcelem {
                                        *bcdualgraph[dualmaxcon as usize]
                                            .wrapping_offset(m as isize) = 0;
                                        m += 1;
                                    }
                                }
                                *bcdualgraph[l as usize].wrapping_offset(bcelem as isize) = bcelem2;
                                totcon += 1;
                            }
                            k += 1;
                        }
                    }
                    j += 1;
                }
                i += 1;
            }
            bc += 1;
        }
        maxbcelem = bcelem;
        step += 1;
    }
    if info != 0 {
        printf(
            c"There are %d connected boundary elements.\n".as_ptr() as *const c_char,
            maxbcelem,
        );
        printf(
            c"There are %d connections altogether in the dual graph.\n".as_ptr() as *const c_char,
            totcon,
        );
    }
    xadj = ivector(0, maxbcelem);
    adjncy = ivector(0, totcon - 1);
    i = 0;
    while i < totcon {
        *adjncy.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    totcon = 0;
    i = 1;
    while i <= maxbcelem {
        *xadj.wrapping_offset((i - 1) as isize) = totcon;
        j = 1;
        while j <= dualmaxcon {
            con = *bcdualgraph[j as usize].wrapping_offset(i as isize);
            if con != 0 {
                *adjncy.wrapping_offset(totcon as isize) = con - 1;
                totcon += 1;
                // NOTE: `if(0)` diagnostic omitted.
            }
            j += 1;
        }
        i += 1;
    }
    *xadj.wrapping_offset(maxbcelem as isize) = totcon;
    free_ivector(neededby, 1, noconnect);
    i = 1;
    while i <= invmaxcon {
        free_ivector(bcinvtopo[i as usize], 1, noconnect);
        i += 1;
    }
    i = 1;
    while i <= dualmaxcon {
        free_ivector(bcdualgraph[i as usize], 1, maxbcelem);
        i += 1;
    }
    // NOTE: C stores numflag = 0 here but never reads it (dropped).
    nn = maxbcelem;
    let npart: *mut idx_t = ivector(0, nn - 1);
    wgtflag = 0;
    vwgt = core::ptr::null_mut();
    adjwgt = core::ptr::null_mut();
    if metisopt == 2 {
        if info != 0 {
            printf(
                c"Starting graph partitioning METIS_PartGraphRecursive.\n".as_ptr()
                    as *const c_char,
            );
        }
        // NOTE: C passes &wgtflag as the vsize slot (quirk preserved).
        METIS_PartGraphRecursive(
            &mut nn,
            &mut ncon,
            xadj,
            adjncy,
            vwgt,
            &mut wgtflag,
            adjwgt,
            &mut nparts,
            core::ptr::null_mut(),
            core::ptr::null_mut(),
            options.as_mut_ptr(),
            &mut edgecut,
            npart,
        );
    } else if metisopt == 3 || metisopt == 4 {
        if info != 0 {
            printf(
                c"Starting graph partitioning METIS_PartGraphKway.\n".as_ptr() as *const c_char,
            );
        }
        ncon = 0;
        wgtflag = 0;
        // NOTE: C passes &wgtflag as the vsize slot (quirk preserved).
        METIS_PartGraphKway(
            &mut nn,
            &mut ncon,
            xadj,
            adjncy,
            vwgt,
            &mut wgtflag,
            adjwgt,
            &mut nparts,
            core::ptr::null_mut(),
            core::ptr::null_mut(),
            options.as_mut_ptr(),
            &mut edgecut,
            npart,
        );
    } else {
        printf(
            c"Unknown Metis option %d\n".as_ptr() as *const c_char,
            metisopt,
        );
    }
    // NOTE: `if(0)` diagnostic omitted.
    free_ivector(xadj, 0, maxbcelem);
    free_ivector(adjncy, 0, totcon - 1);
    (*data).elemconnect = ivector(1, noelements);
    elemconnect = (*data).elemconnect;
    i = 1;
    while i <= noelements {
        *elemconnect.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    bcelem = 0;
    totpartelems = 0;
    bc = 0;
    while bc < MAXBOUNDARIES as c_int {
        if (*bound.wrapping_offset(bc as isize)).created == FALSE {
            bc += 1;
            continue;
        }
        if (*bound.wrapping_offset(bc as isize)).nosides == 0 {
            bc += 1;
            continue;
        }
        i = 1;
        while i <= (*bound.wrapping_offset(bc as isize)).nosides {
            GetBoundaryElement(
                i,
                bound.wrapping_offset(bc as isize),
                data,
                sideind.as_mut_ptr(),
                &mut sideelemtype,
            );
            sidenodes = sideelemtype % 100;
            nohits = 0;
            j = 0;
            while j < sidenodes {
                if *nodeperm.wrapping_offset(sideind[j as usize] as isize) != 0 {
                    nohits += 1;
                }
                j += 1;
            }
            if nohits < sidenodes {
                i += 1;
                continue;
            }
            k = *(*bound.wrapping_offset(bc as isize))
                .parent
                .wrapping_offset(i as isize);
            l = *(*bound.wrapping_offset(bc as isize))
                .parent2
                .wrapping_offset(i as isize);
            // NOTE: `if(0)` diagnostic omitted.
            if k != 0 {
                *elemconnect.wrapping_offset(k as isize) =
                    *npart.wrapping_offset(bcelem as isize) + 1;
                totpartelems += 1;
            }
            if l != 0 {
                *elemconnect.wrapping_offset(l as isize) =
                    *npart.wrapping_offset(bcelem as isize) + 1;
                totpartelems += 1;
            }
            bcelem += 1;
            i += 1;
        }
        bc += 1;
    }
    free_ivector(nodeperm, 1, noknots);
    free_ivector(npart, 0, nn - 1);
    (*data).elemconnectexist = totpartelems;
    if info != 0 {
        printf(
            c"Successfully created boundary partitioning for %d bulk elements\n".as_ptr()
                as *const c_char,
            totpartelems,
        );
    }
    0
}
/// # Safety
///
/// majority-rule extension over valid inverse topology (C contract).
#[no_mangle]
// NOTE: C never uses the bound parameter either (signature only).
pub unsafe extern "C" fn ExtendBoundaryPartitioning(
    data: *mut FemType,
    _bound: *mut BoundaryType,
    elemlayers: c_int,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k, mut l, mut m, mut n): (c_int, c_int, c_int, c_int, c_int, c_int);
    let (mut nonodes, noelements, mut totpartelems, mut nparts, mut set): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    // NOTE: refpart is conditionally assigned but unconditionally read in C
    // (garbage on degenerate 0-node elements); deterministic 0-init here.
    let (mut minpartelems, mut maxpartelems, mut refcount, mut part): (c_int, c_int, c_int, c_int);
    let mut refpart: c_int = 0;
    let (partelems, elemconnect): (*mut c_int, *mut c_int);
    let (invrow, invcol): (*mut c_int, *mut c_int);
    printf(c"Extending boundary partitioning by majority rule\n".as_ptr() as *const c_char);
    CreateInverseTopology(data, info);
    invrow = (*data).invtopo.rows;
    invcol = (*data).invtopo.cols;
    noelements = (*data).noelements;
    // NOTE: C stores noknots here but never reads it (dropped).
    if (*data).elemconnectexist == 0 {
        printf(c"Initial partitioning does not exist!\n".as_ptr() as *const c_char);
        return 1;
    }
    elemconnect = (*data).elemconnect;
    nparts = 0;
    i = 1;
    while i <= (*data).noelements {
        nparts = c_max_int(nparts, *elemconnect.wrapping_offset(i as isize));
        i += 1;
    }
    partelems = ivector(1, nparts);
    i = 1;
    while i <= nparts {
        *partelems.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    totpartelems = 0;
    i = 1;
    while i <= noelements {
        j = *elemconnect.wrapping_offset(i as isize);
        if j != 0 {
            *partelems.wrapping_offset(j as isize) += 1;
            totpartelems += 1;
        }
        i += 1;
    }
    if info != 0 {
        printf(
            c"Initial partitioning with %d elements:\n".as_ptr() as *const c_char,
            totpartelems,
        );
    }
    minpartelems = *partelems.wrapping_offset(1);
    maxpartelems = *partelems.wrapping_offset(1);
    j = 1;
    while j <= nparts {
        minpartelems = c_min_int(minpartelems, *partelems.wrapping_offset(j as isize));
        maxpartelems = c_max_int(maxpartelems, *partelems.wrapping_offset(j as isize));
        j += 1;
    }
    printf(
        c"Initial partitioning has %d to %d elements in partition\n".as_ptr() as *const c_char,
        minpartelems,
        maxpartelems,
    );
    i = 1;
    while i <= nparts {
        *partelems.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    m = 1;
    while m <= elemlayers {
        n = 0;
        j = 1;
        while j <= noelements {
            if *elemconnect.wrapping_offset(j as isize) > 0 {
                j += 1;
                continue;
            }
            set = FALSE;
            nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
            i = 0;
            while i < nonodes {
                k = *(*(*data).topology.wrapping_offset(j as isize)).wrapping_offset(i as isize);
                l = *invrow.wrapping_offset((k - 1) as isize);
                while l < *invrow.wrapping_offset(k as isize) {
                    part = *elemconnect
                        .wrapping_offset((*invcol.wrapping_offset(l as isize) + 1) as isize);
                    if part > 0 {
                        set = TRUE;
                        *partelems.wrapping_offset(part as isize) += 1;
                    }
                    l += 1;
                }
                i += 1;
            }
            if set == 0 {
                j += 1;
                continue;
            }
            refcount = 0;
            i = 0;
            while i < nonodes {
                k = *(*(*data).topology.wrapping_offset(j as isize)).wrapping_offset(i as isize);
                l = *invrow.wrapping_offset((k - 1) as isize);
                while l < *invrow.wrapping_offset(k as isize) {
                    part = *elemconnect
                        .wrapping_offset((*invcol.wrapping_offset(l as isize) + 1) as isize);
                    if part > 0 {
                        if *partelems.wrapping_offset(part as isize) > refcount {
                            refcount = *partelems.wrapping_offset(part as isize);
                            refpart = part;
                        }
                        *partelems.wrapping_offset(part as isize) = 0;
                    }
                    l += 1;
                }
                i += 1;
            }
            // NOTE: `if(0)` diagnostic omitted (runtime-dead).
            *elemconnect.wrapping_offset(j as isize) = -refpart;
            n += 1;
            j += 1;
        }
        j = 1;
        while j <= noelements {
            *elemconnect.wrapping_offset(j as isize) =
                abs(*elemconnect.wrapping_offset(j as isize));
            j += 1;
        }
        if info != 0 {
            printf(
                c"Setting partition for %d elements in layer %d\n".as_ptr() as *const c_char,
                n,
                m,
            );
        }
        if n == 0 {
            printf(
                c"No elements set in the layer, breaking loop %d!\n".as_ptr() as *const c_char,
                m,
            );
            break;
        }
        m += 1;
    }
    i = 1;
    while i <= nparts {
        *partelems.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    totpartelems = 0;
    i = 1;
    while i <= noelements {
        j = *elemconnect.wrapping_offset(i as isize);
        if j != 0 {
            *partelems.wrapping_offset(j as isize) += 1;
            totpartelems += 1;
        }
        i += 1;
    }
    if info != 0 {
        printf(
            c"Extended partitioning of %d elements:\n".as_ptr() as *const c_char,
            totpartelems,
        );
    }
    minpartelems = *partelems.wrapping_offset(1);
    maxpartelems = *partelems.wrapping_offset(1);
    j = 1;
    while j <= nparts {
        minpartelems = c_min_int(minpartelems, *partelems.wrapping_offset(j as isize));
        maxpartelems = c_max_int(maxpartelems, *partelems.wrapping_offset(j as isize));
        // NOTE: `if(0)` diagnostic omitted (runtime-dead).
        j += 1;
    }
    if info != 0 {
        printf(
            c"Extended partitioning has %d to %d elements in partition\n".as_ptr()
                as *const c_char,
            minpartelems,
            maxpartelems,
        );
    }
    (*data).elemconnectexist = totpartelems;
    free_ivector(partelems, 1, nparts);
    if info != 0 {
        printf(
            c"Successfully extended the partitioning by %d layers\n".as_ptr() as *const c_char,
            elemlayers,
        );
    }
    0
}
/// # Safety
///
/// nodal geometric partitioning over valid meshes (C contract).
#[no_mangle]
pub unsafe extern "C" fn PartitionSimpleNodes(
    data: *mut FemType,
    dimpart: *mut c_int,
    dimper: *mut c_int,
    partorder: c_int,
    corder: *mut Real,
    parttol: Real,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k, mut k0, mut kprev, mut l, mut ind): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (mut minpart, mut maxpart): (c_int, c_int);
    let (noknots, noelements, mut elemsinpart, periodic): (c_int, c_int, c_int, c_int);
    let (partitions1, partitions2, mut partitions3, partitions): (c_int, c_int, c_int, c_int);
    let (mut vpartitions1, mut vpartitions2, mut vpartitions3, vpartitions, mut hit): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (indx, nopart, nodepart): (*mut c_int, *mut c_int, *mut c_int);
    
    let (mut x, mut y, mut z, cx, cy, cz): (Real, Real, Real, Real, Real, Real);
    partitions1 = *dimpart.wrapping_offset(0);
    partitions2 = *dimpart.wrapping_offset(1);
    partitions3 = *dimpart.wrapping_offset(2);
    if (*data).dim < 3 {
        partitions3 = 1;
    }
    partitions = partitions1 * partitions2 * partitions3;
    if partitions1 < 2 && partitions2 < 2 && partitions3 < 2 {
        printf(
            c"Some of the divisions must be larger than one: %d %d %d\n".as_ptr()
                as *const c_char,
            partitions1,
            partitions2,
            partitions3,
        );
        bigerror(c"Partitioning not performed".as_ptr() as *const c_char);
    }
    if partitions >= (*data).noelements {
        printf(
            c"There must be fever partitions than elements (%d vs %d)!\n".as_ptr()
                as *const c_char,
            partitions,
            (*data).noelements,
        );
        bigerror(c"Partitioning not performed".as_ptr() as *const c_char);
    }
    if (*data).partitionexist == 0 {
        (*data).partitionexist = TRUE;
        (*data).elempart = ivector(1, (*data).noelements);
        (*data).nodepart = ivector(1, (*data).noknots);
        (*data).nopartitions = partitions;
    }
    nodepart = (*data).nodepart;
    vpartitions1 = partitions1;
    vpartitions2 = partitions2;
    vpartitions3 = partitions3;
    periodic = *dimper.wrapping_offset(0) | *dimper.wrapping_offset(1) | *dimper.wrapping_offset(2);
    if periodic != 0 {
        if *dimper.wrapping_offset(0) != 0 && partitions1 > 1 {
            vpartitions1 *= 2;
        }
        if *dimper.wrapping_offset(1) != 0 && partitions2 > 1 {
            vpartitions2 *= 2;
        }
        if *dimper.wrapping_offset(2) != 0 && partitions3 > 1 {
            vpartitions3 *= 2;
        }
    }
    vpartitions = vpartitions1 * vpartitions2 * vpartitions3;
    nopart = ivector(1, vpartitions);
    noelements = (*data).noelements;
    noknots = (*data).noknots;
    if info != 0 {
        printf(
            c"Making a simple partitioning for %d nodes in %d-dimensions.\n".as_ptr()
                as *const c_char,
            noknots,
            (*data).dim,
        );
    }
    let arrange: *mut Real = dvector(1, noknots);
    indx = ivector(1, noknots);
    if partorder != 0 {
        cx = *corder.wrapping_offset(0);
        cy = *corder.wrapping_offset(1);
        cz = *corder.wrapping_offset(2);
    } else {
        cx = 1.0;
        cy = 0.0001;
        cz = cy * cy;
    }
    // NOTE: C stores z = 0.0 here; dropped (see SimpleElements note).
    i = 1;
    while i <= noknots {
        *nodepart.wrapping_offset(i as isize) = 1;
        i += 1;
    }
    if vpartitions1 > 1 {
        if info != 0 {
            printf(
                c"Ordering 1st direction with (%.3lg*x + %.3lg*y + %.3lg*z)\n".as_ptr()
                    as *const c_char,
                cx,
                cy,
                cz,
            );
        }
        j = 1;
        while j <= noknots {
            x = *(*data).x.wrapping_offset(j as isize);
            y = *(*data).y.wrapping_offset(j as isize);
            z = *(*data).z.wrapping_offset(j as isize);
            *arrange.wrapping_offset(j as isize) = cx * x + cy * y + cz * z;
            j += 1;
        }
        SortIndex(noknots, arrange, indx);
        i = 1;
        while i <= noknots {
            ind = *indx.wrapping_offset(i as isize);
            k = (i * vpartitions1 - 1) / noknots + 1;
            *nodepart.wrapping_offset(ind as isize) = k;
            i += 1;
        }
    }
    if vpartitions2 > 1 {
        if info != 0 {
            printf(c"Ordering in the 2nd direction.\n".as_ptr() as *const c_char);
        }
        j = 1;
        while j <= noknots {
            x = *(*data).x.wrapping_offset(j as isize);
            y = *(*data).y.wrapping_offset(j as isize);
            z = *(*data).z.wrapping_offset(j as isize);
            *arrange.wrapping_offset(j as isize) = -cy * x + cx * y + cz * z;
            j += 1;
        }
        SortIndex(noknots, arrange, indx);
        i = 1;
        while i <= vpartitions {
            *nopart.wrapping_offset(i as isize) = 0;
            i += 1;
        }
        elemsinpart = noknots / (vpartitions1 * vpartitions2);
        j = 1;
        i = 1;
        while i <= noknots {
            ind = *indx.wrapping_offset(i as isize);
            k0 = (*nodepart.wrapping_offset(ind as isize) - 1) * vpartitions2;
            l = 1;
            // NOTE: exact C `for(l=1;l<=vpartitions;l++)` shape (0-trip when
            // vpartitions < 1); the hit-break exits early as in C.
            while l <= vpartitions {
                hit = FALSE;
                if j < vpartitions
                    && *nopart.wrapping_offset((k0 + j) as isize) >= elemsinpart {
                        j += 1;
                        hit = TRUE;
                    }
                if j > 1
                    && *nopart.wrapping_offset((k0 + j - 1) as isize) < elemsinpart {
                        j -= 1;
                        hit = TRUE;
                    }
                if hit == 0 {
                    break;
                }
                l += 1;
            }
            k = k0 + j;
            *nopart.wrapping_offset(k as isize) += 1;
            *nodepart.wrapping_offset(ind as isize) = k;
            i += 1;
        }
    }
    if vpartitions3 > 1 {
        if info != 0 {
            printf(c"Ordering in the 3rd direction.\n".as_ptr() as *const c_char);
        }
        j = 1;
        while j <= noknots {
            x = *(*data).x.wrapping_offset(j as isize);
            y = *(*data).y.wrapping_offset(j as isize);
            z = *(*data).z.wrapping_offset(j as isize);
            *arrange.wrapping_offset(j as isize) = -cz * x - cy * y + cx * z;
            j += 1;
        }
        SortIndex(noknots, arrange, indx);
        i = 1;
        while i <= vpartitions {
            *nopart.wrapping_offset(i as isize) = 0;
            i += 1;
        }
        elemsinpart = noknots / (vpartitions1 * vpartitions2 * vpartitions3);
        j = 1;
        i = 1;
        while i <= noknots {
            ind = *indx.wrapping_offset(i as isize);
            k0 = (*nodepart.wrapping_offset(ind as isize) - 1) * vpartitions3;
            l = 1;
            // NOTE: exact C `for(l=1;l<=vpartitions;l++)` shape (see 2nd dir).
            while l <= vpartitions {
                hit = FALSE;
                if j < vpartitions3
                    && *nopart.wrapping_offset((k0 + j) as isize) >= elemsinpart {
                        j += 1;
                        hit = TRUE;
                    }
                if j > 1
                    && *nopart.wrapping_offset((k0 + j - 1) as isize) < elemsinpart {
                        j -= 1;
                        hit = TRUE;
                    }
                if hit == 0 {
                    break;
                }
                l += 1;
            }
            k = k0 + j;
            *nopart.wrapping_offset(k as isize) += 1;
            *nodepart.wrapping_offset(ind as isize) = k;
            i += 1;
        }
    }
    if parttol > 1.0e-20 {
        // NOTE: `if(0)` dump omitted (runtime-dead).
        j = 1;
        while j <= noknots {
            x = *(*data).x.wrapping_offset(j as isize);
            y = *(*data).y.wrapping_offset(j as isize);
            z = *(*data).z.wrapping_offset(j as isize);
            *arrange.wrapping_offset(j as isize) = cx * x + cy * y + cz * z;
            j += 1;
        }
        SortIndex(noknots, arrange, indx);
        kprev = *nodepart.wrapping_offset(*indx.wrapping_offset(1) as isize);
        l = 0;
        i = 2;
        while i <= noknots {
            ind = *indx.wrapping_offset(i as isize);
            k = *nodepart.wrapping_offset(ind as isize);
            if k != kprev
                && fabs(
                    *arrange.wrapping_offset(i as isize)
                        - *arrange.wrapping_offset((i - 1) as isize),
                ) < parttol
                {
                    *nopart.wrapping_offset(k as isize) -= 1;
                    *nopart.wrapping_offset(kprev as isize) += 1;
                    *nodepart.wrapping_offset(ind as isize) = kprev;
                    k = kprev;
                    // NOTE: `if(0)` diagnostic omitted (runtime-dead).
                    l += 1;
                }
            kprev = k;
            i += 1;
        }
        // NOTE: `if(0)` dump omitted (runtime-dead).
        if info != 0 {
            printf(
                c"Number of partition changes due to tolerance: %d\n".as_ptr() as *const c_char,
                l,
            );
        }
    }
    if periodic != 0 {
        {
            
            let (mut p1, mut p2, mut p3, mut q1, mut q2, mut q3): (
                c_int,
                c_int,
                c_int,
                c_int,
                c_int,
                c_int,
            );
            // NOTE: P/Q match the C locals.
            let (mut P, mut Q): (c_int, c_int);
            let partmap: *mut c_int = ivector(1, vpartitions);
            i = 1;
            while i <= vpartitions {
                *partmap.wrapping_offset(i as isize) = 0;
                i += 1;
            }
            p1 = 1;
            while p1 <= vpartitions1 {
                q1 = p1;
                if *dimper.wrapping_offset(0) != 0 && vpartitions1 > 1 {
                    if q1 == vpartitions1 {
                        q1 = 0;
                    }
                    q1 = q1 / 2 + 1;
                }
                p2 = 1;
                while p2 <= vpartitions2 {
                    q2 = p2;
                    if *dimper.wrapping_offset(1) != 0 && vpartitions2 > 1 {
                        if q2 == vpartitions2 {
                            q2 = 0;
                        }
                        q2 = q2 / 2 + 1;
                    }
                    p3 = 1;
                    while p3 <= vpartitions3 {
                        q3 = p3;
                        if *dimper.wrapping_offset(2) != 0 && vpartitions3 > 1 {
                            if q3 == vpartitions3 {
                                q3 = 0;
                            }
                            q3 = q3 / 2 + 1;
                        }
                        P = vpartitions3 * vpartitions2 * (p1 - 1) + vpartitions3 * (p2 - 1) + p3;
                        Q = partitions3 * partitions2 * (q1 - 1) + partitions3 * (q2 - 1) + q3;
                        *partmap.wrapping_offset(P as isize) = Q;
                        p3 += 1;
                    }
                    p2 += 1;
                }
                p1 += 1;
            }
            i = 1;
            while i <= noknots {
                *nodepart.wrapping_offset(i as isize) =
                    *partmap.wrapping_offset(*nodepart.wrapping_offset(i as isize) as isize);
                i += 1;
            }
            free_ivector(partmap, 1, vpartitions);
        }
    }
    i = 1;
    while i <= partitions {
        *nopart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    i = 1;
    while i <= noknots {
        *nopart.wrapping_offset(*nodepart.wrapping_offset(i as isize) as isize) += 1;
        i += 1;
    }
    minpart = *nopart.wrapping_offset(1);
    maxpart = *nopart.wrapping_offset(1);
    i = 1;
    while i <= partitions {
        minpart = c_min_int(minpart, *nopart.wrapping_offset(i as isize));
        maxpart = c_max_int(maxpart, *nopart.wrapping_offset(i as isize));
        i += 1;
    }
    free_dvector(arrange, 1, noelements);
    free_ivector(nopart, 1, partitions);
    free_ivector(indx, 1, noelements);
    PartitionElementsByNodes(data, info);
    if info != 0 {
        printf(
            c"Successfully made a partitioning with %d to %d nodes.\n".as_ptr() as *const c_char,
            minpart,
            maxpart,
        );
    }
    0
}
// SAFETY: read-only partitioning audit over valid meshes (C contract).
unsafe fn CheckPartitioning(data: *mut FemType, info: c_int) {
    let (mut i, mut j, mut k): (c_int, c_int, c_int);
    let (partitions, mut part, mut part2): (c_int, c_int, c_int);
    let (noknots, noelements, mut mini, mut maxi, mut sumi): (c_int, c_int, c_int, c_int, c_int);
    let (mut hit, mut ind, mut nodesd2, mut elemtype): (c_int, c_int, c_int, c_int);
    let (elempart, nodepart): (*mut c_int, *mut c_int);
    let (elemsinpart, nodesinpart, sharedinpart): (*mut c_int, *mut c_int, *mut c_int);
    noknots = (*data).noknots;
    noelements = (*data).noelements;
    partitions = (*data).nopartitions;
    elemsinpart = ivector(1, partitions);
    nodesinpart = ivector(1, partitions);
    sharedinpart = ivector(1, partitions);
    i = 1;
    while i <= partitions {
        *elemsinpart.wrapping_offset(i as isize) = 0;
        *nodesinpart.wrapping_offset(i as isize) = 0;
        *sharedinpart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    if info != 0 {
        printf(c"Checking for partitioning\n".as_ptr() as *const c_char);
    }
    elempart = (*data).elempart;
    j = 0;
    i = 1;
    while i <= (*data).noelements {
        part = *elempart.wrapping_offset(i as isize);
        if part < 1 || part > partitions {
            j += 1;
        } else {
            *elemsinpart.wrapping_offset(part as isize) += 1;
        }
        i += 1;
    }
    if j != 0 {
        printf(
            c"Bad initial partitioning: %d elements do not belong anywhere!\n".as_ptr()
                as *const c_char,
            j,
        );
        bigerror(c"Can't continue with broken partitioning".as_ptr() as *const c_char);
    }
    nodepart = (*data).nodepart;
    j = 0;
    i = 1;
    while i <= (*data).noknots {
        part = *nodepart.wrapping_offset(i as isize);
        if part < 1 || part > partitions {
            j += 1;
        } else {
            *nodesinpart.wrapping_offset(part as isize) += 1;
        }
        i += 1;
    }
    if j != 0 {
        printf(
            c"Bad initial partitioning: %d nodes do not belong anywhere!\n".as_ptr()
                as *const c_char,
            j,
        );
        bigerror(c"Can't continue with broken partitioning".as_ptr() as *const c_char);
    }
    if (*data).partitiontableexists != 0 {
        i = 1;
        while i <= noknots {
            part = *nodepart.wrapping_offset(i as isize);
            j = 1;
            while j <= (*data).maxpartitiontable {
                part2 = *(*data).partitiontable[j as usize].wrapping_offset(i as isize);
                if part2 == 0 {
                    break;
                }
                if part != part2 {
                    *sharedinpart.wrapping_offset(part2 as isize) += 1;
                }
                j += 1;
            }
            i += 1;
        }
    }
    if info != 0 {
        printf(c"Information on partition bandwidth\n".as_ptr() as *const c_char);
        if partitions <= 10 {
            printf(
                c"Distribution of elements, nodes and shared nodes\n".as_ptr() as *const c_char,
            );
            printf(
                c"     %-10s %-10s %-10s %-10s\n".as_ptr() as *const c_char,
                c"partition".as_ptr() as *const c_char,
                c"elements".as_ptr() as *const c_char,
                c"nodes".as_ptr() as *const c_char,
                c"shared".as_ptr() as *const c_char,
            );
            i = 1;
            while i <= partitions {
                printf(
                    c"     %-10d %-10d %-10d %-10d\n".as_ptr() as *const c_char,
                    i,
                    *elemsinpart.wrapping_offset(i as isize),
                    *nodesinpart.wrapping_offset(i as isize),
                    *sharedinpart.wrapping_offset(i as isize),
                );
                i += 1;
            }
        } else {
            mini = *elemsinpart.wrapping_offset(1);
            maxi = *elemsinpart.wrapping_offset(1);
            i = 1;
            while i <= partitions {
                mini = c_min_int(mini, *elemsinpart.wrapping_offset(i as isize));
                maxi = c_max_int(maxi, *elemsinpart.wrapping_offset(i as isize));
                i += 1;
            }
            printf(
                c"Average %d elements with range %d in partition\n".as_ptr() as *const c_char,
                noelements / partitions,
                maxi - mini,
            );
            mini = *nodesinpart.wrapping_offset(1);
            maxi = *nodesinpart.wrapping_offset(1);
            i = 1;
            while i <= partitions {
                mini = c_min_int(mini, *nodesinpart.wrapping_offset(i as isize));
                maxi = c_max_int(maxi, *nodesinpart.wrapping_offset(i as isize));
                i += 1;
            }
            printf(
                c"Average %d nodes with range %d in partition\n".as_ptr() as *const c_char,
                noknots / partitions,
                maxi - mini,
            );
            sumi = 0;
            mini = *sharedinpart.wrapping_offset(1);
            maxi = *sharedinpart.wrapping_offset(1);
            i = 1;
            while i <= partitions {
                mini = c_min_int(mini, *sharedinpart.wrapping_offset(i as isize));
                maxi = c_max_int(maxi, *sharedinpart.wrapping_offset(i as isize));
                sumi += *sharedinpart.wrapping_offset(i as isize);
                i += 1;
            }
            printf(
                c"Average %d shared nodes with range %d in partition\n".as_ptr() as *const c_char,
                sumi / partitions,
                maxi - mini,
            );
        }
    }
    if (*data).maxpartitiontable == 0 {
        return;
    }
    // NOTE: `if(0)` diagnostics omitted (runtime-dead).
    i = 1;
    while i <= (*data).noelements {
        part = *elempart.wrapping_offset(i as isize);
        elemtype = *(*data).elementtypes.wrapping_offset(i as isize);
        nodesd2 = elemtype % 100;
        j = 0;
        while j < nodesd2 {
            ind = *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize);
            hit = FALSE;
            part2 = 0;
            k = 1;
            while k <= (*data).maxpartitiontable {
                part2 = *(*data).partitiontable[k as usize].wrapping_offset(ind as isize);
                if part == part2 {
                    hit = TRUE;
                }
                if hit != 0 && part == 0 {
                    break;
                }
                k += 1;
            }
            if hit == 0 {
                printf(c"******** Warning *******\n".as_ptr() as *const c_char);
                printf(
                    c"Node %d in element %d does not belong to partition %d (%d)\n".as_ptr()
                        as *const c_char,
                    ind,
                    i,
                    part,
                    part2,
                );
                printf(
                    c"elemtype = %d nodesd2 = %d\n".as_ptr() as *const c_char,
                    elemtype,
                    nodesd2,
                );
                k = 0;
                while k < nodesd2 {
                    printf(
                        c"ind[%d] = %d\n".as_ptr() as *const c_char,
                        k,
                        *(*(*data).topology.wrapping_offset(i as isize))
                            .wrapping_offset(k as isize),
                    );
                    k += 1;
                }
            }
            j += 1;
        }
        i += 1;
    }
    // NOTE: `if(0)` diagnostic omitted (runtime-dead).
    i = 1;
    while i <= (*data).noknots {
        part = *nodepart.wrapping_offset(i as isize);
        hit = FALSE;
        j = 1;
        while j <= (*data).maxpartitiontable {
            part2 = *(*data).partitiontable[j as usize].wrapping_offset(i as isize);
            if part == part2 {
                hit = TRUE;
            }
            if hit != 0 && part == 0 {
                break;
            }
            j += 1;
        }
        if hit == 0 {
            printf(
                c"***** Node %d in partition %d is not in partition list\n".as_ptr()
                    as *const c_char,
                i,
                part,
            );
        }
        i += 1;
    }
}

/// # Safety
///
/// BC-interface ownership repair over valid partitions (C contract).
/// `do {} while` becomes loop+break; `if(0)` blocks omitted.
// Dead `fix1/fix2 = 0` stores mirror C exactly (empty `if(fix1 || fix2)` kept
// verbatim for oscillation-damping intent); lane precedent: setdiv allows.
#[allow(unused_assignments)]
#[no_mangle]
pub unsafe extern "C" fn OptimizePartitioningAtBoundary(
    data: *mut FemType,
    bound: *mut BoundaryType,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k, mut l): (c_int, c_int, c_int, c_int);
    let (mut boundaryelems, mut ind): (c_int, c_int);
    let (mut hit1, mut hit2, mut fix1, mut fix2): (c_int, c_int, c_int, c_int);
    let (mut dompart, mut part1, mut part2, mut newmam, mut mam1, mut mam2, mut nodesd2): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    
    if (*data).partitionexist == 0 {
        printf(
            c"OptimizePartitioningAtBoundary: this should be called only after partitioning\n"
                .as_ptr() as *const c_char,
        );
        bigerror(c"Optimization not performed!".as_ptr() as *const c_char);
    }
    if info != 0 {
        printf(c"Optimizing the partitioning at boundaries.\n".as_ptr() as *const c_char);
    }
    let alteredparent: *mut c_int = ivector(1, (*data).noelements);
    i = 1;
    while i <= (*data).noelements {
        *alteredparent.wrapping_offset(i as isize) = *(*data).elempart.wrapping_offset(i as isize);
        i += 1;
    }
    k = 0;
    loop {
        k += 1;
        boundaryelems = 0;
        j = 0;
        while j < MAXBOUNDARIES as c_int {
            if (*bound.wrapping_offset(j as isize)).created == 0 {
                j += 1;
                continue;
            }
            i = 1;
            while i <= (*bound.wrapping_offset(j as isize)).nosides {
                if (*bound.wrapping_offset(j as isize)).ediscont != 0
                    && *(*bound.wrapping_offset(j as isize))
                        .discont
                        .wrapping_offset(i as isize)
                        != 0
                    {
                        i += 1;
                        continue;
                    }
                mam1 = *(*bound.wrapping_offset(j as isize))
                    .parent
                    .wrapping_offset(i as isize);
                mam2 = abs(*(*bound.wrapping_offset(j as isize))
                    .parent2
                    .wrapping_offset(i as isize));
                if mam1 == 0 || mam2 == 0 {
                    i += 1;
                    continue;
                }
                part1 = *(*data).elempart.wrapping_offset(mam1 as isize);
                part2 = *(*data).elempart.wrapping_offset(mam2 as isize);
                if part1 == part2 {
                    i += 1;
                    continue;
                }
                if (*data).elemconnectexist != 0 {
                    fix1 = (*(*data).elemconnect.wrapping_offset(mam1 as isize) < 0) as c_int;
                    fix2 = (*(*data).elemconnect.wrapping_offset(mam2 as isize) < 0) as c_int;
                    if fix1 != 0 && fix2 != 0 {
                        i += 1;
                        continue;
                    }
                } else {
                    fix1 = 0;
                    fix2 = 0;
                }
                // NOTE: empty `if(fix1 || fix2)` block kept verbatim (C has
                // it empty too; it documents the oscillation-damping intent).
                
                if k < 5 {
                    hit1 = 0;
                    hit2 = 0;
                    nodesd2 = *(*data).elementtypes.wrapping_offset(mam1 as isize) % 100;
                    l = 0;
                    while l < nodesd2 {
                        ind = *(*(*data).topology.wrapping_offset(mam1 as isize))
                            .wrapping_offset(l as isize);
                        if *(*data).nodepart.wrapping_offset(ind as isize) == part1 {
                            hit1 += 1;
                        }
                        if *(*data).nodepart.wrapping_offset(ind as isize) == part2 {
                            hit2 += 1;
                        }
                        l += 1;
                    }
                    nodesd2 = *(*data).elementtypes.wrapping_offset(mam2 as isize) % 100;
                    l = 0;
                    while l < nodesd2 {
                        ind = *(*(*data).topology.wrapping_offset(mam2 as isize))
                            .wrapping_offset(l as isize);
                        if *(*data).nodepart.wrapping_offset(ind as isize) == part1 {
                            hit1 += 1;
                        }
                        if *(*data).nodepart.wrapping_offset(ind as isize) == part2 {
                            hit2 += 1;
                        }
                        l += 1;
                    }
                    fix1 = (hit1 >= hit2) as c_int;
                    // NOTE: C stores fix2 here (and FALSE below), but fix2
                    // is never read after this point (only the earlier empty
                    // if reads it); both stores dropped for warning-freedom.
                } else {
                    fix1 = TRUE;
                }
                if fix1 != 0 {
                    dompart = part1;
                    newmam = mam2;
                } else {
                    dompart = part2;
                    newmam = mam1;
                }
                *(*data).elempart.wrapping_offset(newmam as isize) = dompart;
                boundaryelems += 1;
                // NOTE: `if(0)`/`else if(0)` node-moving blocks omitted.
                i += 1;
            }
            j += 1;
        }
        if info != 0 && boundaryelems != 0 {
            printf(
                c"Round %d: %d bulk elements with BCs removed from interface.\n".as_ptr()
                    as *const c_char,
                k,
                boundaryelems,
            );
        }
        if !(boundaryelems != 0 && k < 10) {
            break;
        }
    }
    j = 0;
    i = 1;
    while i <= (*data).noelements {
        if *alteredparent.wrapping_offset(i as isize)
            == *(*data).elempart.wrapping_offset(i as isize)
        {
            *alteredparent.wrapping_offset(i as isize) = 0;
        } else {
            j += 1;
        }
        i += 1;
    }
    free_ivector(alteredparent, 1, (*data).noelements);
    if info != 0 {
        printf(
            c"Ownership of %d parents was changed at BCs\n".as_ptr() as *const c_char,
            j,
        );
    }
    j = 0;
    while j < MAXBOUNDARIES as c_int {
        if (*bound.wrapping_offset(j as isize)).created == 0 {
            j += 1;
            continue;
        }
        i = 1;
        while i <= (*bound.wrapping_offset(j as isize)).nosides {
            *(*bound.wrapping_offset(j as isize))
                .parent2
                .wrapping_offset(i as isize) = c_max_int(
                0,
                *(*bound.wrapping_offset(j as isize))
                    .parent2
                    .wrapping_offset(i as isize),
            );
            i += 1;
        }
        j += 1;
    }
    // NOTE: `if(0)` diagnostic omitted (runtime-dead).
    0
}
// SAFETY: recursive level assignment over valid CSR graphs (C contract).
unsafe fn Levelize(
    n: c_int,
    level: c_int,
    maxlevel: *mut c_int,
    levels: *mut c_int,
    rows: *mut c_int,
    cols: *mut c_int,
    done: *mut c_int,
) {
    let (mut j, mut k): (c_int, c_int);
    *levels.wrapping_offset(n as isize) = level;
    *done.wrapping_offset(n as isize) = TRUE;
    *maxlevel = c_max_int(*maxlevel, level);
    j = *rows.wrapping_offset(n as isize);
    while j < *rows.wrapping_offset((n + 1) as isize) {
        k = *cols.wrapping_offset(j as isize);
        if *done.wrapping_offset(k as isize) == 0 {
            Levelize(k, level + 1, maxlevel, levels, rows, cols, done);
        }
        j += 1;
    }
}

// SAFETY: Cuthill-McKee renumbering over valid CSR graphs (C contract).
unsafe fn RenumberCuthillMckee(
    nrows: c_int,
    rows: *mut c_int,
    cols: *mut c_int,
    iperm: *mut c_int,
) -> c_int {
    let (mut i, mut j, mut k): (c_int, c_int, c_int);
    let (mut n, mut startn, mut mindegree, mut maxlevel, mut newroot): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (mut bw_bef, mut bw_aft): (c_int, c_int);
    let (level, degree, done): (*mut c_int, *mut c_int, *mut c_int);
    done = ivector(0, nrows - 1);
    level = ivector(0, nrows - 1);
    degree = ivector(0, nrows - 1);
    bw_bef = 0;
    i = 0;
    while i < nrows {
        j = *rows.wrapping_offset(i as isize);
        while j < *rows.wrapping_offset((i + 1) as isize) {
            bw_bef = c_max_int(bw_bef, c_abs_int(*cols.wrapping_offset(j as isize) - i) + 1);
            j += 1;
        }
        *degree.wrapping_offset(i as isize) =
            *rows.wrapping_offset((i + 1) as isize) - *rows.wrapping_offset(i as isize);
        i += 1;
    }
    printf(
        c"RenumberCuthillMckee: Bandwidth before: %d\n".as_ptr() as *const c_char,
        bw_bef,
    );
    startn = 0;
    mindegree = *degree.wrapping_offset(startn as isize);
    i = 0;
    while i < nrows {
        if *degree.wrapping_offset(i as isize) < mindegree {
            startn = i;
            mindegree = *degree.wrapping_offset(i as isize);
        }
        *level.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    maxlevel = 0;
    i = 0;
    while i < nrows {
        *done.wrapping_offset(i as isize) = FALSE;
        i += 1;
    }
    Levelize(startn, 0, &mut maxlevel, level, rows, cols, done);
    newroot = TRUE;
    while newroot != 0 {
        newroot = FALSE;
        mindegree = *degree.wrapping_offset(startn as isize);
        k = startn;
        i = 0;
        while i < nrows {
            if *level.wrapping_offset(i as isize) == maxlevel
                && *degree.wrapping_offset(i as isize) < mindegree {
                    k = i;
                    mindegree = *degree.wrapping_offset(i as isize);
                }
            i += 1;
        }
        if k != startn {
            j = maxlevel;
            maxlevel = 0;
            i = 0;
            while i < nrows {
                *done.wrapping_offset(i as isize) = FALSE;
                i += 1;
            }
            Levelize(k, 0, &mut maxlevel, level, rows, cols, done);
            if j > maxlevel {
                newroot = TRUE;
                startn = j;
            }
        }
    }
    i = 0;
    while i < nrows {
        *done.wrapping_offset(i as isize) = -1;
        *iperm.wrapping_offset(i as isize) = -1;
        i += 1;
    }
    *done.wrapping_offset(0) = startn;
    *iperm.wrapping_offset(startn as isize) = 0;
    i = 1;
    j = 0;
    while j < nrows {
        if *done.wrapping_offset(j as isize) < 0 {
            k = 0;
            while k < nrows {
                if *iperm.wrapping_offset(k as isize) < 0 {
                    *done.wrapping_offset(i as isize) = k;
                    *iperm.wrapping_offset(k as isize) = i;
                    i += 1;
                    break;
                }
                k += 1;
            }
        }
        k = *rows.wrapping_offset(*done.wrapping_offset(j as isize) as isize);
        while k < *rows.wrapping_offset(*done.wrapping_offset(j as isize) as isize + 1) {
            n = *cols.wrapping_offset(k as isize);
            if *iperm.wrapping_offset(n as isize) < 0 {
                *done.wrapping_offset(i as isize) = n;
                *iperm.wrapping_offset(n as isize) = i;
                i += 1;
            }
            k += 1;
        }
        j += 1;
    }
    i = 0;
    while i < nrows {
        *iperm.wrapping_offset(*done.wrapping_offset(i as isize) as isize) = nrows - 1 - i;
        i += 1;
    }
    bw_aft = 0;
    i = 0;
    while i < nrows {
        j = *rows.wrapping_offset(i as isize);
        while j < *rows.wrapping_offset((i + 1) as isize) {
            bw_aft = c_max_int(
                bw_aft,
                c_abs_int(
                    *iperm.wrapping_offset(*cols.wrapping_offset(j as isize) as isize)
                        - *iperm.wrapping_offset(i as isize),
                ) + 1,
            );
            j += 1;
        }
        i += 1;
    }
    printf(
        c"RenumberCuthillMckee: Bandwidth after: %d\n".as_ptr() as *const c_char,
        bw_aft,
    );
    free_ivector(level, 0, nrows - 1);
    free_ivector(done, 0, nrows - 1);
    free_ivector(degree, 0, nrows - 1);
    (bw_aft < bw_bef) as c_int
}

// SAFETY: partition renumbering over valid partition tables (C contract).
unsafe fn RenumberPartitions(data: *mut FemType, info: c_int) {
    let (mut i, mut j, mut k, mut l, mut m): (c_int, c_int, c_int, c_int, c_int);
    let (mut hit, mut con, mut totcon, noelements, noknots, partitions): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (maxneededtimes, mut totneededtimes): (c_int, c_int);
    let (mut part, mut part1, mut part2, bw_reduced): (c_int, c_int, c_int, c_int);
    let (nodepart, elempart): (*mut c_int, *mut c_int);
    
    let (xadj, adjncy): (*mut c_int, *mut c_int);
    let mut partparttable: [*mut c_int; MAXCONNECTIONS] = [core::ptr::null_mut(); MAXCONNECTIONS];
    if info != 0 {
        printf(c"Renumbering partitions to minimize bandwidth.\n".as_ptr() as *const c_char);
    }
    partitions = (*data).nopartitions;
    noelements = (*data).noelements;
    noknots = (*data).noknots;
    maxneededtimes = (*data).maxpartitiontable;
    totneededtimes = 0;
    totcon = 0;
    i = 1;
    while i <= noknots {
        if *(*data).partitiontable[2].wrapping_offset(i as isize) == 0 {
            i += 1;
            continue;
        }
        j = 1;
        while j <= maxneededtimes {
            part1 = *(*data).partitiontable[j as usize].wrapping_offset(i as isize);
            if part1 == 0 {
                break;
            }
            k = 1;
            while k <= maxneededtimes {
                if k == j {
                    k += 1;
                    continue;
                }
                part2 = *(*data).partitiontable[k as usize].wrapping_offset(i as isize);
                if part2 == 0 {
                    break;
                }
                hit = 0;
                l = 1;
                while l <= totneededtimes {
                    if *partparttable[l as usize].wrapping_offset(part1 as isize) == part2 {
                        hit = -1;
                        break;
                    } else if *partparttable[l as usize].wrapping_offset(part1 as isize) == 0 {
                        totcon += 1;
                        *partparttable[l as usize].wrapping_offset(part1 as isize) = part2;
                        hit = 1;
                        break;
                    }
                    l += 1;
                }
                if hit == 0 {
                    totneededtimes += 1;
                    partparttable[totneededtimes as usize] = ivector(1, partitions);
                    m = 1;
                    while m <= partitions {
                        *partparttable[totneededtimes as usize].wrapping_offset(m as isize) = 0;
                        m += 1;
                    }
                    *partparttable[totneededtimes as usize].wrapping_offset(part1 as isize) = part2;
                    totcon += 1;
                }
                k += 1;
            }
            j += 1;
        }
        i += 1;
    }
    if info != 0 {
        printf(
            c"There are %d connections altogether\n".as_ptr() as *const c_char,
            totcon,
        );
        printf(
            c"There are %.3f connections between partitions in average\n".as_ptr()
                as *const c_char,
            1.0 * totcon as Real / partitions as Real,
        );
    }
    xadj = ivector(0, partitions);
    adjncy = ivector(0, totcon - 1);
    i = 0;
    while i < totcon {
        *adjncy.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    totcon = 0;
    i = 1;
    while i <= partitions {
        *xadj.wrapping_offset((i - 1) as isize) = totcon;
        j = 1;
        while j <= totneededtimes {
            con = *partparttable[j as usize].wrapping_offset(i as isize);
            if con == 0 {
                j += 1;
                continue;
            }
            *adjncy.wrapping_offset(totcon as isize) = con - 1;
            totcon += 1;
            j += 1;
        }
        i += 1;
    }
    *xadj.wrapping_offset(partitions as isize) = totcon;
    let perm: *mut c_int = ivector(0, partitions - 1);
    bw_reduced = RenumberCuthillMckee(partitions, xadj, adjncy, perm);
    // NOTE: `if(0 && info)` diagnostic omitted (runtime-dead).
    if bw_reduced != 0 {
        if info != 0 {
            printf(
                c"Successful ordering: moving partitions to new positions\n".as_ptr()
                    as *const c_char,
            );
        }
        nodepart = (*data).nodepart;
        elempart = (*data).elempart;
        i = 1;
        while i <= noelements {
            *elempart.wrapping_offset(i as isize) =
                *perm.wrapping_offset(*elempart.wrapping_offset(i as isize) as isize - 1) + 1;
            i += 1;
        }
        i = 1;
        while i <= noknots {
            *nodepart.wrapping_offset(i as isize) =
                *perm.wrapping_offset(*nodepart.wrapping_offset(i as isize) as isize - 1) + 1;
            i += 1;
        }
        i = 1;
        while i <= noknots {
            j = 1;
            while j <= maxneededtimes {
                part = *(*data).partitiontable[j as usize].wrapping_offset(i as isize);
                if part == 0 {
                    break;
                }
                *(*data).partitiontable[j as usize].wrapping_offset(i as isize) =
                    *perm.wrapping_offset(part as isize - 1) + 1;
                j += 1;
            }
            i += 1;
        }
    }
    i = 1;
    while i <= totneededtimes {
        free_ivector(partparttable[i as usize], 1, partitions);
        i += 1;
    }
    free_ivector(xadj, 0, partitions);
    free_ivector(adjncy, 0, totcon - 1);
}

// SAFETY: pure statistics over the needed vector (C contract).
unsafe fn CheckSharedDeviation(neededvector: *mut c_int, partitions: c_int, info: c_int) -> c_int {
    let mut i: c_int;
    let (mut minshared, mut maxshared, dshared): (c_int, c_int, c_int);
    let (mut sumshared, mut sumshared2, varshared, mut needed, det, ave): (
        Real,
        Real,
        Real,
        Real,
        Real,
        Real,
    );
    sumshared = 0.0;
    sumshared2 = 0.0;
    minshared = *neededvector.wrapping_offset(1);
    maxshared = *neededvector.wrapping_offset(1);
    i = 1;
    while i <= partitions {
        needed = 1.0 * *neededvector.wrapping_offset(i as isize) as Real;
        sumshared += needed;
        sumshared2 += needed * needed;
        maxshared = c_max_int(maxshared, *neededvector.wrapping_offset(i as isize));
        minshared = c_min_int(minshared, *neededvector.wrapping_offset(i as isize));
        i += 1;
    }
    dshared = maxshared - minshared;
    if info != 0 {
        ave = sumshared / partitions as Real;
        det = 1.0 * sumshared2 / partitions as Real
            - 1.0 * (sumshared / partitions as Real) * (sumshared / partitions as Real);
        varshared = sqrt(det);
        printf(
            c"Average number of elements in partition %.3le\n".as_ptr() as *const c_char,
            ave,
        );
        printf(
            c"Maximum deviation in ownership %d\n".as_ptr() as *const c_char,
            dshared,
        );
        printf(
            c"Average deviation in ownership %.3le\n".as_ptr() as *const c_char,
            varshared,
        );
        printf(
            c"Average relative deviation %.2lf %%\n".as_ptr() as *const c_char,
            100.0 * varshared / ave,
        );
    }
    dshared
}
/// # Safety
///
/// iterative load-balance optimization over valid partitions.
/// `goto optimizeownership` becomes loop-continue; block decls become blocks.
#[no_mangle]
// NOTE: C never uses the bound parameter either (signature only).
pub unsafe extern "C" fn OptimizePartitioning(
    data: *mut FemType,
    _bound: *mut BoundaryType,
    noopt: c_int,
    partbw: c_int,
    info: c_int,
) -> c_int {
    // NOTE: n is read at the goto-check on the noopt path where C leaves it
    // indeterminate (short-circuit avoids the read there too); the 0-init is
    // unobservable either way.
    let (mut i, mut j, mut k, mut l, mut m): (c_int, c_int, c_int, c_int, c_int);
    let mut n: c_int = 0;
    let (noelements, partitions, mut ind, mut hit): (c_int, c_int, c_int, c_int);
    let (noknots, mut dshared, mut dshared0): (c_int, c_int, c_int);
    let (elempart, nodepart): (*mut c_int, *mut c_int);
    let (mut sharings, mut maxshared): (c_int, c_int);
    let (mut nodesd2, maxneededtimes): (c_int, c_int);
    // NOTE: NULL-inits in C; assigned on first use paths guarded by the
    // same flags that guard the reads — unobservable inits.
    let mut probnodes: *mut c_int = core::ptr::null_mut();
    let mut optimize: c_int;
    
    let mut somethingdone: c_int = 0;
    let mut elemparts: *mut c_int = core::ptr::null_mut();
    let mut invelemparts: *mut c_int = core::ptr::null_mut();
    let mut knows: *mut *mut c_int = core::ptr::null_mut();
    if (*data).partitionexist == 0 {
        printf(
            c"OptimizePartitioning: this should be called only after partitioning\n".as_ptr()
                as *const c_char,
        );
        bigerror(c"Optimization not performed!".as_ptr() as *const c_char);
    }
    noknots = (*data).noknots;
    noelements = (*data).noelements;
    partitions = (*data).nopartitions;
    if info != 0 {
        printf(
            c"Optimizing for %d partitions\n".as_ptr() as *const c_char,
            partitions,
        );
    }
    elempart = (*data).elempart;
    nodepart = (*data).nodepart;
    maxshared = 0;
    if partitions < 2 {
        printf(
            c"OptimizePartitioning: does not make sense for %d partitions\n".as_ptr()
                as *const c_char,
            partitions,
        );
        bigerror(c"Optimization not performed!".as_ptr() as *const c_char);
    }
    CreatePartitionTable(data, info);
    maxneededtimes = (*data).maxpartitiontable;
    if partbw != 0 {
        RenumberPartitions(data, info);
    }
    printf(c"Checking partitioning before optimization\n".as_ptr() as *const c_char);
    CheckPartitioning(data, info);
    let neededvector: *mut c_int = ivector(1, partitions);
    i = 1;
    while i <= partitions {
        *neededvector.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    i = 1;
    while i <= noknots {
        *neededvector.wrapping_offset(*nodepart.wrapping_offset(i as isize) as isize) += 1;
        i += 1;
    }
    if noopt == 0 {
        optimize = 1;
        probnodes = ivector(1, noknots);
        i = 1;
        while i <= noknots {
            *probnodes.wrapping_offset(i as isize) = 0;
            i += 1;
        }
        printf(
            c"Applying aggressive optimization for load balancing\n".as_ptr() as *const c_char,
        );
    } else {
        optimize = 1;
    }
    // NOTE: `optimizeownership:` backward goto becomes loop-continue.
    loop {
        dshared = CheckSharedDeviation(neededvector, partitions, info);
        if noopt == 0 {
            // NOTE: block-scoped C locals become a Rust block.
            {
                let (mut nochanges, maxrounds, dtarget): (c_int, c_int, c_int);
                maxrounds = 5;
                dtarget = 3;
                // NOTE: C stores nochanges = 0 here; the loop body
                // reassigned it first — dropped (n = 1 kept).
                n = 1;
                while n <= maxrounds {
                    nochanges = 0;
                    i = 1;
                    while i <= noknots {
                        ind = i;
                        k = *(*data).partitiontable[2].wrapping_offset(ind as isize);
                        if k == 0 {
                            i += 1;
                            continue;
                        }
                        if maxneededtimes > 2
                            && *(*data).partitiontable[3].wrapping_offset(ind as isize) != 0 {
                                i += 1;
                                continue;
                            }
                        if *probnodes.wrapping_offset(ind as isize) != 0 {
                            i += 1;
                            continue;
                        }
                        j = *(*data).partitiontable[1].wrapping_offset(ind as isize);
                        if abs(*neededvector.wrapping_offset(j as isize)
                            - *neededvector.wrapping_offset(k as isize))
                            < dtarget
                        {
                            i += 1;
                            continue;
                        }
                        if *neededvector.wrapping_offset(j as isize)
                            < *neededvector.wrapping_offset(k as isize)
                            && *nodepart.wrapping_offset(ind as isize) == k
                        {
                            nochanges += 1;
                            *neededvector.wrapping_offset(j as isize) += 1;
                            *neededvector.wrapping_offset(k as isize) -= 1;
                            *nodepart.wrapping_offset(ind as isize) = j;
                        } else if *neededvector.wrapping_offset(k as isize)
                            < *neededvector.wrapping_offset(j as isize)
                            && *nodepart.wrapping_offset(ind as isize) == j
                        {
                            nochanges += 1;
                            *neededvector.wrapping_offset(k as isize) += 1;
                            *neededvector.wrapping_offset(j as isize) -= 1;
                            *nodepart.wrapping_offset(ind as isize) = k;
                        }
                        i += 1;
                    }
                    if info != 0 && nochanges != 0 {
                        printf(
                            c"Changed the ownership of %d nodes\n".as_ptr() as *const c_char,
                            nochanges,
                        );
                    }
                    dshared0 = dshared;
                    dshared = CheckSharedDeviation(neededvector, partitions, info);
                    if dshared >= dshared0 {
                        break;
                    }
                    somethingdone += nochanges;
                    n += 1;
                }
            }
        }
        if info != 0 {
            printf(c"Checking for problematic sharings\n".as_ptr() as *const c_char);
        }
        m = 0;
        if partitions > 2 {
            if info != 0 {
                printf(
                    c"Optimizing sharing for %d partitions\n".as_ptr() as *const c_char,
                    partitions,
                );
            }
            loop {
                // NOTE: block-scoped C locals become a Rust block.
                {
                    let (mut i1, mut i2, mut e1, mut e2, mut owners): (
                        c_int,
                        c_int,
                        c_int,
                        c_int,
                        c_int,
                    );
                    m += 1;
                    sharings = 0;
                    i1 = 0;
                    i2 = 0;
                    // NOTE: C stores e1 = e2 = 0 here; reassigned before any
                    // read in the hit branch — dropped.
                    maxshared = 27;
                    if m == 1 && optimize == 1 {
                        elemparts = ivector(1, partitions);
                        invelemparts = ivector(1, maxshared);
                        knows = imatrix(1, maxshared, 1, maxshared);
                    }
                    j = 1;
                    while j <= partitions {
                        *elemparts.wrapping_offset(j as isize) = 0;
                        j += 1;
                    }
                    j = 1;
                    while j <= maxshared {
                        *invelemparts.wrapping_offset(j as isize) = 0;
                        j += 1;
                    }
                    j = 1;
                    while j <= maxshared {
                        k = 1;
                        while k <= maxshared {
                            *(*knows.wrapping_offset(j as isize)).wrapping_offset(k as isize) = 0;
                            k += 1;
                        }
                        j += 1;
                    }
                    i = 1;
                    while i <= noelements {
                        let mut ownpart: c_int;
                        owners = 0;
                        nodesd2 = *(*data).elementtypes.wrapping_offset(i as isize) % 100;
                        ownpart = FALSE;
                        j = 0;
                        while j < nodesd2 {
                            ind = *(*(*data).topology.wrapping_offset(i as isize))
                                .wrapping_offset(j as isize);
                            k = *nodepart.wrapping_offset(ind as isize);
                            if k == *elempart.wrapping_offset(i as isize) {
                                ownpart = TRUE;
                            }
                            if *elemparts.wrapping_offset(k as isize) == 0 {
                                owners += 1;
                                *elemparts.wrapping_offset(k as isize) = owners;
                                *invelemparts.wrapping_offset(owners as isize) = k;
                            }
                            j += 1;
                        }
                        if owners - ownpart <= 1 {
                            j = 1;
                            while j <= owners {
                                k = *invelemparts.wrapping_offset(j as isize);
                                *elemparts.wrapping_offset(k as isize) = 0;
                                j += 1;
                            }
                            i += 1;
                            continue;
                        }
                        j = 0;
                        while j < nodesd2 {
                            ind = *(*(*data).topology.wrapping_offset(i as isize))
                                .wrapping_offset(j as isize);
                            l = 1;
                            while l <= maxneededtimes {
                                e1 = *(*data).partitiontable[l as usize]
                                    .wrapping_offset(ind as isize);
                                if e1 == 0 {
                                    break;
                                }
                                e1 = *elemparts.wrapping_offset(e1 as isize);
                                if e1 == 0 {
                                    l += 1;
                                    continue;
                                }
                                k = l + 1;
                                while k <= maxneededtimes {
                                    e2 = *(*data).partitiontable[k as usize]
                                        .wrapping_offset(ind as isize);
                                    if e2 == 0 {
                                        break;
                                    }
                                    e2 = *elemparts.wrapping_offset(e2 as isize);
                                    if e2 == 0 {
                                        k += 1;
                                        continue;
                                    }
                                    *(*knows.wrapping_offset(e1 as isize))
                                        .wrapping_offset(e2 as isize) = TRUE;
                                    *(*knows.wrapping_offset(e2 as isize))
                                        .wrapping_offset(e1 as isize) = TRUE;
                                    k += 1;
                                }
                                l += 1;
                            }
                            j += 1;
                        }
                        hit = FALSE;
                        j = 1;
                        while j <= owners {
                            k = j + 1;
                            while k <= owners {
                                if *(*knows.wrapping_offset(j as isize)).wrapping_offset(k as isize)
                                    == 0
                                {
                                    hit += 1;
                                    i1 = *invelemparts.wrapping_offset(j as isize);
                                    i2 = *invelemparts.wrapping_offset(k as isize);
                                }
                                k += 1;
                            }
                            j += 1;
                        }
                        j = 1;
                        while j <= owners {
                            k = *invelemparts.wrapping_offset(j as isize);
                            *elemparts.wrapping_offset(k as isize) = 0;
                            j += 1;
                        }
                        j = 1;
                        while j <= owners {
                            k = 1;
                            while k <= owners {
                                *(*knows.wrapping_offset(j as isize)).wrapping_offset(k as isize) =
                                    FALSE;
                                k += 1;
                            }
                            j += 1;
                        }
                        if hit != 0 {
                            e1 = 0;
                            e2 = 0;
                            if info != 0
                                && hit + m > 2 {
                                    printf(
                                        c"Partitions %d and %d in element %d (%d owners) oddly related %d times\n"
                                            .as_ptr() as *const c_char,
                                        i1,
                                        i2,
                                        i,
                                        owners,
                                        hit,
                                    );
                                }
                            j = 0;
                            while j < nodesd2 {
                                ind = *(*(*data).topology.wrapping_offset(i as isize))
                                    .wrapping_offset(j as isize);
                                l = 1;
                                while l <= maxneededtimes {
                                    k = *(*data).partitiontable[l as usize]
                                        .wrapping_offset(ind as isize);
                                    if k == 0 {
                                        break;
                                    }
                                    if k == i1 {
                                        e1 += 1;
                                    }
                                    if k == i2 {
                                        e2 += 1;
                                    }
                                    l += 1;
                                }
                                j += 1;
                            }
                            j = 0;
                            while j < nodesd2 {
                                ind = *(*(*data).topology.wrapping_offset(i as isize))
                                    .wrapping_offset(j as isize);
                                k = *nodepart.wrapping_offset(ind as isize);
                                if (k == i1 && e1 < e2) || (k == i2 && e1 >= e2) {
                                    if noopt == 0 {
                                        *probnodes.wrapping_offset(ind as isize) += 1;
                                    }
                                    *nodepart.wrapping_offset(ind as isize) =
                                        *elempart.wrapping_offset(i as isize);
                                    *neededvector.wrapping_offset(
                                        *elempart.wrapping_offset(i as isize) as isize,
                                    ) += 1;
                                    *neededvector.wrapping_offset(k as isize) -= 1;
                                }
                                j += 1;
                            }
                            sharings += 1;
                        }
                        i += 1;
                    }
                    somethingdone += sharings;
                    if info != 0 && sharings != 0 {
                        printf(
                            c"Changed the ownership of %d nodes\n".as_ptr() as *const c_char,
                            sharings,
                        );
                    }
                }
                if !(sharings > 0 && m < 3) {
                    break;
                }
            }
            if info != 0 {
                if sharings != 0 {
                    printf(
                        c"%d problematic sharings may still exist\n".as_ptr() as *const c_char,
                        sharings,
                    );
                } else {
                    printf(
                        c"There shouldn't be any problematic sharings, knock, knock...\n".as_ptr()
                            as *const c_char,
                    );
                }
            }
        }
        // NOTE: `goto optimizeownership` via loop-continue.
        if noopt == 0 && m + n > 10 && optimize < 50 {
            optimize += 1;
            printf(
                c"Performing ownership optimization round %d\n".as_ptr() as *const c_char,
                optimize,
            );
            continue;
        }
        break;
    }
    free_ivector(neededvector, 1, partitions);
    if noopt == 0 {
        free_ivector(probnodes, 1, noknots);
        if partitions > 2 {
            free_ivector(elemparts, 1, partitions);
            free_ivector(invelemparts, 1, maxshared);
            free_imatrix(knows, 1, maxshared, 1, maxshared);
        }
    }
    if somethingdone != 0 {
        if info != 0 {
            printf(
                c"The partitioning was optimized: %d\n".as_ptr() as *const c_char,
                somethingdone,
            );
            printf(c"Checking partitioning after optimization\n".as_ptr() as *const c_char);
        }
        CheckPartitioning(data, info);
    } else {
        if info != 0 {
            printf(c"Partitioning was not altered\n".as_ptr() as *const c_char);
        }
    }
    0
}
/// # Safety
///
/// partitioned mesh export with chdir round-trip and binary paths.
/// `goto next_*_set` becomes loop-continue; the stale-`out` fwrite quirk at
/// the binary-nodes path is preserved exactly (same handle value, same fate).
#[no_mangle]
pub unsafe extern "C" fn SaveElmerInputPartitioned(
    data: *mut FemType,
    bound: *mut BoundaryType,
    prefix: *mut c_char,
    decimals: c_int,
    binary: c_int,
    parthalo: *mut c_int,
    parthypre: c_int,
    subparts: c_int,
    nooverwrite: c_int,
    info: c_int,
) -> c_int {
    let (noknots, noelements, mut sumsides, partitions): (c_int, c_int, c_int, c_int);
    // NOTE: C-dead locals found/l3/ind2/savethis dropped (never used).
    let (mut hit, mut parent, mut parent2, mut maxnosides): (c_int, c_int, c_int, c_int);
    let (mut nodesd2, mut nodesd1, mut discont, maxelemtype, minelemtype): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (mut sidehits, mut elemsides, mut side, mut bctype): (c_int, c_int, c_int, c_int);
    let (mut part, mut otherpart, mut part2, mut part3): (c_int, c_int, c_int, c_int);
    let mut elemtype: c_int;
    // NOTE: &mut-handoff init (see SaveBoundary note in egextra).
    let mut sideelemtype: c_int = 0;
    let (needednodes, neededtwice): (*mut c_int, *mut c_int);
    
    
    let (mut tottypes, mut splitsides): (c_int, c_int);
    let (mut i, mut j, mut k, mut l, mut l2, mut m, n): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let mut ind: c_int;
    let mut sideind: [c_int; MAXNODESD1] = [0; MAXNODESD1];
    // NOTE: C declares elemhit[] but never uses it; dropped.
    let mut filename: [c_char; MAXFILESIZE] = [0; MAXFILESIZE];
    let mut outstyle: [c_char; MAXFILESIZE] = [0; MAXFILESIZE];
    let mut directoryname: [c_char; MAXFILESIZE] = [0; MAXFILESIZE];
    let mut subdirectoryname: [c_char; MAXFILESIZE] = [0; MAXFILESIZE];
    let (neededtimes, elempart, elementsinpart, sidesinpart): (
        *mut c_int,
        *mut c_int,
        *mut c_int,
        *mut c_int,
    );
    let (mut maxneededtimes, mut bcneeded, mut trueparent, mut trueparent2): (
        c_int,
        c_int,
        c_int,
        c_int,
    );
    
    let (sharednodes, ownnodes): (*mut c_int, *mut c_int);
    
    // NOTE: NULL-inits in C; assigned iff reorder, read under the same flag.
    let mut order: *mut c_int = core::ptr::null_mut();
    let mut invorder: *mut c_int = core::ptr::null_mut();
    let (bcnode, bcnodedummy, elementhalo, neededtimes2): (
        *mut c_int,
        *mut c_int,
        *mut c_int,
        *mut c_int,
    );
    let (mut partstart, mut partfin, filesetsize): (c_int, c_int, c_int);
    let (mut nofile, mut nofile2, mut nobcnodes, mut hasbcnodes, mut halomode): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    );
    let (mut halobulkelems, mut halobcs): (c_int, c_int);
    // NOTE: C initializes fail=0 but mkdir overwrites before any read.
    let mut fail: c_int;
    let (mut immersed, mut halocopies, mut anyparthalo): (c_int, c_int, c_int);
    let (mut totsides, mut singleprec): (c_int, c_int);
    // NOTE: out stays NULL until the boundary section assigns it; the
    // binary-nodes fwrite below therefore uses the stale value exactly as C.
    let mut out: *mut FILE = core::ptr::null_mut();
    let mut outfiles: [*mut FILE; MAXPARTITIONS + 1] = [core::ptr::null_mut(); MAXPARTITIONS + 1];
    let (mut sumelementsinpart, mut sumownnodes, mut sumsharednodes, mut sumsidesinpart): (
        c_int,
        c_int,
        c_int,
        c_int,
    );
    if info != 0 {
        printf(c"Saving Elmer mesh in partitioned format\n".as_ptr() as *const c_char);
        if subparts != 0 {
            printf(
                c"There are %d subpartitions\n".as_ptr() as *const c_char,
                subparts,
            );
        }
    }
    if (*data).created == 0 {
        printf(c"You tried to save points that were never created.\n".as_ptr() as *const c_char);
        bigerror(c"No Elmer mesh files saved!".as_ptr() as *const c_char);
    }
    partitions = (*data).nopartitions;
    if partitions == 0 {
        printf(
            c"Tried to save partitioned format without partitions!\n".as_ptr() as *const c_char,
        );
        bigerror(c"No Elmer mesh files saved!".as_ptr() as *const c_char);
    }
    anyparthalo = FALSE;
    i = 0;
    while i < MAXHALOMODES as c_int {
        if *parthalo.wrapping_offset(i as isize) != 0 {
            printf(
                c"Saving halo elements mode active: %d\n".as_ptr() as *const c_char,
                i,
            );
            anyparthalo = TRUE;
        }
        i += 1;
    }
    if subparts < 1 && *parthalo.wrapping_offset(3) != 0 {
        printf(
            c"There can be no layer halo since there are no layers!\n".as_ptr() as *const c_char,
        );
        bigerror(c"No Elmer mesh files saved!".as_ptr() as *const c_char);
    }
    singleprec = FALSE;
    if binary == 2 {
        singleprec = TRUE;
    }
    elempart = (*data).elempart;
    let ownerpart: *mut c_int = (*data).nodepart;
    noelements = (*data).noelements;
    noknots = (*data).noknots;
    minelemtype = 101;
    maxelemtype = GetMaxElementType(data);
    halobulkelems = 0;
    needednodes = ivector(1, partitions);
    neededtwice = ivector(1, partitions);
    sharednodes = ivector(1, partitions);
    ownnodes = ivector(1, partitions);
    let sidetypes: *mut c_int = ivector(minelemtype, maxelemtype);
    let bulktypes: *mut *mut c_int = imatrix(1, partitions, minelemtype, maxelemtype);
    bcnodedummy = ivector(1, noknots);
    let reorder: c_int = parthypre;
    if reorder != 0 {
        order = ivector(1, noknots);
        k = 0;
        j = 1;
        while j <= partitions {
            i = 1;
            while i <= noknots {
                if *ownerpart.wrapping_offset(i as isize) == j {
                    k += 1;
                    *order.wrapping_offset(i as isize) = k;
                }
                i += 1;
            }
            j += 1;
        }
        invorder = ivector(1, noknots);
        i = 1;
        while i <= noknots {
            *invorder.wrapping_offset(*order.wrapping_offset(i as isize) as isize) = i;
            i += 1;
        }
        if info != 0 {
            printf(c"Created continuous parallel ordering\n".as_ptr() as *const c_char);
        }
    }
    bcnode = ivector(1, noknots);
    i = 1;
    while i <= noknots {
        *bcnode.wrapping_offset(i as isize) = FALSE;
        i += 1;
    }
    j = 0;
    while j < MAXBOUNDARIES as c_int {
        if (*bound.wrapping_offset(j as isize)).created == 0 {
            j += 1;
            continue;
        }
        i = 1;
        while i <= (*bound.wrapping_offset(j as isize)).nosides {
            GetBoundaryElement(
                i,
                bound.wrapping_offset(j as isize),
                data,
                sideind.as_mut_ptr(),
                &mut sideelemtype,
            );
            nodesd1 = sideelemtype % 100;
            l = 0;
            while l < nodesd1 {
                *bcnode.wrapping_offset(sideind[l as usize] as isize) = 1;
                l += 1;
            }
            i += 1;
        }
        j += 1;
    }
    nobcnodes = 0;
    i = 1;
    while i <= noknots {
        if *bcnode.wrapping_offset(i as isize) == 1 {
            nobcnodes += 1;
        }
        i += 1;
    }
    if info != 0 {
        printf(
            c"Number of boundary nodes at the boundary: %d\n".as_ptr() as *const c_char,
            nobcnodes,
        );
    }
    i = 1;
    while i <= noknots {
        *bcnodedummy.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    i = 1;
    while i <= noelements {
        k = *(*data).material.wrapping_offset(i as isize);
        elemtype = *(*data).elementtypes.wrapping_offset(i as isize);
        nodesd2 = elemtype % 100;
        j = 0;
        while j < nodesd2 {
            ind = *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize);
            if *bcnode.wrapping_offset(ind as isize) == 1 {
                j += 1;
                continue;
            }
            if *bcnodedummy.wrapping_offset(ind as isize) != 0 {
                if *bcnodedummy.wrapping_offset(ind as isize) != k {
                    *bcnode.wrapping_offset(ind as isize) = 2;
                }
            } else {
                *bcnodedummy.wrapping_offset(ind as isize) = k;
            }
            j += 1;
        }
        i += 1;
    }
    nobcnodes = 0;
    i = 1;
    while i <= noknots {
        if *bcnode.wrapping_offset(i as isize) == 2 {
            nobcnodes += 1;
        }
        i += 1;
    }
    if info != 0 {
        printf(
            c"Number of additional interface nodes: %d\n".as_ptr() as *const c_char,
            nobcnodes,
        );
    }
    sprintf(
        directoryname.as_mut_ptr(),
        c"%s".as_ptr() as *const c_char,
        prefix,
    );
    sprintf(
        subdirectoryname.as_mut_ptr(),
        c"%s.%d".as_ptr() as *const c_char,
        c"partitioning".as_ptr() as *const c_char,
        partitions,
    );
    // NOTE: MINGW32 undefined on this target; mkdir(path,0750) arms taken.
    fail = mkdir(directoryname.as_ptr(), 0o750);
    if info != 0 && fail == 0 {
        printf(
            c"Created mesh directory: %s\n".as_ptr() as *const c_char,
            directoryname.as_ptr(),
        );
    }
    // NOTE: C stores the chdir status into fail, but mkdir overwrites
    // before any read; the call is kept for its directory effect.
    chdir(directoryname.as_ptr());
    fail = mkdir(subdirectoryname.as_ptr(), 0o750);
    if fail != 0 {
        if info != 0 {
            printf(
                c"Reusing existing subdirectory: %s\n".as_ptr() as *const c_char,
                subdirectoryname.as_ptr(),
            );
        }
        if nooverwrite != 0 {
            printf(
                c"Mesh seems to already exist, writing is cancelled!\n".as_ptr() as *const c_char,
            );
            return 0;
        }
    } else {
        if info != 0 {
            printf(
                c"Created subdirectory: %s\n".as_ptr() as *const c_char,
                subdirectoryname.as_ptr(),
            );
        }
    }
    chdir(subdirectoryname.as_ptr());
    if info != 0 {
        printf(
            c"Saving mesh in parallel ElmerSolver format to directory %s/%s.\n".as_ptr()
                as *const c_char,
            directoryname.as_ptr(),
            subdirectoryname.as_ptr(),
        );
    }
    filesetsize = MAXPARTITIONS as c_int;
    if partitions > filesetsize && info != 0 {
        printf(
            c"Saving %d partitions in maximum sets of %d\n".as_ptr() as *const c_char,
            partitions,
            filesetsize,
        );
    }
    elementsinpart = ivector(1, partitions);
    sidesinpart = ivector(1, partitions);
    elementhalo = ivector(1, partitions);
    i = 1;
    while i <= partitions {
        *elementsinpart.wrapping_offset(i as isize) = 0;
        *sidesinpart.wrapping_offset(i as isize) = 0;
        *elementhalo.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    j = 1;
    while j <= partitions {
        i = minelemtype;
        while i <= maxelemtype {
            *(*bulktypes.wrapping_offset(j as isize)).wrapping_offset(i as isize) = 0;
            i += 1;
        }
        j += 1;
    }
    maxneededtimes = (*data).maxpartitiontable;
    neededtimes = ivector(1, noknots);
    i = 1;
    while i <= noknots {
        *neededtimes.wrapping_offset(i as isize) = 0;
        j = 1;
        while j <= maxneededtimes {
            if *(*data).partitiontable[j as usize].wrapping_offset(i as isize) != 0 {
                *neededtimes.wrapping_offset(i as isize) += 1;
            }
            j += 1;
        }
        i += 1;
    }
    if info != 0 {
        printf(
            c"Nodes belong to %d partitions in maximum\n".as_ptr() as *const c_char,
            maxneededtimes,
        );
    }
    partstart = 1;
    partfin = c_min_int(partitions, filesetsize);
    // NOTE: `next_elements_set:` backward goto becomes loop-continue.
    loop {
        part = partstart;
        while part <= partfin {
            if binary != 0 {
                sprintf(
                    filename.as_mut_ptr(),
                    c"%s.%d.%s".as_ptr() as *const c_char,
                    c"part".as_ptr() as *const c_char,
                    part,
                    c"elements.bin".as_ptr() as *const c_char,
                );
            } else {
                sprintf(
                    filename.as_mut_ptr(),
                    c"%s.%d.%s".as_ptr() as *const c_char,
                    c"part".as_ptr() as *const c_char,
                    part,
                    c"elements".as_ptr() as *const c_char,
                );
            }
            nofile = part - partstart + 1;
            outfiles[nofile as usize] = fopen(filename.as_ptr(), c"w".as_ptr() as *const c_char);
            part += 1;
        }
        i = 1;
        while i <= noelements {
            part = *elempart.wrapping_offset(i as isize);
            elemtype = *(*data).elementtypes.wrapping_offset(i as isize);
            nodesd2 = elemtype % 100;
            if part >= partstart && part <= partfin {
                nofile = part - partstart + 1;
                *(*bulktypes.wrapping_offset(part as isize)).wrapping_offset(elemtype as isize) +=
                    1;
                *elementsinpart.wrapping_offset(part as isize) += 1;
                if binary != 0 {
                    fwrite(
                        &i as *const c_int as *const c_void,
                        size_of::<c_int>() as size_t,
                        1,
                        outfiles[nofile as usize],
                    );
                    k = -1;
                    fwrite(
                        &k as *const c_int as *const c_void,
                        size_of::<c_int>() as size_t,
                        1,
                        outfiles[nofile as usize],
                    );
                    fwrite(
                        &mut *(*data).material.wrapping_offset(i as isize) as *mut c_int
                            as *const c_void,
                        size_of::<c_int>() as size_t,
                        1,
                        outfiles[nofile as usize],
                    );
                    fwrite(
                        &elemtype as *const c_int as *const c_void,
                        size_of::<c_int>() as size_t,
                        1,
                        outfiles[nofile as usize],
                    );
                    fwrite(
                        &mut *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(0)
                            as *mut c_int as *const c_void,
                        size_of::<c_int>() as size_t,
                        nodesd2 as size_t,
                        outfiles[nofile as usize],
                    );
                } else {
                    fprintf(
                        outfiles[nofile as usize],
                        c"%d %d %d ".as_ptr() as *const c_char,
                        i,
                        *(*data).material.wrapping_offset(i as isize),
                        elemtype,
                    );
                    j = 0;
                    while j < nodesd2 {
                        ind = *(*(*data).topology.wrapping_offset(i as isize))
                            .wrapping_offset(j as isize);
                        if reorder != 0 {
                            ind = *order.wrapping_offset(ind as isize);
                        }
                        fprintf(
                            outfiles[nofile as usize],
                            c"%d ".as_ptr() as *const c_char,
                            ind,
                        );
                        j += 1;
                    }
                    fprintf(outfiles[nofile as usize], c"\n".as_ptr() as *const c_char);
                }
            }
            if anyparthalo == 0 {
                i += 1;
                continue;
            }
            halocopies = 0;
            halomode = 1;
            while halomode <= 4 {
                if *parthalo.wrapping_offset(halomode as isize) == 0 {
                    halomode += 1;
                    continue;
                }
                if halomode == 3
                    && part <= subparts {
                        // NOTE: block-scoped C local becomes a Rust block.
                        {
                            let mut leftright: c_int;
                            leftright = -1;
                            while leftright <= 1 {
                                part2 = part + leftright;
                                if part2 < 1 || part2 > subparts {
                                    leftright += 2;
                                    continue;
                                }
                                *elementhalo.wrapping_offset(part2 as isize) = i;
                                halocopies += 1;
                                leftright += 2;
                            }
                        }
                    }
                otherpart = 0;
                hasbcnodes = FALSE;
                j = 0;
                while j < nodesd2 {
                    ind = *(*(*data).topology.wrapping_offset(i as isize))
                        .wrapping_offset(j as isize);
                    if *neededtimes.wrapping_offset(ind as isize) > 1 {
                        otherpart += 1;
                    }
                    if *bcnode.wrapping_offset(ind as isize) != 0 {
                        hasbcnodes = TRUE;
                    }
                    j += 1;
                }
                if halomode == 1 && otherpart != 0 {
                    elemsides = elemtype / 100;
                    if elemsides == 8 {
                        immersed = (otherpart >= 4) as c_int;
                        elemsides = 6;
                    } else if elemsides == 7 || elemsides == 6 {
                        immersed = (otherpart >= 3) as c_int;
                        elemsides = 5;
                    } else if elemsides == 5 {
                        immersed = (otherpart >= 3) as c_int;
                        elemsides = 4;
                    } else {
                        immersed = (otherpart >= 2) as c_int;
                    }
                    if immersed != 0 {
                        side = 0;
                        while side < elemsides {
                            GetElementSide(
                                i,
                                side,
                                1,
                                data,
                                sideind.as_mut_ptr(),
                                &mut sideelemtype,
                            );
                            l = 1;
                            while l <= *neededtimes.wrapping_offset(sideind[0] as isize) {
                                part2 = *(*data).partitiontable[l as usize]
                                    .wrapping_offset(sideind[0] as isize);
                                if part2 == part {
                                    l += 1;
                                    continue;
                                }
                                if *elementhalo.wrapping_offset(part2 as isize) == i {
                                    l += 1;
                                    continue;
                                }
                                if part2 < partstart || part2 > partfin {
                                    l += 1;
                                    continue;
                                }
                                sidehits = 1;
                                k = 1;
                                while k < sideelemtype % 100 {
                                    l2 = 1;
                                    while l2
                                        <= *neededtimes
                                            .wrapping_offset(sideind[k as usize] as isize)
                                    {
                                        part3 = *(*data).partitiontable[l2 as usize]
                                            .wrapping_offset(sideind[k as usize] as isize);
                                        if part2 == part3 {
                                            sidehits += 1;
                                        }
                                        l2 += 1;
                                    }
                                    k += 1;
                                }
                                if sidehits < sideelemtype % 100 {
                                    l += 1;
                                    continue;
                                }
                                // NOTE: `if(0)` diagnostic omitted.
                                *elementhalo.wrapping_offset(part2 as isize) = i;
                                halocopies += 1;
                                l += 1;
                            }
                            side += 1;
                        }
                    }
                }
                if halomode == 2 && hasbcnodes != 0 {
                    j = 0;
                    while j < nodesd2 {
                        ind = *(*(*data).topology.wrapping_offset(i as isize))
                            .wrapping_offset(j as isize);
                        if *bcnode.wrapping_offset(ind as isize) == 0 {
                            j += 1;
                            continue;
                        }
                        k = 1;
                        while k <= *neededtimes.wrapping_offset(ind as isize) {
                            part2 =
                                *(*data).partitiontable[k as usize].wrapping_offset(ind as isize);
                            if part == part2 {
                                k += 1;
                                continue;
                            }
                            if *elementhalo.wrapping_offset(part2 as isize) == i {
                                k += 1;
                                continue;
                            }
                            if part2 < partstart || part2 > partfin {
                                k += 1;
                                continue;
                            }
                            // NOTE: `if(0)` diagnostic omitted.
                            *elementhalo.wrapping_offset(part2 as isize) = i;
                            halocopies += 1;
                            k += 1;
                        }
                        j += 1;
                    }
                }
                if halomode == 4 && otherpart != 0 {
                    j = 0;
                    while j < nodesd2 {
                        ind = *(*(*data).topology.wrapping_offset(i as isize))
                            .wrapping_offset(j as isize);
                        if *neededtimes.wrapping_offset(ind as isize) == 1 {
                            j += 1;
                            continue;
                        }
                        l = 1;
                        while l <= *neededtimes.wrapping_offset(ind as isize) {
                            part2 =
                                *(*data).partitiontable[l as usize].wrapping_offset(ind as isize);
                            if part2 == part {
                                l += 1;
                                continue;
                            }
                            if *elementhalo.wrapping_offset(part2 as isize) == i {
                                l += 1;
                                continue;
                            }
                            if part2 < partstart || part2 > partfin {
                                l += 1;
                                continue;
                            }
                            *elementhalo.wrapping_offset(part2 as isize) = i;
                            halocopies += 1;
                            l += 1;
                        }
                        j += 1;
                    }
                }
                halomode += 1;
            }
            if halocopies != 0 {
                halobulkelems += halocopies;
                part2 = partstart;
                while part2 <= partfin {
                    if part == part2 {
                        part2 += 1;
                        continue;
                    }
                    if *elementhalo.wrapping_offset(part2 as isize) != i {
                        part2 += 1;
                        continue;
                    }
                    nofile2 = part2 - partstart + 1;
                    if binary != 0 {
                        fwrite(
                            &i as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            outfiles[nofile2 as usize],
                        );
                        fwrite(
                            &part as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            outfiles[nofile2 as usize],
                        );
                        fwrite(
                            &mut *(*data).material.wrapping_offset(i as isize) as *mut c_int
                                as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            outfiles[nofile2 as usize],
                        );
                        fwrite(
                            &elemtype as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            outfiles[nofile2 as usize],
                        );
                        j = 0;
                        while j < nodesd2 {
                            ind = *(*(*data).topology.wrapping_offset(i as isize))
                                .wrapping_offset(j as isize);
                            if reorder != 0 {
                                ind = *order.wrapping_offset(ind as isize);
                            }
                            fwrite(
                                &ind as *const c_int as *const c_void,
                                size_of::<c_int>() as size_t,
                                1,
                                outfiles[nofile2 as usize],
                            );
                            j += 1;
                        }
                    } else {
                        fprintf(
                            outfiles[nofile2 as usize],
                            c"%d/%d %d %d ".as_ptr() as *const c_char,
                            i,
                            part,
                            *(*data).material.wrapping_offset(i as isize),
                            elemtype,
                        );
                        j = 0;
                        while j < nodesd2 {
                            ind = *(*(*data).topology.wrapping_offset(i as isize))
                                .wrapping_offset(j as isize);
                            if reorder != 0 {
                                ind = *order.wrapping_offset(ind as isize);
                            }
                            fprintf(
                                outfiles[nofile2 as usize],
                                c"%d ".as_ptr() as *const c_char,
                                ind,
                            );
                            j += 1;
                        }
                        fprintf(
                            outfiles[nofile2 as usize],
                            c"\n".as_ptr() as *const c_char,
                        );
                    }
                    *(*bulktypes.wrapping_offset(part2 as isize))
                        .wrapping_offset(elemtype as isize) += 1;
                    *elementsinpart.wrapping_offset(part2 as isize) += 1;
                    j = 0;
                    while j < nodesd2 {
                        ind = *(*(*data).topology.wrapping_offset(i as isize))
                            .wrapping_offset(j as isize);
                        hit = FALSE;
                        k = 1;
                        while k <= maxneededtimes {
                            part3 =
                                *(*data).partitiontable[k as usize].wrapping_offset(ind as isize);
                            if part3 == part2 {
                                hit = TRUE;
                            }
                            if part3 == 0 {
                                break;
                            }
                            k += 1;
                        }
                        if hit == 0 {
                            if k <= maxneededtimes {
                                *(*data).partitiontable[k as usize].wrapping_offset(ind as isize) =
                                    part2;
                            } else {
                                maxneededtimes += 1;
                                // NOTE: `if(0)` diagnostic omitted.
                                (*data).partitiontable[k as usize] = ivector(1, noknots);
                                m = 1;
                                while m <= noknots {
                                    *(*data).partitiontable[k as usize]
                                        .wrapping_offset(m as isize) = 0;
                                    m += 1;
                                }
                                *(*data).partitiontable[k as usize].wrapping_offset(ind as isize) =
                                    part2;
                            }
                        }
                        j += 1;
                    }
                    part2 += 1;
                }
            }
            i += 1;
        }
        if anyparthalo != 0 && info != 0 {
            printf(
                c"There are %d bulk elements in the halo.\n".as_ptr() as *const c_char,
                halobulkelems,
            );
        }
        part = partstart;
        while part <= partfin {
            nofile = part - partstart + 1;
            fclose(outfiles[nofile as usize]);
            part += 1;
        }
        if partfin < partitions {
            partstart = partfin + 1;
            partfin = c_min_int(partfin + filesetsize, partitions);
            // NOTE: `goto next_elements_set` via loop-continue.
            continue;
        }
        break;
    }
    // NOTE: `part.n.elements saved` comment position.
    if anyparthalo != 0 {
        // NOTE: block-scoped C locals become a Rust block.
        {
            let mut halonodes: c_int;
            neededtimes2 = ivector(1, noknots);
            halonodes = 0;
            // NOTE: C accumulates k/l here but never reads them back (final
            // values discarded; later sections reassign first). The counts
            // and stores are dropped; the neededtimes2 work stays.
            i = 1;
            while i <= noknots {
                *neededtimes2.wrapping_offset(i as isize) = 0;
                j = 1;
                while j <= maxneededtimes {
                    if *(*data).partitiontable[j as usize].wrapping_offset(i as isize) != 0 {
                        *neededtimes2.wrapping_offset(i as isize) += 1;
                    }
                    j += 1;
                }
                halonodes += *neededtimes2.wrapping_offset(i as isize)
                    - *neededtimes.wrapping_offset(i as isize);
                i += 1;
            }
            if (*data).maxpartitiontable < maxneededtimes {
                (*data).maxpartitiontable = maxneededtimes;
                if info != 0 {
                    printf(
                        c"With the halos nodes belong to %d partitions in maximum\n".as_ptr()
                            as *const c_char,
                        maxneededtimes,
                    );
                }
            }
            if info != 0 {
                printf(
                    c"There are %d additional shared nodes resulting from the halo.\n".as_ptr()
                        as *const c_char,
                    halonodes,
                );
            }
        }
    } else {
        neededtimes2 = neededtimes;
    }
    if (*data).dim == 2 {
        sprintf(
            outstyle.as_mut_ptr(),
            c"%%d %%d %%.%dg %%.%dg 0.0\n".as_ptr() as *const c_char,
            decimals,
            decimals,
        );
    } else {
        sprintf(
            outstyle.as_mut_ptr(),
            c"%%d %%d %%.%dg %%.%dg %%.%dg\n".as_ptr() as *const c_char,
            decimals,
            decimals,
            decimals,
        );
    }
    if info != 0 {
        printf(
            c"Saving mesh for %d partitions\n".as_ptr() as *const c_char,
            partitions,
        );
    }
    partstart = 1;
    partfin = c_min_int(partitions, filesetsize);
    i = 1;
    while i <= partitions {
        *needednodes.wrapping_offset(i as isize) = 0;
        *neededtwice.wrapping_offset(i as isize) = 0;
        *sharednodes.wrapping_offset(i as isize) = 0;
        *ownnodes.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    // NOTE: `next_nodes_set:` backward goto becomes loop-continue.
    loop {
        part = partstart;
        while part <= partfin {
            if binary != 0 {
                if singleprec != 0 {
                    sprintf(
                        filename.as_mut_ptr(),
                        c"%s.%d.%s".as_ptr() as *const c_char,
                        c"part".as_ptr() as *const c_char,
                        part,
                        c"nodes.sbin".as_ptr() as *const c_char,
                    );
                } else {
                    sprintf(
                        filename.as_mut_ptr(),
                        c"%s.%d.%s".as_ptr() as *const c_char,
                        c"part".as_ptr() as *const c_char,
                        part,
                        c"nodes.bin".as_ptr() as *const c_char,
                    );
                }
            } else {
                sprintf(
                    filename.as_mut_ptr(),
                    c"%s.%d.%s".as_ptr() as *const c_char,
                    c"part".as_ptr() as *const c_char,
                    part,
                    c"nodes".as_ptr() as *const c_char,
                );
            }
            nofile = part - partstart + 1;
            outfiles[nofile as usize] = fopen(filename.as_ptr(), c"w".as_ptr() as *const c_char);
            part += 1;
        }
        l = 1;
        while l <= noknots {
            i = l;
            if reorder != 0 {
                i = *invorder.wrapping_offset(l as isize);
            }
            // NOTE: commented-out neededtimes2 loop kept as maxneededtimes
            // loop exactly as compiled in C.
            j = 1;
            while j <= maxneededtimes {
                k = *(*data).partitiontable[j as usize].wrapping_offset(i as isize);
                if k == 0 {
                    break;
                }
                if k < partstart || k > partfin {
                    j += 1;
                    continue;
                }
                nofile = k - partstart + 1;
                ind = i;
                if reorder != 0 {
                    ind = *order.wrapping_offset(i as isize);
                }
                if binary != 0 {
                    // NOTE: block-scoped C arrays become a Rust block.
                    {
                        let mut coords: [c_double; 3] = [0.0; 3];
                        let mut scoords: [c_float; 3] = [0.0; 3];
                        // NOTE: C writes the index to `out` (stale handle
                        // from an earlier section, possibly NULL) — the
                        // pristine bug is preserved exactly.
                        fwrite(
                            &ind as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        if singleprec != 0 {
                            scoords[0] = *(*data).x.wrapping_offset(i as isize) as c_float;
                            scoords[1] = *(*data).y.wrapping_offset(i as isize) as c_float;
                            scoords[2] = *(*data).z.wrapping_offset(i as isize) as c_float;
                            fwrite(
                                scoords.as_ptr() as *const c_void,
                                size_of::<c_float>() as size_t,
                                3,
                                outfiles[nofile as usize],
                            );
                        } else {
                            coords[0] = *(*data).x.wrapping_offset(i as isize);
                            coords[1] = *(*data).y.wrapping_offset(i as isize);
                            coords[2] = *(*data).z.wrapping_offset(i as isize);
                            fwrite(
                                coords.as_ptr() as *const c_void,
                                size_of::<c_double>() as size_t,
                                3,
                                outfiles[nofile as usize],
                            );
                        }
                    }
                } else {
                    if (*data).dim == 2 {
                        fprintf(
                            outfiles[nofile as usize],
                            outstyle.as_ptr(),
                            ind,
                            -1,
                            *(*data).x.wrapping_offset(i as isize),
                            *(*data).y.wrapping_offset(i as isize),
                        );
                    } else if (*data).dim == 3 {
                        fprintf(
                            outfiles[nofile as usize],
                            outstyle.as_ptr(),
                            ind,
                            -1,
                            *(*data).x.wrapping_offset(i as isize),
                            *(*data).y.wrapping_offset(i as isize),
                            *(*data).z.wrapping_offset(i as isize),
                        );
                    }
                }
                *needednodes.wrapping_offset(k as isize) += 1;
                if k == *ownerpart.wrapping_offset(i as isize) {
                    *ownnodes.wrapping_offset(k as isize) += 1;
                } else {
                    *sharednodes.wrapping_offset(k as isize) += 1;
                }
                j += 1;
            }
            l += 1;
        }
        part = partstart;
        while part <= partfin {
            nofile = part - partstart + 1;
            fclose(outfiles[nofile as usize]);
            part += 1;
        }
        if partfin < partitions {
            partstart = partfin + 1;
            partfin = c_min_int(partfin + filesetsize, partitions);
            // NOTE: `goto next_nodes_set` via loop-continue.
            continue;
        }
        break;
    }
    // NOTE: `part.n.nodes saved` comment position.
    partstart = 1;
    partfin = c_min_int(partitions, filesetsize);
    // NOTE: `next_shared_set:` backward goto becomes loop-continue.
    loop {
        part = partstart;
        while part <= partfin {
            sprintf(
                filename.as_mut_ptr(),
                c"%s.%d.%s".as_ptr() as *const c_char,
                c"part".as_ptr() as *const c_char,
                part,
                c"shared".as_ptr() as *const c_char,
            );
            nofile = part - partstart + 1;
            outfiles[nofile as usize] = fopen(filename.as_ptr(), c"w".as_ptr() as *const c_char);
            part += 1;
        }
        l = 1;
        while l <= noknots {
            i = l;
            if reorder != 0 {
                i = *invorder.wrapping_offset(l as isize);
            }
            if *neededtimes2.wrapping_offset(i as isize) <= 1 {
                l += 1;
                continue;
            }
            j = 1;
            while j <= *neededtimes2.wrapping_offset(i as isize) {
                k = *(*data).partitiontable[j as usize].wrapping_offset(i as isize);
                // NOTE: commented-out parthalo[2] arm omitted (never compiled).
                if k < partstart || k > partfin {
                    j += 1;
                    continue;
                }
                nofile = k - partstart + 1;
                ind = i;
                if reorder != 0 {
                    ind = *order.wrapping_offset(i as isize);
                }
                *neededtwice.wrapping_offset(k as isize) += 1;
                // NOTE: commented-out parthalo[2] arm omitted (never compiled).
                fprintf(
                    outfiles[nofile as usize],
                    c"%d %d %d".as_ptr() as *const c_char,
                    ind,
                    *neededtimes2.wrapping_offset(i as isize),
                    *ownerpart.wrapping_offset(i as isize),
                );
                m = 1;
                while m <= *neededtimes2.wrapping_offset(i as isize) {
                    if *(*data).partitiontable[m as usize].wrapping_offset(i as isize)
                        != *ownerpart.wrapping_offset(i as isize)
                    {
                        fprintf(
                            outfiles[nofile as usize],
                            c" %d".as_ptr() as *const c_char,
                            *(*data).partitiontable[m as usize].wrapping_offset(i as isize),
                        );
                    }
                    m += 1;
                }
                fprintf(outfiles[nofile as usize], c"\n".as_ptr() as *const c_char);
                j += 1;
            }
            l += 1;
        }
        part = partstart;
        while part <= partfin {
            nofile = part - partstart + 1;
            fclose(outfiles[nofile as usize]);
            part += 1;
        }
        if partfin < partitions {
            partstart = partfin + 1;
            partfin = c_min_int(partfin + filesetsize, partitions);
            // NOTE: `goto next_shared_set` via loop-continue.
            continue;
        }
        break;
    }
    // NOTE: `part.n.shared saved` comment position.
    // NOTE: C stores discont = FALSE here; reassigned before any read.
    splitsides = 0;
    maxnosides = 0;
    j = 0;
    while j < MAXBOUNDARIES as c_int {
        maxnosides = c_max_int(maxnosides, (*bound.wrapping_offset(j as isize)).nosides);
        j += 1;
    }
    halobcs = 0;
    part = 1;
    while part <= partitions {
        // NOTE: block-scoped C locals become a Rust block; the C-dead
        // bcneeded2/closeparent/closeparent2/haloelem are dropped.
        {
            let mut saveelem: c_int;
            if binary != 0 {
                sprintf(
                    filename.as_mut_ptr(),
                    c"%s.%d.%s".as_ptr() as *const c_char,
                    c"part".as_ptr() as *const c_char,
                    part,
                    c"boundary.bin".as_ptr() as *const c_char,
                );
            } else {
                sprintf(
                    filename.as_mut_ptr(),
                    c"%s.%d.%s".as_ptr() as *const c_char,
                    c"part".as_ptr() as *const c_char,
                    part,
                    c"boundary".as_ptr() as *const c_char,
                );
            }
            out = fopen(filename.as_ptr(), c"w".as_ptr() as *const c_char);
            i = minelemtype;
            while i <= maxelemtype {
                *sidetypes.wrapping_offset(i as isize) = 0;
                i += 1;
            }
            sumsides = 0;
            totsides = 0;
            j = 0;
            while j < MAXBOUNDARIES as c_int {
                if (*bound.wrapping_offset(j as isize)).created == 0 {
                    j += 1;
                    continue;
                }
                if (*bound.wrapping_offset(j as isize)).nosides == 0 {
                    j += 1;
                    continue;
                }
                i = 1;
                while i <= (*bound.wrapping_offset(j as isize)).nosides {
                    GetBoundaryElement(
                        i,
                        bound.wrapping_offset(j as isize),
                        data,
                        sideind.as_mut_ptr(),
                        &mut sideelemtype,
                    );
                    bctype = *(*bound.wrapping_offset(j as isize))
                        .types
                        .wrapping_offset(i as isize);
                    nodesd1 = sideelemtype % 100;
                    parent = *(*bound.wrapping_offset(j as isize))
                        .parent
                        .wrapping_offset(i as isize);
                    parent2 = *(*bound.wrapping_offset(j as isize))
                        .parent2
                        .wrapping_offset(i as isize);
                    totsides += 1;
                    trueparent = FALSE;
                    trueparent2 = FALSE;
                    if parent != 0 {
                        trueparent = (*elempart.wrapping_offset(parent as isize) == part) as c_int;
                    }
                    if parent2 != 0 {
                        trueparent2 =
                            (*elempart.wrapping_offset(parent2 as isize) == part) as c_int;
                    }
                    saveelem = FALSE;
                    part2 = 0;
                    if trueparent != 0 || trueparent2 != 0 {
                        if parent != 0 && trueparent == 0 {
                            splitsides += 1;
                            if *parthalo.wrapping_offset(1) == 0
                                && *parthalo.wrapping_offset(2) == 0
                                && *parthalo.wrapping_offset(4) == 0
                            {
                                parent = 0;
                            }
                        } else if parent2 != 0 && trueparent2 == 0 {
                            splitsides += 1;
                            if *parthalo.wrapping_offset(1) == 0
                                && *parthalo.wrapping_offset(2) == 0
                                && *parthalo.wrapping_offset(4) == 0
                            {
                                parent2 = 0;
                            }
                        }
                        saveelem = TRUE;
                    }
                    if saveelem == 0 && *parthalo.wrapping_offset(3) != 0 && part <= subparts {
                        // NOTE: block-scoped C local becomes a Rust block.
                        {
                            let mut leftright: c_int;
                            leftright = -1;
                            while leftright <= 1 {
                                part2 = part + leftright;
                                if part2 < 1 || part2 > subparts {
                                    leftright += 2;
                                    continue;
                                }
                                if parent != 0 {
                                    trueparent = (*elempart.wrapping_offset(parent as isize)
                                        == part2)
                                        as c_int;
                                }
                                if parent2 != 0 {
                                    trueparent2 = (*elempart.wrapping_offset(parent2 as isize)
                                        == part2)
                                        as c_int;
                                }
                                if trueparent != 0 || trueparent2 != 0 {
                                    if parent != 0 && trueparent == 0 {
                                        parent = 0;
                                    } else if parent2 != 0 && trueparent2 == 0 {
                                        parent2 = 0;
                                    }
                                    saveelem = TRUE;
                                    break;
                                }
                                leftright += 2;
                            }
                        }
                        if saveelem != 0 {
                            halobcs += 1;
                        }
                    }
                    if saveelem == 0
                        && (*parthalo.wrapping_offset(2) != 0 || *parthalo.wrapping_offset(4) != 0)
                    {
                        bcneeded = 0;
                        l = 0;
                        while l < nodesd1 {
                            ind = sideind[l as usize];
                            k = 1;
                            while k <= *neededtimes.wrapping_offset(ind as isize) {
                                if part
                                    == *(*data).partitiontable[k as usize]
                                        .wrapping_offset(ind as isize)
                                {
                                    bcneeded += 1;
                                }
                                k += 1;
                            }
                            l += 1;
                        }
                        saveelem = bcneeded;
                        if saveelem != 0 {
                            halobcs += 1;
                        }
                    }
                    if saveelem == 0 {
                        i += 1;
                        continue;
                    }
                    sumsides += 1;
                    *sidetypes.wrapping_offset(sideelemtype as isize) += 1;
                    if binary != 0 {
                        k = part2;
                        if k == 0 {
                            k = -1;
                        }
                        fwrite(
                            &totsides as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &k as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &bctype as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &parent as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &parent2 as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &sideelemtype as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        l = 0;
                        while l < nodesd1 {
                            k = sideind[l as usize];
                            if reorder != 0 {
                                k = *order.wrapping_offset(k as isize);
                            }
                            fwrite(
                                &k as *const c_int as *const c_void,
                                size_of::<c_int>() as size_t,
                                1,
                                out,
                            );
                            l += 1;
                        }
                    } else {
                        if part2 != 0 {
                            fprintf(
                                out,
                                c"%d/%d %d %d %d %d".as_ptr() as *const c_char,
                                totsides,
                                part2,
                                bctype,
                                parent,
                                parent2,
                                sideelemtype,
                            );
                        } else {
                            fprintf(
                                out,
                                c"%d %d %d %d %d".as_ptr() as *const c_char,
                                totsides,
                                bctype,
                                parent,
                                parent2,
                                sideelemtype,
                            );
                        }
                        l = 0;
                        while l < nodesd1 {
                            k = sideind[l as usize];
                            if reorder != 0 {
                                k = *order.wrapping_offset(k as isize);
                            }
                            fprintf(out, c" %d".as_ptr() as *const c_char, k);
                            l += 1;
                        }
                        fprintf(out, c"\n".as_ptr() as *const c_char);
                    }
                    i += 1;
                }
                // NOTE: C `for` increment (fall-through path; the guard
                // arms increment before their `continue`s).
                j += 1;
            }
            // NOTE: C `for(j=0;...)` re-init for the discont pass.
            j = 0;
            while j < MAXBOUNDARIES as c_int {
                if (*bound.wrapping_offset(j as isize)).created == 0 {
                    j += 1;
                    continue;
                }
                i = 1;
                while i <= (*bound.wrapping_offset(j as isize)).nosides {
                    if (*bound.wrapping_offset(j as isize)).ediscont != 0 {
                        discont = *(*bound.wrapping_offset(j as isize))
                            .discont
                            .wrapping_offset(i as isize);
                    } else {
                        discont = FALSE;
                    }
                    if *(*bound.wrapping_offset(j as isize))
                        .parent2
                        .wrapping_offset(i as isize)
                        == 0
                        || discont == 0
                    {
                        i += 1;
                        continue;
                    }
                    GetElementSide(
                        *(*bound.wrapping_offset(j as isize))
                            .parent2
                            .wrapping_offset(i as isize),
                        *(*bound.wrapping_offset(j as isize))
                            .side2
                            .wrapping_offset(i as isize),
                        -*(*bound.wrapping_offset(j as isize))
                            .normal
                            .wrapping_offset(i as isize),
                        data,
                        sideind.as_mut_ptr(),
                        &mut sideelemtype,
                    );
                    nodesd1 = sideelemtype % 100;
                    totsides += 1;
                    bcneeded = 0;
                    l = 0;
                    while l < nodesd1 {
                        ind = sideind[l as usize];
                        k = 1;
                        while k <= *neededtimes.wrapping_offset(ind as isize) {
                            if part
                                == *(*data).partitiontable[k as usize].wrapping_offset(ind as isize)
                            {
                                bcneeded += 1;
                            }
                            k += 1;
                        }
                        l += 1;
                    }
                    if bcneeded < nodesd1 {
                        i += 1;
                        continue;
                    }
                    trueparent = (*elempart.wrapping_offset(
                        *(*bound.wrapping_offset(j as isize))
                            .parent2
                            .wrapping_offset(i as isize) as isize,
                    ) == part) as c_int;
                    if trueparent == 0 {
                        i += 1;
                        continue;
                    }
                    sumsides += 1;
                    *sidetypes.wrapping_offset(sideelemtype as isize) += 1;
                    if binary != 0 {
                        k = -1;
                        fwrite(
                            &totsides as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &k as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &mut *(*bound.wrapping_offset(j as isize))
                                .types
                                .wrapping_offset(i as isize)
                                as *mut c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &mut *(*bound.wrapping_offset(j as isize))
                                .parent2
                                .wrapping_offset(i as isize)
                                as *mut c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &mut *(*bound.wrapping_offset(j as isize))
                                .parent
                                .wrapping_offset(i as isize)
                                as *mut c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        fwrite(
                            &sideelemtype as *const c_int as *const c_void,
                            size_of::<c_int>() as size_t,
                            1,
                            out,
                        );
                        l = 0;
                        while l < nodesd1 {
                            k = sideind[l as usize];
                            if reorder != 0 {
                                k = *order.wrapping_offset(k as isize);
                            }
                            fwrite(
                                &k as *const c_int as *const c_void,
                                size_of::<c_int>() as size_t,
                                1,
                                out,
                            );
                            l += 1;
                        }
                    } else {
                        fprintf(
                            out,
                            c"%d %d %d %d %d".as_ptr() as *const c_char,
                            totsides,
                            *(*bound.wrapping_offset(j as isize))
                                .types
                                .wrapping_offset(i as isize),
                            *(*bound.wrapping_offset(j as isize))
                                .parent2
                                .wrapping_offset(i as isize),
                            *(*bound.wrapping_offset(j as isize))
                                .parent
                                .wrapping_offset(i as isize),
                            sideelemtype,
                        );
                        l = 0;
                        while l < nodesd1 {
                            k = sideind[l as usize];
                            if reorder != 0 {
                                k = *order.wrapping_offset(k as isize);
                            }
                            fprintf(out, c" %d".as_ptr() as *const c_char, k);
                            l += 1;
                        }
                        fprintf(out, c"\n".as_ptr() as *const c_char);
                    }
                    i += 1;
                }
                // NOTE: C `for` increment (fall-through path).
                j += 1;
            }
            *sidesinpart.wrapping_offset(part as isize) = sumsides;
            fclose(out);
            // NOTE: `part.n.boundary` end comment position.
            tottypes = 0;
            i = minelemtype;
            while i <= maxelemtype {
                if *(*bulktypes.wrapping_offset(part as isize)).wrapping_offset(i as isize) != 0 {
                    tottypes += 1;
                }
                if *sidetypes.wrapping_offset(i as isize) != 0 {
                    tottypes += 1;
                }
                i += 1;
            }
            sprintf(
                filename.as_mut_ptr(),
                c"%s.%d.%s".as_ptr() as *const c_char,
                c"part".as_ptr() as *const c_char,
                part,
                c"header".as_ptr() as *const c_char,
            );
            out = fopen(filename.as_ptr(), c"w".as_ptr() as *const c_char);
            fprintf(
                out,
                c"%-6d %-6d %-6d\n".as_ptr() as *const c_char,
                *needednodes.wrapping_offset(part as isize),
                *elementsinpart.wrapping_offset(part as isize),
                sumsides,
            );
            fprintf(out, c"%-6d\n".as_ptr() as *const c_char, tottypes);
            i = minelemtype;
            while i <= maxelemtype {
                if *(*bulktypes.wrapping_offset(part as isize)).wrapping_offset(i as isize) != 0 {
                    fprintf(
                        out,
                        c"%-6d %-6d\n".as_ptr() as *const c_char,
                        i,
                        *(*bulktypes.wrapping_offset(part as isize)).wrapping_offset(i as isize),
                    );
                }
                i += 1;
            }
            i = minelemtype;
            while i <= maxelemtype {
                if *sidetypes.wrapping_offset(i as isize) != 0 {
                    fprintf(
                        out,
                        c"%-6d %-6d\n".as_ptr() as *const c_char,
                        i,
                        *sidetypes.wrapping_offset(i as isize),
                    );
                }
                i += 1;
            }
            fprintf(
                out,
                c"%-6d %-6d\n".as_ptr() as *const c_char,
                *neededtwice.wrapping_offset(part as isize),
                0,
            );
            fclose(out);
            if info != 0 {
                if part == 1 {
                    printf(
                        c"   %-5s %-10s %-10s %-8s %-8s\n".as_ptr() as *const c_char,
                        c"part".as_ptr() as *const c_char,
                        c"elements".as_ptr() as *const c_char,
                        c"nodes".as_ptr() as *const c_char,
                        c"shared".as_ptr() as *const c_char,
                        c"bc elems".as_ptr() as *const c_char,
                    );
                }
                printf(
                    c"   %-5d %-10d %-10d %-8d %-8d\n".as_ptr() as *const c_char,
                    part,
                    *elementsinpart.wrapping_offset(part as isize),
                    *ownnodes.wrapping_offset(part as isize),
                    *sharednodes.wrapping_offset(part as isize),
                    *sidesinpart.wrapping_offset(part as isize),
                );
            }
            part += 1;
        }
    }
    // NOTE: `part.n.header` end comment position.
    sumelementsinpart = 0;
    sumownnodes = 0;
    sumsharednodes = 0;
    sumsidesinpart = 0;
    i = 1;
    while i <= partitions {
        sumelementsinpart += *elementsinpart.wrapping_offset(i as isize);
        sumownnodes += *ownnodes.wrapping_offset(i as isize);
        sumsharednodes += *sharednodes.wrapping_offset(i as isize);
        sumsidesinpart += *sidesinpart.wrapping_offset(i as isize);
        i += 1;
    }
    n = partitions;
    printf(
        c"----------------------------------------------------------------------------------------------\n"
            .as_ptr() as *const c_char,
    );
    printf(
        c"   ave   %-10.1f %-10.1f %-8.1f %-8.1f\n".as_ptr() as *const c_char,
        1.0 * sumelementsinpart as Real / n as Real,
        1.0 * sumownnodes as Real / n as Real,
        1.0 * sumsharednodes as Real / n as Real,
        1.0 * sumsidesinpart as Real / n as Real,
    );
    if info != 0 && halobcs != 0 {
        printf(
            c"Number of boundary elements associated with halo: %d\n".as_ptr() as *const c_char,
            halobcs,
        );
    }
    if splitsides != 0 && anyparthalo == 0 {
        printf(
            c"************************* Warning ****************************\n".as_ptr()
                as *const c_char,
        );
        printf(
            c"Number or boundary elements split at between parents: %d\n".as_ptr()
                as *const c_char,
            splitsides,
        );
        printf(
            c"This could be a problem for internal jump conditions\n".as_ptr() as *const c_char,
        );
        printf(
            c"You could try to use '-halobc' flag as remedy with ElmerSolver.\n".as_ptr()
                as *const c_char,
        );
        printf(
            c"**************************************************************\n".as_ptr()
                as *const c_char,
        );
    }
    if anyparthalo != 0 {
        free_ivector(neededtimes2, 1, noknots);
    }
    chdir(c"..".as_ptr() as *const c_char);
    chdir(c"..".as_ptr() as *const c_char);
    if reorder != 0 {
        free_ivector(order, 1, noknots);
    }
    free_ivector(needednodes, 1, partitions);
    free_ivector(neededtwice, 1, partitions);
    free_ivector(sharednodes, 1, partitions);
    free_ivector(ownnodes, 1, partitions);
    free_ivector(sidetypes, minelemtype, maxelemtype);
    free_imatrix(bulktypes, 1, partitions, minelemtype, maxelemtype);
    if info != 0 {
        printf(c"Writing of partitioned mesh finished\n".as_ptr() as *const c_char);
    }
    0
}
// SAFETY: pure element-family mapping table.
unsafe fn LinearNodes(elemtype: c_int) -> c_int {
    let fam2nodemap: [c_int; 9] = [0, 1, 2, 3, 4, 4, 5, 6, 8];
    let elemfamily: c_int = elemtype / 100;
    // NOTE: C indexes fam2nodemap[elemfamily] unguarded; family codes are
    // small non-negative type prefixes on every real mesh (fail-loud here).
    fam2nodemap[elemfamily as usize]
}

/// # Safety
///
/// METIS mesh partitioning over valid single-type meshes.
/// NOTE: minnodes/mintype/maxtype carry 0-inits for the empty-mesh proof
/// (assigned in-loop on every real mesh; unobservable).
#[no_mangle]
pub unsafe extern "C" fn PartitionMetisMesh(
    data: *mut FemType,
    eg: *mut ElmergridType,
    partitions: c_int,
    dual: c_int,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k): (c_int, c_int, c_int);
    let (periodic, noelements, noknots): (c_int, c_int, c_int);
    let mut nodesd2: c_int;
    // NOTE: etype declared but never used in C.
    // NOTE: 0-inits for the empty-mesh proof (assigned in-loop on every
    // real mesh; unobservable). elemtype is assigned before every read.
    let mut elemtype: c_int;
    let (mut mintype, mut maxtype, mut minnodes): (c_int, c_int, c_int) = (0, 0, 0);
    // NOTE: highorder is initialized 0 and never reassigned in C (the
    // PartitionNodesByElements arm is dead but kept for shape fidelity).
    let highorder: c_int = 0;
    let mut ninv: c_int;
    let (neededby, inpart): (*mut c_int, *mut c_int);
    // NOTE: assigned iff periodic; reads are periodic-guarded (flag never
    // changes mid-function) — null init unobservable.
    let mut indxper: *mut c_int = core::ptr::null_mut();
    let (metistopo, eptr, npart, epart): (*mut idx_t, *mut idx_t, *mut idx_t, *mut idx_t);
    let (mut ne, mut nn, mut ncommon, mut nparts): (idx_t, idx_t, idx_t, idx_t);
    // NOTE: METIS-written, never read by C — 0-init for the &mut handoff.
    let mut edgecut: idx_t = 0;
    let mut options: [idx_t; METIS_NOPTIONS] = [0; METIS_NOPTIONS];
    SetMetisOptions(options.as_mut_ptr(), eg, info);
    if info != 0 {
        printf(
            c"Making a Metis partitioning for %d elements in %d-dimensions.\n".as_ptr()
                as *const c_char,
            (*data).noelements,
            (*data).dim,
        );
    }
    noelements = (*data).noelements;
    noknots = (*data).noknots;
    i = 1;
    while i <= noelements {
        elemtype = *(*data).elementtypes.wrapping_offset(i as isize);
        nodesd2 = LinearNodes(elemtype);
        if i == 1 {
            mintype = elemtype;
            maxtype = elemtype;
            minnodes = nodesd2;
        } else {
            mintype = c_min_int(mintype, elemtype);
            maxtype = c_max_int(maxtype, elemtype);
            minnodes = c_min_int(minnodes, nodesd2);
        }
        i += 1;
    }
    if info != 0 {
        if mintype == maxtype {
            printf(
                c"All elements are of type %d\n".as_ptr() as *const c_char,
                mintype,
            );
        } else {
            printf(
                c"Minimum element type is %d\n".as_ptr() as *const c_char,
                mintype,
            );
            printf(
                c"Maximum element type is %d\n".as_ptr() as *const c_char,
                maxtype,
            );
        }
    }
    if minnodes >= 8 {
        ncommon = 4;
    } else if minnodes > 4 {
        ncommon = 3;
    } else {
        ncommon = 2;
    }
    if info != 0 {
        printf(
            c"Minimum number of linear nodes in elements %d\n".as_ptr() as *const c_char,
            minnodes,
        );
        printf(
            c"Requiring number of nodes in dual graph %d\n".as_ptr() as *const c_char,
            ncommon,
        );
    }
    if (*data).partitionexist == 0 {
        (*data).partitionexist = TRUE;
        (*data).elempart = ivector(1, (*data).noelements);
        (*data).nodepart = ivector(1, (*data).noknots);
        (*data).nopartitions = partitions;
    }
    inpart = (*data).elempart;
    periodic = (*data).periodicexist;
    if periodic != 0 {
        if info != 0 {
            printf(c"There seems to be periodic boundaries\n".as_ptr() as *const c_char);
        }
        indxper = (*data).periodic;
    }
    ne = noelements;
    // NOTE: C stores nn = noknots here, but `nn = j` below overwrites
    // before any read; dropped.
    neededby = ivector(1, noknots);
    // NOTE: C stores numflag = 0 here but never reads it (dropped).
    nparts = partitions;
    k = 0;
    i = 1;
    while i <= noknots {
        *neededby.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    if periodic != 0 {
        i = 1;
        while i <= noelements {
            nodesd2 = LinearNodes(*(*data).elementtypes.wrapping_offset(i as isize));
            j = 0;
            while j < nodesd2 {
                *neededby.wrapping_offset(*indxper.wrapping_offset(
                    *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize)
                        as isize,
                ) as isize) = 1;
                k += 1;
                j += 1;
            }
            i += 1;
        }
    } else {
        i = 1;
        while i <= noelements {
            nodesd2 = LinearNodes(*(*data).elementtypes.wrapping_offset(i as isize));
            j = 0;
            while j < nodesd2 {
                *neededby.wrapping_offset(
                    *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize)
                        as isize,
                ) = 1;
                k += 1;
                j += 1;
            }
            i += 1;
        }
    }
    j = 0;
    i = 1;
    while i <= noknots {
        if *neededby.wrapping_offset(i as isize) != 0 {
            j += 1;
            *neededby.wrapping_offset(i as isize) = j;
        }
        i += 1;
    }
    nn = j;
    if info != 0 {
        if nn == noknots {
            printf(
                c"Using all %d possible nodes in the Metis graph\n".as_ptr() as *const c_char,
                nn,
            );
        } else {
            printf(
                c"Using %d nodes of %d possible nodes in the Metis graph\n".as_ptr()
                    as *const c_char,
                nn,
                noknots,
            );
        }
        printf(
            c"Allocating mesh topology of size %d\n".as_ptr() as *const c_char,
            k,
        );
    }
    eptr = ivector(0, noelements);
    metistopo = ivector(0, k - 1);
    k = 0;
    *eptr.wrapping_offset(0) = k;
    if periodic != 0 {
        i = 1;
        while i <= noelements {
            nodesd2 = LinearNodes(*(*data).elementtypes.wrapping_offset(i as isize));
            j = 0;
            while j < nodesd2 {
                *metistopo.wrapping_offset(k as isize) = *neededby.wrapping_offset(
                    *indxper.wrapping_offset(
                        *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize)
                            as isize,
                    ) as isize,
                ) - 1;
                k += 1;
                j += 1;
            }
            *eptr.wrapping_offset(i as isize) = k;
            i += 1;
        }
    } else if nn < noknots {
        i = 1;
        while i <= noelements {
            nodesd2 = LinearNodes(*(*data).elementtypes.wrapping_offset(i as isize));
            j = 0;
            while j < nodesd2 {
                *metistopo.wrapping_offset(k as isize) = *neededby.wrapping_offset(
                    *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize)
                        as isize,
                ) - 1;
                k += 1;
                j += 1;
            }
            *eptr.wrapping_offset(i as isize) = k;
            i += 1;
        }
    } else {
        i = 1;
        while i <= noelements {
            nodesd2 = LinearNodes(*(*data).elementtypes.wrapping_offset(i as isize));
            j = 0;
            while j < nodesd2 {
                *metistopo.wrapping_offset(k as isize) =
                    *(*(*data).topology.wrapping_offset(i as isize)).wrapping_offset(j as isize)
                        - 1;
                k += 1;
                j += 1;
            }
            *eptr.wrapping_offset(i as isize) = k;
            i += 1;
        }
    }
    npart = ivector(0, nn - 1);
    i = 0;
    while i < nn {
        *npart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    epart = ivector(0, ne - 1);
    i = 0;
    while i < ne {
        *epart.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    if dual != 0 {
        if info != 0 {
            printf(
                c"Starting graph partitioning METIS_PartMeshDual.\n".as_ptr() as *const c_char,
            );
        }
        METIS_PartMeshDual(
            &mut ne,
            &mut nn,
            eptr,
            metistopo,
            core::ptr::null_mut(),
            core::ptr::null_mut(),
            &mut ncommon,
            &mut nparts,
            core::ptr::null_mut(),
            options.as_mut_ptr(),
            &mut edgecut,
            epart,
            npart,
        );
        if info != 0 {
            printf(
                c"Finished graph partitioning METIS_PartMeshDual.\n".as_ptr() as *const c_char,
            );
        }
    } else {
        if info != 0 {
            printf(
                c"Starting graph partitioning METIS_PartMeshNodal.\n".as_ptr() as *const c_char,
            );
        }
        METIS_PartMeshNodal(
            &mut ne,
            &mut nn,
            eptr,
            metistopo,
            core::ptr::null_mut(),
            core::ptr::null_mut(),
            &mut nparts,
            core::ptr::null_mut(),
            options.as_mut_ptr(),
            &mut edgecut,
            epart,
            npart,
        );
        if info != 0 {
            printf(
                c"Finished graph partitioning METIS_PartMeshNodal.\n".as_ptr() as *const c_char,
            );
        }
    }
    ninv = 0;
    i = 1;
    while i <= noelements {
        *inpart.wrapping_offset(i as isize) = *epart.wrapping_offset((i - 1) as isize) + 1;
        if *inpart.wrapping_offset(i as isize) < 1
            || *inpart.wrapping_offset(i as isize) > partitions
        {
            if ninv < 10 {
                printf(
                    c"Invalid partition %d for element %d\n".as_ptr() as *const c_char,
                    *inpart.wrapping_offset(i as isize),
                    i,
                );
            }
            ninv += 1;
        }
        i += 1;
    }
    if ninv != 0 {
        printf(
            c"Number of invalid element partitions by Metis %d\n".as_ptr() as *const c_char,
            ninv,
        );
        bigerror(c"Cannot continue with invalid partitioning!".as_ptr() as *const c_char);
    }
    if highorder != 0 {
        PartitionNodesByElements(data, info);
    } else {
        if info != 0 {
            printf(c"Set the partition given by Metis for each node\n".as_ptr() as *const c_char);
        }
        i = 1;
        while i <= noknots {
            if periodic != 0 {
                j = *neededby.wrapping_offset(*indxper.wrapping_offset(i as isize) as isize);
            } else if nn < noknots {
                j = *neededby.wrapping_offset(i as isize);
            } else {
                j = i;
            }
            if j == 0 {
                printf(
                    c"Cannot set partitioning for node %d\n".as_ptr() as *const c_char,
                    i,
                );
            }
            *(*data).nodepart.wrapping_offset(i as isize) =
                *npart.wrapping_offset((j - 1) as isize) + 1;
            if *(*data).nodepart.wrapping_offset(i as isize) < 1
                || *(*data).nodepart.wrapping_offset(i as isize) > partitions
            {
                if ninv < 10 {
                    printf(
                        c"Invalid partition %d for node %d\n".as_ptr() as *const c_char,
                        *(*data).nodepart.wrapping_offset(i as isize),
                        i,
                    );
                }
                ninv += 1;
            }
            i += 1;
        }
        if ninv != 0 {
            printf(
                c"Number of invalid node partitions by Metis %d\n".as_ptr() as *const c_char,
                ninv,
            );
            bigerror(c"Cannot continue with invalid partitioning!".as_ptr() as *const c_char);
        }
    }
    free_ivector(neededby, 1, noknots);
    free_ivector(metistopo, 0, k - 1);
    free_ivector(epart, 0, noelements - 1);
    free_ivector(npart, 0, nn - 1);
    free_ivector(eptr, 0, noelements + 1);
    if info != 0 {
        printf(
            c"Successfully made a Metis partition using the element mesh.\n".as_ptr()
                as *const c_char,
        );
    }
    0
}
/// # Safety
///
/// METIS graph partitioning with skip-path (nparts==1) preserved.
/// `goto skippart` becomes a skip flag (frees skipped on that path, as in C).
#[no_mangle]
pub unsafe extern "C" fn PartitionMetisGraph(
    data: *mut FemType,
    bound: *mut BoundaryType,
    eg: *mut ElmergridType,
    partitions: c_int,
    metisopt: c_int,
    dual: c_int,
    info: c_int,
) -> c_int {
    let (mut i, mut j, mut k): (c_int, c_int, c_int);
    let (noelements, noknots): (c_int, c_int);
    // NOTE: errstat is written but only read via `if(errstat)` right after.
    let errstat: c_int;
    // NOTE: nn/totcon/maxcon/xadj/adjncy are assigned in the dual/nodal arms
    // and read in the shared section iff !skip; the skip correlation is
    // invisible to the compiler, so deterministic inits stand in (dead on
    // every path that reads — the arms overwrite first — and unread on the
    // skip path, exactly as C's indeterminate values).
    #[allow(unused_assignments)]
    let (mut nn, mut ncon, mut con, mut maxcon, mut totcon): (
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    ) = (0, 0, 0, 0, 0);
    #[allow(unused_assignments)]
    let (mut xadj, mut adjncy, mut vwgt, mut adjwgt): (
        *mut c_int,
        *mut c_int,
        *mut c_int,
        *mut c_int,
    ) = (
        core::ptr::null_mut(),
        core::ptr::null_mut(),
        core::ptr::null_mut(),
        core::ptr::null_mut(),
    );
    // NOTE: numflag is assigned 0 and never read in C (dropped).
    // npart is null on the skip path where nparts == 1 skips its only
    // post-skippart read; overwritten first on all other paths.
    let mut npart: *mut c_int = core::ptr::null_mut();
    let (mut nparts, mut maxconset): (c_int, c_int);
    // NOTE: METIS-written, never read by C — 0-init for the &mut handoff.
    let mut edgecut: c_int = 0;
    // NOTE: assigned 0/1 before METIS reads (C shape).
    let mut wgtflag: c_int;
    // NOTE: C reassigns the metisopt parameter (Kway enforcement); mirrored
    // with a mutable shadow so later reads see the enforced value.
    let mut metisopt = metisopt;
    let dualgraph: *mut CRSType;
    let mut options: [idx_t; METIS_NOPTIONS] = [0; METIS_NOPTIONS];
    let mut sideind: [c_int; MAXNODESD1] = [0; MAXNODESD1];
    let mut sideelemtype: c_int = 0;
    // NOTE: con2 declared but never used in C (dropped).
    let (mut sidenodes, mut bc, mut bctype): (c_int, c_int, c_int);
    let (mut j2, mut ind, mut ind2): (c_int, c_int, c_int);
    // NOTE: fn-top maxweight/sidehits/sidepartitions/par are C-dead (live
    // only as block-local shadows in the weights/integrity blocks below).
    if info != 0 {
        printf(
            c"Making a Metis partitioning for %d nodes in %d-dimensions.\n".as_ptr()
                as *const c_char,
            (*data).noknots,
            (*data).dim,
        );
    }
    if partitions < 2 {
        bigerror(
            c"There should be at least two partitions for partitioning!".as_ptr()
                as *const c_char,
        );
    }
    noknots = (*data).noknots;
    noelements = (*data).noelements;
    maxconset = 0;
    if (*data).periodicexist != 0 && dual != 0 {
        printf(c"Dual graph not implemented for periodic system!\n".as_ptr() as *const c_char);
        printf(c"Enforcing nodal graph for partitioning\n".as_ptr() as *const c_char);
    }
    // NOTE: C reassigns the dual parameter; mirrored with a shadow binding
    // so later reads see the enforced value.
    let dual = if (*data).periodicexist != 0 && dual != 0 {
        FALSE
    } else {
        dual
    };
    SetMetisOptions(options.as_mut_ptr(), eg, info);
    nparts = partitions;
    // NOTE: skip flag for `goto skippart` (nparts == 1 in the dual arm).
    let mut skip = false;
    if dual != 0 {
        if (*eg).partbcz > 1 || (*eg).partbcr != 0 {
            PartitionConnectedElements1D(data, bound, eg, info);
        } else if (*eg).partbcmetis > 1 {
            PartitionConnectedElementsMetis(data, bound, eg, (*eg).partbcmetis, metisopt, info);
        } else if (*eg).connect != 0 {
            PartitionConnectedElementsStraight(data, bound, eg, info);
        }
        if (*data).nodeconnectexist != 0 {
            errstat = ExtendBoundaryPartitioning(data, bound, (*eg).partbclayers, info);
            if errstat != 0 {
                bigerror(
                    c"Extend boundary partitioning returned error!".as_ptr() as *const c_char
                );
            }
        }
        if (*data).elemconnectexist != 0 {
            maxconset = 0;
            i = 1;
            while i <= noelements {
                maxconset = c_max_int(maxconset, *(*data).elemconnect.wrapping_offset(i as isize));
                i += 1;
            }
            nparts -= maxconset;
        }
        if nparts == 1 {
            if info != 0 {
                printf(
                    c"Just one partition left, skipping 2nd part of hybrid partitioning!\n"
                        .as_ptr() as *const c_char,
                );
            }
            skip = true;
        }
        if !skip {
            CreateDualGraph(data, TRUE, info);
            dualgraph = &mut (*data).dualgraph;
            nn = (*dualgraph).rowsize;
            xadj = (*dualgraph).rows;
            adjncy = (*dualgraph).cols;
            totcon = (*dualgraph).colsize;
        }
    }
    if !skip {
        if dual == 0 {
            CreateNodalGraph(data, TRUE, info);
            maxcon = (*data).nodalmaxconnections;
            nn = noknots;
            totcon = 0;
            i = 1;
            while i <= nn {
                j = 0;
                while j < maxcon {
                    con = *(*(*data).nodalgraph.as_mut_ptr().wrapping_offset(j as isize))
                        .wrapping_offset(i as isize);
                    if con != 0 {
                        totcon += 1;
                    }
                    j += 1;
                }
                i += 1;
            }
            xadj = ivector(0, nn);
            adjncy = ivector(0, totcon - 1);
            i = 0;
            while i < totcon {
                *adjncy.wrapping_offset(i as isize) = 0;
                i += 1;
            }
            totcon = 0;
            i = 1;
            while i <= nn {
                *xadj.wrapping_offset((i - 1) as isize) = totcon;
                j = 0;
                while j < maxcon {
                    con = *(*(*data).nodalgraph.as_mut_ptr().wrapping_offset(j as isize))
                        .wrapping_offset(i as isize);
                    if con != 0 {
                        *adjncy.wrapping_offset(totcon as isize) = con - 1;
                        totcon += 1;
                    }
                    j += 1;
                }
                i += 1;
            }
            *xadj.wrapping_offset(nn as isize) = totcon;
        }
        if info != 0 {
            printf(
                c"There are %d connections altogether in the graph.\n".as_ptr() as *const c_char,
                totcon,
            );
        }
        // NOTE: C stores numflag = 0 here but never reads it (dropped).
        npart = ivector(0, nn - 1);
        wgtflag = 0;
        ncon = 1;
        vwgt = core::ptr::null_mut();
        adjwgt = core::ptr::null_mut();
        if dual == 0 {
            if (*data).periodicexist != 0 || (*eg).connect != 0 {
                if info != 0 {
                    printf(
                        c"Creating weight for %d connections.\n".as_ptr() as *const c_char,
                        totcon,
                    );
                }
                wgtflag = 1;
                adjwgt = ivector(0, totcon - 1);
                i = 0;
                while i < totcon {
                    *adjwgt.wrapping_offset(i as isize) = 1;
                    i += 1;
                }
                if metisopt != 3 {
                    printf(
                        c"For weighted partitioning Metis subroutine METIS_PartGraphKway is enforced\n"
                            .as_ptr() as *const c_char,
                    );
                    metisopt = 3;
                }
            }
            if (*data).periodicexist != 0 {
                if info != 0 {
                    printf(
                        c"Setting periodic connections to dominate %d\n".as_ptr()
                            as *const c_char,
                        totcon,
                    );
                }
                i = 0;
                while i < noknots {
                    j = *(*data).periodic.wrapping_offset((i + 1) as isize) - 1;
                    if j == i {
                        i += 1;
                        continue;
                    }
                    k = *xadj.wrapping_offset(i as isize);
                    while k < *xadj.wrapping_offset((i + 1) as isize) {
                        if *adjncy.wrapping_offset(k as isize) == j {
                            *adjwgt.wrapping_offset(k as isize) = maxcon;
                        }
                        k += 1;
                    }
                    i += 1;
                }
            }
            if (*eg).connect != 0 {
                // NOTE: block-scoped C locals become a Rust block.
                {
                    
                    let maxweight: c_int = noknots + noelements;
                    printf(
                        c"Adding weight of %d for constrained nodes\n".as_ptr() as *const c_char,
                        maxweight,
                    );
                    con = 1;
                    while con <= (*eg).connect {
                        bctype = *(*eg)
                            .connectbounds
                            .as_mut_ptr()
                            .wrapping_offset((con - 1) as isize);
                        bc = 0;
                        while bc < MAXBOUNDARIES as c_int {
                            if (*bound.wrapping_offset(bc as isize)).created == FALSE {
                                bc += 1;
                                continue;
                            }
                            if (*bound.wrapping_offset(bc as isize)).nosides == 0 {
                                bc += 1;
                                continue;
                            }
                            i = 1;
                            while i <= (*bound.wrapping_offset(bc as isize)).nosides {
                                if *(*bound.wrapping_offset(bc as isize))
                                    .types
                                    .wrapping_offset(i as isize)
                                    != bctype
                                {
                                    i += 1;
                                    continue;
                                }
                                GetBoundaryElement(
                                    i,
                                    bound.wrapping_offset(bc as isize),
                                    data,
                                    sideind.as_mut_ptr(),
                                    &mut sideelemtype,
                                );
                                sidenodes = sideelemtype % 100;
                                j = 0;
                                while j < sidenodes {
                                    j2 = 0;
                                    while j2 < sidenodes {
                                        if j == j2 {
                                            j2 += 1;
                                            continue;
                                        }
                                        ind = sideind[j as usize] - 1;
                                        ind2 = sideind[j2 as usize] - 1;
                                        k = *xadj.wrapping_offset(ind as isize);
                                        while k < *xadj.wrapping_offset((ind + 1) as isize) {
                                            if *adjncy.wrapping_offset(k as isize) == ind2 {
                                                *adjwgt.wrapping_offset(k as isize) = maxweight;
                                            }
                                            k += 1;
                                        }
                                        j2 += 1;
                                    }
                                    j += 1;
                                }
                                i += 1;
                            }
                            bc += 1;
                        }
                        con += 1;
                    }
                }
            }
        }
        if metisopt == 2 {
            if info != 0 {
                printf(
                    c"Starting graph partitioning METIS_PartGraphRecursive.\n".as_ptr()
                        as *const c_char,
                );
            }
            // NOTE: C passes &wgtflag as the vsize slot (quirk preserved).
            METIS_PartGraphRecursive(
                &mut nn,
                &mut ncon,
                xadj,
                adjncy,
                vwgt,
                &mut wgtflag,
                adjwgt,
                &mut nparts,
                core::ptr::null_mut(),
                core::ptr::null_mut(),
                options.as_mut_ptr(),
                &mut edgecut,
                npart,
            );
        } else if metisopt == 3 || metisopt == 4 {
            ncon = 1;
            wgtflag = 0;
            if info != 0 {
                printf(
                    c"Starting graph partitioning METIS_PartGraphKway.\n".as_ptr()
                        as *const c_char,
                );
            }
            // NOTE: C passes &wgtflag as the vsize slot (quirk preserved).
            METIS_PartGraphKway(
                &mut nn,
                &mut ncon,
                xadj,
                adjncy,
                vwgt,
                &mut wgtflag,
                adjwgt,
                &mut nparts,
                core::ptr::null_mut(),
                core::ptr::null_mut(),
                options.as_mut_ptr(),
                &mut edgecut,
                npart,
            );
        } else {
            printf(
                c"Unknown Metis option %d\n".as_ptr() as *const c_char,
                metisopt,
            );
        }
        if info != 0 {
            printf(c"Finished Metis graph partitioning call.\n".as_ptr() as *const c_char);
        }
        free_ivector(adjncy, 0, totcon - 1);
        free_ivector(xadj, 0, nn);
        if wgtflag == 1 {
            free_ivector(adjwgt, 0, totcon - 1);
        }
    }
    // NOTE: `skippart:` label target.
    if (*data).partitionexist == 0 {
        (*data).partitionexist = TRUE;
        (*data).elempart = ivector(1, noelements);
        (*data).nodepart = ivector(1, noknots);
        (*data).nopartitions = partitions;
    }
    if dual != 0 {
        if (*data).elemconnectexist != 0 {
            i = 1;
            while i <= noelements {
                j = *(*data).elemconnect.wrapping_offset(i as isize);
                if nparts == 1 {
                    if j != 0 {
                        *(*data).elempart.wrapping_offset(i as isize) = j;
                    } else {
                        *(*data).elempart.wrapping_offset(i as isize) = maxconset + 1;
                    }
                } else {
                    if j < 0 {
                        *(*data).elempart.wrapping_offset(i as isize) = -j;
                    } else {
                        *(*data).elempart.wrapping_offset(i as isize) =
                            *npart.wrapping_offset((j - 1) as isize) + 1 + maxconset;
                    }
                }
                i += 1;
            }
        } else {
            i = 1;
            while i <= nn {
                *(*data).elempart.wrapping_offset(i as isize) =
                    *npart.wrapping_offset((i - 1) as isize) + 1;
                i += 1;
            }
        }
        PartitionNodesByElements(data, info);
    } else {
        i = 1;
        while i <= nn {
            *(*data).nodepart.wrapping_offset(i as isize) =
                *npart.wrapping_offset((i - 1) as isize) + 1;
            i += 1;
        }
        PartitionElementsByNodes(data, info);
        if (*eg).connect == 0 {
            // NOTE: C block-declares these locally (shadowing fn-top ints);
            // mirrored with block bindings (sidehits is a pointer here).
            {
                let (mut con, mut bc, mut bctype, mut sidenodes): (c_int, c_int, c_int, c_int);
                // NOTE: &mut-handoff init (callee overwrites before any read).
                let mut sideelemtype: c_int = 0;
                let (mut par, mut ind): (c_int, c_int);
                let mut sideind: [c_int; MAXNODESD1] = [0; MAXNODESD1];
                
                let mut sidepartitions: c_int;
                printf(c"Checking connection integrity\n".as_ptr() as *const c_char);
                let sidehits: *mut c_int = ivector(1, partitions);
                con = 1;
                while con <= (*eg).connect {
                    bctype = *(*eg)
                        .connectbounds
                        .as_mut_ptr()
                        .wrapping_offset((con - 1) as isize);
                    i = 1;
                    while i <= partitions {
                        *sidehits.wrapping_offset(i as isize) = 0;
                        i += 1;
                    }
                    bc = 0;
                    while bc < MAXBOUNDARIES as c_int {
                        if (*bound.wrapping_offset(bc as isize)).created == FALSE {
                            bc += 1;
                            continue;
                        }
                        if (*bound.wrapping_offset(bc as isize)).nosides == 0 {
                            bc += 1;
                            continue;
                        }
                        i = 1;
                        while i <= (*bound.wrapping_offset(bc as isize)).nosides {
                            if *(*bound.wrapping_offset(bc as isize))
                                .types
                                .wrapping_offset(i as isize)
                                != bctype
                            {
                                i += 1;
                                continue;
                            }
                            // NOTE: `if(1)` arm taken (else GetElementSide
                            // compiled but dead); transcribed directly.
                            GetBoundaryElement(
                                i,
                                bound.wrapping_offset(bc as isize),
                                data,
                                sideind.as_mut_ptr(),
                                &mut sideelemtype,
                            );
                            sidenodes = sideelemtype % 100;
                            j = 0;
                            while j < sidenodes {
                                ind = sideind[j as usize];
                                par = *(*data).nodepart.wrapping_offset(ind as isize);
                                *sidehits.wrapping_offset(par as isize) += 1;
                                j += 1;
                            }
                            i += 1;
                        }
                        bc += 1;
                    }
                    sidepartitions = 0;
                    i = 1;
                    while i <= partitions {
                        if *sidehits.wrapping_offset(i as isize) != 0 {
                            sidepartitions += 1;
                        }
                        i += 1;
                    }
                    if sidepartitions != 1 {
                        printf(
                            c"PartitionMetisGraph: side %d belongs to %d partitions\n".as_ptr()
                                as *const c_char,
                            bctype,
                            sidepartitions,
                        );
                        bigerror(
                            c"Parallel constraints might not be set!".as_ptr() as *const c_char
                        );
                    }
                    con += 1;
                }
            }
        }
    }
    if nparts > 1 {
        free_ivector(npart, 0, nn - 1);
    }
    if info != 0 {
        if dual != 0 {
            printf(
                c"Successfully made a Metis partition using the dual graph.\n".as_ptr()
                    as *const c_char,
            );
        } else {
            printf(
                c"Successfully made a Metis partition using the nodal graph.\n".as_ptr()
                    as *const c_char,
            );
        }
    }
    0
}
/// # Safety
///
/// METIS fill-reducing reorder over the valid nodal graph.
#[no_mangle]
pub unsafe extern "C" fn ReorderElementsMetis(data: *mut FemType, info: c_int) -> c_int {
    let (mut i, mut j, mut k): (c_int, c_int, c_int);
    let (mut nn, mut totcon, maxcon, mut con): (c_int, c_int, c_int, c_int);
    let (noelements, noknots, mut nonodes): (c_int, c_int, c_int);
    // NOTE: numflag doubles as the vwgt slot (C passes &numflag); edgecut
    // does not exist in this function (unlike the partitioners).
    let mut numflag: c_int;
    let (xadj, adjncy): (*mut c_int, *mut c_int);
    // NOTE: NULL-inits in C; assigned before METIS reads.
    let (iperm, mut perm): (*mut c_int, *mut c_int);
    
    let (newx, newy, newz): (*mut Real, *mut Real, *mut Real);
    let mut options: [c_int; 8] = [0; 8];
    noelements = (*data).noelements;
    noknots = (*data).noknots;
    if info != 0 {
        printf(
            c"Reordering %d knots and %d elements using Metis reordering routine.\n".as_ptr()
                as *const c_char,
            noknots,
            noelements,
        );
    }
    perm = core::ptr::null_mut();
    i = CalculateIndexwidth(data, FALSE, perm);
    if info != 0 {
        printf(
            c"Indexwidth of the original node order is %d.\n".as_ptr() as *const c_char,
            i,
        );
    }
    CreateNodalGraph(data, TRUE, info);
    maxcon = (*data).nodalmaxconnections;
    totcon = 0;
    i = 1;
    while i <= noknots {
        j = 0;
        while j < maxcon {
            con = *(*(*data).nodalgraph.as_mut_ptr().wrapping_offset(j as isize))
                .wrapping_offset(i as isize);
            if con != 0 {
                totcon += 1;
            }
            j += 1;
        }
        i += 1;
    }
    if info != 0 {
        printf(
            c"There are %d connections altogether\n".as_ptr() as *const c_char,
            totcon,
        );
    }
    xadj = ivector(0, noknots);
    adjncy = ivector(0, totcon - 1);
    i = 0;
    while i < totcon {
        *adjncy.wrapping_offset(i as isize) = 0;
        i += 1;
    }
    totcon = 0;
    i = 1;
    while i <= noknots {
        *xadj.wrapping_offset((i - 1) as isize) = totcon;
        j = 0;
        while j < maxcon {
            con = *(*(*data).nodalgraph.as_mut_ptr().wrapping_offset(j as isize))
                .wrapping_offset(i as isize);
            if con != 0 {
                *adjncy.wrapping_offset(totcon as isize) = con - 1;
                totcon += 1;
            }
            j += 1;
        }
        i += 1;
    }
    *xadj.wrapping_offset(noknots as isize) = totcon;
    nn = noknots;
    numflag = 0;
    i = 0;
    while i < 8 {
        options[i as usize] = 0;
        i += 1;
    }
    perm = ivector(0, noknots - 1);
    iperm = ivector(0, noknots - 1);
    if info != 0 {
        printf(c"Starting Metis reordering routine.\n".as_ptr() as *const c_char);
    }
    METIS_NodeND(
        &mut nn,
        xadj,
        adjncy,
        &mut numflag,
        options.as_mut_ptr(),
        perm,
        iperm,
    );
    if info != 0 {
        printf(c"Finished Metis reordering routine.\n".as_ptr() as *const c_char);
    }
    if info != 0 {
        printf(c"Moving knots to new positions\n".as_ptr() as *const c_char);
    }
    newx = dvector(1, (*data).noknots);
    newy = dvector(1, (*data).noknots);
    newz = dvector(1, (*data).noknots);
    i = 1;
    while i <= (*data).noknots {
        // NOTE: C x[perm[i-1]+1] (0-based METIS output shifted to 1-based).
        *newx.wrapping_offset(i as isize) = *(*data)
            .x
            .wrapping_offset((*perm.wrapping_offset((i - 1) as isize) + 1) as isize);
        *newy.wrapping_offset(i as isize) = *(*data)
            .y
            .wrapping_offset((*perm.wrapping_offset((i - 1) as isize) + 1) as isize);
        *newz.wrapping_offset(i as isize) = *(*data)
            .z
            .wrapping_offset((*perm.wrapping_offset((i - 1) as isize) + 1) as isize);
        i += 1;
    }
    free_dvector((*data).x, 1, (*data).noknots);
    free_dvector((*data).y, 1, (*data).noknots);
    free_dvector((*data).z, 1, (*data).noknots);
    (*data).x = newx;
    (*data).y = newy;
    (*data).z = newz;
    if info != 0 {
        printf(c"Changing the element topology\n".as_ptr() as *const c_char);
    }
    let newtopology: *mut *mut c_int = imatrix(1, noelements, 0, (*data).maxnodes - 1);
    j = 1;
    while j <= noelements {
        nonodes = *(*data).elementtypes.wrapping_offset(j as isize) % 100;
        i = 0;
        while i < nonodes {
            k = *(*(*data).topology.wrapping_offset(j as isize)).wrapping_offset(i as isize);
            *(*newtopology.wrapping_offset(j as isize)).wrapping_offset(i as isize) =
                *iperm.wrapping_offset((k - 1) as isize) + 1;
            i += 1;
        }
        j += 1;
    }
    free_imatrix((*data).topology, 1, noelements, 0, (*data).maxnodes - 1);
    (*data).topology = newtopology;
    i = CalculateIndexwidth(data, FALSE, perm);
    if info != 0 {
        printf(
            c"Indexwidth of the new node order is %d.\n".as_ptr() as *const c_char,
            i,
        );
    }
    // NOTE: `if(0)` diagnostic omitted (runtime-dead).
    free_ivector(iperm, 0, noknots - 1);
    free_ivector(perm, 0, noknots - 1);
    0
}
