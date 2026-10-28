//! Port of `fem/src/mpif_stubs.F90` (Elmer @ 9f6af2f85) for lane-femcore-b,
//! through the gfortran ABI shim (ADR-027, T017 shape A).
//!
//! Unit: `fortran:fem/src/mpif_stubs.F90` — file-scope external
//! subprograms plus one module variable. The emitted gfortran symbols are
//! the 39 `mpi_*_`/`pd*upd_`/`pz*upd_` entry points and the module
//! variable `__mpi_stub_MOD_mpi_stub_init`; every dummy travels by
//! reference (a pointer to the caller's argument).
//!
//! The unit is Elmer's single-rank MPI fallback: the seven subroutines
//! with bodies implement the degenerate communicator (size 1, rank 0,
//! split is identity, the requested thread level passes through and the
//! init flag tracks `mpi_init`/`mpi_finalize`); the other 32 are empty
//! bodies that only return. What the port owns is that exact behavior —
//! no numerics, no allocation, no branching beyond it.
//!
//! Faithfulness notes:
//! - `mpi_initialized` always reports `.FALSE.` (the stub never consults
//!   the flag); the port does the same.
//! - gfortran `LOGICAL` is 4 bytes (`.FALSE.` = 0, `.TRUE.` = 1); the
//!   module variable lives in BSS as a zero-initialized `c_int`, exactly
//!   like the Fortran module storage.
//! - The no-argument stubs take no parameters at all, matching the
//!   Fortran call sequence (nothing is passed, nothing is read).

use std::os::raw::c_int;

/// `mpi_stub::mpi_stub_init`: `.FALSE.` (0) until `mpi_init`, `.FALSE.`
/// again after `mpi_finalize`. BSS storage, like the Fortran module
/// variable.
#[export_name = "__mpi_stub_MOD_mpi_stub_init"]
pub static mut MPI_STUB_INIT: c_int = 0;

fn set_flag(value: c_int) {
    unsafe {
        *std::ptr::addr_of_mut!(MPI_STUB_INIT) = value;
    }
}

/// `MPI_INIT`: serial fallback — success, flag set.
/// # Safety
///
/// `ierr` must point to a writable `c_int`.
#[export_name = "mpi_init_"]
pub unsafe extern "C" fn port_mpi_init(ierr: *mut c_int) {
    unsafe {
        *ierr = 0;
    }
    set_flag(1);
}

/// `MPI_INITIALIZED`: the stub always reports `.FALSE.`, flag or not.
/// # Safety
///
/// `init` and `ierr` must point to writable `c_int`s.
#[export_name = "mpi_initialized_"]
pub unsafe extern "C" fn port_mpi_initialized(init: *mut c_int, ierr: *mut c_int) {
    unsafe {
        *ierr = 0;
        *init = 0;
    }
}

/// `MPI_INIT_THREAD`: the requested level passes straight through.
/// # Safety
///
/// `req` must point to a valid `c_int`; `prov` and `ierr` must point to
/// writable `c_int`s.
#[export_name = "mpi_init_thread_"]
pub unsafe extern "C" fn port_mpi_init_thread(
    req: *const c_int,
    prov: *mut c_int,
    ierr: *mut c_int,
) {
    unsafe {
        *prov = *req;
        *ierr = 0;
    }
}

/// `MPI_FINALIZE`: serial fallback — success, flag cleared.
/// # Safety
///
/// `ierr` must point to a writable `c_int`.
#[export_name = "mpi_finalize_"]
pub unsafe extern "C" fn port_mpi_finalize(ierr: *mut c_int) {
    unsafe {
        *ierr = 0;
    }
    set_flag(0);
}

/// `MPI_COMM_SIZE`: the degenerate communicator has size 1.
/// # Safety
///
/// `csize` and `ierr` must point to writable `c_int`s (`comm` is read but
/// ignored, matching the stub).
#[export_name = "mpi_comm_size_"]
pub unsafe extern "C" fn port_mpi_comm_size(
    _comm: *const c_int,
    csize: *mut c_int,
    ierr: *mut c_int,
) {
    unsafe {
        *ierr = 0;
        *csize = 1;
    }
}

/// `MPI_COMM_RANK`: the degenerate communicator has rank 0.
/// # Safety
///
/// `rank` and `ierr` must point to writable `c_int`s (`comm` is read but
/// ignored, matching the stub).
#[export_name = "mpi_comm_rank_"]
pub unsafe extern "C" fn port_mpi_comm_rank(
    _comm: *const c_int,
    rank: *mut c_int,
    ierr: *mut c_int,
) {
    unsafe {
        *rank = 0;
        *ierr = 0;
    }
}

/// `MPI_COMM_SPLIT`: the split of a singleton is itself.
/// # Safety
///
/// `comm` must point to a valid `c_int`; `newcomm` and `ierr` must point
/// to writable `c_int`s.
#[export_name = "mpi_comm_split_"]
pub unsafe extern "C" fn port_mpi_comm_split(
    comm: *const c_int,
    _color: *const c_int,
    _key: *const c_int,
    newcomm: *mut c_int,
    ierr: *mut c_int,
) {
    unsafe {
        *newcomm = *comm;
        *ierr = 0;
    }
}

macro_rules! empty_stub {
    ($export:literal, $fn:ident) => {
        #[export_name = $export]
        pub extern "C" fn $fn() {}
    };
}

empty_stub!("mpi_scan_", port_mpi_scan);
empty_stub!("mpi_allreduce_", port_mpi_allreduce);
empty_stub!("mpi_buffer_detach_", port_mpi_buffer_detach);
empty_stub!("mpi_recv_", port_mpi_recv);
empty_stub!("mpi_buffer_attach_", port_mpi_buffer_attach);
empty_stub!("mpi_send_", port_mpi_send);
empty_stub!("mpi_barrier_", port_mpi_barrier);
empty_stub!("mpi_wait_", port_mpi_wait);
empty_stub!("mpi_waitany_", port_mpi_waitany);
empty_stub!("mpi_bsend_", port_mpi_bsend);
empty_stub!("mpi_comm_free_", port_mpi_comm_free);
empty_stub!("mpi_waitall_", port_mpi_waitall);
empty_stub!("mpi_comm_group_", port_mpi_comm_group);
empty_stub!("mpi_group_incl_", port_mpi_group_incl);
empty_stub!("mpi_comm_create_", port_mpi_comm_create);
empty_stub!("mpi_irecv_", port_mpi_irecv);
empty_stub!("mpi_isend_", port_mpi_isend);
empty_stub!("mpi_bcast_", port_mpi_bcast);
empty_stub!("mpi_allgather_", port_mpi_allgather);
empty_stub!("mpi_allgatherv_", port_mpi_allgatherv);
empty_stub!("mpi_alltoallv_", port_mpi_alltoallv);
empty_stub!("mpi_alltoall_", port_mpi_alltoall);
empty_stub!("mpi_gatherv_", port_mpi_gatherv);
empty_stub!("mpi_gather_", port_mpi_gather);
empty_stub!("mpi_reduce_", port_mpi_reduce);
empty_stub!("mpi_get_count_", port_mpi_get_count);
empty_stub!("pdseupd_", port_pdseupd);
empty_stub!("pdneupd_", port_pdneupd);
empty_stub!("pzneupd_", port_pzneupd);
empty_stub!("pdsaupd_", port_pdsaupd);
empty_stub!("pdnaupd_", port_pdnaupd);
empty_stub!("pznaupd_", port_pznaupd);
