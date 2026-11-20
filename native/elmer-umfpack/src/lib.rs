//! Bit-exact Rust port of the UMFPACK 4.4 sparse LU and the AMD 1.1 ordering
//! that pristine Elmer bundles (`umfpack/src/{umfpack,amd}`), real `double`
//! with `int` indices (the `di` / `DINT` build Elmer's `umf4*` wrappers use),
//! compiled as Elmer compiles it (`NBLAS`, `NDEBUG`, `-O3` without FMA).
//!
//! The port follows the C sources statement by statement so that orderings,
//! pivots and every floating-point operation happen in the same order as the
//! C library; the result is bitwise identical, not merely close.
//!
//! Licence of the original code (retained as its licence requires):
//!
//! ```text
//! UMFPACK Version 4.4, Copyright 1995-2005 by Timothy A. Davis.
//! All Rights Reserved.
//! AMD Version 1.1 (Jan. 21, 2004), Copyright (c) 2004 by Timothy A. Davis,
//! Patrick R. Amestoy, and Iain S. Duff.
//!
//! UMFPACK License:
//!
//!     Your use or distribution of UMFPACK or any modified version of
//!     UMFPACK implies that you agree to this License.
//!
//!     This library is free software; you can redistribute it and/or
//!     modify it under the terms of the GNU Lesser General Public
//!     License as published by the Free Software Foundation; either
//!     version 2.1 of the License, or (at your option) any later version.
//!
//!     This library is distributed in the hope that it will be useful,
//!     but WITHOUT ANY WARRANTY; without even the implied warranty of
//!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//!     Lesser General Public License for more details.
//!
//!     You should have received a copy of the GNU Lesser General Public
//!     License along with this library; if not, write to the Free Software
//!     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
//!     USA
//!
//!     Permission is hereby granted to use or copy this program under the
//!     terms of the GNU LGPL, provided that the Copyright, this License,
//!     and the Availability of the original version is retained on all copies.
//!     User documentation of any code that uses this code or any modified
//!     version of this code must cite the Copyright, this License, the
//!     Availability note, and "Used by permission." Permission to modify
//!     the code and to distribute modified code is granted, provided the
//!     Copyright, this License, and the Availability note are retained,
//!     and a notice that the code was modified is included.
//!
//! Availability:
//!
//!     http://www.cise.ufl.edu/research/sparse/umfpack
//!     http://www.cise.ufl.edu/research/sparse/amd
//! ```
//!
//! Modification notice: this is a translation of the C sources to Rust for
//! the Elmer Rust rewrite (2026). Debugging, reporting, MATLAB and complex
//! code paths are omitted; the computational paths are translated unchanged.
//! Used by permission.

// The port keeps the C loop bounds, argument lists and array names so that it
// can be checked line by line against the original.
#![allow(
    clippy::int_plus_one,
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::needless_late_init,
    clippy::if_same_then_else,
    clippy::explicit_counter_loop,
    clippy::neg_cmp_op_on_partial_ord,
    non_snake_case,
    non_upper_case_globals
)]

pub mod amd;
pub mod ptr;
pub mod umf;
