//! A C-style pointer for the translation of UMFPACK's pointer arithmetic.
//!
//! UMFPACK carves several arrays out of one workspace and lets them overlap
//! (`Ci + nz`, `Rp = Fr_npivcol`, the Numeric object's `Memory` of `Unit`s
//! holding integers, doubles, element headers and tuples). `P<T>` keeps that
//! structure: a copyable pointer with an offset, indexed with C's signed
//! integers. Debug builds check every access against the allocation it was
//! derived from, so a translation slip panics instead of corrupting memory.

use core::marker::PhantomData;
use core::ops::{Add, Index, IndexMut, Sub};

pub struct P<T> {
    base: *mut T,
    /// Number of `T` in the allocation (checked in debug builds).
    len: isize,
    off: isize,
    _t: PhantomData<T>,
}

impl<T> Clone for P<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for P<T> {}

impl<T> P<T> {
    pub fn null() -> Self {
        P {
            base: core::ptr::null_mut(),
            len: 0,
            off: 0,
            _t: PhantomData,
        }
    }

    /// Pointer to the first element of `v`. The vector must outlive every
    /// use of the pointer and must not be reallocated meanwhile.
    pub fn new(v: &mut [T]) -> Self {
        P {
            base: v.as_mut_ptr(),
            len: v.len() as isize,
            off: 0,
            _t: PhantomData,
        }
    }

    pub fn is_null(&self) -> bool {
        self.base.is_null()
    }

    #[inline(always)]
    fn at(&self, i: isize) -> *mut T {
        let k = self.off + i;
        debug_assert!(
            !self.base.is_null() && k >= 0 && k < self.len,
            "P: index {i} (offset {}) outside allocation of {}",
            self.off,
            self.len
        );
        // SAFETY: in-bounds by construction (checked in debug builds); the
        // allocation outlives the pointer (caller contract of `new`).
        unsafe { self.base.offset(k) }
    }

    /// Read element `i` (C `p[i]`).
    #[inline(always)]
    pub fn get(&self, i: i32) -> T
    where
        T: Copy,
    {
        unsafe { *self.at(i as isize) }
    }

    /// Write element `i`.
    #[inline(always)]
    pub fn set(&self, i: i32, v: T) {
        unsafe { *self.at(i as isize) = v }
    }

    /// Mutable reference to element `i` (for `+=` and friends). The
    /// reference must not outlive the statement it is used in.
    #[allow(clippy::mut_from_ref)]
    #[inline(always)]
    pub fn r(&self, i: i32) -> &mut T {
        unsafe { &mut *self.at(i as isize) }
    }

    /// `p - q` for two pointers into the same allocation.
    pub fn diff(&self, q: P<T>) -> i32 {
        debug_assert!(self.base == q.base);
        (self.off - q.off) as i32
    }

    /// Reinterpret as a pointer to `U` (C cast). The offset must be a whole
    /// number of `U`s.
    pub fn cast<U>(&self) -> P<U> {
        let st = core::mem::size_of::<T>() as isize;
        let su = core::mem::size_of::<U>() as isize;
        debug_assert!((self.off * st) % su == 0, "P::cast: misaligned offset");
        P {
            base: self.base as *mut U,
            len: self.len * st / su,
            off: self.off * st / su,
            _t: PhantomData,
        }
    }

    /// Raw address, for pointer comparisons.
    pub fn addr(&self) -> usize {
        (self.base as usize)
            .wrapping_add((self.off as usize).wrapping_mul(core::mem::size_of::<T>()))
    }
}

impl<T> Add<i32> for P<T> {
    type Output = P<T>;
    #[inline(always)]
    fn add(self, k: i32) -> P<T> {
        P {
            off: self.off + k as isize,
            ..self
        }
    }
}

impl<T> Sub<i32> for P<T> {
    type Output = P<T>;
    #[inline(always)]
    fn sub(self, k: i32) -> P<T> {
        P {
            off: self.off - k as isize,
            ..self
        }
    }
}

impl<T> Index<i32> for P<T> {
    type Output = T;
    #[inline(always)]
    fn index(&self, i: i32) -> &T {
        unsafe { &*self.at(i as isize) }
    }
}

impl<T> IndexMut<i32> for P<T> {
    #[inline(always)]
    fn index_mut(&mut self, i: i32) -> &mut T {
        unsafe { &mut *self.at(i as isize) }
    }
}

impl<T> PartialEq for P<T> {
    fn eq(&self, o: &Self) -> bool {
        self.addr() == o.addr()
    }
}

impl<T> PartialOrd for P<T> {
    fn partial_cmp(&self, o: &Self) -> Option<core::cmp::Ordering> {
        self.addr().partial_cmp(&o.addr())
    }
}
