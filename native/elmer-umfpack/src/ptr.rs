//! A C-style pointer for the translation of UMFPACK's pointer arithmetic.
//!
//! UMFPACK carves several arrays out of one workspace and lets them overlap
//! (`Ci + nz`, `Rp = Fr_npivcol`, the Numeric object's `Memory` of `Unit`s
//! holding integers, doubles, element headers and tuples). `P<T>` keeps that
//! structure: a copyable pointer (allocation base plus byte offset), indexed
//! with C's signed integers and freely cast between element types. Debug
//! builds check every access against the allocation it was derived from, so
//! a translation slip panics instead of corrupting memory.

use core::marker::PhantomData;
use core::ops::{Add, Index, IndexMut, Sub};

pub struct P<T> {
    base: *mut u8,
    /// Size of the allocation in bytes (checked in debug builds).
    len: isize,
    /// Byte offset of the pointer from `base`.
    off: isize,
    _t: PhantomData<*mut T>,
}

impl<T> Clone for P<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for P<T> {}

impl<T> core::fmt::Debug for P<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "P({:p}+{})", self.base, self.off)
    }
}

const fn sz<T>() -> isize {
    core::mem::size_of::<T>() as isize
}

impl<T> P<T> {
    pub fn null() -> Self {
        P {
            base: core::ptr::null_mut(),
            len: 0,
            off: 0,
            _t: PhantomData,
        }
    }

    /// Pointer to the first element of `v`. The storage must outlive every
    /// use of the pointer and must not be reallocated meanwhile.
    pub fn new(v: &mut [T]) -> Self {
        P {
            base: v.as_mut_ptr() as *mut u8,
            len: core::mem::size_of_val(v) as isize,
            off: 0,
            _t: PhantomData,
        }
    }

    pub fn is_null(&self) -> bool {
        self.base.is_null()
    }

    #[inline(always)]
    fn at(&self, i: isize) -> *mut T {
        let k = self.off + i * sz::<T>();
        debug_assert!(
            !self.base.is_null() && k >= 0 && k + sz::<T>() <= self.len,
            "P: index {i} (byte offset {}) outside allocation of {} bytes",
            self.off,
            self.len
        );
        // SAFETY: in bounds by construction (checked in debug builds); the
        // allocation outlives the pointer (caller contract of `new`).
        unsafe { self.base.offset(k) as *mut T }
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

    /// Mutable reference to element `i` (C `p[i]` / `p->field`). Like a C
    /// pointer dereference it is not tied to a borrow: the caller keeps it
    /// only while the allocation lives and no other reference to the same
    /// element is used.
    #[allow(clippy::mut_from_ref)]
    #[inline(always)]
    pub fn r<'a>(&self, i: i32) -> &'a mut T {
        unsafe { &mut *self.at(i as isize) }
    }

    /// `p - q` in elements, for two pointers into the same allocation.
    pub fn diff(&self, q: P<T>) -> i32 {
        debug_assert!(self.base == q.base);
        ((self.off - q.off) / sz::<T>()) as i32
    }

    /// Reinterpret as a pointer to `U` (C cast).
    #[inline(always)]
    pub fn cast<U>(&self) -> P<U> {
        P {
            base: self.base,
            len: self.len,
            off: self.off,
            _t: PhantomData,
        }
    }

    /// Raw address, for pointer comparisons.
    pub fn addr(&self) -> usize {
        (self.base as usize).wrapping_add(self.off as usize)
    }
}

impl<T> Add<i32> for P<T> {
    type Output = P<T>;
    #[inline(always)]
    fn add(self, k: i32) -> P<T> {
        P {
            off: self.off + k as isize * sz::<T>(),
            ..self
        }
    }
}

impl<T> Sub<i32> for P<T> {
    type Output = P<T>;
    #[inline(always)]
    fn sub(self, k: i32) -> P<T> {
        P {
            off: self.off - k as isize * sz::<T>(),
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
