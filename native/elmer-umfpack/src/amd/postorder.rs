//! `AMD_postorder` and `AMD_post_tree`: depth-first postorder of the
//! assembly tree, largest child last.
//!
//! AMD Version 1.1, Copyright (c) 2004 by Timothy A. Davis, Patrick R.
//! Amestoy, and Iain S. Duff.  Translated to Rust.

use super::*;

#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn amd_postorder(
    nn: Int,
    Parent: &[Int],
    Nv: &[Int],
    Fsize: &[Int],
    Order: &mut [Int],
    Child: &mut [Int],
    Sibling: &mut [Int],
    Stack: &mut [Int],
) {
    let n = nn as usize;
    for j in 0..n {
        Child[j] = EMPTY;
        Sibling[j] = EMPTY;
    }
    for j in (0..n).rev() {
        if Nv[j] > 0 {
            let parent = Parent[j];
            if parent != EMPTY {
                Sibling[j] = Child[parent as usize];
                Child[parent as usize] = j as Int;
            }
        }
    }
    for i in 0..n {
        if Nv[i] > 0 && Child[i] != EMPTY {
            let mut fprev = EMPTY;
            let mut maxfrsize = EMPTY;
            let mut bigfprev = EMPTY;
            let mut bigf = EMPTY;
            let mut f = Child[i];
            while f != EMPTY {
                let frsize = Fsize[f as usize];
                if frsize >= maxfrsize {
                    maxfrsize = frsize;
                    bigfprev = fprev;
                    bigf = f;
                }
                fprev = f;
                f = Sibling[f as usize];
            }
            let fnext = Sibling[bigf as usize];
            if fnext != EMPTY {
                if bigfprev == EMPTY {
                    Child[i] = fnext;
                } else {
                    Sibling[bigfprev as usize] = fnext;
                }
                Sibling[bigf as usize] = EMPTY;
                Sibling[fprev as usize] = bigf;
            }
        }
    }
    for o in Order.iter_mut().take(n) {
        *o = EMPTY;
    }
    let mut k: Int = 0;
    for i in 0..n {
        if Parent[i] == EMPTY && Nv[i] > 0 {
            k = amd_post_tree(i as Int, k, Child, Sibling, Order, Stack);
        }
    }
}

#[allow(non_snake_case)]
pub fn amd_post_tree(
    root: Int,
    mut k: Int,
    Child: &mut [Int],
    Sibling: &[Int],
    Order: &mut [Int],
    Stack: &mut [Int],
) -> Int {
    let mut head: Int = 0;
    Stack[0] = root;
    while head >= 0 {
        let i = Stack[head as usize];
        let iu = i as usize;
        if Child[iu] != EMPTY {
            let mut f = Child[iu];
            while f != EMPTY {
                head += 1;
                f = Sibling[f as usize];
            }
            let mut h = head;
            let mut f = Child[iu];
            while f != EMPTY {
                Stack[h as usize] = f;
                h -= 1;
                f = Sibling[f as usize];
            }
            Child[iu] = EMPTY;
        } else {
            head -= 1;
            Order[iu] = k;
            k += 1;
        }
    }
    k
}
