#![allow(dead_code)]

use crate::params::N;

#[repr(align(32))]
#[derive(Debug)]
pub struct Poly{
    pub coeffs: [i32; N],
}

impl Poly {
    pub fn new() -> Self {
        Self {
            coeffs: [0i32; N],
        }
    }

}
