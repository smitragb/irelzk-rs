#![allow(dead_code)]
use crate::poly_arith::poly::Poly;
use crate::params::{K, L, M};

pub struct PolyVec<const S: usize> {
    pub vec: [Poly; S]
}

impl<const S: usize> PolyVec<S> {
    pub fn new() -> Self {
        Self {
            vec: std::array::from_fn(|_| Poly::new())
        }
    }

    pub fn freeze(&mut self) {
        for i in 0..S {
            self.vec[i].freeze();
        }
    }

    pub fn add(&mut self, other: &PolyVec<S>) {
        for i in 0..S {
            self.vec[i].add(&other.vec[i]);
        }
    }

    pub fn add_other(w: &mut PolyVec<S>, u: &PolyVec<S>, v: &PolyVec<S>) {
        for i in 0..S {
            Poly::add_other(&mut w.vec[i], &u.vec[i], &v.vec[i]);
        }
    }

    pub fn sub(&mut self, other: &PolyVec<S>) {
        for i in 0..S {
            self.vec[i].sub(&other.vec[i]);
        }
    }

    pub fn sub_other(w: &mut PolyVec<S>, u: &PolyVec<S>, v: &PolyVec<S>) {
        for i in 0..S {
            Poly::sub_other(&mut w.vec[i], &u.vec[i], &v.vec[i]);
        }
    }

    pub fn vec_ntt(&mut self) {
        for i in 0..S {
            self.vec[i].ntt();
        }
    }

    pub fn vec_inverse_ntt(&mut self) {
        for i in 0..S {
            self.vec[i].inverse_ntt();
        }
    }

    pub fn vec_inverse_ntt_tomont(&mut self) {
        for i in 0..S {
            self.vec[i].inverse_ntt_tomont();
        }
    }

    pub fn pointwise_acc_montgomery(u: &PolyVec<S>, v: &PolyVec<S>) -> Poly {
        let mut r = Poly::new();
        let mut t = Poly::new();
        Poly::pointwise_montgomery_other(&mut r, &u.vec[0], &v.vec[0]);
        for i in 1..S {
            Poly::pointwise_montgomery_other(&mut t, &u.vec[i], &v.vec[i]);
            r.add(&t);
        }
        r
    }

    pub fn scale_montgomery(&mut self, s: i32) {
        for i in 0..S {
            self.vec[i].scale_montgomery(s);
        }
    }

    pub fn scale_montgomery_other(v: &mut PolyVec<S>, u: &PolyVec<S>, s: i32) {
        for i in 0..S {
            Poly::scale_montgomery_other(&mut v.vec[i], &u.vec[i], s);
        }
    }

    pub fn vec_check_norm(v: &PolyVec<S>, b: u32) -> bool {
        for i in 0..S {
            if Poly::check_norm(&v.vec[i], b) {
                return true;
            }
        }
        return false;
    }

    pub fn vec_power2round (v1: &mut PolyVec<S>, v0: &mut PolyVec<S>, v: &mut PolyVec<S>) {
        for i in 0..S {
            Poly::power2round(&mut v1.vec[i], &mut v0.vec[i], &mut v.vec[i]);
        }
    }

    pub fn vec_decompose (v1: &mut PolyVec<S>, v0: &mut PolyVec<S>, v: &mut PolyVec<S>) {
        for i in 0..S {
            Poly::decompose(&mut v1.vec[i], &mut v0.vec[i], &mut v.vec[i]);
        }
    }

    pub fn vec_makehint (h: &mut PolyVec<S>, v1: &PolyVec<S>, v0: &mut PolyVec<S>) {
        for i in 0..S {
            Poly::makehint(&mut h.vec[i], &v1.vec[i], &mut v0.vec[i]);
        }
    }

    pub fn vec_usehint (v1: &mut PolyVec<S>, v: &mut PolyVec<S>, h: &PolyVec<S>) {
        for i in 0..S {
            Poly::usehint(&mut v1.vec[i], &mut v.vec[i], &h.vec[i]);
        }
    }
}

pub type PolyVecK = PolyVec<K>;
pub type PolyVecL = PolyVec<L>;
pub type PolyVecM = PolyVec<M>;
