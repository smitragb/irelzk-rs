#![allow(dead_code)]
#![allow(unused_imports)]
use array_init::array_init;

use crate::{
    crypto::aes256::Aes256Ctx, 
    params::{K, M, SYMBYTES}, 
    poly_arith::{
        poly::Poly, 
        polyvec::{PolyVecK, PolyVecL, PolyVecM}
    }
};

pub struct Comm {
   t0: PolyVecK,
   tm: PolyVecM,
}

pub struct CommRnd {
    s: PolyVecL,
    e: PolyVecK,
    m: PolyVecM,
}

pub struct CommKey {
    b0: [PolyVecL; K],
    bt: [PolyVecM; K],
    bm: [PolyVecL; M],
}

