#![allow(dead_code)]
#![allow(unused_assignments)]

use crate::{
    comm::commitment::{Comm, CommKey, CommRnd}, 
    crypto::aes256::Aes256Ctx, 
    params::{M, N, R, SYMBYTES}, 
    poly_arith::{
        consts::MONTSQ, 
        poly::Poly, 
        polyvec::{PolyVecL, PolyVecM}
    }
};

use super::product::poly_shift;

pub fn proof (
    vprime: &mut [Poly; R],
       msg: &PolyVecM,
     chash: &[u8; SYMBYTES],
         g: &[PolyVecM; R]
) -> Poly {
    let mut nonce = 0;
    let mut state = Aes256Ctx::init(&chash, nonce);
    let mut h = Poly::new();
    let vtmp: [Poly; R] = std::array::from_fn(|j| {
        let mut tmp = Poly::new();
        Poly::add_other(&mut tmp, &g[j].vec[0], &g[j].vec[1]);
        tmp.sub(&g[j].vec[2]);
        tmp
    });
    let mut mpr = Poly::new();
    Poly::add_other(&mut mpr, &msg.vec[0], &msg.vec[1]);
    mpr.sub(&msg.vec[2]);
    h = msg.vec[M-1].clone(); 
    for j in 0..R {
        vprime[j] = g[j].vec[M-1].clone();
    }
    for i in 0..R {
        state.select(nonce);
        nonce += 1;
        let mut gamma   = Poly::new();
        let mut atgamma = Poly::new();
        let mut tmp     = Poly::new();
        let mut tmp1    = Poly::new();
        
        Poly::uniform_preinit(&mut gamma, &mut state);
        for j in 0..(N-1) {
            atgamma.coeffs[j] = 2*gamma.coeffs[j] - gamma.coeffs[j+1];
        }
        atgamma.coeffs[N-1] = 2*gamma.coeffs[N-1];

        for j in 0..4 {
            Poly::pointwise_montgomery_other(&mut tmp, &gamma, &vtmp[j]);
            Poly::pointwise_montgomery_other(&mut tmp1, &atgamma, &g[j].vec[3]);
            tmp.sub(&tmp1);
            tmp.trace65_ntt();
            let tmp2 = tmp.clone();
            poly_shift(&mut tmp, &tmp2, i);
            vprime[j].add(&tmp);
        }
        
        Poly::pointwise_montgomery_other(&mut tmp, &gamma, &mpr);
        Poly::pointwise_montgomery_other(&mut tmp1, &atgamma, &msg.vec[3]);
        tmp.sub(&tmp1);
        tmp.trace65_ntt();
        let tmp2 = tmp.clone();
        poly_shift(&mut tmp, &tmp2, i);
        h.add(&tmp);
    }
    
    h.inverse_ntt();
    h.reduce();
    h.freeze();
    for i in 0..R {
        vprime[i].freeze();
    }
    h
}

pub fn verify (
    vprime: &mut [Poly; R],
     chash: &[u8; SYMBYTES],
         h: &Poly,
         c: &[Poly; R],
         z: &[CommRnd; R],
        tp: &Comm,
       ckp: &CommKey
) -> bool {
    if h.coeffs[0] != 0 || h.coeffs[1] != 0 || h.coeffs[2] != 0 || h.coeffs[3] != 0 {
        return true;
    }
    
    let mut nonce = 0;
    let mut state = Aes256Ctx::init(&chash, nonce);
    let gamma: [Poly; R] = std::array::from_fn(|_| {
        state.select(nonce);
        nonce += 1;
        let mut a = Poly::new();
        Poly::uniform_preinit(&mut a, &mut state);
        a
    });
    let atgamma: [Poly; R] = std::array::from_fn(|i| {
        let mut a = Poly::new();
        for j in 0..(N-1) {
            a.coeffs[j] = 2*gamma[i].coeffs[j] - gamma[i].coeffs[j+1];
        }
        a.coeffs[N-1] = 2*gamma[i].coeffs[N-1];
        a
    });
    let mut bpr = ckp.bm[0].clone();
    bpr.add(&ckp.bm[1]);
    bpr.sub(&ckp.bm[3]);

    let mut tpr = tp.tm.vec[0].clone();
    tpr.add(&tp.tm.vec[1]);
    tpr.sub(&tp.tm.vec[3]);

    let mut hhat = h.clone();
    hhat.ntt();
    for i in 0..4 {
        let mut zshat = z[i].s.clone();
        zshat.vec_ntt();
        
        let mut vtmp1 = PolyVecL::pointwise_acc_montgomery(&bpr, &zshat);
        let mut vtmp2 = PolyVecL::pointwise_acc_montgomery(&ckp.bm[3], &zshat);
        vprime[i] = PolyVecL::pointwise_acc_montgomery(&ckp.bm[M-1], &zshat);

        let mut chat = c[i].clone();
        chat.ntt();
        let mut tmp = Poly::new();
        
        Poly::pointwise_montgomery_other(&mut tmp, &chat, &tpr);
        vtmp1.sub(&tmp);
        vtmp1.scale_montgomery(MONTSQ as i32);
        
        Poly::pointwise_montgomery_other(&mut tmp, &chat, &tp.tm.vec[3]);
        vtmp2.sub(&tmp);
        vtmp2.scale_montgomery(MONTSQ as i32);
        
        Poly::sub_other(&mut tmp, &tp.tm.vec[5], &hhat);
        tmp.pointwise_montgomery(&chat);
        vprime[i].sub(&tmp);
        vprime[i].scale_montgomery(MONTSQ as i32);

        Poly::add_other(&mut tmp, &z[i].em.vec[0], &z[i].em.vec[1]);
        tmp.sub(&z[i].em.vec[3]);
        tmp.ntt();
        vtmp1.add(&tmp);

        tmp = z[i].em.vec[3];
        tmp.ntt();
        vtmp2.add(&tmp);

        tmp = z[i].em.vec[5];
        tmp.ntt();
        vprime[i].add(&tmp);

        for j in 0..R {
            Poly::pointwise_montgomery_other(&mut tmp, &gamma[j], &vtmp1);
            Poly::pointwise_montgomery_other(&mut chat, &atgamma[j], &vtmp2);
            tmp.sub(&chat);
            tmp.trace65_ntt();
            let tmp2 = tmp.clone();
            poly_shift(&mut tmp, &tmp2, j);
            vprime[i].add(&tmp);
        }

        vprime[i].freeze();
    }

    false
}
