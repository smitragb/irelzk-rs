#![allow(dead_code)]
#![allow(unused_imports)]
use crate::{
    crypto::{aes256::Aes256Ctx, shake::Shake128}, 
    params::{BETA, GAMMA1, GAMMA2, K, L, M, N, R, SYMBYTES}, 
    poly_arith::{
        consts::MONTSQ, poly::Poly, polyvec::{PolyVecK, PolyVecL, PolyVecM}
    }
};

use super::commitment::{self, Comm, CommKey, CommRnd};

pub fn challenge_prehash (c: &mut [Poly; R], chash: &[u8; N/4]) {
    let lut: Vec<i32> = vec![0, 0, 1, -1];
    for i in 0..R {
        for j in 0..(N/(4*R)) {
            let idx = (N/(4*R))*i + j;
            let b = chash[idx];
            c[i].coeffs[16 * j + 0 * R] = lut[((b >> 0) & 3) as usize];
            c[i].coeffs[16 * j + 1 * R] = lut[((b >> 2) & 3) as usize];
            c[i].coeffs[16 * j + 2 * R] = lut[((b >> 4) & 3) as usize];
            c[i].coeffs[16 * j + 3 * R] = lut[((b >> 6) & 3) as usize];
        }
    }
}

pub fn challenge (w: &[PolyVecK; R]) {
    let mut chash = [0u8; N/4];
    let bytes = bytemuck::cast_slice(w);
    Shake128::hash(&mut chash, &bytes);
    let mut c: [Poly; R] = std::array::from_fn(|_| Poly::new());
    challenge_prehash(&mut c, &chash);
}

pub fn generate_y (seed: &[u8; SYMBYTES], nonce: u16) -> [CommRnd; R] {
    let mut n = nonce as u64;
    let mut state = Aes256Ctx::init(&seed, n);

    std::array::from_fn(|_| {
        let y = CommRnd::generate_y(&mut state, n);
        n += (L+M) as u64;
        y
    })
}

pub fn first (
    w1: &mut [PolyVecK; R],
     g: &mut [PolyVecM; R],
     y: &mut [CommRnd; R],
    ck: &CommKey
) {
    for i in 0..R {
        
        y[i].s.vec_ntt();
        y[i].em.vec_ntt();
    
        for j in 0..K {
            w1[i].vec[j] = PolyVecL::pointwise_acc_montgomery(&ck.b0[j], &y[i].s);
        } 

        for j in 0..K {
            let tmp = PolyVecM::pointwise_acc_montgomery(&ck.bt[j], &y[i].em);
            w1[i].vec[j].add(&tmp);
        }

        w1[i].vec_inverse_ntt_tomont();
        w1[i].vec_decompose(&mut y[i].e);

        for j in 0..M {
            g[i].vec[j] = PolyVecL::pointwise_acc_montgomery(&ck.bm[j], &y[i].s);
        }
        g[i].scale_montgomery(MONTSQ as i32);
        g[i].add(&y[i].em);
    }
}

pub fn last (
     z: &mut [CommRnd; R],
     y: &mut [CommRnd; R],
     r: &CommRnd,
     c: &[Poly; R],
    w1: &[PolyVecK; R],
 t0low: &PolyVecK
) -> bool {
    let chat: [Poly; R] = std::array::from_fn(|i| {
        let mut chat_i = c[i].clone();
        chat_i.ntt();
        chat_i.scale_montgomery(MONTSQ as i32);
        chat_i
    });

    for i in 0..R {
        for j in 0..L {
            Poly::pointwise_montgomery_other(&mut z[i].s.vec[j], &chat[i], &r.s.vec[j]);
        }
        z[i].s.add(&y[i].s);
        z[i].s.vec_inverse_ntt();
        z[i].s.reduce();
        let mut bound = GAMMA1 as u32 - BETA as u32;
        if PolyVecL::vec_check_norm(&z[i].s, bound) {
            return true;
        }

        for j in 0..M {
            Poly::pointwise_montgomery_other(&mut z[i].em.vec[j], &chat[i], &r.em.vec[j]);
        }
        z[i].em.add(&y[i].em);
        z[i].em.vec_inverse_ntt();
        z[i].em.reduce();
        if PolyVecM::vec_check_norm(&z[i].em, bound) {
            return true;
        }

        for j in 0..K {
            Poly::pointwise_montgomery_other(&mut z[i].e.vec[j], &chat[i], &r.e.vec[j]);
        }
        z[i].e.vec_inverse_ntt();
        y[i].e.sub(&z[i].e);
        bound = GAMMA2 as u32 - BETA as u32;
        if PolyVecK::vec_check_norm(&y[i].e, bound) {
            return true;
        }
    }

    for i in 0..R {
        for j in 0..K {
            Poly::pointwise_montgomery_other(&mut z[i].e.vec[j], &chat[i], &t0low.vec[j]);
        }
        z[i].e.vec_inverse_ntt();
        if PolyVecK::vec_check_norm(&z[i].e, GAMMA2 as u32) {
            return true;
        }
        y[i].e.add(&z[i].e);
        PolyVecK::vec_makehint(&mut z[i].e, &w1[i], &mut y[i].e);
    }

    return false;
}

pub fn verify_first (
    w1: &mut [PolyVecK; R],
     c: &[Poly; R],
     z: &[CommRnd; R],
    tp: &Comm,
   ckp: &CommKey
) -> bool {
    for i in 0..R {
        let bound = GAMMA1 as u32 - BETA as u32;
        return PolyVecL::vec_check_norm(&z[i].s, bound) || PolyVecM::vec_check_norm(&z[i].em, bound);
    }

    for i in 0..R {
        let mut zshat = z[i].s;
        let mut zmhat = z[i].em;
        zshat.vec_ntt();
        zmhat.vec_ntt();
        
        for j in 0..K {
            w1[i].vec[j] = PolyVecL::pointwise_acc_montgomery(&ckp.b0[j], &zshat);
        }
        for j in 0..K {
            let tmp = PolyVecM::pointwise_acc_montgomery(&ckp.bt[j], &zmhat);
            w1[i].vec[j].add(&tmp);
        }

        let mut chat = c[i];
        chat.ntt();
        for j in 0..K {
            let mut tmp = Poly::new();
            Poly::scale_montgomery_other(&mut tmp, &tp.t0.vec[j], 4128752);
            tmp.pointwise_montgomery(&chat);
            w1[i].vec[j].sub(&tmp);
        }
        w1[i].vec_inverse_ntt_tomont();
        w1[i].vec_usehint(&z[i].e); 
    } 
    return false;
}

pub fn verify_last (
        c: &[Poly; R],
    chash: &[u8; N/4]
) -> bool {
    let mut c2: [Poly; R] = std::array::from_fn(|_| Poly::new());
    challenge_prehash(&mut c2, chash);
    for i in 0..R {
        for j in 0..N {
            if c[i].coeffs[j] != c2[i].coeffs[j] {
                return true;
            }
        }
    }
    return false;
}
