#![allow(dead_code)]
#![allow(unused_imports)]

use rand::{rngs::OsRng, RngCore};

use crate::{
    crypto::aes256::Aes256Ctx, 
    params::{K, M, SYMBYTES}, 
    poly_arith::{
        consts::MONTSQ, poly::Poly, polyvec::{PolyVecK, PolyVecL, PolyVecM}
    }
};

pub struct Comm {
   t0: PolyVecK,
   tm: PolyVecM,
}

pub struct CommRnd {
    s: PolyVecL,
    e: PolyVecK,
    em: PolyVecM,
}

pub struct CommKey {
    b0: [PolyVecL; K],
    bt: [PolyVecM; K],
    bm: [PolyVecL; M],
}

impl CommKey {
    pub fn expand(rho: &[u8; SYMBYTES]) -> Self {
        let mut state = Aes256Ctx::init(rho, 0);
        let b0: [PolyVecL; K] = std::array::from_fn(|i| {
            PolyVecL {
                vec: std::array::from_fn(|j| {
                    let nonce = ((i as u64)<<16) + (j as u64);
                    state.select(nonce);
                    let mut a = Poly::new();
                    Poly::uniform_preinit(&mut a, &mut state);
                    a
                })
            }
        });

        let bt: [PolyVecM; K] = std::array::from_fn(|i| {
            PolyVecM {
                vec: std::array::from_fn(|j| {
                    let nonce = (((K+i) as u64) << 16) + (j as u64);
                    state.select(nonce);
                    let mut a = Poly::new();
                    Poly::uniform_preinit(&mut a, &mut state);
                    a
                })
            } 
        });
        
        let bm: [PolyVecL; M] = std::array::from_fn(|i| {
            PolyVecL {
                vec: std::array::from_fn(|j| {
                    let nonce = (((2*K+i) as u64) <<16) + (j as u64);
                    state.select(nonce);
                    let mut a = Poly::new();
                    Poly::uniform_preinit(&mut a, &mut state);
                    a
                })
            }
        });

        Self { b0, bt, bm }
    }
}

impl CommRnd {
    pub fn generate() -> Self {
        let mut buf = [0u8; SYMBYTES];
        OsRng.fill_bytes(&mut buf);
        let mut nonce = 0;
        let mut state = Aes256Ctx::init(&buf, nonce);
        let mut make_vec = || {
            nonce += 1;
            let mut a = Poly::new();
            state.select(nonce);
            Poly::trinary_preinit(&mut a, &mut state);
            a
        };
        Self {
             s: PolyVecL { vec: std::array::from_fn(|_| make_vec()) },
             e: PolyVecK { vec: std::array::from_fn(|_| make_vec()) },
            em: PolyVecM { vec: std::array::from_fn(|_| make_vec()) },
        }
    }    
}

impl Comm {
    pub fn commit(ck: &CommKey, r: &mut CommRnd, msg: &PolyVecM) -> Self {
        r.s.vec_ntt();
        r.e.vec_ntt();
        r.em.vec_ntt();
        
        let mut t0: PolyVecK = PolyVecK {
            vec: std::array::from_fn(|i| {
                PolyVecL::pointwise_acc_montgomery(&ck.b0[i], &r.s)
            })
        };
        let mut tm: PolyVecM = PolyVecM {
            vec: std::array::from_fn(|i| {
                PolyVecL::pointwise_acc_montgomery(&ck.bm[i], &r.s)
            })
        };
        let tag: PolyVecK = PolyVecK {
            vec: std::array::from_fn(|i| {
                PolyVecM::pointwise_acc_montgomery(&ck.bt[i], &r.em)
            })
        };
        t0.add(&tag);
        t0.scale_montgomery(MONTSQ as i32);
        t0.add(&r.e);
        tm.scale_montgomery(MONTSQ as i32);
        tm.add(&r.em);
        tm.add(msg);
        Self { t0, tm }
    }
}
