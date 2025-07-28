#![allow(dead_code)]
#![allow(unused_assignments)]

use bytemuck::{bytes_of, cast_slice};
use rand::{rngs::OsRng, RngCore};

use crate::{
    add,
    comm::{
        commitment::{Comm, CommKey, CommRnd}, 
        opening
    }, 
    crypto::shake::{Shake128, SHAKE128_RATE}, 
    params::{K, L, M, N, R, SYMBYTES}, 
    poly_arith::{
        poly::Poly, 
        polyvec::{PolyVecK, PolyVecM}
    }
};

pub struct Proof {
    h: Poly,
    c: [Poly; R],
    z: [CommRnd; R],
}

impl Proof {
    pub fn prove (
        rho: &[u8; SYMBYTES],
          a: &[u64; 2],
          b: &[u64; 2],
    ) -> (Proof, Comm) 
    {
        let mut nonce = 0;
        let mut seed  = [0u8; SYMBYTES];
        let mut thash = [0u8; SYMBYTES];
        let mut chash = [0u8; SHAKE128_RATE];

        OsRng.fill_bytes(&mut seed);
        let mut x = 0;
        
        let mut msg = PolyVecM::new();
        for i in 0..64 {
            let f = ((a[0] >> i) & 1) as i32;
            let g = ((b[0] >> i) & 1) as i32;
            x = x + f + g;
            msg.vec[0].coeffs[i] = f;
            msg.vec[1].coeffs[i] = g;
            msg.vec[2].coeffs[i] = x & 1;
            x >>= 1;
            msg.vec[3].coeffs[i] = x;
        }
        for i in 0..64 {
            let f = ((a[1] >> i) & 1) as i32;
            let g = ((b[1] >> i) & 1) as i32;
            x = x + f + g;
            msg.vec[0].coeffs[64 + i] = f;
            msg.vec[1].coeffs[64 + i] = g;
            msg.vec[2].coeffs[64 + i] = x & 1;
            x >>= 1;
            msg.vec[3].coeffs[64 + i] = x;
        }
        
        Poly::uniform_random(&mut msg.vec[M-1], &seed, nonce);
        nonce += 1;
        for i in 0..R {
            msg.vec[M-1].coeffs[i] = 0;
        }
        msg.vec[M-1].ntt();

        let ck = CommKey::expand(&rho);
        let mut r = CommRnd::generate();
        let mut t = Comm::commit(&ck, &mut r, &msg);

        t.t0.vec_inverse_ntt();
        let mut t0low = PolyVecK::new();
        PolyVecK::vec_power2round(&mut t.t0, &mut t0low);
        t.t0.vec_ntt();
        t0low.vec_ntt();
       
        let tmslice: &[Poly] = &t.tm.vec[..(M-1)]; 
        let t0bytes: &[u8] = bytes_of(&t.t0);
        let tmbytes: &[u8] = cast_slice(tmslice);
        let lastbytes: &[u8] = bytes_of(&t.tm.vec[M-1]);

        let mut shake128_state = Shake128::init();
        shake128_state.absorb(rho);
        shake128_state.absorb(t0bytes);
        shake128_state.absorb(tmbytes);
        shake128_state.absorb(lastbytes);
        shake128_state.finalize();
        shake128_state.squeeze(&mut thash);

        let mut tmp = Poly::new();
        let mut c = [Poly::new(); R];
        let mut z: [CommRnd; R] = std::array::from_fn(|_| CommRnd::new());
        let mut h = Poly::new();

        loop {
            let mut  w1: [PolyVecK; R] = std::array::from_fn(|_| PolyVecK::new());
            let mut   g: [PolyVecM; R] = std::array::from_fn(|_| PolyVecM::new());
            let mut vpr = [Poly::new(); R];
            let mut   y = opening::generate_y(&seed, nonce); 
            opening::first(&mut w1, &mut g, &mut y, &ck);
            nonce += (R*(K+L+M)) as u16;
            
            let w1bytes: &[u8] = cast_slice(&w1);
            shake128_state = Shake128::init();
            shake128_state.absorb(&thash);
            shake128_state.absorb(w1bytes);
            shake128_state.finalize();
            shake128_state.squeezeblocks(&mut chash, 1);

            let prod_slice: &[u8; SYMBYTES] = (&chash[..SYMBYTES]).try_into().unwrap();
            let lin_slice : &[u8; SYMBYTES] = (&chash[SYMBYTES..(2*SYMBYTES)]).try_into().unwrap(); 
            let v = add::product::proof(&mut msg, &g, prod_slice);
            Poly::add_other(&mut tmp, &t.tm.vec[M-2], &msg.vec[M-2]); 
            h = add::linear::proof(&mut vpr, &msg, lin_slice, &g); 

            let vprbytes: &[u8] = cast_slice(&vpr);
            shake128_state = Shake128::init();
            shake128_state.absorb(prod_slice);
            shake128_state.absorb(lin_slice);
            shake128_state.absorb(bytes_of(&tmp));
            shake128_state.absorb(bytes_of(&v));
            shake128_state.absorb(bytes_of(&h));
            shake128_state.absorb(vprbytes);
            shake128_state.finalize();
            shake128_state.squeezeblocks(&mut chash, 1);

            let cslice: &[u8; N/4] = (&chash[..(N/4)]).try_into().unwrap();
            opening::challenge_prehash(&mut c, cslice);


            if !opening::last(&mut z, &mut y, &r, &c, &w1, &t0low) {
                break;
            }
        }
        t.tm.vec[M-2] = tmp.clone();
        ( Proof { h, c, z }, t )
    }

    pub fn verify (
          p: &Proof,
          t: &Comm,
        rho: &[u8; SYMBYTES]
    ) -> bool {
        let mut thash = [0u8; SYMBYTES];
        let mut chash = [0u8; SHAKE128_RATE];
        let mut shake128_state = Shake128::init();

        let tmslice: &[Poly] = &t.tm.vec[..(M-1)]; 
        let t0bytes: &[u8] = bytes_of(&t.t0);
        let tmbytes: &[u8] = cast_slice(tmslice);
        let lastbytes: &[u8] = bytes_of(&t.tm.vec[M-1]);
        shake128_state.absorb(rho);
        shake128_state.absorb(t0bytes);
        shake128_state.absorb(tmbytes);
        shake128_state.absorb(lastbytes);
        shake128_state.finalize();
        shake128_state.squeeze(&mut thash);

        let ck = CommKey::expand(rho);
        let mut w1: [PolyVecK; R] = std::array::from_fn(|_| PolyVecK::new());
        if opening::verify_first(&mut w1, &p.c, &p.z, t, &ck) {
            return true
        }

        let w1bytes: &[u8] = cast_slice(&w1);
        shake128_state = Shake128::init();
        shake128_state.absorb(&thash);
        shake128_state.absorb(&w1bytes);
        shake128_state.finalize();
        shake128_state.squeezeblocks(&mut chash, 1);

        let mut v = Poly::new();
        let mut vpr = [Poly::new(); R];
        let prod_slice: &[u8; SYMBYTES] = (&chash[..SYMBYTES]).try_into().unwrap();
        let lin_slice : &[u8; SYMBYTES] = (&chash[SYMBYTES..(2*SYMBYTES)]).try_into().unwrap();
        if add::product::verify(&mut v, prod_slice, &p.c, &p.z, t, &ck) {
            return true
        }
        if add::linear::verify(&mut vpr, lin_slice, &p.h, &p.c, &p.z, t, &ck) {
            return true
        }


        let vprbytes: &[u8] = cast_slice(&vpr);
        shake128_state = Shake128::init();
        shake128_state.absorb(prod_slice);
        shake128_state.absorb(lin_slice);
        shake128_state.absorb(bytes_of(&t.tm.vec[M-2]));
        shake128_state.absorb(bytes_of(&v));
        shake128_state.absorb(bytes_of(&p.h));
        shake128_state.absorb(vprbytes);
        shake128_state.finalize();
        shake128_state.squeezeblocks(&mut chash, 1);

        let cslice: &[u8; N/4] = (&chash[..(N/4)]).try_into().unwrap();
        if opening::verify_last(&p.c, cslice) {
            return true
        }

        false
    }
}
