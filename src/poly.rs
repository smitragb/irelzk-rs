#![allow(dead_code)]
use std::{arch::x86_64::*, ops::{Deref, DerefMut}};

use crate::{
    crypto::aes256::{Aes256Ctx, AES256CTR_BLOCKBYTES}, 
    params::{SYMBYTES, N, Q},
    consts::{QDATA, QINV, REJIDX, _8XQ}, 
    ntt::{forward_ntt, inverse_ntt}
};

#[repr(align(32))]
#[derive(Debug, Clone)]
pub struct Poly{
    pub coeffs: [i32; N],
}

#[repr(align(32))]
pub struct AlignedBuf<const S: usize>([u8; S]);

impl<const S: usize> Deref for AlignedBuf<S> {
    type Target = [u8; S]; 
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const S:usize> DerefMut for AlignedBuf<S> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

pub const REJ_UNIFORM_BUFLEN: usize = ((512+AES256CTR_BLOCKBYTES - 1)/AES256CTR_BLOCKBYTES)*AES256CTR_BLOCKBYTES;
pub const POLY_UNIFORM_NBLOCKS: usize = (512+AES256CTR_BLOCKBYTES - 1)/AES256CTR_BLOCKBYTES;

impl Poly {
    pub fn new() -> Self {
        Self {
            coeffs: [0i32; N],
        }
    }

    pub fn reduce (&mut self) {
        unsafe {
            let mask = _mm256_set1_epi32((1<<30) - 1);
            let coeffs_ptr = self.coeffs.as_mut_ptr();
            for i in 0..(N/8) {
                let mut f = _mm256_load_si256(coeffs_ptr.add(8*i) as *const __m256i);
                let mut t = _mm256_srai_epi32(f, 30);
                f = _mm256_and_si256(f, mask);
                f = _mm256_sub_epi32(f, t);
                t = _mm256_slli_epi32(t, 18);
                f = _mm256_add_epi32(f, t);
                _mm256_store_si256(coeffs_ptr.add(8*i) as *mut __m256i, f);
            }
        }    
    }

    pub fn ntt (&mut self) {
        forward_ntt(&mut self.coeffs, &QDATA.0); 
        self.reduce();
    } 

    pub fn scale_montgomery(&mut self, s: i32) {
        unsafe {
            let qdata_ptr = QDATA.0.as_ptr();
            let self_ptr = self.coeffs.as_mut_ptr();
            let prod = ((s as i64) * (QINV as i64)) as i32;

            let q = _mm256_load_si256(qdata_ptr.add(_8XQ) as *const __m256i);
            let lo = _mm256_set1_epi32(prod);
            let hi = _mm256_set1_epi32(s);

            for i in (0..N).step_by(8) {
                let mut f0 = _mm256_load_si256(self_ptr.add(i) as *const __m256i);
                let mut f1 = _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(f0)));
                let mut g0 = _mm256_mul_epi32(f0, lo);
                let mut g1 = _mm256_mul_epi32(f1, lo);
                f0 = _mm256_mul_epi32(f0, hi);
                f1 = _mm256_mul_epi32(f1, hi);
                g0 = _mm256_mul_epi32(g0, q);
                g1 = _mm256_mul_epi32(g1, q);
                f0 = _mm256_sub_epi32(f0, g0);
                f1 = _mm256_sub_epi32(f1, g1);
                f0 = _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(f0)));
                f0 = _mm256_blend_epi32(f0, f1, 0xAA);
                _mm256_store_si256(self_ptr.add(i) as *mut __m256i, f0);
            }

        }
    }

    pub fn scale_montgomery_other(r: &mut Poly, a: &Poly, s: i32) {
        unsafe {
            let qdata_ptr = QDATA.0.as_ptr();
            let a_ptr = a.coeffs.as_ptr();
            let r_ptr = r.coeffs.as_mut_ptr();
            let prod = ((s as i64) * (QINV as i64)) as i32;

            let q = _mm256_load_si256(qdata_ptr.add(_8XQ) as *const __m256i);
            let lo = _mm256_set1_epi32(prod);
            let hi = _mm256_set1_epi32(s);

            for i in (0..N).step_by(8) {
                let mut f0 = _mm256_load_si256(a_ptr.add(i) as *const __m256i);
                let mut f1 = _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(f0)));
                let mut g0 = _mm256_mul_epi32(f0, lo);
                let mut g1 = _mm256_mul_epi32(f1, lo);
                f0 = _mm256_mul_epi32(f0, hi);
                f1 = _mm256_mul_epi32(f1, hi);
                g0 = _mm256_mul_epi32(g0, q);
                g1 = _mm256_mul_epi32(g1, q);
                f0 = _mm256_sub_epi32(f0, g0);
                f1 = _mm256_sub_epi32(f1, g1);
                f0 = _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(f0)));
                f0 = _mm256_blend_epi32(f0, f1, 0xAA);
                _mm256_store_si256(r_ptr.add(i) as *mut __m256i, f0);
            }
        }
    }

    pub fn inverse_ntt(&mut self) {
        self.scale_montgomery(33554432);
        inverse_ntt(&mut self.coeffs, &QDATA.0);
    }

    pub fn inverse_ntt_tomont(&mut self) {
        self.scale_montgomery(-132153352);
        inverse_ntt(&mut self.coeffs, &QDATA.0);
    }

    pub fn rej_uniform (r: &mut [i32], buf: &[u8]) -> usize {
        let mut ctr = 0;
        let mut pos = 0;
        let len = r.len();
        let buflen = buf.len();

        while ctr < len && pos+4 <= buflen {
            let mut t = buf[pos] as i32;
            t |= (buf[pos+1] as i32) << 8;
            t |= (buf[pos+2] as i32) << 16;
            t |= (buf[pos+3] as i32) << 24;
            t &= (1<<30) - 1;

            pos += 4;
            if t < Q {
                r[ctr] = t;
                ctr += 1;
            }
        }
        ctr
    }
    
    pub fn rej_uniform_avx (r: &mut [i32; 128], buf: &[u8; REJ_UNIFORM_BUFLEN]) -> usize {
        let mut ctr = 0;
        let mut pos = 0;

        let qdata_ptr = QDATA.0.as_ptr();
        let buf_ptr = buf.as_ptr();
        let r_ptr = r.as_mut_ptr();
        unsafe {
            let bound = _mm256_load_si256(qdata_ptr.add(_8XQ) as *const __m256i);
            let mask = _mm256_set1_epi32((1<<30) - 1);

            while pos <= REJ_UNIFORM_BUFLEN - 32 && ctr <= 120 {
                let mut d = _mm256_load_si256(buf_ptr.add(pos) as *const __m256i);
                d = _mm256_and_si256(d, mask);
                pos += 32;

                let mut tmp = _mm256_sub_epi32(d, bound);
                let good = _mm256_movemask_ps(_mm256_castsi256_ps(tmp));
                if good == 255 {
                    _mm256_storeu_si256(r_ptr.add(ctr) as *mut __m256i, d);
                    ctr += 8;
                    continue;
                }

                let idx_ptr = REJIDX[good as usize].as_ptr() as *const __m128i;
                tmp = _mm256_cvtepu8_epi32(_mm_loadl_epi64(idx_ptr));
                d = _mm256_permutevar8x32_epi32(d, tmp);
                _mm256_storeu_si256(r_ptr.add(ctr) as *mut __m256i, d);
                ctr += (good as usize).count_ones() as usize;
            }
        }
        ctr
    }

    pub fn uniform_preinit(r: &mut Poly, state: &mut Aes256Ctx) {
        let mut ctr = 0;
        const BUFSIZE: usize = POLY_UNIFORM_NBLOCKS * AES256CTR_BLOCKBYTES;
        let mut buf = AlignedBuf::<BUFSIZE>([0u8; BUFSIZE]);
        unsafe { 
            state.squeezeblocks(&mut buf.0, POLY_UNIFORM_NBLOCKS);
            ctr += Self::rej_uniform_avx(&mut r.coeffs, &buf.0);

            while ctr < N {
                state.squeezeblocks(&mut buf.0, 1);
                ctr += Self::rej_uniform(&mut r.coeffs[ctr..], &buf.0)
            }
        }
    }

    pub fn uniform_random(seed: &[u8; SYMBYTES], nonce: u16) -> Poly {
        let mut r = Poly::new();
        unsafe{
            let mut state = Aes256Ctx::init(seed, nonce as u64);
            Self::uniform_preinit(&mut r, &mut state);
        }
        r
    }
}
