#![allow(dead_code)]
use std::{arch::x86_64::*, ops::{Deref, DerefMut}};

use crate::{
    consts::{QDATA, QINV, REJIDX, _8XQ, _8XQINV}, crypto::aes256::{Aes256Ctx, AES256CTR_BLOCKBYTES}, ntt::{forward_ntt, inverse_ntt}, params::{N, Q, SYMBYTES}, rounding::{decompose_avx, makehint_avx, power2round_avx, usehint_avx}
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
pub const POLY_UNIFORM_GAMMA_NBLOCKS: usize = (304+AES256CTR_BLOCKBYTES - 1)/AES256CTR_BLOCKBYTES;
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

    pub fn uniform_random(r: &mut Poly, seed: &[u8; SYMBYTES], nonce: u16) {
        unsafe{
            let mut state = Aes256Ctx::init(seed, nonce as u64);
            Self::uniform_preinit(r, &mut state);
        }
    }

    pub fn trinary_preinit(r: &mut Poly, state: &mut Aes256Ctx) {
        const BUFSIZE: usize = N/2;
        let mut buf = AlignedBuf::<BUFSIZE>([0u8; BUFSIZE]); 
        unsafe {
            let lut = _mm256_set1_epi32(0xA815);
            let mask32 = _mm256_cmpeq_epi32(lut, lut);
            let mask4 = _mm256_srli_epi32(mask32, 28);
            let mask2 = _mm256_srli_epi32(mask32, 30);
            let r_ptr = r.coeffs.as_mut_ptr();

            state.squeezeblocks(&mut buf.0, 1);

            let buf_ptr = buf.0.as_ptr();
            for i in 0..(N/16) {
                let mut f = _mm256_cvtepu8_epi32(_mm_loadl_epi64(buf_ptr.add(8*i) as *const __m128i)); 
                let mut g = _mm256_srli_epi32(f, 4);
                f = _mm256_and_si256(f, mask4);
                f = _mm256_srlv_epi32(lut, f);
                g = _mm256_srlv_epi32(lut, g);
                let h = _mm256_unpacklo_epi32(f, g);
                g = _mm256_unpackhi_epi32(f, g);
                f = _mm256_permute2x128_si256(h, g, 0x20);
                g = _mm256_permute2x128_si256(h, g, 0x31);
                f = _mm256_and_si256(f, mask2);
                g = _mm256_and_si256(g, mask2);
                f = _mm256_add_epi32(f, mask32);
                g = _mm256_add_epi32(g, mask32);
                _mm256_store_si256(r_ptr.add(16*i + 0) as *mut __m256i, f);
                _mm256_store_si256(r_ptr.add(16*i + 8) as *mut __m256i, g);
            }
        }
    }

    pub fn trinary(r: &mut Poly, seed: &[u8; SYMBYTES], nonce: u16) {
        unsafe { 
            let mut state = Aes256Ctx::init(seed, nonce as u64); 
            Self::trinary_preinit(r, &mut state);
        }
    }

    pub fn uniform_gamma_preinit(r: &mut Poly, state: &mut Aes256Ctx) {
        const BUFSIZE:usize = POLY_UNIFORM_GAMMA_NBLOCKS * AES256CTR_BLOCKBYTES;
        let mut buf = AlignedBuf::<BUFSIZE>([0u8; BUFSIZE]);
        let mut pos = 0;
        unsafe {
            let mask  = _mm256_set1_epi32(0x7FFFF);
            let min   = _mm256_set1_epi32(-(1 << 18));
            let idx32 = _mm256_set_epi32(5,2,7,4,1,6,3,0);
            let idx8  = _mm256_set_epi8(-1,10, 9, 8,-1, 8, 7, 6,
                                         6, 5, 4, 3,-1, 3, 2, 1,
                                        -1, 9, 8, 7, 7, 6, 5, 4,
                                        -1, 4, 3, 2,-1, 2, 1, 0);

            state.squeezeblocks(&mut buf.0, POLY_UNIFORM_GAMMA_NBLOCKS);
            let buf_ptr = buf.0.as_ptr();
            let r_ptr = r.coeffs.as_mut_ptr();

            for i in 0..(N/8) {
                let mut f = _mm256_loadu_si256(buf_ptr.add(pos) as *const __m256i);    
                pos += 19;
                f = _mm256_permute4x64_epi64(f, 0x94);
                f = _mm256_shuffle_epi8(f, idx8);
                f = _mm256_srlv_epi32(f, idx32);
                f = _mm256_and_si256(f, mask);
                f = _mm256_add_epi32(f, min);
                _mm256_store_si256(r_ptr.add(8*i) as *mut __m256i, f);
            }
        }
    }

    pub fn uniform_gamma(r: &mut Poly, seed: &[u8; SYMBYTES], nonce: u16) {
        unsafe {
            let mut state = Aes256Ctx::init(seed, nonce as u64);
            Self::uniform_gamma_preinit(r, &mut state);
        }
    }
    
    pub fn freeze(&mut self) {
        unsafe {
            let qdata_ptr = QDATA.0.as_ptr();
            let q = _mm256_load_si256(qdata_ptr.add(_8XQ) as *const __m256i);
            let hq = _mm256_srli_epi32(q, 1);

            let coeffs_ptr = self.coeffs.as_mut_ptr();
            
            for i in (0..N).step_by(8) {
                let mut f = _mm256_load_si256(coeffs_ptr.add(i) as *const __m256i);
                let mut t = _mm256_cmpgt_epi32(f, hq);
                t = _mm256_and_si256(t, q);
                f = _mm256_sub_epi32(f, t);
                _mm256_store_si256(coeffs_ptr.add(i) as *mut __m256i, f);
            }

        }
    } 

    pub fn add(&mut self, other: &Poly) {
        let coeffs_ptr = self.coeffs.as_mut_ptr();
        let other_coeffs_ptr = other.coeffs.as_ptr();
        unsafe {
            for i in (0..N).step_by(8) {
                let f = _mm256_load_si256(coeffs_ptr.add(i) as *const __m256i);
                let g = _mm256_load_si256(other_coeffs_ptr.add(8) as *const __m256i);
                let h = _mm256_add_epi32(f, g);
                _mm256_store_si256(coeffs_ptr.add(i) as *mut __m256i, h);
            }
            self.reduce();
        }
    }

    pub fn add_other(r: &mut Poly, a: &Poly, b: &Poly) {
        let a_ptr = a.coeffs.as_ptr();
        let b_ptr = b.coeffs.as_ptr();
        let r_ptr = r.coeffs.as_mut_ptr();
        unsafe {
            for i in (0..N).step_by(8) {
                let f = _mm256_load_si256(a_ptr.add(i) as *const __m256i);
                let g = _mm256_load_si256(b_ptr.add(i) as *const __m256i);
                let h = _mm256_add_epi32(f, g);
                _mm256_store_si256(r_ptr.add(i) as *mut __m256i, h);
            }
            r.reduce();
        }
    }

    pub fn sub(&mut self, other: &Poly) {
        let coeffs_ptr = self.coeffs.as_mut_ptr();
        let other_coeffs_ptr = other.coeffs.as_ptr();
        unsafe {
            for i in (0..N).step_by(8) {
                let f = _mm256_load_si256(coeffs_ptr.add(i) as *const __m256i);
                let g = _mm256_load_si256(other_coeffs_ptr.add(8) as *const __m256i);
                let h = _mm256_sub_epi32(f, g);
                _mm256_store_si256(coeffs_ptr.add(i) as *mut __m256i, h);
            }
            self.reduce();
        }
    }

    pub fn sub_other(r: &mut Poly, a: &Poly, b: &Poly) {
        let a_ptr = a.coeffs.as_ptr();
        let b_ptr = b.coeffs.as_ptr();
        let r_ptr = r.coeffs.as_mut_ptr();
        unsafe {
            for i in (0..N).step_by(8) {
                let f = _mm256_load_si256(a_ptr.add(i) as *const __m256i);
                let g = _mm256_load_si256(b_ptr.add(i) as *const __m256i);
                let h = _mm256_sub_epi32(f, g);
                _mm256_store_si256(r_ptr.add(i) as *mut __m256i, h);
            }
            r.reduce();
        }
    }

    pub fn pointwise_montgomery(&mut self, other: &Poly) {
        let qdata_ptr = QDATA.0.as_ptr();
        let coeffs_ptr = self.coeffs.as_mut_ptr();
        let other_coeffs_ptr = other.coeffs.as_ptr();
        unsafe {
            let q    = _mm256_load_si256(qdata_ptr.add(_8XQ) as *const __m256i);
            let qinv = _mm256_load_si256(qdata_ptr.add(_8XQINV) as *const __m256i);
            for i in (0..N).step_by(8) {
                let mut f0 = _mm256_load_si256(coeffs_ptr.add(i) as *const __m256i);
                let mut g0 = _mm256_load_si256(other_coeffs_ptr.add(i) as *const __m256i);
                let mut f1 = _mm256_mul_epi32(f0, g0);
                        f0 = _mm256_srli_epi64(f0, 32);
                        g0 = _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(g0)));
                        f0 = _mm256_mul_epi32(f0, g0);
                let mut g1 = _mm256_mul_epi32(f1, qinv);
                        g0 = _mm256_mul_epi32(f0, qinv);
                        g1 = _mm256_mul_epi32(g1, q);
                        g0 = _mm256_mul_epi32(g0, q);
                        f1 = _mm256_sub_epi32(f1, g1);
                        f0 = _mm256_sub_epi32(f0, g0);
                        f1 = _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(f1)));
                        f0 = _mm256_blend_epi32(f1, f0, 0xaa);
                        _mm256_store_si256(coeffs_ptr.add(i) as *mut __m256i, f0);
            }
        }
    }

    pub fn pointwise_montgomery_other(r: &mut Poly, a: &Poly, b: &Poly) {
        let qdata_ptr = QDATA.0.as_ptr();
        let r_ptr = r.coeffs.as_mut_ptr();
        let a_ptr = a.coeffs.as_ptr();
        let b_ptr = b.coeffs.as_ptr();
        unsafe {
            let q    = _mm256_load_si256(qdata_ptr.add(_8XQ) as *const __m256i);
            let qinv = _mm256_load_si256(qdata_ptr.add(_8XQINV) as *const __m256i);
            for i in (0..N).step_by(8) {
                let mut f0 = _mm256_load_si256(a_ptr.add(i) as *const __m256i);
                let mut g0 = _mm256_load_si256(b_ptr.add(i) as *const __m256i);
                let mut f1 = _mm256_mul_epi32(f0, g0);
                        f0 = _mm256_srli_epi64(f0, 32);
                        g0 = _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(g0)));
                        f0 = _mm256_mul_epi32(f0, g0);
                let mut g1 = _mm256_mul_epi32(f1, qinv);
                        g0 = _mm256_mul_epi32(f0, qinv);
                        g1 = _mm256_mul_epi32(g1, q);
                        g0 = _mm256_mul_epi32(g0, q);
                        f1 = _mm256_sub_epi32(f1, g1);
                        f0 = _mm256_sub_epi32(f0, g0);
                        f1 = _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(f1)));
                        f0 = _mm256_blend_epi32(f1, f0, 0xaa);
                        _mm256_store_si256(r_ptr.add(i) as *mut __m256i, f0);
            }
        }
    }

    pub fn check_norm(a: &Poly, b: u32) -> bool {
        let qdata_ptr = QDATA.0.as_ptr();
        let a_ptr = a.coeffs.as_ptr();
        unsafe {
            let q = _mm256_load_si256(qdata_ptr.add(_8XQ) as *const __m256i);
            let hq = _mm256_srli_epi32(q, 1);
            let bound = _mm256_set1_epi32((b-1) as i32);
            let mut t = _mm256_setzero_si256();
            for i in 0..(N/8) {
                let mut f = _mm256_load_si256(a_ptr.add(8*i) as *const __m256i);
                let mut g = _mm256_cmpgt_epi32(f, hq);
                        g = _mm256_and_si256(g, q);
                        f = _mm256_sub_epi32(f,g);
                        g = _mm256_srai_epi32(f,31);
                        f = _mm256_xor_si256(f,g);
                        f = _mm256_sub_epi32(f,g);
                        f = _mm256_cmpgt_epi32(f,bound);
                        t = _mm256_or_si256(t,f);
            }
            _mm256_testz_si256(t, t) == 0
        }
    }
    
    pub fn sigma(a: &Poly, k: isize) -> Poly {
        let mut t = Poly::new();
        let mut j: usize = 0; 
        for i in 0..N {
            let mut x = a.coeffs[i];
            let branch_check = -((j & N) as i32) >> 31;
            let bitflip = x ^ -x;
            let new_idx = j & (N-1);
            x ^= branch_check & bitflip;
            t.coeffs[new_idx] = x;
            j += k as usize;
        }
        t
    } 
    
    pub fn sigma65_ntt (&mut self) {
        let coeffs_ptr = self.coeffs.as_mut_ptr();
        unsafe {
            for i in 0..(N/64) {
                let f0 = _mm256_load_si256(coeffs_ptr.add(32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(coeffs_ptr.add(32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(coeffs_ptr.add(32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(coeffs_ptr.add(32 * i + 24) as *const __m256i);

                _mm256_store_si256(coeffs_ptr.add(32 * i +  0) as *mut __m256i, f2);
                _mm256_store_si256(coeffs_ptr.add(32 * i +  8) as *mut __m256i, f3);
                _mm256_store_si256(coeffs_ptr.add(32 * i + 16) as *mut __m256i, f1);
                _mm256_store_si256(coeffs_ptr.add(32 * i + 24) as *mut __m256i, f0);

                let f0 = _mm256_load_si256(coeffs_ptr.add(N/2 + 32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(coeffs_ptr.add(N/2 + 32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(coeffs_ptr.add(N/2 + 32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(coeffs_ptr.add(N/2 + 32 * i + 24) as *const __m256i);
                
                _mm256_store_si256(coeffs_ptr.add(N/2 + 32 * i +  0) as *mut __m256i, f3);
                _mm256_store_si256(coeffs_ptr.add(N/2 + 32 * i +  8) as *mut __m256i, f2);
                _mm256_store_si256(coeffs_ptr.add(N/2 + 32 * i + 16) as *mut __m256i, f0);
                _mm256_store_si256(coeffs_ptr.add(N/2 + 32 * i + 24) as *mut __m256i, f1);
            }
        }
    } 
   
    pub fn sigma65_ntt_other(r: &mut Poly, a: &Poly) {
        let r_ptr = r.coeffs.as_mut_ptr();
        let a_ptr = a.coeffs.as_ptr();
        unsafe {
            for i in 0..(N/64) {
                let f0 = _mm256_load_si256(a_ptr.add(32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(a_ptr.add(32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(a_ptr.add(32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(a_ptr.add(32 * i + 24) as *const __m256i);

                _mm256_store_si256(r_ptr.add(32 * i +  0) as *mut __m256i, f2);
                _mm256_store_si256(r_ptr.add(32 * i +  8) as *mut __m256i, f3);
                _mm256_store_si256(r_ptr.add(32 * i + 16) as *mut __m256i, f1);
                _mm256_store_si256(r_ptr.add(32 * i + 24) as *mut __m256i, f0);

                let f0 = _mm256_load_si256(a_ptr.add(N/2 + 32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(a_ptr.add(N/2 + 32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(a_ptr.add(N/2 + 32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(a_ptr.add(N/2 + 32 * i + 24) as *const __m256i);
                
                _mm256_store_si256(r_ptr.add(N/2 + 32 * i +  0) as *mut __m256i, f3);
                _mm256_store_si256(r_ptr.add(N/2 + 32 * i +  8) as *mut __m256i, f2);
                _mm256_store_si256(r_ptr.add(N/2 + 32 * i + 16) as *mut __m256i, f0);
                _mm256_store_si256(r_ptr.add(N/2 + 32 * i + 24) as *mut __m256i, f1);
            }
        }
    }

    pub fn sigma129_ntt(&mut self) {
        let coeffs_ptr = self.coeffs.as_mut_ptr();
        unsafe {
            for i in 0..(N/32) {
                let f0 = _mm256_load_si256(coeffs_ptr.add(32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(coeffs_ptr.add(32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(coeffs_ptr.add(32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(coeffs_ptr.add(32 * i + 24) as *const __m256i);

                _mm256_store_si256(coeffs_ptr.add(32 * i +  0) as *mut __m256i, f1);
                _mm256_store_si256(coeffs_ptr.add(32 * i +  8) as *mut __m256i, f0);
                _mm256_store_si256(coeffs_ptr.add(32 * i + 16) as *mut __m256i, f3);
                _mm256_store_si256(coeffs_ptr.add(32 * i + 24) as *mut __m256i, f2);
            }
        }
    }

    pub fn sigma129_ntt_other(r: &mut Poly, a: &Poly) {
        let r_ptr = r.coeffs.as_mut_ptr();
        let a_ptr = a.coeffs.as_ptr();
        unsafe {
            for i in 0..(N/32) {
                let f0 = _mm256_load_si256(a_ptr.add(32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(a_ptr.add(32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(a_ptr.add(32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(a_ptr.add(32 * i + 24) as *const __m256i);

                _mm256_store_si256(r_ptr.add(32 * i +  0) as *mut __m256i, f1);
                _mm256_store_si256(r_ptr.add(32 * i +  8) as *mut __m256i, f0);
                _mm256_store_si256(r_ptr.add(32 * i + 16) as *mut __m256i, f3);
                _mm256_store_si256(r_ptr.add(32 * i + 24) as *mut __m256i, f2);
            }
        } 
    }
    
    pub fn sigma193_ntt(&mut self) {
        let coeffs_ptr = self.coeffs.as_mut_ptr();
        unsafe {
            for i in 0..(N/64) {
                let f0 = _mm256_load_si256(coeffs_ptr.add(32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(coeffs_ptr.add(32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(coeffs_ptr.add(32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(coeffs_ptr.add(32 * i + 24) as *const __m256i);

                _mm256_store_si256(coeffs_ptr.add(32 * i +  0) as *mut __m256i, f3);
                _mm256_store_si256(coeffs_ptr.add(32 * i +  8) as *mut __m256i, f2);
                _mm256_store_si256(coeffs_ptr.add(32 * i + 16) as *mut __m256i, f0);
                _mm256_store_si256(coeffs_ptr.add(32 * i + 24) as *mut __m256i, f1);

                let f0 = _mm256_load_si256(coeffs_ptr.add(N/2 + 32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(coeffs_ptr.add(N/2 + 32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(coeffs_ptr.add(N/2 + 32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(coeffs_ptr.add(N/2 + 32 * i + 24) as *const __m256i);
                
                _mm256_store_si256(coeffs_ptr.add(N/2 + 32 * i +  0) as *mut __m256i, f2);
                _mm256_store_si256(coeffs_ptr.add(N/2 + 32 * i +  8) as *mut __m256i, f3);
                _mm256_store_si256(coeffs_ptr.add(N/2 + 32 * i + 16) as *mut __m256i, f1);
                _mm256_store_si256(coeffs_ptr.add(N/2 + 32 * i + 24) as *mut __m256i, f0);
            }
        }
    }

    pub fn sigma193_ntt_other(r: &mut Poly, a: &Poly) {
        let r_ptr = r.coeffs.as_mut_ptr();
        let a_ptr = a.coeffs.as_ptr();
        unsafe {
            for i in 0..(N/64) {
                let f0 = _mm256_load_si256(a_ptr.add(32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(a_ptr.add(32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(a_ptr.add(32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(a_ptr.add(32 * i + 24) as *const __m256i);

                _mm256_store_si256(r_ptr.add(32 * i +  0) as *mut __m256i, f3);
                _mm256_store_si256(r_ptr.add(32 * i +  8) as *mut __m256i, f2);
                _mm256_store_si256(r_ptr.add(32 * i + 16) as *mut __m256i, f0);
                _mm256_store_si256(r_ptr.add(32 * i + 24) as *mut __m256i, f1);

                let f0 = _mm256_load_si256(a_ptr.add(N/2 + 32 * i +  0) as *const __m256i);
                let f1 = _mm256_load_si256(a_ptr.add(N/2 + 32 * i +  8) as *const __m256i);
                let f2 = _mm256_load_si256(a_ptr.add(N/2 + 32 * i + 16) as *const __m256i);
                let f3 = _mm256_load_si256(a_ptr.add(N/2 + 32 * i + 24) as *const __m256i);
                
                _mm256_store_si256(r_ptr.add(N/2 + 32 * i +  0) as *mut __m256i, f2);
                _mm256_store_si256(r_ptr.add(N/2 + 32 * i +  8) as *mut __m256i, f3);
                _mm256_store_si256(r_ptr.add(N/2 + 32 * i + 16) as *mut __m256i, f1);
                _mm256_store_si256(r_ptr.add(N/2 + 32 * i + 24) as *mut __m256i, f0);
            }
        }
    }
    
    pub fn trace65_ntt(&mut self) {
        let coeffs_ptr = self.coeffs.as_mut_ptr();
        unsafe {
            let mask = _mm256_set1_epi32((1 << 30) - 1);
            for i in 0..(N/32) {
                let mut f0 = _mm256_load_si256(coeffs_ptr.add(32 * i +  0) as *const __m256i);
                let mut f1 = _mm256_load_si256(coeffs_ptr.add(32 * i +  8) as *const __m256i);
                let mut f2 = _mm256_load_si256(coeffs_ptr.add(32 * i + 16) as *const __m256i);
                let mut f3 = _mm256_load_si256(coeffs_ptr.add(32 * i + 24) as *const __m256i);
                
                f0 = _mm256_add_epi32(f0,f1);
                f2 = _mm256_add_epi32(f2,f3);
                f1 = _mm256_srai_epi32(f0,30);
                f3 = _mm256_srai_epi32(f2,30);
                f0 = _mm256_and_si256(f0,mask);
                f2 = _mm256_and_si256(f2,mask);
                f0 = _mm256_sub_epi32(f0,f1);
                f2 = _mm256_sub_epi32(f2,f3);
                f1 = _mm256_slli_epi32(f1,18);
                f3 = _mm256_slli_epi32(f3,18);
                f0 = _mm256_add_epi32(f0,f1);
                f2 = _mm256_add_epi32(f2,f3);

                f0 = _mm256_add_epi32(f0,f2);
                f1 = _mm256_srai_epi32(f0,30);
                f0 = _mm256_and_si256(f0,mask);
                f0 = _mm256_sub_epi32(f0,f1);
                f1 = _mm256_slli_epi32(f1,18);
                f0 = _mm256_add_epi32(f0,f1);
                
                _mm256_store_si256(coeffs_ptr.add(32 * i +  0) as *mut __m256i, f0);
                _mm256_store_si256(coeffs_ptr.add(32 * i +  8) as *mut __m256i, f0);
                _mm256_store_si256(coeffs_ptr.add(32 * i + 16) as *mut __m256i, f0);
                _mm256_store_si256(coeffs_ptr.add(32 * i + 24) as *mut __m256i, f0);
            }
        }
    }

    pub fn trace65_ntt_other (r: &mut Poly, a: &Poly) {
        let r_ptr = r.coeffs.as_mut_ptr();
        let a_ptr = a.coeffs.as_ptr();
        unsafe {
            let mask = _mm256_set1_epi32((1 << 30) - 1);
            for i in 0..(N/32) {
                let mut f0 = _mm256_load_si256(a_ptr.add(32 * i +  0) as *const __m256i);
                let mut f1 = _mm256_load_si256(a_ptr.add(32 * i +  8) as *const __m256i);
                let mut f2 = _mm256_load_si256(a_ptr.add(32 * i + 16) as *const __m256i);
                let mut f3 = _mm256_load_si256(a_ptr.add(32 * i + 24) as *const __m256i);
                
                f0 = _mm256_add_epi32(f0,f1);
                f2 = _mm256_add_epi32(f2,f3);
                f1 = _mm256_srai_epi32(f0,30);
                f3 = _mm256_srai_epi32(f2,30);
                f0 = _mm256_and_si256(f0,mask);
                f2 = _mm256_and_si256(f2,mask);
                f0 = _mm256_sub_epi32(f0,f1);
                f2 = _mm256_sub_epi32(f2,f3);
                f1 = _mm256_slli_epi32(f1,18);
                f3 = _mm256_slli_epi32(f3,18);
                f0 = _mm256_add_epi32(f0,f1);
                f2 = _mm256_add_epi32(f2,f3);

                f0 = _mm256_add_epi32(f0,f2);
                f1 = _mm256_srai_epi32(f0,30);
                f0 = _mm256_and_si256(f0,mask);
                f0 = _mm256_sub_epi32(f0,f1);
                f1 = _mm256_slli_epi32(f1,18);
                f0 = _mm256_add_epi32(f0,f1);
                
                _mm256_store_si256(r_ptr.add(32 * i +  0) as *mut __m256i, f0);
                _mm256_store_si256(r_ptr.add(32 * i +  8) as *mut __m256i, f0);
                _mm256_store_si256(r_ptr.add(32 * i + 16) as *mut __m256i, f0);
                _mm256_store_si256(r_ptr.add(32 * i + 24) as *mut __m256i, f0);
            }
        }
    }

    pub fn power2round (a1: &mut Poly, a0: &mut Poly, a: &mut Poly) {
        a.reduce();
        a.freeze();
        power2round_avx(&mut a1.coeffs, &mut a0.coeffs, &a.coeffs);
    }

    pub fn decompose (a1: &mut Poly, a0: &mut Poly, a: &mut Poly) {
        a.reduce();
        a.freeze();
        decompose_avx(&mut a1.coeffs, &mut a0.coeffs, &a.coeffs);
    }

    pub fn makehint (h: &mut Poly, a1: &Poly, a0: &mut Poly) {
        a0.freeze();
        makehint_avx(&mut h.coeffs, &a1.coeffs, &a0.coeffs);
    }

    pub fn usehint (b1: &mut Poly, a: &mut Poly, h: &Poly) {
        a.reduce();
        a.freeze();
        usehint_avx(&mut b1.coeffs, &a.coeffs, &h.coeffs);
    }
}
