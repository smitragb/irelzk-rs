#![allow(dead_code)]
use std::arch::x86_64::*;

use crate::params::{D, GAMMA2, N};

pub fn power2round_avx(a1: &mut [i32; N], a0: &mut [i32; N], a: &[i32; N]) {
    unsafe {
        let mask = _mm256_set1_epi32(-(1 << D));
        let half = _mm256_set1_epi32((1 << (D-1)) - 1);

        let a_ptr = a.as_ptr();
        let a0_ptr = a0.as_mut_ptr();
        let a1_ptr = a1.as_mut_ptr();

        for i in 0..(N/8) {
            let f = _mm256_load_si256(a_ptr.add(8*i) as *const __m256i);
            let mut f1 = _mm256_add_epi32(f, half);
            let mut f0 = _mm256_and_si256(f1, mask);
            f1 = _mm256_srai_epi32(f1, D as i32);
            f0 = _mm256_sub_epi32(f, f0);
            _mm256_store_si256(a1_ptr.add(8*i) as *mut __m256i, f1);
            _mm256_store_si256(a0_ptr.add(8*i) as *mut __m256i, f0);
        }
    }
}

pub fn power2round_avx_self (a1: &mut [i32; N], a0: &mut [i32; N]) {
    unsafe {
        let mask = _mm256_set1_epi32(-(1 << D));
        let half = _mm256_set1_epi32((1 << (D-1)) - 1);

        let a0_ptr = a0.as_mut_ptr();
        let a1_ptr = a1.as_mut_ptr();

        for i in 0..(N/8) {
            let f = _mm256_load_si256(a1_ptr.add(8*i) as *const __m256i);
            let mut f1 = _mm256_add_epi32(f, half);
            let mut f0 = _mm256_and_si256(f1, mask);
            f1 = _mm256_srai_epi32(f1, D as i32);
            f0 = _mm256_sub_epi32(f, f0);
            _mm256_store_si256(a1_ptr.add(8*i) as *mut __m256i, f1);
            _mm256_store_si256(a0_ptr.add(8*i) as *mut __m256i, f0);
        }
    }
}

pub fn decompose_avx(a1: &mut [i32; N], a0: &mut [i32; N], a: &[i32; N]) {
    unsafe {
        let a_ptr = a.as_ptr();
        let a0_ptr = a0.as_mut_ptr();
        let a1_ptr = a1.as_mut_ptr();

        let half = _mm256_set1_epi32(1 << 17);
        let mask = _mm256_set1_epi32(-(1<<18));
        let max = _mm256_set1_epi32(2048);

        for i in 0..(N/8) {
            let f = _mm256_load_si256(a_ptr.add(8*i) as *const __m256i);
            let mut t0 = _mm256_srai_epi32(f, 12);
            let mut t1 = _mm256_srai_epi32(f, 24);
            t0 = _mm256_add_epi32(t0, f);
            t1 = _mm256_add_epi32(t1, half);
            t0 = _mm256_add_epi32(t0, t1);
            let mut f1 = _mm256_srai_epi32(t0, 18);

            t0 = _mm256_and_si256(t0, mask);
            t1 = _mm256_slli_epi32(f1, 6);
            let mut f0 = _mm256_sub_epi32(f, t0);
            f0 = _mm256_add_epi32(f0, t1);

            t0 = _mm256_cmpeq_epi32(f1, max);
            f1 = _mm256_xor_si256(f1, t0);
            f1 = _mm256_sub_epi32(f1, t0);
            f0 = _mm256_add_epi32(f0, t0);

            _mm256_store_si256(a1_ptr.add(8*i) as *mut __m256i, f1);
            _mm256_store_si256(a0_ptr.add(8*i) as *mut __m256i, f0);
        }
    }
}

pub fn decompose_avx_self(a1: &mut [i32; N], a0: &mut [i32; N]) {
    unsafe {
        let a0_ptr = a0.as_mut_ptr();
        let a1_ptr = a1.as_mut_ptr();

        let half = _mm256_set1_epi32(1 << 17);
        let mask = _mm256_set1_epi32(-(1<<18));
        let max = _mm256_set1_epi32(2048);

        for i in 0..(N/8) {
            let f = _mm256_load_si256(a1_ptr.add(8*i) as *const __m256i);
            let mut t0 = _mm256_srai_epi32(f, 12);
            let mut t1 = _mm256_srai_epi32(f, 24);
            t0 = _mm256_add_epi32(t0, f);
            t1 = _mm256_add_epi32(t1, half);
            t0 = _mm256_add_epi32(t0, t1);
            let mut f1 = _mm256_srai_epi32(t0, 18);

            t0 = _mm256_and_si256(t0, mask);
            t1 = _mm256_slli_epi32(f1, 6);
            let mut f0 = _mm256_sub_epi32(f, t0);
            f0 = _mm256_add_epi32(f0, t1);

            t0 = _mm256_cmpeq_epi32(f1, max);
            f1 = _mm256_xor_si256(f1, t0);
            f1 = _mm256_sub_epi32(f1, t0);
            f0 = _mm256_add_epi32(f0, t0);

            _mm256_store_si256(a1_ptr.add(8*i) as *mut __m256i, f1);
            _mm256_store_si256(a0_ptr.add(8*i) as *mut __m256i, f0);
        }
    }
}

pub fn makehint_avx(h: &mut [i32; N], a1: &[i32; N], a0: &[i32; N]) {
    unsafe {
        let h_ptr = h.as_mut_ptr();
        let a1_ptr = a1.as_ptr();
        let a0_ptr = a0.as_ptr();

        let blo = _mm256_set1_epi32(-GAMMA2);
        let bhi = _mm256_set1_epi32(GAMMA2);
        let min = _mm256_set1_epi32(-2048);

        for i in 0..(N/8) {
            let f0 = _mm256_load_si256(a0_ptr.add(8*i) as *const __m256i);
            let f1 = _mm256_load_si256(a1_ptr.add(8*i) as *const __m256i);

            let mut g0 = _mm256_cmpgt_epi32(blo, f0);
            let mut g1 = _mm256_cmpgt_epi32(f0, bhi);
            g0 = _mm256_or_si256(g0, g1);
            g1 = _mm256_cmpeq_epi32(blo, f0);
            let g2 = _mm256_cmpeq_epi32(f1, min);
            g1 = _mm256_andnot_si256(g2, g1);
            g0 = _mm256_or_si256(g0, g1);
            g0 = _mm256_sign_epi32(g0, g0);
            _mm256_store_si256(h_ptr.add(8*i) as *mut __m256i, g0);
        }
    }
}

#[repr(align(32))]
struct Aligned([i32; N]);

pub fn usehint_avx(b1: &mut [i32; N], a: &[i32; N], hint: &[i32; N]) {
    let mut a0 = Aligned([0; N]);
    unsafe {
        let off = _mm256_set1_epi32(2048);
        let mask = _mm256_set1_epi32(4095);
        decompose_avx(b1, &mut a0.0, a);

        let a0_ptr = a0.0.as_ptr();
        let b1_ptr = b1.as_mut_ptr();
        let hint_ptr = hint.as_ptr();

        for i in 0..(N/8) {
            let f = _mm256_load_si256(a0_ptr.add(8*i) as *const __m256i);
            let mut g = _mm256_load_si256(b1_ptr.add(8*i) as *const __m256i);
            let mut h = _mm256_load_si256(hint_ptr.add(8*i) as *const __m256i);
            h = _mm256_sign_epi32(h, f);
            g = _mm256_add_epi32(g, h);
            g = _mm256_add_epi32(g, off);
            g = _mm256_and_si256(g, mask);
            g = _mm256_sub_epi32(g, off);
            _mm256_store_si256(b1_ptr.add(8*i) as *mut __m256i, g);
        }
    }
}

pub fn usehint_avx_self(b1: &mut [i32; N], hint: &[i32; N]) {
    let mut a0 = Aligned([0; N]);
    unsafe {
        let off = _mm256_set1_epi32(2048);
        let mask = _mm256_set1_epi32(4095);
        decompose_avx_self(b1, &mut a0.0);

        let a0_ptr = a0.0.as_ptr();
        let b1_ptr = b1.as_mut_ptr();
        let hint_ptr = hint.as_ptr();

        for i in 0..(N/8) {
            let f = _mm256_load_si256(a0_ptr.add(8*i) as *const __m256i);
            let mut g = _mm256_load_si256(b1_ptr.add(8*i) as *const __m256i);
            let mut h = _mm256_load_si256(hint_ptr.add(8*i) as *const __m256i);
            h = _mm256_sign_epi32(h, f);
            g = _mm256_add_epi32(g, h);
            g = _mm256_add_epi32(g, off);
            g = _mm256_and_si256(g, mask);
            g = _mm256_sub_epi32(g, off);
            _mm256_store_si256(b1_ptr.add(8*i) as *mut __m256i, g);
        }
    }

}
