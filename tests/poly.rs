#![allow(dead_code)]
use irelzk_rs::{params::{N, SYMBYTES}, poly::Poly};
use rand::{rngs::OsRng, RngCore};

fn bitrev7(a: u8) -> u8 {
    let mut t = 
         (a &  1) << 6;
    t |= (a &  2) << 4;
    t |= (a &  4) << 2;
    t |= (a &  8) << 0;
    t |= (a & 16) >> 2;
    t |= (a & 32) >> 4;
    t |= (a & 64) >> 6;
    t
}

fn idx (i: isize) -> isize {
    let mut r = i/32;
    let i_mod = i%32;
    r *= 32;
    r += 8*(i_mod % 4);
    r += i_mod/4;
    r
}

fn idxinv(i: isize) -> isize {
    let mut r = i/32;
    let i_mod = i%32;
    r *= 32;
    r += 4 * (i_mod % 8);
    r += i_mod / 8;
    r
}

fn lut (i: isize, k: isize) -> isize {
    let mut r = bitrev7(idxinv(i) as u8) as isize;
    r = (r*k + k/2) % 128;
    r = idx(bitrev7(r as u8) as isize);
    r
}

#[test]
fn test_uniform_random() {
    let mut seed = [0u8; SYMBYTES];
    OsRng.fill_bytes(&mut seed);
    let mut a = Poly::new();
    Poly::uniform_random(&mut a, &seed, 0);
    let mut b = a.clone();
    b.freeze();
    b.sub(&a);
    b.freeze();
    for i in 0..N {
        assert_eq!(b.coeffs[i], 0, "Failing at index {}", i);
    }
}

#[test] 
fn test_trinary_output() {
    let mut seed = [0u8; SYMBYTES];
    OsRng.fill_bytes(&mut seed);
    let mut a = Poly::new();
    Poly::trinary(&mut a, &seed, 0);
    for i in 0..N {
        let x = a.coeffs[i];
        assert!(x == 0 || x == 1 || x == -1, "Failing at index {}", i);
    }
}

#[test]
fn test_gamma_output() {
    let mut seed = [0u8; SYMBYTES];
    OsRng.fill_bytes(&mut seed);
    let mut a = Poly::new();
    Poly::uniform_gamma(&mut a, &seed, 0);
    assert!(Poly::check_norm(&a, 1 << 17));
}

#[test]
fn test_sigma() {
    let mut seed = [0u8; SYMBYTES];
    OsRng.fill_bytes(&mut seed);
    let mut f = Poly::new();
    Poly::uniform_random(&mut f, &seed, 0);
    let mut g = Poly::sigma(&f, 65);
    f.ntt(); g.ntt();
    f.freeze(); g.freeze();
    let mut h = Poly::new();
    for i in 0..N {
        let idx = lut(i as isize,65) as usize;
        h.coeffs[i] = f.coeffs[idx];
    }
    Poly::sub_other(&mut f, &g, &h);
    for i in 0..N {
        let index = idx(i as isize) as usize;
        assert_eq!(f.coeffs[index], 0, "Failing at index {}", i);
    }
}

#[test]
fn test_trace65() {
    let mut seed = [0u8; SYMBYTES];
    OsRng.fill_bytes(&mut seed);
    let mut f = Poly::new();
    Poly::uniform_random(&mut f, &seed, 0);
    let mut h = Poly::sigma(&f, 65);
    let mut g = Poly::new();
    Poly::add_other(&mut g, &f, &h);
    h = Poly::sigma(&f, 129);
    g.add(&h);
    h = Poly::sigma(&f, 193);
    g.add(&h);
    f.ntt(); g.ntt();
    f.trace65_ntt();
    f.sub(&g);
    f.freeze();
    for i in 0..N {
        let index = idx(i as isize) as usize;
        assert_eq!(f.coeffs[index], 0, "Failing at index {}", i);
    }
}
