#![allow(dead_code)]
#![allow(unused_imports)]

use irelzk_rs::{params::{N, Q, SYMBYTES}, poly::Poly};
use rand::{rngs::OsRng, RngCore};

fn pow (base: i32, exp: usize) -> i32 {
    if exp == 0 {
        return 1;
    } else if exp == 1 {
        return base; 
    }

    let mut  r = pow(base, exp / 2);
    r = ((r as i64) * (r as i64) % (Q as i64)) as i32;
    
    if exp % 2 != 0 {
        r = ((r as i64) * (base as i64) % (Q as i64)) as i32
    }
    r
}

fn idx (i: usize) -> usize {
    let mut r = i/32;
    let i_mod = i % 32;
    r *= 32;
    r += 8*(i_mod % 4);
    r += i_mod/4;
    r
} 

#[test]
fn test_ntt_x() {
    let mut p = Poly::new();
    p.coeffs[1] = 1;
    let mut zeta = [0i32; N];

    p.ntt();
    for i in 0..N {
        zeta[i] = p.coeffs[idx(i)] % Q;
        let zeta_n = pow(zeta[i], N) + 1;
        assert_eq!(zeta_n % Q, 0);
        for j in 0..i {
            let diff = zeta[i] - zeta[j] % Q;
            assert_ne!(diff, 0, "Failing for ({},{})", i, j);
        }
    } 
}

#[test]
fn test_ntt_monomial() {
    let mut p = Poly::new();
    p.coeffs[1] = 1;
    let mut zeta = [0i32; N];
    let mut zetapow = [[0i32; N]; N];
    p.ntt();
    for i in 0..N {
        zeta[i] = p.coeffs[idx(i)] % Q; 
    }
    
    for i in 0..N {
        for j in 0..N {
            zetapow[i][j] = pow(zeta[i], j);
        }
    }

    for i in 0..N {
        let mut b = Poly::new();
        b.coeffs[i] = 1;
        b.ntt();
        for j in 0..N {
            b.coeffs[idx(j)] %= Q;
            let diff = (b.coeffs[idx(j)] - zetapow[j][i]) % Q;
            assert_eq!(diff, 0, "Failing for ({}, {})", i, j);
        }
    }
}

#[test] 
fn test_invntt() {
    let mut seed = [0u8; SYMBYTES];
    OsRng.fill_bytes(&mut seed);
    let mut a = Poly::uniform_random(&seed, 0);
    let b = a.clone();
    a.ntt();
    a.inverse_ntt();
    for i in 0..N {
        let diff = (a.coeffs[i] - b.coeffs[i]) % Q;
        assert_eq!(diff, 0, "Failed at {}", i);
    }
}
