#![allow(dead_code)]

use crate::{
    comm::commitment::{Comm, CommKey, CommRnd}, 
    crypto::aes256::Aes256Ctx, 
    params::{M, N, R, SYMBYTES}, 
    poly_arith::{
        consts::{MONTSQ, NTTX, NTTX2, NTTX3, NTTX64}, 
        poly::Poly, 
        polyvec::{PolyVecL, PolyVecM}
    }
};

#[inline(always)]
fn poly_sigmainv_ntt(r: &mut Poly, a: &Poly, i: usize) {
    if i == 0 {
        *r = *a;
    } else if i == 1 {
        Poly::sigma193_ntt_other(r, a);
    } else if i == 2 {
        Poly::sigma129_ntt_other(r, a);
    } else if i == 3 {
        Poly::sigma65_ntt_other(r, a);
    }
}

#[inline(always)]
pub fn poly_shift(r: &mut Poly, a: &Poly, i: usize) {
    if i == 0 {
        *r = *a;
    } else if i == 1 {
        Poly::pointwise_montgomery_other(r, a, &NTTX);
    } else if i == 2 {
        Poly::pointwise_montgomery_other(r, a, &NTTX2);
    } else if i == 3 {
        Poly::pointwise_montgomery_other(r, a, &NTTX3);
    }
}

fn autobase_proof (r: &mut [Poly; R], a: &[PolyVecM; R], idx: usize) {
    let mut b: [Poly; R] = std::array::from_fn(|_| Poly::new());
    let mut rclone = r.clone();

    b[0] = a[0].vec[idx].clone();
    Poly::pointwise_montgomery_other(&mut b[1], &a[1].vec[idx], &NTTX);
    Poly::pointwise_montgomery_other(&mut b[2], &a[2].vec[idx], &NTTX2);
    Poly::pointwise_montgomery_other(&mut b[3], &a[3].vec[idx], &NTTX3);

    let btmp = b.clone();
    Poly::add_other(&mut rclone[0], &btmp[0], &btmp[2]);
    Poly::add_other(&mut rclone[1], &btmp[1], &btmp[3]);
    Poly::sub_other(&mut r[2], &rclone[0], &rclone[1]);
    Poly::add_other(&mut r[0], &rclone[0], &rclone[1]);

    Poly::sub_other(&mut b[1], &btmp[1], &btmp[3]);
    b[1].pointwise_montgomery(&NTTX64);
    Poly::sub_other(&mut b[0], &btmp[0], &btmp[2]);
    Poly::add_other(&mut r[1], &b[0], &b[1]);
    Poly::sub_other(&mut r[3], &b[0], &b[1]);

    r[1].sigma193_ntt();
    r[2].sigma129_ntt();
    r[3].sigma65_ntt();
}

fn autobase_verify (f: &mut [Poly; R]) {
    let mut b = [Poly::new(); R];
    let mut fclone = f.clone();
    b[0] = fclone[0].clone();
    Poly::pointwise_montgomery_other(&mut b[1], &fclone[1], &NTTX);
    Poly::pointwise_montgomery_other(&mut b[2], &fclone[2], &NTTX2);
    Poly::pointwise_montgomery_other(&mut b[3], &fclone[3], &NTTX3);

    let btmp = b.clone();
    Poly::add_other(&mut fclone[0], &btmp[0], &btmp[2]);
    Poly::add_other(&mut fclone[1], &btmp[1], &btmp[3]);
    Poly::sub_other(&mut f[2], &fclone[0], &fclone[1]);
    Poly::add_other(&mut f[0], &fclone[0], &fclone[1]);

    Poly::sub_other(&mut b[1], &btmp[1], &btmp[3]);
    b[1].pointwise_montgomery(&NTTX64);
    Poly::sub_other(&mut b[0], &btmp[0], &btmp[2]);
    Poly::add_other(&mut f[2], &b[0], &b[1]);
    Poly::sub_other(&mut f[3], &b[0], &b[1]);

    f[1].sigma193_ntt();
    f[2].sigma129_ntt();
    f[3].sigma65_ntt();
}


pub fn proof (
      msg: &mut PolyVecM, 
        g: &[PolyVecM; R],
    chash: &[u8; SYMBYTES]
) -> Poly {
    let mut v = Poly::new();
    let mut nonce = 0 as u64;
    let mut state = Aes256Ctx::init(chash, nonce);
    let mut a: [Poly; R] = std::array::from_fn(|_| Poly::new()); 
    let alpha: [Poly; 4] = std::array::from_fn(|_| {
        state.select(nonce);
        nonce += 1;
        let mut a = Poly::new();
        Poly::uniform_preinit(&mut a, &mut state);
        a
    });
    let beta: [Poly; R] = std::array::from_fn(|_| {
        state.select(nonce);
        nonce += 1;
        let mut a = Poly::new();
        Poly::uniform_preinit(&mut a, &mut state);
        a
    });

    msg.vec[M-2] = Poly::new();
    let mut tmp = Poly::new();
    
    for i in 0..4 {
        autobase_proof(&mut a, g, i);
        for j in 0..R {
            Poly::pointwise_montgomery_other(&mut tmp, &a[j], &a[j]);
            tmp.pointwise_montgomery(&alpha[i]);
            tmp.pointwise_montgomery(&beta[j]);
            v.add(&tmp);
        }

        let mut mprime = Poly::new();
        for j in 0..N {
            mprime.coeffs[j] = 1 - (msg.vec[i].coeffs[j] << 1);
        }
        for j in 0..R {
            poly_sigmainv_ntt(&mut tmp, &mprime, j);

            for k in 0..N {
                tmp.coeffs[k] *= a[j].coeffs[k];
            }
            tmp.pointwise_montgomery(&alpha[i]);
            tmp.pointwise_montgomery(&beta[j]);
            msg.vec[M-2].add(&tmp);
        }
    }
    v.scale_montgomery(MONTSQ as i32);
    for j in 0..R {
        poly_shift(&mut tmp, &g[j].vec[M-2], j);
        v.add(&tmp);
    }
    v.freeze();
    v
}

pub fn verify (
        v: &mut Poly,
    chash: &[u8; SYMBYTES],
        c: &[Poly; R],
        z: &[CommRnd; R],
       tp: &Comm,
      ckp: &CommKey,
) -> bool {
    let mut nonce = 0;
    let mut state = Aes256Ctx::init(&chash, nonce);
    let alpha: [Poly; 4] = std::array::from_fn(|_| {
        state.select(nonce);
        nonce += 1;
        let mut a = Poly::new();
        Poly::uniform_preinit(&mut a, &mut state);
        a
    });
    let beta: [Poly; R] = std::array::from_fn(|_| {
        state.select(nonce);
        nonce += 1;
        let mut a = Poly::new();
        Poly::uniform_preinit(&mut a, &mut state);
        a
    });
    let mut zshat: [PolyVecL; R] = std::array::from_fn(|i| { z[i].s  });
    let mut zmhat: [PolyVecM; R] = std::array::from_fn(|i| { z[i].em });
    let mut chat : [Poly; R]     = std::array::from_fn(|i| { c[i] });
    for j in 0..R {
        zshat[j].vec_ntt();
        zmhat[j].vec_ntt();
        chat[j].ntt();
    }

    let mut tmp = Poly::new();
    let mut cfull = Poly::new();
    Poly::pointwise_montgomery_other(&mut tmp, &chat[1], &NTTX);
    Poly::add_other(&mut cfull, &chat[0], &tmp);
    Poly::pointwise_montgomery_other(&mut tmp, &chat[2], &NTTX2);
    cfull.add(&tmp);
    Poly::pointwise_montgomery_other(&mut tmp, &chat[3], &NTTX3);
    cfull.add(&tmp);

    for i in 0..4 {
        let mut f: [Poly; R] = std::array::from_fn(|j| {
            let mut a = PolyVecL::pointwise_acc_montgomery(&ckp.bm[i], &zshat[j]);
            let mut tmp1 = Poly::new();
            Poly::pointwise_montgomery_other(&mut tmp1, &chat[j], &tp.tm.vec[i]);
            a.sub(&tmp1);
            a.scale_montgomery(MONTSQ as i32);
            a.add(&zmhat[j].vec[i]);
            a
        });
        autobase_verify(&mut f);
        for j in 0..R {
            Poly::add_other(&mut tmp, &f[j], &cfull);
            f[j].pointwise_montgomery(&tmp);
            f[j].pointwise_montgomery(&alpha[i]);
            f[j].pointwise_montgomery(&beta[j]);
            v.add(&f[j]);
        }
    }
    v.scale_montgomery(MONTSQ as i32);

    for j in 0..R {
        tmp = PolyVecL::pointwise_acc_montgomery(&ckp.bm[M-2], &zshat[j]);
        tmp.scale_montgomery(MONTSQ as i32);
        tmp.add(&zmhat[j].vec[M-2]);
        let tmp1 = tmp.clone();
        poly_shift(&mut tmp, &tmp1, j);
        v.add(&tmp);
    }
    Poly::pointwise_montgomery_other(&mut tmp, &cfull, &tp.tm.vec[M-2]);
    tmp.scale_montgomery(MONTSQ as i32);
    v.sub(&tmp);
    v.freeze();

    false
}
