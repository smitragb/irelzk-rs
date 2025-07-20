use irelzk_rs::crypto::aes256::Aes256Ctx;

#[test]
fn test_squeezeblocks() {
    if !is_x86_feature_detected!("aes") {
        eprintln!("skipping test: aes-ni not available on this cpu");
        return;
    }

    let key = [0u8; 32];
    let nonce = 43;
    let mut ctx1 = unsafe { Aes256Ctx::init(&key, nonce) };
    let mut ctx2 = unsafe { Aes256Ctx::init(&key, nonce) };
    
    let mut out1 = vec![0u8; 256];
    let mut out2 = vec![0u8; 256];
    let nblocks = 256/64;

    unsafe {
        ctx1.squeezeblocks(&mut out1, nblocks);
        ctx2.squeezeblocks(&mut out2, nblocks);
    }

    assert_eq!(out1, out2);
}

#[test]
fn test_prf() {
    if !is_x86_feature_detected!("aes") {
        eprintln!("skipping test: aes-ni not available on this cpu");
        return;
    }

    let key = [0u8; 32];
    let nonce = 1231;
    
    let mut out1 = vec![0u8; 256];
    let mut out2 = vec![0u8; 256];

    unsafe {
        Aes256Ctx::prf(&mut out1, &key, nonce);
        Aes256Ctx::prf(&mut out2, &key, nonce);
    }

    assert_eq!(out1, out2);
}
