use std::collections::HashMap;

use irelzk_rs::crypto::shake::Shake128;

#[test]
fn test_shake() {
    let mut shake_handle = Shake128::init();
    let input1 = [1u8; 10];
    let input2 = [1u8; 10];
    let mut out1 = [0u8; 10];
    let mut out2 = [0u8; 10];
    shake_handle.absorb(&input1);
    shake_handle.finalize();
    shake_handle.squeeze(&mut out1);
    
    shake_handle = Shake128::init();
    shake_handle.absorb(&input2);
    shake_handle.finalize();
    shake_handle.squeeze(&mut out2);
    assert_eq!(out1, out2);
}

#[test]
fn test_shake_againt_known_vectors() {
    let mut map: HashMap<&'static str, &'static str>  = HashMap::new();
    map.insert("", "7f9c2ba4e88f827d616045507605853e");
    map.insert("0e", "fa996dafaa208d72287c23bc4ed4bfd5");
    map.insert("d9e8", "c7211512340734235bb8d3c4651495aa");
    map.insert("1b3b6e", "d7335497e4cd3666885edbb0824d7a75");

    for (k,v) in map {
        let input = hex::decode(k).unwrap();
        let expected_output:[u8; 16] = hex::decode(v)
            .unwrap()
            .try_into()
            .expect("Expected output must be 16 bytes");
        let mut out = [0u8; 16];
        let mut shake_handle = Shake128::init();
        shake_handle.absorb(&input);
        shake_handle.finalize();
        shake_handle.squeeze(&mut out);
        assert_eq!(out, expected_output, "Failed for ({}, {}) -- Output: {}", k, v, hex::encode(out));
    }
}
