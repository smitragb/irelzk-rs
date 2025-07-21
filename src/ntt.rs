extern "C" {
    pub fn ntt_avx(r: *mut i32, qdata: *const i32);
    pub fn invntt_avx(r: *mut i32, qdata: *const i32);
}

pub fn forward_ntt(r: &mut [i32], qdata: [i32]) {
    unsafe {
        ntt_avx(r.as_mut_ptr(), qdata.as_ptr());
    }
}

pub fn inverse_ntt(r: &mut [i32], qdata: [i32]) {
    unsafe {
        invntt_avx(r.as_mut_ptr(), qdata.as_ptr());
    }
}
