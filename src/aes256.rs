#![cfg(target_arch = "x86_64")]

use std::arch::x86_64::*;

#[repr(C, align(16))]
#[allow(dead_code)]
pub struct Aes256Ctx {
    pub rkeys: [__m128i; 16],
    pub n: __m128i,
}

#[allow(dead_code)]
pub const AES256CTR_BLOCKBYTES: u8 = 64;

#[allow(dead_code)]
impl Aes256Ctx {

    #[target_feature(enable = "aes")]
    #[target_feature(enable = "sse2")]
    #[target_feature(enable = "ssse3")]
    pub unsafe fn init(key: &[u8; 32], nonce: u64) -> Self {
        let mut rkeys = [_mm_setzero_si128(); 16];
        let mut idx = 0;
        let n: __m128i = _mm_loadl_epi64(&nonce as *const u64 as *const __m128i);

        let key0 = _mm_loadu_si128(key.as_ptr() as *const __m128i);
        let key1 = _mm_loadu_si128(key.as_ptr().add(16) as *const __m128i);

        rkeys[idx] = key0;
        idx += 1;
        
        let mut temp0 = key0;
        let mut temp2 = key1;
        let mut temp4 = _mm_setzero_si128();

        macro_rules! BLOCK1 {
            ($imm: expr) => {
                let mut temp1 = _mm_aeskeygenassist_si128(temp2, $imm);
                rkeys[idx] = temp2;
                idx += 1;
                    
                temp4 = _mm_castps_si128(_mm_shuffle_ps(
                        _mm_castsi128_ps(temp4), _mm_castsi128_ps(temp0), 0x10
                ));
                temp0 = _mm_xor_si128(temp0, temp4);

                temp4 = _mm_castps_si128(_mm_shuffle_ps(
                        _mm_castsi128_ps(temp4), _mm_castsi128_ps(temp0), 0x8c
                ));
                temp0 = _mm_xor_si128(temp0, temp4);
                    
                temp1 = _mm_castps_si128(_mm_shuffle_ps(
                        _mm_castsi128_ps(temp1), _mm_castsi128_ps(temp1), 0xff
                ));
                temp0 = _mm_xor_si128(temp0, temp1);
            };
        }

        macro_rules! BLOCK2 {
            ($imm: expr) => {
                let mut temp1 = _mm_aeskeygenassist_si128(temp0, $imm);
                rkeys[idx] = temp0;
                idx += 1;

                temp4 = _mm_castps_si128(_mm_shuffle_ps(
                        _mm_castsi128_ps(temp4), _mm_castsi128_ps(temp2), 0x10
                ));
                temp2 = _mm_xor_si128(temp2, temp4);

                temp4 = _mm_castps_si128(_mm_shuffle_ps(
                        _mm_castsi128_ps(temp4), _mm_castsi128_ps(temp2), 0x8c
                ));
                temp2 = _mm_xor_si128(temp2, temp4);
                    
                temp1 = _mm_castps_si128(_mm_shuffle_ps(
                        _mm_castsi128_ps(temp1), _mm_castsi128_ps(temp1), 0xaa
                ));
                temp2 = _mm_xor_si128(temp2, temp1);
            };
        }

        BLOCK1!(0x01);
        BLOCK2!(0x01);
            
        BLOCK1!(0x02);
        BLOCK2!(0x02);

        BLOCK1!(0x04);
        BLOCK2!(0x04);

        BLOCK1!(0x08);
        BLOCK2!(0x08);
            
        BLOCK1!(0x10);
        BLOCK2!(0x10);
            
        BLOCK1!(0x20);
        BLOCK2!(0x20);

        BLOCK1!(0x40);

        rkeys[idx] = temp0;

        Aes256Ctx { rkeys, n }
    } 


    #[target_feature(enable = "aes")]
    #[target_feature(enable = "sse2")]
    #[target_feature(enable = "ssse3")]
    #[allow(dead_code)]
    pub unsafe fn prf (out: &mut [u8], seed: &[u8; 32], nonce: u64) {
        let mut ctx = Aes256Ctx::init(seed, nonce);
        let chunk_size = 64;
        let mut chunks = out.chunks_exact_mut(chunk_size);
        for chunk in &mut chunks {
            ctx.encrypt4(chunk.try_into().unwrap());
        }
        let final_chunk = chunks.into_remainder();
        if !final_chunk.is_empty() {
            let mut buf = [0u8; 64];
            ctx.encrypt4(&mut buf);
            final_chunk.copy_from_slice(&buf[..final_chunk.len()]);
        }
    } 
    
    pub fn select(&mut self, nonce: u64) {
        unsafe {self.n = _mm_loadl_epi64(&nonce as *const u64 as *const __m128i);}
    }

    #[target_feature(enable = "aes")]
    #[target_feature(enable = "sse2")]
    #[target_feature(enable = "ssse3")]
    unsafe fn encrypt4 (&mut self, out: &mut [u8; 64]) {
        /* Load current counter value */
        let f = _mm_load_si128(&self.n);

        /* Increase counter in 4 consecutive blocks */
        let mut t = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
        let mut f0 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(0, 0)), t);
        let mut f1 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(1, 0)), t);
        let mut f2 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(2, 0)), t);
        let mut f3 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(3, 0)), t);

        /* Write counter for next iteration, increased by 4 */
        let inc_n = _mm_add_epi64(f, _mm_set_epi64x(4, 0));
        _mm_store_si128(&mut self.n, inc_n); 

        /* Actual AES encryption, 4x interleaved */
        t  = _mm_load_si128(&self.rkeys[0]);
        f0 = _mm_xor_si128(f0, t);
        f1 = _mm_xor_si128(f1, t);
        f2 = _mm_xor_si128(f2, t);
        f3 = _mm_xor_si128(f3, t);

        for i in 1..14 {
            t = _mm_load_si128(&self.rkeys[i]);
            f0 = _mm_aesenc_si128(f0, t);
            f1 = _mm_aesenc_si128(f1, t);
            f2 = _mm_aesenc_si128(f2, t);
            f3 = _mm_aesenc_si128(f3, t);
        }
        
        t = _mm_load_si128(&self.rkeys[14]);
        f0 = _mm_aesenclast_si128(f0, t);
        f1 = _mm_aesenclast_si128(f1, t);
        f2 = _mm_aesenclast_si128(f2, t);
        f3 = _mm_aesenclast_si128(f3, t);

        /* Write results */
        _mm_storeu_si128(out.as_mut_ptr().add( 0) as *mut __m128i, f0);
        _mm_storeu_si128(out.as_mut_ptr().add(16) as *mut __m128i, f1);
        _mm_storeu_si128(out.as_mut_ptr().add(32) as *mut __m128i, f2);
        _mm_storeu_si128(out.as_mut_ptr().add(48) as *mut __m128i, f3);

    }

    #[target_feature(enable = "aes")]
    #[target_feature(enable = "sse2")]
    #[target_feature(enable = "ssse3")]
    pub unsafe fn squeezeblocks(&mut self, out: &mut [u8], nblocks: usize) {
        assert!(out.len() >= 64*nblocks);
        for chunk in out.chunks_exact_mut(64) {
            unsafe {
                self.encrypt4(chunk.try_into().unwrap());
            }
        }
    }
}

