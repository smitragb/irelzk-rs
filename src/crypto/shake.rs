#![allow(dead_code)]
use crate::crypto::keccak::KeccakState;

pub const SHAKE128_RATE: usize = 168;
pub struct Shake128 {
    state: KeccakState,
}

impl Shake128 {
    pub fn init() -> Self {
        Self {
            state: KeccakState::init(),
        }
    }

    pub fn absorb(&mut self, input: &[u8]) {
        self.state.pos = KeccakState::absorb(
            &mut self.state.s, 
            SHAKE128_RATE,
            self.state.pos,
            input
        );
    }

    pub fn finalize(&mut self) {
        KeccakState::finalize(
            &mut self.state.s, 
            SHAKE128_RATE, 
            self.state.pos,
            0x1F
        );
        self.state.pos = 0;
    }

    pub fn squeezeblocks(&mut self, out: &mut [u8], nblocks: usize) {
        KeccakState::squeezeblocks(
            out, 
            nblocks, 
            &mut self.state.s, 
            SHAKE128_RATE
        );
    }

    pub fn squeeze(&mut self, out: &mut [u8]) {
        self.state.pos = KeccakState::squeeze(
            out, 
            &mut self.state.s, 
            SHAKE128_RATE, 
            self.state.pos
        );
    }
}
