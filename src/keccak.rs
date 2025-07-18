#![allow(dead_code)]
#![allow(non_snake_case)]
pub struct KeccakState {
    pub s: [u64; 25],
    pub pos: usize,
}

const NROUNDS: usize = 24;
pub const KECCAK_F_ROUND_CONSTANTS: [u64; NROUNDS] = [
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808a,
    0x8000000080008000,
    0x000000000000808b,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008a,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000a,
    0x000000008000808b,
    0x800000000000008b,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800a,
    0x800000008000000a,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008,
];

impl KeccakState {
    fn load64 (x: &[u8; 8]) -> u64 {
        u64::from_le_bytes(*x)
    }

    fn store64 (x: u64) -> [u8; 8] {
        x.to_le_bytes()
    }

    pub fn init() -> Self {
        Self {
            s: [0u64; 25],
            pos: 0,
        }
    }

    fn keccak_f1600_state_permute(state: &mut [u64; 25]) { 
        let mut Aba = state[ 0]; let mut Abe = state[ 1]; let mut Abi = state[ 2]; let mut Abo = state[ 3]; let mut Abu = state[ 4];
        let mut Aga = state[ 5]; let mut Age = state[ 6]; let mut Agi = state[ 7]; let mut Ago = state[ 8]; let mut Agu = state[ 9];
        let mut Aka = state[10]; let mut Ake = state[11]; let mut Aki = state[12]; let mut Ako = state[13]; let mut Aku = state[14];
        let mut Ama = state[15]; let mut Ame = state[16]; let mut Ami = state[17]; let mut Amo = state[18]; let mut Amu = state[19];
        let mut Asa = state[20]; let mut Ase = state[21]; let mut Asi = state[22]; let mut Aso = state[23]; let mut Asu = state[24];

        for round in (0..NROUNDS).step_by(2) 
        {
            let mut BCa = Aba ^ Aga ^ Aka ^ Ama ^ Asa;
            let mut BCe = Abe ^ Age ^ Ake ^ Ame ^ Ase;
            let mut BCi = Abi ^ Agi ^ Aki ^ Ami ^ Asi;
            let mut BCo = Abo ^ Ago ^ Ako ^ Amo ^ Aso;
            let mut BCu = Abu ^ Agu ^ Aku ^ Amu ^ Asu;

            let mut Da = BCu ^ BCe.rotate_left(1);
            let mut De = BCa ^ BCi.rotate_left(1);
            let mut Di = BCe ^ BCo.rotate_left(1);
            let mut Do = BCi ^ BCu.rotate_left(1);
            let mut Du = BCo ^ BCa.rotate_left(1);

            Aba ^= Da;
            BCa = Aba;
            Age ^= De;
            BCe = Age.rotate_left(44);
            Aki ^= Di;
            BCi = Aki.rotate_left(43);
            Amo ^= Do;
            BCo = Amo.rotate_left(21);
            Asu ^= Du;
            BCu = Asu.rotate_left(14);

            let mut Eba = BCa ^ ((!BCe) & BCi);
            Eba ^= KECCAK_F_ROUND_CONSTANTS[round];
            let mut Ebe = BCe ^ ((!BCi) & BCo);
            let mut Ebi = BCi ^ ((!BCo) & BCu); 
            let mut Ebo = BCo ^ ((!BCu) & BCa);
            let mut Ebu = BCu ^ ((!BCa) & BCe);

            Abo ^= Do;
            BCa = Abo.rotate_left(28);
            Agu ^= Du;
            BCe = Agu.rotate_left(20);
            Aka ^= Da;
            BCi = Aka.rotate_left(3);
            Ame ^= De;
            BCo = Ame.rotate_left(45);
            Asi ^= Di;
            BCu = Asi.rotate_left(61);

            let mut Ega = BCa ^ ((!BCe) & BCi);
            let mut Ege = BCe ^ ((!BCi) & BCo);
            let mut Egi = BCi ^ ((!BCo) & BCu);
            let mut Ego = BCo ^ ((!BCu) & BCa); 
            let mut Egu = BCu ^ ((!BCa) & BCe);

            Abe ^= De;
            BCa = Abe.rotate_left(1);
            Agi ^= Di;
            BCe = Agi.rotate_left(6);
            Ako ^= Do;
            BCi = Ako.rotate_left(25);
            Amu ^= Du;
            BCo = Amu.rotate_left(8);
            Asa ^= Da; 
            BCu = Asa.rotate_left(18);

            let mut Eka = BCa ^ ((!BCe) & BCi);
            let mut Eke = BCe ^ ((!BCi) & BCo);
            let mut Eki = BCi ^ ((!BCo) & BCu);
            let mut Eko = BCo ^ ((!BCu) & BCa); 
            let mut Eku = BCu ^ ((!BCa) & BCe);
            
            Abu ^= Du;
            BCa = Abu.rotate_left(27);
            Aga ^= Da;
            BCe = Aga.rotate_left(36);
            Ake ^= De;
            BCi = Ake.rotate_left(10);
            Ami ^= Di;
            BCo = Ami.rotate_left(15);
            Aso ^= Do;
            BCu = Aso.rotate_left(56);

            let mut Ema = BCa ^ ((!BCe) & BCi);
            let mut Eme = BCe ^ ((!BCi) & BCo);
            let mut Emi = BCi ^ ((!BCo) & BCu);
            let mut Emo = BCo ^ ((!BCu) & BCa); 
            let mut Emu = BCu ^ ((!BCa) & BCe);

            Abi ^= Di;
            BCa = Abi.rotate_left(62);
            Ago ^= Do;
            BCe = Ago.rotate_left(55);
            Aku ^= Du;
            BCi = Aku.rotate_left(39);
            Ama ^= Da;
            BCo = Ama.rotate_left(41);
            Ase ^= De;
            BCu = Ase.rotate_left(2);
            
            let mut Esa = BCa ^ ((!BCe) & BCi);
            let mut Ese = BCe ^ ((!BCi) & BCo);
            let mut Esi = BCi ^ ((!BCo) & BCu);
            let mut Eso = BCo ^ ((!BCu) & BCa); 
            let mut Esu = BCu ^ ((!BCa) & BCe);

            // Prepare theta
            BCa = Eba ^ Ega ^ Eka ^ Ema ^ Esa;
            BCe = Ebe ^ Ege ^ Eke ^ Eme ^ Ese;
            BCi = Ebi ^ Egi ^ Eki ^ Emi ^ Esi;
            BCo = Ebo ^ Ego ^ Eko ^ Emo ^ Eso;
            BCu = Ebu ^ Egu ^ Eku ^ Emu ^ Esu;

            Da = BCu ^ BCe.rotate_left(1);
            De = BCa ^ BCi.rotate_left(1);
            Di = BCe ^ BCo.rotate_left(1);
            Do = BCi ^ BCu.rotate_left(1);
            Du = BCo ^ BCa.rotate_left(1);

            Eba ^= Da;
            BCa = Eba;
            Ege ^= De;
            BCe = Ege.rotate_left(44);
            Eki ^= Di;
            BCi = Eki.rotate_left(43);
            Emo ^= Do;
            BCo = Emo.rotate_left(21);
            Esu ^= Du;
            BCu = Esu.rotate_left(14);

            Aba = BCa ^ ((!BCe) & BCi);
            Aba ^= KECCAK_F_ROUND_CONSTANTS[round+1];
            Abe = BCe ^ ((!BCi) & BCo);
            Abi = BCi ^ ((!BCo) & BCu); 
            Abo = BCo ^ ((!BCu) & BCa);
            Abu = BCu ^ ((!BCa) & BCe);

            Ebo ^= Do;
            BCa = Ebo.rotate_left(28);
            Egu ^= Du;
            BCe = Egu.rotate_left(20);
            Eka ^= Da;
            BCi = Eka.rotate_left(3);
            Eme ^= De;
            BCo = Eme.rotate_left(45);
            Esi ^= Di;
            BCu = Esi.rotate_left(61);

            Aga = BCa ^ ((!BCe) & BCi);
            Age = BCe ^ ((!BCi) & BCo);
            Agi = BCi ^ ((!BCo) & BCu);
            Ago = BCo ^ ((!BCu) & BCa); 
            Agu = BCu ^ ((!BCa) & BCe);

            Ebe ^= De;
            BCa = Ebe.rotate_left(1);
            Egi ^= Di;
            BCe = Egi.rotate_left(6);
            Eko ^= Do;
            BCi = Eko.rotate_left(25);
            Emu ^= Du;
            BCo = Emu.rotate_left(8);
            Esa ^= Da; 
            BCu = Esa.rotate_left(18);

            Aka = BCa ^ ((!BCe) & BCi);
            Ake = BCe ^ ((!BCi) & BCo);
            Aki = BCi ^ ((!BCo) & BCu);
            Ako = BCo ^ ((!BCu) & BCa); 
            Aku = BCu ^ ((!BCa) & BCe);
            
            Ebu ^= Du;
            BCa = Ebu.rotate_left(27);
            Ega ^= Da;
            BCe = Ega.rotate_left(36);
            Eke ^= De;
            BCi = Eke.rotate_left(10);
            Emi ^= Di;
            BCo = Emi.rotate_left(15);
            Eso ^= Do;
            BCu = Eso.rotate_left(56);

            Ama = BCa ^ ((!BCe) & BCi);
            Ame = BCe ^ ((!BCi) & BCo);
            Ami = BCi ^ ((!BCo) & BCu);
            Amo = BCo ^ ((!BCu) & BCa); 
            Amu = BCu ^ ((!BCa) & BCe);

            Ebi ^= Di;
            BCa = Ebi.rotate_left(62);
            Ego ^= Do;
            BCe = Ego.rotate_left(55);
            Eku ^= Du;
            BCi = Eku.rotate_left(39);
            Ema ^= Da;
            BCo = Ema.rotate_left(41);
            Ese ^= De;
            BCu = Ese.rotate_left(2);
            
            Asa = BCa ^ ((!BCe) & BCi);
            Ase = BCe ^ ((!BCi) & BCo);
            Asi = BCi ^ ((!BCo) & BCu);
            Aso = BCo ^ ((!BCu) & BCa); 
            Asu = BCu ^ ((!BCa) & BCe);
        }

        state[ 0] = Aba; state[ 1] = Abe; state[ 2] = Abi; state[ 3] = Abo; state[ 4] = Abu; 
        state[ 5] = Aga; state[ 6] = Age; state[ 7] = Agi; state[ 8] = Ago; state[ 9] = Agu; 
        state[10] = Aka; state[11] = Ake; state[12] = Aki; state[13] = Ako; state[14] = Aku; 
        state[15] = Ama; state[16] = Ame; state[17] = Ami; state[18] = Amo; state[19] = Amu; 
        state[20] = Asa; state[21] = Ase; state[22] = Asi; state[23] = Aso; state[24] = Asu; 
    }

    pub fn absorb(
        s: &mut [u64; 25], 
        r: usize, 
        mut pos: usize, 
        mut m: &[u8]
    ) -> usize {
        let mut t = [0u8; 8];

        if (pos & 7) != 0 {
            let mut i = pos & 7;
            while i < 8 && !m.is_empty() {
                t[i] = m[0];
                m = &m[1..];
                i += 1;
                pos += 1;
            }
            s[(pos-i)/8] ^= Self::load64(&t);
        }

        if pos != 0 && m.len() >= r - pos {
            let block_bytes = r - pos;
            for i in 0..(block_bytes/8) {
                let mut chunk = [0u8; 8];
                chunk.copy_from_slice(&m[8*i..8*i+8]);
                s[pos/8 + i] ^= Self::load64(&chunk);
            }
            m = &m[block_bytes..];
            pos = 0;
            Self::keccak_f1600_state_permute(s);
        }

        while m.len() >= r {
            for i in 0..(r/8) {
                let mut chunk = [0u8; 8];
                chunk.copy_from_slice(&m[8*i..8*i+8]);
                s[i] ^= Self::load64(&chunk);
            }
            m = &m[r..];
            Self::keccak_f1600_state_permute(s);
        }

        let block_count = m.len() / 8;
        for i in 0..block_count {
            let mut chunk = [0u8; 8];
            chunk.copy_from_slice(&m[8*i..8*i+8]);
            s[i] ^= Self::load64(&chunk);
        }
        m = &m[block_count*8..];
        pos += 8*block_count;

        if !m.is_empty() {
            t.fill(0);
            for i in 0..m.len() {
                t[i] = m[i];
            }
            s[pos/8] ^= Self::load64(&t);
            pos += m.len();
        }

        pos
    }

    pub fn finalize(
        s: &mut [u64; 25],
        r: usize,
        pos: usize,
        p: u8
    ) {
        let i = pos >> 3;
        let j = pos & 7;
        s[i] ^= (p as u64) << 8*j;
        s[r/8-1] ^= (1u64) << 63;
    }

    pub fn squeezeblocks (
        out: &mut [u8],
        nblocks: usize,
        s: &mut [u64; 25],
        r: usize
    ) {
        let mut nblks = nblocks;
        let mut offset = 0;
        while nblks > 0 {
            Self::keccak_f1600_state_permute(s);
            for i in 0..(r/8) {
                let bytes = Self::store64(s[i]);
                out[offset + 8*i..offset + 8*i + 8].copy_from_slice(&bytes);
            } 
            offset += r;
            nblks -= 1;
        }
    }

    pub fn squeeze (
        out: &mut [u8],
        s: &mut [u64; 25],
        r: usize,
        mut pos: usize,
    ) -> usize {
        let mut outlen = out.len();
        let mut offset = 0;

        if (pos & 7) != 0 {
            let t = Self::store64(s[pos/8]);
            let mut i = pos & 7;
            while i < 8 && outlen > 0 {
                out[offset] = t[i];
                offset += 1; i += 1;
                pos += 1; outlen -= 1;
            }
        }
        
        if pos != 0 && outlen >= r - pos {
            for i in 0..(r-pos)/8 {
                let t =  Self::store64(s[pos/8+i]);
                out[offset + 8*i..offset + 8*i + 8].copy_from_slice(&t);
            }
            outlen -= r-pos;
            offset += r-pos;
            pos = 0;
        }

        while outlen >= r {
            Self::keccak_f1600_state_permute(s);
            for i in 0..(r/8) {
                let t = Self::store64(s[i]);
                out[offset + 8*i..offset + 8*i + 8].copy_from_slice(&t);
            }
            outlen -= r;
            offset += r;
        }

        if outlen == 0 {
            return pos;
        }
        
        if pos == 0 {
            Self::keccak_f1600_state_permute(s);
        }
        
        let remaining_blocks = outlen/8;
        for i in 0..remaining_blocks {
            let t = Self::store64(s[pos/8+i]);
            out[offset + 8*i..offset+ 8*i +8].copy_from_slice(&t);
        }
        outlen -= 8 * remaining_blocks;
        offset += 8 * remaining_blocks;
        pos += 8 * remaining_blocks;

        if outlen > 0 {
            let t = Self::store64(s[pos/8]);
            out[offset..offset+outlen].copy_from_slice(&t[..outlen]);
            pos += outlen;
        }

        pos
    }
}
