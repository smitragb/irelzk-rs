mod params;
/*
mod keccak;
mod shake;
mod aes256;
*/
use params::{GAMMA2, Q};

fn main() {
    println!("Hello, world!, {}, {}", Q, GAMMA2);
}
