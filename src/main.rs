mod aes256;
mod params;
mod keccak;
mod shake;
use params::{GAMMA2, Q};

fn main() {
    println!("Hello, world!, {}, {}", Q, GAMMA2);
}
