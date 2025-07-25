pub mod params;
pub mod poly;
pub mod rounding;
pub mod ntt;
pub mod consts;
pub mod crypto {
    pub mod aes256;
    pub mod shake;
    pub mod keccak;
}
pub mod polyvec;
use consts::{MONT, QDATA, _8XDIV, _8XQ, _8XQINV, _PMASK, _ZETAS, _ZETAS_QINV};
use params::Q;

fn main() {
    println!("Hello, world!, {}, {}", Q, MONT);
    println!("Checking the consts QDATA array: ");
    println!("QDATA[{}] = {}", _8XQ, QDATA[_8XQ]);
    println!("QDATA[{}] = {}", _8XQINV, QDATA[_8XQINV]);
    println!("QDATA[{}] = {}", _8XDIV, QDATA[_8XDIV]);
    println!("QDATA[{}] = {}", _PMASK, QDATA[_PMASK]);
    println!("QDATA[{}] = {}", _ZETAS, QDATA[_ZETAS]);
    println!("QDATA[{}] = {}", _ZETAS_QINV, QDATA[_ZETAS_QINV]);
    println!("QDATA[{}] = {}", 287, QDATA[287]);
}
