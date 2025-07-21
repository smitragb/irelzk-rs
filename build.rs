fn main() {
    cc::Build::new()
        .file("src/asm/ntt.S")
        .include("src/asm/")
        .flag("-mavx2")
        .compile("ntt_avx");

    cc::Build::new()
        .file("src/asm/invntt.S")
        .include("src/asm/")
        .flag("-mavx2")
        .compile("invntt_avx");

    // Rebuild if changed
    println!("cargo:rerun-if-changed=src/asm/ntt.S");
    println!("cargo:rerun-if-changed=src/asm/invntt.S");
    println!("cargo:rerun-if-changed=src/asm/shuffle.inc");

}
