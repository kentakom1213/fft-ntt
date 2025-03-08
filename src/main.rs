use fft::{ntt::FFT, num::Fp};

fn main() {
    let arr = vec![314, 159, 265, 358, 979, 323, 846, 264];

    let fft = FFT(Fp::new(12289).unwrap());

    let encoded = fft.fft(&arr).unwrap();
    let decoded = fft.ifft(&encoded).unwrap();

    eprintln!("arr: {:?}", arr);
    eprintln!("encoded: {:?}", encoded);
    eprintln!("decoded: {:?}", decoded);

    assert_eq!(arr, decoded);
    println!("Success!");
}
