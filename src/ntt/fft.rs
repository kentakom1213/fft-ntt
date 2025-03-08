//! 高速フーリエ変換の実装

use crate::num::Fp;

/// 高速フーリエ変換の実装
pub struct FFT(Fp);

impl FFT {
    /// 入力された配列をフーリエ変換する
    pub fn fft(&self, X: &[u64]) -> Result<Vec<u64>, &'static str> {
        let (i, X) = self.extend_array(X)?;
        let w = self.0.root_pow2m(i)?;

        Ok(self.fft_core(X, w))
    }

    /// 入力された配列をフーリエ逆変換する
    pub fn ifft(&self, F: &[u64]) -> Result<Vec<u64>, &'static str> {
        let (i, F) = self.extend_array(F)?;
        let w = self.0.root_pow2m(i)?;
        let winv = self.0.inv(w);

        let mut res = self.fft_core(F, winv);
        let n = res.len();

        // 逆変換後の配列を正規化
        let inv_n = self.0.inv(n as u64);
        res.iter_mut().for_each(|v| *v = self.0.mul(*v, inv_n));

        Ok(res)
    }

    /// フーリエ変換，フーリエ逆変換の共通部分
    ///
    /// - `w`: 回転演算子
    fn fft_core(&self, X: Vec<u64>, w: u64) -> Vec<u64> {
        let n = X.len();

        if n == 1 {
            return X.to_vec();
        }

        let (X_even, X_odd): (Vec<u64>, Vec<u64>) = (0..n / 2)
            .map(|i| {
                let l = X[i];
                let r = X[i + n / 2];
                (
                    self.0.add(l, r),
                    self.0.mul(self.0.sub(l, r), self.0.pow(w, i)),
                )
            })
            .collect();

        // 再帰的にFFT
        let new_w = self.0.pow(w, 2);

        let Y_even = self.fft_core(X_even, new_w);
        let Y_odd = self.fft_core(X_odd, new_w);

        // マージ
        Y_even
            .into_iter()
            .zip(Y_odd.into_iter())
            .flat_map(|(e, o)| [e, o])
            .collect()
    }

    /// 長さが 2 べきになるように配列を生成する
    ///
    /// **Arguments**
    /// - `array`: 配列
    ///
    /// **Returns**
    /// - `(i, res)`: 配列の長さを 2^i に拡張した結果
    fn extend_array(&self, array: &[u64]) -> Result<(usize, Vec<u64>), &'static str> {
        let n = array.len();
        // 2^i >= n となるような最小の i
        let mut i = 0;
        let mut n_ = 1;
        while n_ < n {
            i += 1;
            n_ *= 2;
        }
        if i > self.0.k {
            return Err("The prime p does not have enough factors of 2 in (p - 1).");
        }
        // 配列を生成
        let mut res = array.to_vec();
        // 残りをゼロ埋め
        res.extend(std::iter::repeat_n(0, n_ - n));

        Ok((i, res))
    }
}

// ===== テスト =====
#[cfg(test)]
mod test {
    use rand::{rng, Rng};
    use rstest::rstest;

    use crate::num::Fp;

    use super::FFT;

    #[test]
    fn test_extend_array() {
        let arr_1 = vec![1, 2, 3];
        let arr_2 = vec![1, 2, 3, 4];
        let arr_3 = vec![1, 2, 3, 4, 5];

        let fp = Fp::new(5).unwrap();
        let fft = FFT(fp);

        assert_eq!(fft.extend_array(&arr_1), Ok((2, vec![1, 2, 3, 0])));
        assert_eq!(fft.extend_array(&arr_2), Ok((2, vec![1, 2, 3, 4])));
        assert!(fft.extend_array(&arr_3).is_err());
    }

    #[test]
    fn test_fft() {
        {
            let arr = vec![1, 2, 3, 4];
            let fp = Fp::new(5).unwrap();
            eprintln!("\nfp = {:?}", fp);

            let fft = FFT(fp);

            let res = fft.fft(&arr).unwrap();
            eprintln!("fft({:?}) = {:?}", arr, res);

            let res2 = fft.ifft(&res).unwrap();
            eprintln!("ifft({:?}) = {:?}", res, res2);

            assert_eq!(res2, arr);
        }

        {
            let arr = vec![3, 1, 4, 1, 5, 9];
            let fp = Fp::new(17).unwrap();

            eprintln!("\nfp = {:?}", fp);

            let fft = FFT(fp);

            let res = fft.fft(&arr).unwrap();
            eprintln!("fft({:?}) = {:?}", arr, res);

            let res2 = fft.ifft(&res).unwrap();
            eprintln!("ifft({:?}) = {:?}", res, res2);

            let arr_ext = vec![3, 1, 4, 1, 5, 9, 0, 0];
            assert_eq!(res2, arr_ext);
        }

        {
            let arr = vec![31415, 92653, 58979, 32384, 62643, 38327, 95028];
            let fp = Fp::new(5767169).unwrap();

            eprintln!("\nfp = {:?}", fp);

            let fft = FFT(fp);

            let res = fft.fft(&arr).unwrap();
            eprintln!("fft({:?}) = {:?}", arr, res);

            let res2 = fft.ifft(&res).unwrap();
            eprintln!("ifft({:?}) = {:?}", res, res2);

            let arr_ext = vec![31415, 92653, 58979, 32384, 62643, 38327, 95028, 0];
            assert_eq!(res2, arr_ext);
        }

        {
            let arr = vec![
                31415926, 53589793, 23846264, 33832795, 02884197, 16939937, 51058209, 74944592,
            ];

            let fp = Fp::new(998244353).unwrap();

            eprintln!("\nfp = {:?}", fp);

            let fft = FFT(fp);

            let res = fft.fft(&arr).unwrap();
            eprintln!("fft({:?}) = {:?}", arr, res);

            let res2 = fft.ifft(&res).unwrap();
            eprintln!("ifft({:?}) = {:?}", res, res2);

            assert_eq!(res2, arr);
        }
    }

    #[rstest(
        size,
        p,
        case(500, 5767169),
        case(500, 5767169),
        case(3000, 5767169),
        case(500, 998244353),
        case(500, 998244353),
        case(3000, 998244353),
        // too large case
        case(200000, 998244353),
        case(200000, 998244353),
        case(200000, 998244353),
    )]
    fn test_fft_large(size: usize, p: u64) {
        let mut rng = rng();

        let arr: Vec<u64> = (0..size).map(|_| rng.random_range(0..p)).collect();

        let dft = FFT(Fp::new(p).unwrap());

        let res = dft.fft(&arr).unwrap();
        let res2 = dft.ifft(&res).unwrap();

        assert_eq!(&res2[..size], arr);
    }
}
