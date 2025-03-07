//! 離散フーリエ変換の実装

use crate::num::Fp;

/// 離散フーリエ変換の実装
pub struct DFT(Fp);

impl DFT {
    /// 入力された配列をフーリエ変換する
    pub fn dft(&self, X: &[u64]) -> Result<Vec<u64>, &'static str> {
        let (i, X) = self.extend_array(X)?;
        let w = self.0.root_pow2m(i)?;

        Ok(self.dft_inner(&X, w))
    }

    /// 入力された配列をフーリエ逆変換する
    pub fn idft(&self, F: &[u64]) -> Result<Vec<u64>, &'static str> {
        let (i, F) = self.extend_array(F)?;
        let w = self.0.root_pow2m(i)?;
        let winv = self.0.inv(w);

        let mut res = self.dft_inner(&F, winv);
        let n = res.len();

        // 逆変換後の配列を正規化
        let inv_n = self.0.inv(n as u64);
        res.iter_mut().for_each(|v| *v = self.0.mul(*v, inv_n));

        Ok(res)
    }

    /// フーリエ変換，フーリエ逆変換の共通部分
    ///
    /// - `w`: 回転演算子
    fn dft_inner(&self, X: &Vec<u64>, w: u64) -> Vec<u64> {
        let n = X.len();

        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let rot = self.0.pow(w, i * j);
                        self.0.mul(X[j], rot)
                    })
                    .fold(0, |acc, v| self.0.add(acc, v))
            })
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
    use crate::num::Fp;

    use super::DFT;

    #[test]
    fn test_extend_array() {
        let arr_1 = vec![1, 2, 3];
        let arr_2 = vec![1, 2, 3, 4];
        let arr_3 = vec![1, 2, 3, 4, 5];

        let fp = Fp::new(5).unwrap();
        let dft = DFT(fp);

        assert_eq!(dft.extend_array(&arr_1), Ok((2, vec![1, 2, 3, 0])));
        assert_eq!(dft.extend_array(&arr_2), Ok((2, vec![1, 2, 3, 4])));
        assert!(dft.extend_array(&arr_3).is_err());
    }

    #[test]
    fn test_dft() {
        {
            let arr = vec![1, 2, 3, 4];
            let fp = Fp::new(5).unwrap();
            eprintln!("fp = {:?}", fp);

            let dft = DFT(fp);

            let res = dft.dft(&arr).unwrap();
            eprintln!("dft({:?}) = {:?}", arr, res);

            let res2 = dft.idft(&res).unwrap();
            eprintln!("idft({:?}) = {:?}", res, res2);

            assert_eq!(res2, arr);
        }

        {
            let arr = vec![3, 1, 4, 1, 5, 9];
            let fp = Fp::new(17).unwrap();

            eprintln!("fp = {:?}", fp);

            let dft = DFT(fp);

            let res = dft.dft(&arr).unwrap();
            eprintln!("dft({:?}) = {:?}", arr, res);

            let res2 = dft.idft(&res).unwrap();
            eprintln!("idft({:?}) = {:?}", res, res2);

            let arr_ext = vec![3, 1, 4, 1, 5, 9, 0, 0];
            assert_eq!(res2, arr_ext);
        }

        {
            let arr = vec![31415, 92653, 58979, 32384, 62643, 38327, 95028];
            let fp = Fp::new(5767169).unwrap();

            eprintln!("fp = {:?}", fp);

            let dft = DFT(fp);

            let res = dft.dft(&arr).unwrap();
            eprintln!("dft({:?}) = {:?}", arr, res);

            let res2 = dft.idft(&res).unwrap();
            eprintln!("idft({:?}) = {:?}", res, res2);

            let arr_ext = vec![31415, 92653, 58979, 32384, 62643, 38327, 95028, 0];
            assert_eq!(res2, arr_ext);
        }
    }
}
