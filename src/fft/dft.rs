//! 離散フーリエ変換の実装

use crate::num::Fp;

/// 離散フーリエ変換の実装
pub struct DFT(Fp);

impl DFT {
    /// 入力された配列をフーリエ変換する
    pub fn dft(&self, array: &[u64]) -> Result<Vec<u64>, &'static str> {
        todo!()
    }

    /// 入力された配列をフーリエ逆変換する
    pub fn idft(&self, array: &[u64]) -> Result<Vec<u64>, &'static str> {
        todo!()
    }

    /// フーリエ変換，フーリエ逆変換の共通部分
    ///
    /// - `W`: 回転演算子
    fn dft_inner(&self, array: &[u64], W: u64) -> Result<Vec<u64>, &'static str> {
        todo!()
    }

    /// 長さが 2 べきになるように配列を生成する
    fn extend_array(&self, array: &[u64]) -> Result<Vec<u64>, &'static str> {
        let n = array.len();
        // 2^i >= n となるような最小の i
        let mut i = 0;
        let mut n_ = 1;
        while n_ < n {
            i += 1;
            n_ *= 2;
        }
        if i > self.0.count_2 {
            return Err("The prime p does not have enough factors of 2 in (p - 1).");
        }
        // 配列を生成
        let mut res = array.to_vec();
        // 残りをゼロ埋め
        res.extend(std::iter::repeat_n(0, n_ - n));

        Ok(res)
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

        assert_eq!(dft.extend_array(&arr_1), Ok(vec![1, 2, 3, 0]));
        assert_eq!(dft.extend_array(&arr_2), Ok(vec![1, 2, 3, 4]));
        assert!(dft.extend_array(&arr_3).is_err());
    }
}
