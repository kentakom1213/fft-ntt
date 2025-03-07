//! 有限体の実装

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// 有限体の実装
#[derive(Debug)]
pub struct Fp {
    /// mod p
    pub p: u64,
    /// p の原始根
    pub root: u64,
    /// root の逆元
    pub rinv: u64,
}

impl Fp {
    /// 初期化する
    fn new(p: u64) -> Result<Self, &'static str> {
        Ok(Self {
            p,
            root: 0,
            rinv: 0,
        })
    }

    // /// Fpの原始根を探索する
    // fn find_root(p: u64, factors: &Vec<(u64, u64)>) -> u64 {
    //     // x が Fp の原始根であるか判定する
    //     let is_ok = |x: u64| {
    //         factors.iter().all(|&(pi, _)| {
    //             x.pow((p - 1) / pi)
    //         })
    //     };
    // }

    /// 素因数分解
    fn factorize(mut x: u64) -> Vec<(u64, u64)> {
        let mut res = vec![];

        for p in std::iter::once(2).chain((1..).map(|x| 2 * x + 1)) {
            if p * p > x {
                break;
            }

            let mut cnt = 0;
            while x % p == 0 {
                cnt += 1;
                x /= p;
            }

            if cnt > 0 {
                res.push((p, cnt));
            }
        }

        if x > 1 {
            res.push((x, 1));
        }

        res
    }

    // ===== 基本的な演算の実装 =====
    /// 0 <= a < p となるように正規化
    fn normalize(&self, a: u64) -> u64 {
        if a < self.p {
            return a;
        }
        a % self.p
    }

    /// a + b (mod p)
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let a = self.normalize(a);
        let b = self.normalize(b);

        let mut res = a + b;
        if res >= self.p {
            res -= self.p;
        }
        res
    }

    /// - a (mod p)
    pub fn neg(&self, a: u64) -> u64 {
        let a = self.normalize(a);

        self.p - a
    }

    /// a - b (mod p)
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        let a = self.normalize(a);
        let b = self.normalize(b);

        self.add(a, self.neg(b))
    }

    /// a * b (mod p)
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        let a = self.normalize(a);
        let b = self.normalize(b);

        a * b % self.p
    }

    /// a ^ b (mod p)
    pub fn pow(&self, a: u64, mut b: u64) -> u64 {
        let mut a = self.normalize(a);
        let mut res = 1;
        while b > 0 {
            if b & 1 == 1 {
                res = self.mul(res, a);
            }
            a = self.mul(a, a);
            b >>= 1;
        }
        res
    }

    /// a^(-1) mod p
    pub fn inv(&self, a: u64) -> u64 {
        self.pow(a, self.p - 2)
    }
}

// ===== テスト =====
#[cfg(test)]
mod test {
    use super::Fp;

    const P: u64 = 998244353;

    #[test]
    fn test_factorize() {
        assert_eq!(Fp::factorize(0), vec![]);
        assert_eq!(Fp::factorize(1), vec![]);
        assert_eq!(Fp::factorize(2), vec![(2, 1)]);
        assert_eq!(Fp::factorize(12), vec![(2, 2), (3, 1)]);
        assert_eq!(Fp::factorize(1024), vec![(2, 10)]);
        assert_eq!(Fp::factorize(3628800), vec![(2, 8), (3, 4), (5, 2), (7, 1)]);
        assert_eq!(Fp::factorize(1048576), vec![(2, 20)]);
        assert_eq!(Fp::factorize(998244353), vec![(998244353, 1)]);
        assert_eq!(Fp::factorize(1000000007), vec![(1000000007, 1)]);
    }

    #[test]
    fn test_add() {
        let fp = Fp::new(P).unwrap();

        assert_eq!(fp.add(2, 10), 12);
        assert_eq!(fp.add(2, P - 1), 1);
    }

    #[test]
    fn test_pow() {
        let fp = Fp::new(P).unwrap();

        assert_eq!(fp.pow(2, 1), 2);
        assert_eq!(fp.pow(2, 2), 4);
        assert_eq!(fp.pow(2, 3), 8);
        assert_eq!(fp.pow(2, 4), 16);
        assert_eq!(fp.pow(2, 5), 32);
        assert_eq!(fp.pow(2, 6), 64);
        assert_eq!(fp.pow(2, 7), 128);
    }

    #[test]
    fn test_inv() {
        let fp = Fp::new(P).unwrap();

        for x in (1..=10).chain(100000..=100010).chain(100000000..=100000010) {
            let xinv = fp.inv(x);

            eprintln!("inverse of {x:?} mod {P} is {xinv:?}");
            assert_eq!(fp.mul(x, xinv), 1);
        }
    }
}
