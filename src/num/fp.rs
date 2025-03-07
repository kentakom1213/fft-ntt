//! 有限体の実装

/// 有限体の実装
#[derive(Debug)]
pub struct Fp {
    /// mod p
    pub p: u64,
    /// p の原始根
    pub root: u64,
    /// root の逆元
    pub rinv: u64,
    /// (p-1) に素因数として含まれる 2 の個数
    pub count_2: usize,
}

impl Fp {
    /// 初期化する
    pub fn new(p: u64) -> Result<Self, &'static str> {
        // p は素数である必要がある
        if Fp::factorize(p).len() > 1 {
            return Err("`p` should be prime number.");
        }

        // (p - 1) を素因数分解
        let factors = Fp::factorize(p - 1);

        // 原始根を探索
        let root = Self::find_root(p, &factors);

        // (p-1) に素因数として含まれる 2 の個数
        let count_2 = factors[0].1 as usize;

        Ok(Self {
            p,
            root,
            rinv: Self::_inv(p, root),
            count_2,
        })
    }

    /// Fpの原始根を探索する
    fn find_root(p: u64, factors: &Vec<(u64, u64)>) -> u64 {
        // x が Fp の原始根であるか判定する
        let is_ok = |x: u64| {
            factors
                .iter()
                .all(|&(pi, _)| Self::_pow(p, x, (p - 1) / pi) != 0)
        };

        (1..p)
            .find(|x| is_ok(*x))
            // 原始根の存在性より
            .unwrap()
    }

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
    fn normalize(p: u64, a: u64) -> u64 {
        if a < p {
            return a;
        }
        a % p
    }

    /// a + b (mod p)
    fn _add(p: u64, a: u64, b: u64) -> u64 {
        let a = Self::normalize(p, a);
        let b = Self::normalize(p, b);

        let mut res = a + b;
        if res >= p {
            res -= p;
        }
        res
    }

    /// - a (mod p)
    fn _neg(p: u64, a: u64) -> u64 {
        let a = Self::normalize(p, a);

        p - a
    }

    /// a - b (mod p)
    fn _sub(p: u64, a: u64, b: u64) -> u64 {
        let a = Self::normalize(p, a);
        let b = Self::normalize(p, b);

        Self::_add(p, a, Self::_neg(p, b))
    }

    /// a * b (mod p)
    fn _mul(p: u64, a: u64, b: u64) -> u64 {
        let a = Self::normalize(p, a);
        let b = Self::normalize(p, b);

        a * b % p
    }

    /// a ^ b (mod p)
    fn _pow(p: u64, a: u64, mut b: u64) -> u64 {
        let mut a = Self::normalize(p, a);
        let mut res = 1;
        while b > 0 {
            if b & 1 == 1 {
                res = Self::_mul(p, res, a);
            }
            a = Self::_mul(p, a, a);
            b >>= 1;
        }
        res
    }

    /// a^(-1) mod p
    fn _inv(p: u64, a: u64) -> u64 {
        Self::_pow(p, a, p - 2)
    }

    // ===== 公開する演算 =====
    /// a + b (mod p)
    pub fn add(&self, a: u64, b: u64) -> u64 {
        Self::_add(self.p, a, b)
    }
    /// -a (mod p)
    pub fn neg(&self, a: u64) -> u64 {
        Self::_neg(self.p, a)
    }
    /// a - b (mod p)
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        Self::_sub(self.p, a, b)
    }
    /// a * b (mod p)
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        Self::_mul(self.p, a, b)
    }
    /// a ^ b (mod p)
    pub fn pow(&self, a: u64, b: u64) -> u64 {
        Self::_pow(self.p, a, b)
    }
    /// a^(-1) (mod p)
    pub fn inv(&self, a: u64) -> u64 {
        Self::_inv(self.p, a)
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
