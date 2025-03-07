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
    /// p = 2^k * m + 1 となるような m
    pub k: usize,
    /// p = 2^k * m + 1 となるような m
    pub m: u64,
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
        let k = factors[0].1 as usize;

        Ok(Self {
            p,
            root,
            rinv: Self::_inv(p, root),
            k,
            m: (p - 1) >> k,
        })
    }

    /// Fpの原始根を探索する
    fn find_root(p: u64, factors: &Vec<(u64, u64)>) -> u64 {
        // x が Fp の原始根であるか判定する
        let is_ok = |x: u64| {
            factors
                .iter()
                .all(|&(pi, _)| Self::_pow(p, x, (p - 1) / pi) != 1)
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
    pub fn pow(&self, a: u64, b: usize) -> u64 {
        Self::_pow(self.p, a, b as u64)
    }
    /// a^(-1) (mod p)
    pub fn inv(&self, a: u64) -> u64 {
        Self::_inv(self.p, a)
    }
    /// 2^(1 / 2^a) (mod p)
    pub fn root_pow2m(&self, a: usize) -> Result<u64, &'static str> {
        if a > self.k {
            return Err("The prime p does not have enough factors of 2 in (p - 1).");
        }

        Ok(Self::_pow(self.p, self.root, self.m << (self.k - a)))
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
    fn test_mul() {
        let fp = Fp::new(P).unwrap();

        assert_eq!(fp.mul(2, 10), 20);
        assert_eq!(fp.mul(2, P - 1), P - 2);
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

    #[test]
    fn test_find_root() {
        let fp5 = Fp::new(5).unwrap();
        assert_eq!(fp5.root, 2);

        let fp7 = Fp::new(7).unwrap();
        assert_eq!(fp7.root, 3);

        let fp11 = Fp::new(11).unwrap();
        assert_eq!(fp11.root, 2);

        let fpP = Fp::new(P).unwrap();
        assert_eq!(fpP.root, 3);
    }

    #[test]
    fn test_pow2m_root() {
        {
            let fp5 = Fp::new(5).unwrap();

            let fp5_0 = fp5.root_pow2m(0).unwrap();
            eprintln!("2^(1 / 2^0) mod 5 = {:?}", fp5_0);
            assert_eq!(fp5.pow(fp5_0, 1 << 0), 1);

            let fp5_1 = fp5.root_pow2m(1).unwrap();
            eprintln!("2^(1 / 2^1) mod 5 = {:?}", fp5_1);
            assert_eq!(fp5.pow(fp5_1, 1 << 1), 1);

            let fp5_2 = fp5.root_pow2m(2).unwrap();
            eprintln!("2^(1 / 2^2) mod 5 = {:?}", fp5_2);
            assert_eq!(fp5.pow(fp5_2, 1 << 2), 1);

            assert!(fp5.root_pow2m(3).is_err());
            eprintln!("2^(1 / 2^3) mod 5 = {:?}", fp5.root_pow2m(3));
            assert_eq!(
                fp5.root_pow2m(3),
                Err("The prime p does not have enough factors of 2 in (p - 1).")
            );
        }

        {
            let fpP = Fp::new(P).unwrap();

            let fpP_0 = fpP.root_pow2m(0).unwrap();
            eprintln!("2^(1 / 2^0) mod {P} = {:?}", fpP_0);
            assert_eq!(fpP.pow(fpP_0, 1 << 0), 1);

            let fpP_1 = fpP.root_pow2m(1).unwrap();
            eprintln!("2^(1 / 2^1) mod {P} = {:?}", fpP_1);
            assert_eq!(fpP.pow(fpP_1, 1 << 1), 1);

            let fpP_2 = fpP.root_pow2m(2).unwrap();
            eprintln!("2^(1 / 2^2) mod {P} = {:?}", fpP_2);
            assert_eq!(fpP.pow(fpP_2, 1 << 2), 1);

            let fpP_3 = fpP.root_pow2m(3).unwrap();
            eprintln!("2^(1 / 2^3) mod {P} = {:?}", fpP_3);
            assert_eq!(fpP.pow(fpP_3, 1 << 3), 1);

            let fpP_4 = fpP.root_pow2m(4).unwrap();
            eprintln!("2^(1 / 2^4) mod {P} = {:?}", fpP_4);
            assert_eq!(fpP.pow(fpP_4, 1 << 4), 1);

            let fpP_5 = fpP.root_pow2m(5).unwrap();
            eprintln!("2^(1 / 2^5) mod {P} = {:?}", fpP_5);
            assert_eq!(fpP.pow(fpP_5, 1 << 5), 1);

            let fpP_6 = fpP.root_pow2m(6).unwrap();
            eprintln!("2^(1 / 2^6) mod {P} = {:?}", fpP_6);
            assert_eq!(fpP.pow(fpP_6, 1 << 6), 1);

            let fpP_7 = fpP.root_pow2m(7).unwrap();
            eprintln!("2^(1 / 2^7) mod {P} = {:?}", fpP_7);
            assert_eq!(fpP.pow(fpP_7, 1 << 7), 1);

            let fpP_8 = fpP.root_pow2m(8).unwrap();
            eprintln!("2^(1 / 2^8) mod {P} = {:?}", fpP_8);
            assert_eq!(fpP.pow(fpP_8, 1 << 8), 1);

            let fpP_9 = fpP.root_pow2m(9).unwrap();
            eprintln!("2^(1 / 2^9) mod {P} = {:?}", fpP_9);
            assert_eq!(fpP.pow(fpP_9, 1 << 9), 1);

            let fpP_10 = fpP.root_pow2m(10).unwrap();
            eprintln!("2^(1 / 2^10) mod {P} = {:?}", fpP_10);
            assert_eq!(fpP.pow(fpP_10, 1 << 10), 1);

            let fpP_20 = fpP.root_pow2m(20).unwrap();
            eprintln!("2^(1 / 2^20) mod {P} = {:?}", fpP_20);
            assert_eq!(fpP.pow(fpP_20, 1 << 20), 1);

            let fpP_21 = fpP.root_pow2m(21).unwrap();
            eprintln!("2^(1 / 2^21) mod {P} = {:?}", fpP_21);
            assert_eq!(fpP.pow(fpP_21, 1 << 21), 1);

            let fpP_22 = fpP.root_pow2m(22).unwrap();
            eprintln!("2^(1 / 2^22) mod {P} = {:?}", fpP_22);
            assert_eq!(fpP.pow(fpP_22, 1 << 22), 1);

            let fpP_23 = fpP.root_pow2m(23).unwrap();
            eprintln!("2^(1 / 2^23) mod {P} = {:?}", fpP_23);
            assert_eq!(fpP.pow(fpP_23, 1 << 23), 1);

            assert!(fpP.root_pow2m(24).is_err());
            eprintln!("2^(1 / 2^24) mod {P} = {:?}", fpP.root_pow2m(24));
            assert_eq!(
                fpP.root_pow2m(24),
                Err("The prime p does not have enough factors of 2 in (p - 1).")
            );
        }
    }
}
