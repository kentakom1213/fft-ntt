//! 有限体の実装

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// 有限体の実装
///
/// `P`: 基数（素数のみ）
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Fp<const P: u64>(pub u64);

// ===== 基本的な演算の実装 =====
impl<const P: u64> Add for Fp<P> {
    type Output = Fp<P>;
    fn add(self, rhs: Self) -> Self::Output {
        let mut res = self.0 + rhs.0;
        if res >= P {
            res -= P;
        }
        Fp(res)
    }
}

impl<const P: u64> Neg for Fp<P> {
    type Output = Fp<P>;
    fn neg(self) -> Self::Output {
        Fp(P - self.0)
    }
}

impl<const P: u64> Sub for Fp<P> {
    type Output = Fp<P>;
    fn sub(self, rhs: Self) -> Self::Output {
        self + rhs.neg()
    }
}

impl<const P: u64> Mul for Fp<P> {
    type Output = Fp<P>;
    fn mul(self, rhs: Self) -> Self::Output {
        Fp(self.0 * rhs.0 % P)
    }
}

// ===== 代入演算 =====
impl<const P: u64> AddAssign for Fp<P> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const P: u64> SubAssign for Fp<P> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const P: u64> MulAssign for Fp<P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

// ===== べき乗 =====

/// あまりを取るべき乗の実装
pub trait ModPow {
    /// calc pow(a, x) mod P
    fn pow(self, x: u64) -> Self;
    /// calc inverse of x mod P
    fn inv(self) -> Self;
}

impl<const P: u64> ModPow for Fp<P> {
    fn pow(self, mut x: u64) -> Self {
        let mut tmp = self;
        let mut res = Fp(1);
        while x > 0 {
            if x & 1 == 1 {
                res *= tmp;
            }
            tmp = tmp * tmp;
            x >>= 1;
        }
        res
    }
    fn inv(self) -> Self {
        self.pow(P - 2)
    }
}

// ===== テスト =====
#[cfg(test)]
mod test {
    use crate::num::fp::ModPow;

    use super::Fp;

    const P: u64 = 998244353;

    #[test]
    fn test_add() {
        let mut a: Fp<P> = Fp(2);
        let b = Fp(P - 1);

        assert_eq!(a + Fp(10), Fp(12));
        assert_eq!(a + b, Fp(1));

        a += Fp(1);

        assert_eq!(a, Fp(3));

        a += b;

        assert_eq!(a, Fp(2));
    }

    #[test]
    fn test_pow() {
        let a: Fp<P> = Fp(2);

        assert_eq!(a.pow(1), Fp(2));
        assert_eq!(a.pow(2), Fp(4));
        assert_eq!(a.pow(3), Fp(8));
        assert_eq!(a.pow(4), Fp(16));
        assert_eq!(a.pow(5), Fp(32));
        assert_eq!(a.pow(6), Fp(64));
        assert_eq!(a.pow(7), Fp(128));
    }

    #[test]
    fn test_inv() {
        for i in (1..=10).chain(100000..=100010).chain(100000000..=100000010) {
            let x: Fp<P> = Fp(i);
            let xinv = x.inv();

            eprintln!("inverse of {x:?} mod {P} is {xinv:?}");
            assert_eq!(x * xinv, Fp(1));
        }
    }
}
