use core::f64;
use std::{
    borrow::Borrow,
    cmp::Ordering,
    fmt::Display,
    hash::Hash,
    iter::Sum,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
    sync::Arc,
};

use anyhow::{Error, Result, anyhow};

use crate::{
    ebi_number::{EbiNumber, Infinite, Normal, Round},
    exact::MaybeExact,
    fraction::EPSILON,
};

use super::ebi_number::{One, Signed, Zero};

#[derive(Debug, Clone, Copy)]
pub struct FractionF64(pub f64);

impl FractionF64 {
    pub fn two() -> Self {
        Self(2.0)
    }

    pub fn one_minus(self) -> Self {
        Self(1.0 - self.0)
    }

    pub fn is_sign_negative(&self) -> bool {
        self.0.is_sign_negative()
    }

    pub fn is_sign_positive(&self) -> bool {
        self.0.is_sign_positive()
    }

    /// Returns true if the value is Infinity (does not matter positive or negative)
    pub fn is_infinite(&self) -> bool {
        self.0.is_infinite()
    }

    pub fn is_nan(&self) -> bool {
        self.0.is_nan()
    }

    pub fn infinity() -> Self {
        Self(f64::INFINITY)
    }

    pub fn neg_infinity() -> Self {
        Self(f64::NEG_INFINITY)
    }

    pub fn nan() -> Self {
        Self(f64::NAN)
    }

    pub fn sqrt_abs(&self, _decimal_places: u32) -> FractionF64 {
        Self(self.0.abs().sqrt())
    }

    /**
     * 1/self
     */
    pub fn recip(&self) -> Self {
        Self(self.0.recip())
    }
}

impl EbiNumber for FractionF64 {}

impl MaybeExact for FractionF64 {
    type Approximate = f64;
    type Exact = fraction::BigFraction;

    fn is_exact(&self) -> bool {
        false
    }

    fn extract_approx(&self) -> Result<f64> {
        Ok(self.0)
    }

    fn extract_exact(&self) -> Result<&fraction::BigFraction> {
        Err(anyhow!("cannot extract a fraction from a float"))
    }
}

impl One for FractionF64 {
    fn one() -> Self {
        Self(1.0)
    }

    fn is_one(&self) -> bool {
        (self.0 - 1.0).abs() - &EPSILON < 0.0
    }
}

impl Zero for FractionF64 {
    fn zero() -> Self {
        Self(0.0)
    }

    fn is_zero(&self) -> bool {
        self.0.abs() - &EPSILON < 0.0
    }
}

impl Signed for FractionF64 {
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }

    fn is_positive(&self) -> bool {
        self.0 != 0f64 && self.0 > EPSILON
    }

    fn is_negative(&self) -> bool {
        self.0 != 0f64 && self.0 < -EPSILON
    }

    fn is_not_negative(&self) -> bool {
        self.0.is_not_negative()
    }

    fn is_not_positive(&self) -> bool {
        self.0.is_not_positive()
    }
}

impl Infinite for FractionF64 {
    fn is_infinite(&self) -> bool {
        self.0.is_infinite()
    }
}

impl Normal for FractionF64 {
    fn is_nan(&self) -> bool {
        self.0.is_nan()
    }
}

impl Round for FractionF64 {
    fn floor(self) -> Self {
        FractionF64(self.0.floor())
    }

    fn ceil(self) -> Self {
        FractionF64(self.0.ceil())
    }
}

impl Default for FractionF64 {
    fn default() -> Self {
        Self(Default::default())
    }
}

impl PartialEq for FractionF64 {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (FractionF64(l0), FractionF64(r0)) => l0 - EPSILON <= *r0 && *r0 <= l0 + EPSILON,
        }
    }
}

impl Display for FractionF64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl From<f64> for FractionF64 {
    fn from(value: f64) -> Self {
        Self(value)
    }
}

impl From<&FractionF64> for FractionF64 {
    fn from(value: &FractionF64) -> Self {
        value.clone()
    }
}

impl From<Arc<FractionF64>> for FractionF64 {
    fn from(value: Arc<FractionF64>) -> Self {
        value.as_ref().clone()
    }
}

impl From<&Arc<FractionF64>> for FractionF64 {
    fn from(value: &Arc<FractionF64>) -> Self {
        value.as_ref().clone()
    }
}

impl FromStr for FractionF64 {
    type Err = Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match f64::from_str(s) {
            Ok(f) => Ok(Self(f)),
            Err(_) => match fraction::Fraction::from_str(s) {
                Ok(f) => Ok(Self(format!("{:.20}", f).parse::<f64>()?)),
                Err(e) => Err(e.into()),
            },
        }
    }
}

impl Eq for FractionF64 {}

impl PartialOrd for FractionF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Hash for FractionF64 {
    /**
     * For good reasons, Rust does not support hashing of doubles. However, we need it to store distributions in a hashmap.
     * Approximate arithmetic is discouraged
     */
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        f64::to_bits(self.0).hash(state)
    }
}

impl Ord for FractionF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.0.is_nan() && other.0.is_nan() {
            Ordering::Equal
        } else if self.0.is_nan() {
            Ordering::Less
        } else if other.0.is_nan() {
            Ordering::Greater
        } else if self.0 == f64::INFINITY {
            if other.0 == f64::INFINITY {
                Ordering::Equal
            } else {
                Ordering::Greater
            }
        } else if other.0 == f64::INFINITY {
            Ordering::Less
        } else if self.0 == f64::NEG_INFINITY {
            if other.0 == f64::NEG_INFINITY {
                Ordering::Equal
            } else {
                Ordering::Less
            }
        } else if other.0 == f64::NEG_INFINITY {
            Ordering::Greater
        } else {
            self.0.partial_cmp(&other.0).unwrap()
        }
    }
}

impl Sum for FractionF64 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |sum, f| &sum + &f)
    }
}

impl<'a> Sum<&'a FractionF64> for FractionF64 {
    fn sum<I: Iterator<Item = &'a FractionF64>>(iter: I) -> Self {
        iter.fold(FractionF64::zero(), |sum, f| &sum + f)
    }
}

impl Neg for FractionF64 {
    type Output = FractionF64;

    fn neg(self) -> Self::Output {
        Self(self.0.neg())
    }
}

impl<'a> Neg for &'a FractionF64 {
    type Output = FractionF64;

    fn neg(self) -> Self::Output {
        FractionF64(self.0.neg())
    }
}

impl Add<&FractionF64> for &FractionF64 {
    type Output = FractionF64;

    fn add(self, rhs: &FractionF64) -> Self::Output {
        FractionF64(self.0.add(rhs.0))
    }
}

impl Add<FractionF64> for FractionF64 {
    type Output = FractionF64;

    fn add(self, rhs: FractionF64) -> Self::Output {
        FractionF64(self.0.add(rhs.0))
    }
}

impl<T> AddAssign<T> for FractionF64
where
    T: Borrow<FractionF64>,
{
    fn add_assign(&mut self, rhs: T) {
        let rhs = rhs.borrow();
        self.0.add_assign(rhs.0)
    }
}

impl Sub<&FractionF64> for &FractionF64 {
    type Output = FractionF64;

    fn sub(self, rhs: &FractionF64) -> Self::Output {
        FractionF64(self.0.sub(rhs.0))
    }
}

impl Sub<FractionF64> for FractionF64 {
    type Output = FractionF64;

    fn sub(self, rhs: FractionF64) -> Self::Output {
        FractionF64(self.0.sub(rhs.0))
    }
}

impl<T> SubAssign<T> for FractionF64
where
    T: Borrow<FractionF64>,
{
    fn sub_assign(&mut self, rhs: T) {
        let rhs = rhs.borrow();
        self.0.sub_assign(rhs.0)
    }
}

impl Mul<&FractionF64> for &FractionF64 {
    type Output = FractionF64;

    fn mul(self, rhs: &FractionF64) -> Self::Output {
        FractionF64(self.0.mul(rhs.0))
    }
}

impl Mul<FractionF64> for FractionF64 {
    type Output = FractionF64;

    fn mul(self, rhs: FractionF64) -> Self::Output {
        FractionF64(self.0.mul(rhs.0))
    }
}

impl<T> MulAssign<T> for FractionF64
where
    T: Borrow<FractionF64>,
{
    fn mul_assign(&mut self, rhs: T) {
        let rhs = rhs.borrow();
        self.0.mul_assign(rhs.0)
    }
}

impl Div<&FractionF64> for &FractionF64 {
    type Output = FractionF64;

    fn div(self, rhs: &FractionF64) -> Self::Output {
        FractionF64(self.0.div(rhs.0))
    }
}

impl Div<FractionF64> for FractionF64 {
    type Output = FractionF64;

    fn div(self, rhs: FractionF64) -> Self::Output {
        FractionF64(self.0.div(rhs.0))
    }
}

impl<T> DivAssign<T> for FractionF64
where
    T: Borrow<FractionF64>,
{
    fn div_assign(&mut self, rhs: T) {
        let rhs = rhs.borrow();
        self.0.div_assign(rhs.0)
    }
}

//======================== primitive types ========================//

impl Mul<f64> for FractionF64 {
    type Output = FractionF64;

    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl Div<f64> for FractionF64 {
    type Output = FractionF64;

    fn div(self, rhs: f64) -> Self::Output {
        Self(self.0 / rhs)
    }
}

impl Add<f64> for FractionF64 {
    type Output = FractionF64;

    fn add(self, rhs: f64) -> Self::Output {
        Self(self.0 + rhs)
    }
}

impl Sub<f64> for FractionF64 {
    type Output = FractionF64;

    fn sub(self, rhs: f64) -> Self::Output {
        Self(self.0 - rhs)
    }
}

macro_rules! from {
    ($t:ident) => {
        impl From<$t> for FractionF64 {
            fn from(value: $t) -> Self {
                Self(value as f64)
            }
        }
    };
}

macro_rules! from_signed {
    ($t:ident) => {
        impl From<$t> for FractionF64 {
            fn from(value: $t) -> Self {
                Self(value as f64)
            }
        }
    };
}

macro_rules! from_tuple_u_u {
    ($t:ident,$tt:ident) => {
        impl From<($t, $tt)> for FractionF64 {
            fn from(value: ($t, $tt)) -> Self {
                Self(value.0 as f64 / value.1 as f64)
            }
        }
    };
}

macro_rules! from_tuple_u_i {
    ($t:ident,$tt:ident) => {
        impl From<($t, $tt)> for FractionF64 {
            fn from(value: ($t, $tt)) -> Self {
                Self(value.0 as f64 / value.1 as f64)
            }
        }
    };
}

macro_rules! from_tuple_i_u {
    ($t:ident,$tt:ident) => {
        impl From<($t, $tt)> for FractionF64 {
            fn from(value: ($t, $tt)) -> Self {
                Self(value.0 as f64 / value.1 as f64)
            }
        }
    };
}

macro_rules! from_tuple_i_i {
    ($t:ident,$tt:ident) => {
        impl From<($t, $tt)> for FractionF64 {
            fn from(value: ($t, $tt)) -> Self {
                Self(value.0 as f64 / value.1 as f64)
            }
        }
    };
}

macro_rules! add {
    ($t:ident) => {
        impl<'a> Add<$t> for &'a FractionF64 {
            type Output = FractionF64;

            fn add(self, rhs: $t) -> Self::Output {
                let rhs: FractionF64 = rhs.into();
                self.add(&rhs)
            }
        }
    };
}

macro_rules! add_assign {
    ($t:ident) => {
        impl AddAssign<$t> for FractionF64 {
            fn add_assign(&mut self, rhs: $t) {
                let rhs: FractionF64 = rhs.into();
                self.add_assign(rhs)
            }
        }
    };
}

macro_rules! sub {
    ($t:ident) => {
        impl<'a> Sub<$t> for &'a FractionF64 {
            type Output = FractionF64;

            fn sub(self, rhs: $t) -> Self::Output {
                let rhs: FractionF64 = rhs.into();
                self.sub(&rhs)
            }
        }
    };
}

macro_rules! sub_assign {
    ($t:ident) => {
        impl SubAssign<$t> for FractionF64 {
            fn sub_assign(&mut self, rhs: $t) {
                let rhs: FractionF64 = rhs.into();
                self.sub_assign(rhs)
            }
        }
    };
}

macro_rules! mul {
    ($t:ident) => {
        impl<'a> Mul<$t> for &'a FractionF64 {
            type Output = FractionF64;

            fn mul(self, rhs: $t) -> Self::Output {
                let rhs: FractionF64 = rhs.into();
                self.mul(&rhs)
            }
        }
    };
}

macro_rules! mul_assign {
    ($t:ident) => {
        impl MulAssign<$t> for FractionF64 {
            fn mul_assign(&mut self, rhs: $t) {
                let rhs: FractionF64 = rhs.into();
                self.mul_assign(rhs)
            }
        }
    };
}

macro_rules! div {
    ($t:ident) => {
        impl<'a> Div<$t> for &'a FractionF64 {
            type Output = FractionF64;

            fn div(self, rhs: $t) -> Self::Output {
                let rhs: FractionF64 = rhs.into();
                self.div(&rhs)
            }
        }
    };
}

macro_rules! div_assign {
    ($t:ident) => {
        impl DivAssign<$t> for FractionF64 {
            fn div_assign(&mut self, rhs: $t) {
                let rhs: FractionF64 = rhs.into();
                self.div_assign(rhs)
            }
        }
    };
}

macro_rules! ttype_tuple {
    ($t:ident) => {
        from_tuple_u_u!($t, usize);
        from_tuple_u_u!($t, u128);
        from_tuple_u_u!($t, u64);
        from_tuple_u_u!($t, u32);
        from_tuple_u_u!($t, u16);
        from_tuple_u_u!($t, u8);
        from_tuple_u_i!($t, i128);
        from_tuple_u_i!($t, i64);
        from_tuple_u_i!($t, i32);
        from_tuple_u_i!($t, i16);
        from_tuple_u_i!($t, i8);
    };
}

macro_rules! ttype_tuple_signed {
    ($t:ident) => {
        from_tuple_i_u!($t, usize);
        from_tuple_i_u!($t, u128);
        from_tuple_i_u!($t, u64);
        from_tuple_i_u!($t, u32);
        from_tuple_i_u!($t, u16);
        from_tuple_i_u!($t, u8);
        from_tuple_i_i!($t, i64);
        from_tuple_i_i!($t, i32);
        from_tuple_i_i!($t, i16);
        from_tuple_i_i!($t, i8);
    };
}

macro_rules! ttype {
    ($t:ident) => {
        from!($t);
        ttype_tuple!($t);
        add!($t);
        add_assign!($t);
        sub!($t);
        sub_assign!($t);
        mul!($t);
        mul_assign!($t);
        div!($t);
        div_assign!($t);
    };
}

macro_rules! ttype_signed {
    ($t:ident) => {
        from_signed!($t);
        ttype_tuple_signed!($t);
        add!($t);
        add_assign!($t);
        sub!($t);
        sub_assign!($t);
        mul!($t);
        mul_assign!($t);
        div!($t);
        div_assign!($t);
    };
}

ttype!(usize);
ttype!(u128);
ttype!(u64);
ttype!(u32);
ttype!(u16);
ttype!(u8);
ttype_signed!(i128);
ttype_signed!(i64);
ttype_signed!(i32);
ttype_signed!(i16);
ttype_signed!(i8);

#[cfg(test)]
mod tests {
    use std::ops::Neg;

    use crate::{
        ebi_number::{One, Signed, Zero},
        fraction_f64::FractionF64,
    };

    #[test]
    fn fraction_neg() {
        let one = FractionF64::one();
        assert!(one.is_positive());
        let one = one.neg();
        assert!(one.is_negative());
    }

    #[test]
    fn fraction_exact() {
        let zero = FractionF64::one().one_minus();

        assert!(zero.is_zero());
    }
}
