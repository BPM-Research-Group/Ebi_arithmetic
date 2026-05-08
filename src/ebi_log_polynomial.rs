use crate::{
    MaybeExact, One, Zero, fraction::{
        approximate::Approximate, fraction_enum::FractionEnum, fraction_exact::FractionExact,
        fraction_f64::FractionF64,
    }, log::LogOf, log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    }
};
use std::{
    fmt::{Debug, Display},
    ops::{AddAssign, DivAssign, MulAssign, Neg, SubAssign},
};

pub trait EbiLogPolynomial<T>:
    MaybeExact
    + From<T>
    + From<usize>
    + From<u128>
    + From<u64>
    + From<u32>
    + From<u16>
    + From<u8>
    + From<i128>
    + From<i64>
    + From<i32>
    + From<i16>
    + From<i8>
    + LogOf<T>
    + for<'a> LogOf<&'a T>
    + MulAssign<T>
    + for<'a> MulAssign<&'a T>
    + MulAssign<usize>
    + MulAssign<u128>
    + MulAssign<u64>
    + MulAssign<u32>
    + MulAssign<u16>
    + MulAssign<u8>
    + MulAssign<i128>
    + MulAssign<i64>
    + MulAssign<i32>
    + MulAssign<i16>
    + MulAssign<i8>
    + DivAssign<T>
    + for<'a> DivAssign<&'a T>
    + DivAssign<usize>
    + DivAssign<u128>
    + DivAssign<u64>
    + DivAssign<u32>
    + DivAssign<u16>
    + DivAssign<u8>
    + DivAssign<i128>
    + DivAssign<i64>
    + DivAssign<i32>
    + DivAssign<i16>
    + DivAssign<i8>
    + AddAssign
    + for<'a> AddAssign<&'a Self>
    + AddAssign<T>
    + for<'a> AddAssign<&'a T>
    + AddAssign<usize>
    + AddAssign<u128>
    + AddAssign<u64>
    + AddAssign<u32>
    + AddAssign<u16>
    + AddAssign<u8>
    + AddAssign<i128>
    + AddAssign<i64>
    + AddAssign<i32>
    + AddAssign<i16>
    + AddAssign<i8>
    + SubAssign
    + for<'a> SubAssign<&'a Self>
    + SubAssign<T>
    + for<'a> SubAssign<&'a T>
    + SubAssign<usize>
    + SubAssign<u128>
    + SubAssign<u64>
    + SubAssign<u32>
    + SubAssign<u16>
    + SubAssign<u8>
    + SubAssign<i128>
    + SubAssign<i64>
    + SubAssign<i32>
    + SubAssign<i16>
    + SubAssign<i8>
    + Sized
    + Clone
    + Debug
    + Display
    + Eq
    + PartialEq
    + Approximate
    + Zero
    + One
    + Neg
{
}

impl EbiLogPolynomial<FractionEnum> for LogPolynomialEnum {}
impl EbiLogPolynomial<FractionExact> for LogPolynomialExact {}
impl EbiLogPolynomial<FractionF64> for LogPolynomialF64 {}
