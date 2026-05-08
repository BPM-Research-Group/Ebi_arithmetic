use crate::{
    LogOf, MaybeExact, One, Zero,
    fraction::{
        approximate::Approximate, fraction_enum::FractionEnum, fraction_exact::FractionExact,
        fraction_f64::FractionF64,
    },
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

pub trait EbiLogPolynomial<T>:
    MaybeExact
    + From<T>
    + LogOf<T>
    + MulAssign<T>
    + for<'a> MulAssign<&'a T>
    + DivAssign<T>
    + for<'a> DivAssign<&'a T>
    + AddAssign
    + AddAssign<T>
    + for<'a> AddAssign<&'a T>
    + SubAssign
    + SubAssign<T>
    + for<'a> SubAssign<&'a T>
    + Sized
    + Clone
    + Approximate
    + Zero
    + One
{
}

impl EbiLogPolynomial<FractionEnum> for LogPolynomialEnum {}
impl EbiLogPolynomial<FractionExact> for LogPolynomialExact {}
impl EbiLogPolynomial<FractionF64> for LogPolynomialF64 {}
