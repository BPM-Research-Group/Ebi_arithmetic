use crate::fraction::{
    fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
};

pub trait EbiNumber: Zero + One + Round + Clone {}

impl EbiNumber for FractionEnum {}
impl EbiNumber for FractionF64 {}
impl EbiNumber for FractionExact {}
impl EbiNumber for f32 {}
impl EbiNumber for f64 {}
impl EbiNumber for usize {}
impl EbiNumber for u128 {}
impl EbiNumber for u64 {}
impl EbiNumber for u32 {}
impl EbiNumber for u16 {}
impl EbiNumber for u8 {}
impl EbiNumber for i128 {}
impl EbiNumber for i64 {}
impl EbiNumber for i32 {}
impl EbiNumber for i16 {}
impl EbiNumber for i8 {}

pub trait Zero: Sized {
    fn zero() -> Self;

    fn set_zero(&mut self) {
        *self = Zero::zero();
    }

    fn is_zero(&self) -> bool;
}

pub trait One: Sized {
    fn one() -> Self;

    fn set_one(&mut self) {
        *self = One::one();
    }

    fn is_one(&self) -> bool;
}

pub trait Signed: Sized {
    fn abs(self) -> Self;

    /// Returns true if the number is positive and false if the number is zero or negative.
    fn is_positive(&self) -> bool;

    /// Returns true if the number is negative and false if the number is zero or positive.
    fn is_negative(&self) -> bool;

    /// For exact arithmetic: Returns true if the number is positive or zero.
    /// For approximate arithmetic: returns true if the number is larger than -epsilon
    fn is_not_negative(&self) -> bool {
        !self.is_negative()
    }

    /// For exact arithmetic: Returns true if the number is negative or zero.
    /// For approximate arithmetic: returns true if the number is smaller than epsilon
    fn is_not_positive(&self) -> bool {
        !self.is_positive()
    }
}

pub trait Round: Sized {
    /// Returns the largest integer less than or equal to `self`.
    fn floor(self) -> Self;

    /// Returns the smallest integer greater than or equal to `self`.
    fn ceil(self) -> Self;
}

pub trait Recip: Sized {
    /// Takes the reciprocal (inverse) of a number, `1/x`.
    fn recip(&self) -> Self;
}
