use malachite::Rational;

use crate::{
    Zero,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use std::ops::MulAssign;

impl MulAssign<FractionExact> for LogPolynomialExact {
    fn mul_assign(&mut self, rhs: FractionExact) {
        if rhs.is_zero() {
            self.argument2coefficient.clear();
        } else {
            self.argument2coefficient
                .iter_mut()
                .for_each(|(_, coefficient)| *coefficient *= &rhs.0);
        }
    }
}

impl MulAssign<&FractionExact> for LogPolynomialExact {
    fn mul_assign(&mut self, rhs: &FractionExact) {
        if rhs.is_zero() {
            self.argument2coefficient.clear();
        } else {
            self.argument2coefficient
                .iter_mut()
                .for_each(|(_, coefficient)| *coefficient *= &rhs.0);
        }
    }
}

impl MulAssign<Rational> for LogPolynomialExact {
    fn mul_assign(&mut self, rhs: Rational) {
        if rhs.is_zero() {
            self.argument2coefficient.clear();
        } else {
            self.argument2coefficient
                .iter_mut()
                .for_each(|(_, coefficient)| *coefficient *= &rhs);
        }
    }
}

impl MulAssign<&Rational> for LogPolynomialExact {
    fn mul_assign(&mut self, rhs: &Rational) {
        if rhs.is_zero() {
            self.argument2coefficient.clear();
        } else {
            self.argument2coefficient
                .iter_mut()
                .for_each(|(_, coefficient)| *coefficient *= rhs);
        }
    }
}

impl MulAssign<FractionF64> for LogPolynomialF64 {
    fn mul_assign(&mut self, rhs: FractionF64) {
        self.0 *= rhs.0
    }
}

impl MulAssign<&FractionF64> for LogPolynomialF64 {
    fn mul_assign(&mut self, rhs: &FractionF64) {
        self.0 *= rhs.0
    }
}

impl MulAssign<f64> for LogPolynomialF64 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0 *= rhs
    }
}

impl MulAssign<&f64> for LogPolynomialF64 {
    fn mul_assign(&mut self, rhs: &f64) {
        self.0 *= rhs
    }
}

impl MulAssign<FractionEnum> for LogPolynomialEnum {
    fn mul_assign(&mut self, rhs: FractionEnum) {
        if self.matches_fraction(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), FractionEnum::Approx(f)) => lp.mul_assign(f),
                (LogPolynomialEnum::Exact(lp), FractionEnum::Exact(f)) => lp.mul_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

impl MulAssign<&FractionEnum> for LogPolynomialEnum {
    fn mul_assign(&mut self, rhs: &FractionEnum) {
        if self.matches_fraction(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), FractionEnum::Approx(f)) => lp.mul_assign(f),
                (LogPolynomialEnum::Exact(lp), FractionEnum::Exact(f)) => lp.mul_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

macro_rules! mulassign_primitive {
    ($t:ty) => {
        impl MulAssign<$t> for LogPolynomialExact {
            fn mul_assign(&mut self, rhs: $t) {
                self.mul_assign(Rational::from(rhs))
            }
        }

        impl MulAssign<$t> for LogPolynomialF64 {
            fn mul_assign(&mut self, rhs: $t) {
                self.mul_assign(rhs as f64)
            }
        }

        impl MulAssign<$t> for LogPolynomialEnum {
            fn mul_assign(&mut self, rhs: $t) {
                match self {
                    LogPolynomialEnum::Exact(lp) => lp.mul_assign(rhs),
                    LogPolynomialEnum::Approx(lp) => lp.mul_assign(rhs),
                    LogPolynomialEnum::CannotCombineExactAndApprox => {}
                }
            }
        }
    };
}
mulassign_primitive!(usize);
mulassign_primitive!(u128);
mulassign_primitive!(u64);
mulassign_primitive!(u32);
mulassign_primitive!(u16);
mulassign_primitive!(u8);
mulassign_primitive!(i128);
mulassign_primitive!(i64);
mulassign_primitive!(i32);
mulassign_primitive!(i16);
mulassign_primitive!(i8);