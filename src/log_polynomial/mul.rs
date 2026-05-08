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