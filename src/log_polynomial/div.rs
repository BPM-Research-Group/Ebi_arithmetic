use malachite::Rational;

use crate::{
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use std::ops::DivAssign;

impl DivAssign<FractionExact> for LogPolynomialExact {
    fn div_assign(&mut self, rhs: FractionExact) {
        self.argument2coefficient
            .iter_mut()
            .for_each(|(_, coefficient)| *coefficient /= &rhs.0);
    }
}

impl DivAssign<&FractionExact> for LogPolynomialExact {
    fn div_assign(&mut self, rhs: &FractionExact) {
        self.argument2coefficient
            .iter_mut()
            .for_each(|(_, coefficient)| *coefficient /= &rhs.0);
    }
}

impl DivAssign<Rational> for LogPolynomialExact {
    fn div_assign(&mut self, rhs: Rational) {
        self.argument2coefficient
            .iter_mut()
            .for_each(|(_, coefficient)| *coefficient /= &rhs);
    }
}

impl DivAssign<&Rational> for LogPolynomialExact {
    fn div_assign(&mut self, rhs: &Rational) {
        self.argument2coefficient
            .iter_mut()
            .for_each(|(_, coefficient)| *coefficient /= rhs);
    }
}

impl DivAssign<FractionF64> for LogPolynomialF64 {
    fn div_assign(&mut self, rhs: FractionF64) {
        self.0 /= rhs.0
    }
}

impl DivAssign<&FractionF64> for LogPolynomialF64 {
    fn div_assign(&mut self, rhs: &FractionF64) {
        self.0 /= rhs.0
    }
}

impl DivAssign<f64> for LogPolynomialF64 {
    fn div_assign(&mut self, rhs: f64) {
        self.0 /= rhs
    }
}

impl DivAssign<&f64> for LogPolynomialF64 {
    fn div_assign(&mut self, rhs: &f64) {
        self.0 /= rhs
    }
}

impl DivAssign<FractionEnum> for LogPolynomialEnum {
    fn div_assign(&mut self, rhs: FractionEnum) {
        if self.matches_fraction(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), FractionEnum::Approx(f)) => lp.div_assign(f),
                (LogPolynomialEnum::Exact(lp), FractionEnum::Exact(f)) => lp.div_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

impl DivAssign<&FractionEnum> for LogPolynomialEnum {
    fn div_assign(&mut self, rhs: &FractionEnum) {
        if self.matches_fraction(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), FractionEnum::Approx(f)) => lp.div_assign(f),
                (LogPolynomialEnum::Exact(lp), FractionEnum::Exact(f)) => lp.div_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}
