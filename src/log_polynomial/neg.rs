use crate::log_polynomial::{
    log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
    log_polynomial_f64::LogPolynomialF64,
};
use malachite::base::num::arithmetic::traits::NegAssign;
use std::ops::Neg;

impl Neg for LogPolynomialExact {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.argument2coefficient
            .iter_mut()
            .for_each(|(_, y)| y.neg_assign());
        self
    }
}

impl Neg for LogPolynomialF64 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.neg())
    }
}

impl Neg for LogPolynomialEnum {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            LogPolynomialEnum::Approx(f) => LogPolynomialEnum::Approx(f.neg()),
            LogPolynomialEnum::Exact(f) => LogPolynomialEnum::Exact(f.neg()),
            LogPolynomialEnum::CannotCombineExactAndApprox => {
                LogPolynomialEnum::CannotCombineExactAndApprox
            }
        }
    }
}
