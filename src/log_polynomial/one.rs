use crate::{
    One, is_exact_globally,
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use fnv::FnvBuildHasher;
use malachite::{
    Natural, Rational,
    base::num::basic::traits::{One as _, Two},
};
use std::collections::HashMap;

impl One for LogPolynomialExact {
    fn one() -> Self {
        let mut argument2coefficient = HashMap::<_, _, FnvBuildHasher>::default();
        argument2coefficient.insert(Natural::TWO, Rational::ONE);
        Self {
            argument2coefficient,
        }
    }

    fn is_one(&self) -> bool {
        if self.argument2coefficient.len() == 1 {
            if let Some((argument, coefficient)) = self.argument2coefficient.iter().next() {
                return coefficient.is_one() && argument == &Natural::TWO;
            }
        }
        false
    }
}

impl One for LogPolynomialF64 {
    fn one() -> Self {
        Self(1.0)
    }

    fn is_one(&self) -> bool {
        self.0.is_one()
    }
}

impl One for LogPolynomialEnum {
    fn one() -> Self {
        if is_exact_globally() {
            Self::Exact(LogPolynomialExact::one())
        } else {
            Self::Approx(LogPolynomialF64::one())
        }
    }

    fn is_one(&self) -> bool {
        match self {
            LogPolynomialEnum::Approx(log_polynomial_f64) => log_polynomial_f64.is_one(),
            LogPolynomialEnum::Exact(log_polynomial_exact) => log_polynomial_exact.is_one(),
            LogPolynomialEnum::CannotCombineExactAndApprox => false,
        }
    }
}
