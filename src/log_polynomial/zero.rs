use crate::{
    Zero, is_exact_globally,
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use fnv::FnvBuildHasher;
use std::collections::HashMap;

impl Zero for LogPolynomialExact {
    fn zero() -> Self {
        Self {
            argument2coefficient: HashMap::<_, _, FnvBuildHasher>::default(),
        }
    }

    fn is_zero(&self) -> bool {
        self.argument2coefficient.is_empty()
    }
}

impl Zero for LogPolynomialF64 {
    fn zero() -> Self {
        Self(0.0)
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl Zero for LogPolynomialEnum {
    fn zero() -> Self {
        if is_exact_globally() {
            Self::Exact(LogPolynomialExact::zero())
        } else {
            Self::Approx(LogPolynomialF64::zero())
        }
    }

    fn is_zero(&self) -> bool {
        match self {
            LogPolynomialEnum::Approx(log_polynomial_f64) => log_polynomial_f64.is_zero(),
            LogPolynomialEnum::Exact(log_polynomial_exact) => log_polynomial_exact.is_zero(),
            LogPolynomialEnum::CannotCombineExactAndApprox => false
        }
    }
}
