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
use fnv::FnvBuildHasher;
use malachite::{Natural, Rational, base::num::basic::traits::Two};
use std::collections::HashMap;

impl From<FractionExact> for LogPolynomialExact {
    fn from(value: FractionExact) -> Self {
        if value.is_zero() {
            Self::zero()
        } else {
            let mut argument2coefficient = HashMap::<_, _, FnvBuildHasher>::default();
            argument2coefficient.insert(Natural::TWO, value.0);
            Self {
                argument2coefficient,
            }
        }
    }
}

impl From<&FractionExact> for LogPolynomialExact {
    fn from(value: &FractionExact) -> Self {
        if value.is_zero() {
            Self::zero()
        } else {
            let mut argument2coefficient = HashMap::<_, _, FnvBuildHasher>::default();
            argument2coefficient.insert(Natural::TWO, value.0.clone());
            Self {
                argument2coefficient,
            }
        }
    }
}

impl From<&Rational> for LogPolynomialExact {
    fn from(value: &Rational) -> Self {
        if value.is_zero() {
            Self::zero()
        } else {
            let mut argument2coefficient = HashMap::<_, _, FnvBuildHasher>::default();
            argument2coefficient.insert(Natural::TWO, value.clone());
            Self {
                argument2coefficient,
            }
        }
    }
}

macro_rules! from_rational {
    ($t:ty) => {
        impl From<$t> for LogPolynomialExact {
            fn from(value: $t) -> Self {
                if value.is_zero() {
                    Self::zero()
                } else {
                    let mut argument2coefficient = HashMap::<_, _, FnvBuildHasher>::default();
                    argument2coefficient.insert(Natural::TWO, value.into());
                    Self {
                        argument2coefficient,
                    }
                }
            }
        }
    };
}
from_rational!(Rational);
from_rational!(usize);
from_rational!(u128);
from_rational!(u64);
from_rational!(u32);
from_rational!(u16);
from_rational!(u8);
from_rational!(i128);
from_rational!(i64);
from_rational!(i32);
from_rational!(i16);
from_rational!(i8);

impl From<FractionF64> for LogPolynomialF64 {
    fn from(value: FractionF64) -> Self {
        Self(value.0)
    }
}

impl From<f64> for LogPolynomialF64 {
    fn from(value: f64) -> Self {
        Self(value)
    }
}

impl From<&f64> for LogPolynomialF64 {
    fn from(value: &f64) -> Self {
        Self(*value)
    }
}

impl From<FractionEnum> for LogPolynomialEnum {
    fn from(value: FractionEnum) -> Self {
        match value {
            FractionEnum::Exact(f) => LogPolynomialEnum::Exact(LogPolynomialExact::from(f)),
            FractionEnum::Approx(f) => LogPolynomialEnum::Approx(LogPolynomialF64::from(f)),
            FractionEnum::CannotCombineExactAndApprox => todo!(),
        }
    }
}
