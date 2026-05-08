use crate::{
    Zero,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
    is_exact_globally,
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use fnv::FnvBuildHasher;
use malachite::{Natural, Rational, base::num::basic::traits::Two};
use std::collections::HashMap;

impl From<Rational> for LogPolynomialExact {
    fn from(value: Rational) -> Self {
        if value.is_zero() {
            Self::zero()
        } else {
            let mut argument2coefficient = HashMap::<_, _, FnvBuildHasher>::default();
            argument2coefficient.insert(Natural::TWO, value);
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

impl From<FractionExact> for LogPolynomialExact {
    fn from(value: FractionExact) -> Self {
        value.0.into()
    }
}

impl From<&FractionExact> for LogPolynomialExact {
    fn from(value: &FractionExact) -> Self {
        (&value.0).into()
    }
}

macro_rules! from_primitive {
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

        impl From<$t> for LogPolynomialF64 {
            fn from(value: $t) -> Self {
                Self(value as f64)
            }
        }

        impl From<$t> for LogPolynomialEnum {
            fn from(value: $t) -> Self {
                if is_exact_globally() {
                    Self::Exact(value.into())
                } else {
                    Self::Approx(value.into())
                }
            }
        }
    };
}
from_primitive!(usize);
from_primitive!(u128);
from_primitive!(u64);
from_primitive!(u32);
from_primitive!(u16);
from_primitive!(u8);
from_primitive!(i128);
from_primitive!(i64);
from_primitive!(i32);
from_primitive!(i16);
from_primitive!(i8);

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
