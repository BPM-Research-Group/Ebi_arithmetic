use crate::{Zero, log_polynomial::log_polynomial_exact::LogPolynomialExact};
use fnv::FnvBuildHasher;
use malachite::{Natural, Rational, base::num::basic::traits::Two};
use std::collections::HashMap;

impl<T> From<T> for LogPolynomialExact
where
    T: Into<Rational>,
{
    fn from(value: T) -> Self {
        let value = value.into();
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
