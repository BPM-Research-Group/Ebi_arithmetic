use crate::fraction::fraction_f64::FractionF64;
use std::fmt::{Debug, Display};

#[derive(Clone)]
pub struct LogPolynomialF64(pub(crate) f64);

impl Display for LogPolynomialF64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}

impl Debug for LogPolynomialF64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(&self.0, f)
    }
}

impl PartialEq for LogPolynomialF64 {
    fn eq(&self, other: &Self) -> bool {
        FractionF64::from(self.0) == FractionF64::from(other.0)
    }
}

impl Eq for LogPolynomialF64 {}
