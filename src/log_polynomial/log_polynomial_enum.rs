use crate::{
    fraction::fraction_enum::FractionEnum,
    log_polynomial::{
        log_polynomial_exact::LogPolynomialExact, log_polynomial_f64::LogPolynomialF64,
    },
};
use anyhow::{Result, anyhow};
use std::{
    fmt::{Debug, Display},
    io::Write,
};

#[derive(Clone, Eq, PartialEq)]
pub enum LogPolynomialEnum {
    Approx(LogPolynomialF64),
    Exact(LogPolynomialExact),
    CannotCombineExactAndApprox,
}

impl LogPolynomialEnum {
    pub fn export(&self, f: &mut dyn Write) -> Result<()> {
        match self {
            LogPolynomialEnum::Approx(lp) => lp.export(f),
            LogPolynomialEnum::Exact(lp) => lp.export(f),
            LogPolynomialEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }

    /**
     * Returns whether the two given objects are either both exact or both approximate
     */
    pub(crate) fn matches_fraction(&self, rhs: &FractionEnum) -> bool {
        match (self, rhs) {
            (LogPolynomialEnum::Exact(_), FractionEnum::Exact(_)) => true,
            (LogPolynomialEnum::Approx(_), FractionEnum::Approx(_)) => true,
            _ => false,
        }
    }

    /**
     * Returns whether the two given objects are either both exact or both approximate
     */
    pub(crate) fn matches(&self, rhs: &LogPolynomialEnum) -> bool {
        match (self, rhs) {
            (LogPolynomialEnum::Exact(_), LogPolynomialEnum::Exact(_)) => true,
            (LogPolynomialEnum::Approx(_), LogPolynomialEnum::Approx(_)) => true,
            _ => false,
        }
    }
}

impl Display for LogPolynomialEnum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LogPolynomialEnum::Approx(lp) => std::fmt::Display::fmt(&lp, f),
            LogPolynomialEnum::Exact(lp) => std::fmt::Display::fmt(&lp, f),
            LogPolynomialEnum::CannotCombineExactAndApprox => {
                writeln!(f, "Cannot combine exact and approximate arithmetic.")
            }
        }
    }
}

impl Debug for LogPolynomialEnum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Approx(arg0) => f.debug_tuple("Approx").field(arg0).finish(),
            Self::Exact(arg0) => f.debug_tuple("Exact").field(arg0).finish(),
            Self::CannotCombineExactAndApprox => write!(f, "CannotCombineExactAndApprox"),
        }
    }
}
