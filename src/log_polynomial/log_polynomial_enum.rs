use crate::{
    fraction::fraction_enum::FractionEnum,
    log_polynomial::{
        log_polynomial_exact::LogPolynomialExact, log_polynomial_f64::LogPolynomialF64,
    },
};

#[derive(Clone)]
pub enum LogPolynomialEnum {
    Approx(LogPolynomialF64),
    Exact(LogPolynomialExact),
    CannotCombineExactAndApprox,
}

impl LogPolynomialEnum {
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
