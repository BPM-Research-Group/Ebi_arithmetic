use crate::{
    Zero, fraction::fraction_exact::FractionExact, log_polynomial::log_polynomial_exact::LogPolynomialExact
};
use std::ops::MulAssign;

impl MulAssign<FractionExact> for LogPolynomialExact {
    fn mul_assign(&mut self, rhs: FractionExact) {
        if rhs.is_zero() {
            self.argument2coefficient.clear();
        } else {
            self.argument2coefficient
                .iter_mut()
                .for_each(|(_, coefficient)| *coefficient *= &rhs.0);
        }
    }
}

impl MulAssign<&FractionExact> for LogPolynomialExact {
    fn mul_assign(&mut self, rhs: &FractionExact) {
        if rhs.is_zero() {
            self.argument2coefficient.clear();
        } else {
            self.argument2coefficient
                .iter_mut()
                .for_each(|(_, coefficient)| *coefficient *= &rhs.0);
        }
    }
}