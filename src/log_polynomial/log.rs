use crate::{Log, Zero, log_polynomial::log_polynomial_exact::LogPolynomialExact};
use anyhow::{Result, anyhow};
use malachite::{Natural, Rational, base::num::{basic::traits::One, factorization::traits::ExpressAsPower}};

impl<T> Log<T> for LogPolynomialExact where T: Into<Natural> {
    fn log(argument: T) -> Result<Self>
    where
        Self: Sized,
    {
        let mut arg = argument.into();
        if arg.is_zero() {
            return Err(anyhow!("Cannot take a logarithm of 0."));
        }
        let mut result = Self::zero();

        let mut factor = Rational::ONE;

        //if argument = x^n, then use that log(x^n) = n log (x)
        if let Some((root, exponent)) = arg.express_as_power() {
            factor *= Rational::from(exponent);
            arg = root;
        }

        todo!("prime factorisation");

        result.add_raw(factor, arg);

        Ok(result)
    }
}
