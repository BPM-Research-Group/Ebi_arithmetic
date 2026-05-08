use crate::{Log, Zero, log_polynomial::log_polynomial_exact::LogPolynomialExact};
use anyhow::{Result, anyhow};
use malachite::{Natural, Rational, base::num::basic::traits::One};
use prime_factorization::Factorization;

impl<T> Log<T> for LogPolynomialExact
where
    T: Into<Natural>,
{
    fn log(argument: T) -> Result<Self>
    where
        Self: Sized,
    {
        let arg = argument.into();
        if arg.is_zero() {
            return Err(anyhow!("Cannot take a logarithm of 0."));
        }

        //Prime factorisation library does not support larger numbers than u128.
        //For now, give an error.
        let arg = match u128::try_from(&arg) {
            Ok(arg) => arg,
            Err(_) => return Err(anyhow!("Number {arg} is too large to take logarithm of.")),
        };

        let mut result = Self::zero();

        //perform the prime factorisation
        for argument in Factorization::run(arg).factors {
            result.add_raw(Rational::ONE, argument.into());
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use crate::{Log, fraction::approximate::Approximate, log_polynomial::log_polynomial_exact::LogPolynomialExact};

    #[test]
    fn log() {
        let a = LogPolynomialExact::log(4usize).unwrap();
        let b = LogPolynomialExact::from(2);

        println!("{:?}", a);
        assert_eq!(a, b);

        println!("{:?}", a.approximate());
    }
}
