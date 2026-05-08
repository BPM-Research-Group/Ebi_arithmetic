use crate::{
    LogOf, Zero,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use anyhow::{Result, anyhow};
use malachite::{Natural, Rational, base::num::basic::traits::One};
use prime_factorization::Factorization;

impl LogOf<Natural> for LogPolynomialExact {
    fn log_of(argument: &Natural) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        if argument.is_zero() {
            return Err(anyhow!("Cannot take a logarithm of 0."));
        }

        //Prime factorisation library does not support larger numbers than u128.
        //For now, give an error.
        let arg = match u128::try_from(argument) {
            Ok(arg) => arg,
            Err(_) => {
                return Err(anyhow!(
                    "Number {argument} is too large to take logarithm of."
                ));
            }
        };

        let mut result = LogPolynomialExact::zero();

        //perform the prime factorisation
        for argument in Factorization::run(arg).factors {
            result.add_raw(Rational::ONE, argument.into());
        }

        Ok(result)
    }

    fn n_log_n_of(argument: &Natural) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        if argument.is_zero() {
            return Err(anyhow!("Cannot take a logarithm of 0."));
        }

        let coefficient: Rational = argument.clone().into();

        //Prime factorisation library does not support larger numbers than u128.
        //For now, give an error.
        let arg = match u128::try_from(argument) {
            Ok(arg) => arg,
            Err(_) => {
                return Err(anyhow!(
                    "Number {argument} is too large to take logarithm of."
                ));
            }
        };

        let mut result = LogPolynomialExact::zero();

        //perform the prime factorisation
        for argument in Factorization::run(arg).factors {
            result.add_raw(coefficient.clone(), argument.into());
        }

        Ok(result)
    }
}

impl LogOf<FractionExact> for LogPolynomialExact {
    fn log_of(argument: &FractionExact) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        let (x, y) = argument.0.numerator_and_denominator_ref();
        let mut result = LogPolynomialExact::log_of(x)?;
        result -= LogPolynomialExact::log_of(y)?;
        Ok(result)
    }

    fn n_log_n_of(argument: &FractionExact) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        let (x, y) = argument.0.numerator_and_denominator_ref();
        let mut result = LogPolynomialExact::log_of(x)?;
        result -= LogPolynomialExact::log_of(y)?;
        result *= argument;

        Ok(result)
    }
}

impl LogOf<FractionF64> for LogPolynomialF64 {
    fn log_of(argument: &FractionF64) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self(argument.0.log2()))
    }

    fn n_log_n_of(argument: &FractionF64) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self(argument.0 * argument.0.log2()))
    }
}

impl LogOf<FractionEnum> for LogPolynomialEnum {
    fn log_of(argument: &FractionEnum) -> Result<Self>
    where
        Self: Sized,
    {
        match argument {
            FractionEnum::Exact(f) => LogPolynomialExact::log_of(f),
            FractionEnum::Approx(f) => Ok(LogPolynomialEnum::Approx(LogPolynomialF64::from(f))),
            FractionEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }

    fn n_log_n_of(argument: &FractionEnum) -> Result<Self>
    where
        Self: Sized,
    {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use malachite::Natural;

    use crate::{
        LogOf,
        fraction::{approximate::Approximate, fraction_exact::FractionExact},
        log_polynomial::log_polynomial_exact::LogPolynomialExact,
    };

    #[test]
    fn log() {
        let a = LogPolynomialExact::log_of(&Natural::from(4usize)).unwrap();
        let b = LogPolynomialExact::from(2);

        println!("{:?}", a);
        assert_eq!(a, b);

        println!("{:?}", a.approximate());
    }

    #[test]
    fn nlogn() {
        let a = LogPolynomialExact::n_log_n_of(&Natural::from(4usize)).unwrap();
        let b = LogPolynomialExact::from(8);

        assert_eq!(a, b);
    }

    #[test]
    fn log_fraction() {
        let half = FractionExact::from((1, 2));
        let nlogn = LogPolynomialExact::n_log_n_of(&half).unwrap();
        let result = LogPolynomialExact::from(-FractionExact::from((1, 2)));

        assert_eq!(nlogn, result);

        assert_eq!(nlogn.approximate().unwrap(), -0.5);
    }
}
