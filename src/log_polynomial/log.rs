use crate::{
    Signed, Zero,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
    log::LogOf,
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use anyhow::{Result, anyhow};
use malachite::{Natural, Rational, base::num::basic::traits::One};
use prime_factorization::Factorization;

// ===== exact =====

impl LogOf<u128> for LogPolynomialExact {
    fn log_of(argument: u128) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        if argument.is_not_positive() {
            return Err(anyhow!("Can only take a logarithm of a positive number."));
        }

        let mut result = LogPolynomialExact::zero();

        //perform the prime factorisation
        for factor in Factorization::run(argument).factors {
            result.add_raw(Rational::ONE, factor.into());
        }

        Ok(result)
    }

    fn n_log_n_of(argument: u128) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        if argument.is_not_positive() {
            return Err(anyhow!("Can only take a logarithm of a positive number."));
        }

        let coefficient: Rational = argument.clone().into();

        let mut result = LogPolynomialExact::zero();

        //perform the prime factorisation
        for factor in Factorization::run(argument).factors {
            result.add_raw(coefficient.clone(), factor.into());
        }

        Ok(result)
    }
}

impl LogOf<&Natural> for LogPolynomialExact {
    fn log_of(argument: &Natural) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        //Prime factorisation library does not support larger numbers than u128.
        //For now, give an error.
        match u128::try_from(argument) {
            Ok(arg) => LogPolynomialExact::log_of(arg),
            Err(_) => Err(anyhow!(
                "Number {argument} is too large to take logarithm of."
            )),
        }
    }

    fn n_log_n_of(argument: &Natural) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        //Prime factorisation library does not support larger numbers than u128.
        //For now, give an error.
        match u128::try_from(argument) {
            Ok(arg) => LogPolynomialExact::n_log_n_of(arg),
            Err(_) => Err(anyhow!(
                "Number {argument} is too large to take logarithm of."
            )),
        }
    }
}

impl LogOf<&Rational> for LogPolynomialExact {
    fn log_of(argument: &Rational) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        if argument.is_not_positive() {
            return Err(anyhow!("Can only take a logarithm of a positive number."));
        }
        let (x, y) = argument.numerator_and_denominator_ref();
        let mut result = LogPolynomialExact::log_of(x)?;
        result -= LogPolynomialExact::log_of(y)?;
        Ok(result)
    }

    fn n_log_n_of(argument: &Rational) -> Result<LogPolynomialExact>
    where
        Self: Sized,
    {
        if argument.is_not_positive() {
            return Err(anyhow!("Can only take a logarithm of a positive number."));
        }
        let (x, y) = argument.numerator_and_denominator_ref();
        let mut result = LogPolynomialExact::log_of(x)?;
        result -= LogPolynomialExact::log_of(y)?;
        result *= argument;

        Ok(result)
    }
}

impl LogOf<FractionExact> for LogPolynomialExact {
    fn log_of(argument: FractionExact) -> Result<LogPolynomialExact> {
        LogPolynomialExact::log_of(&argument.0)
    }

    fn n_log_n_of(argument: FractionExact) -> Result<LogPolynomialExact> {
        LogPolynomialExact::n_log_n_of(&argument.0)
    }
}

impl LogOf<&FractionExact> for LogPolynomialExact {
    fn log_of(argument: &FractionExact) -> Result<LogPolynomialExact> {
        LogPolynomialExact::log_of(&argument.0)
    }

    fn n_log_n_of(argument: &FractionExact) -> Result<LogPolynomialExact> {
        LogPolynomialExact::n_log_n_of(&argument.0)
    }
}

// ===== approximate =====

impl LogOf<f64> for LogPolynomialF64 {
    fn log_of(argument: f64) -> Result<Self>
    where
        Self: Sized,
    {
        if Signed::is_not_positive(&argument) {
            return Err(anyhow!("Can only take a logarithm of a positive number."));
        }
        Ok(Self(argument.log2()))
    }

    fn n_log_n_of(argument: f64) -> Result<Self>
    where
        Self: Sized,
    {
        if argument.is_not_positive() {
            return Err(anyhow!("Can only take a logarithm of a positive number."));
        }
        Ok(Self(argument * argument.log2()))
    }
}

impl LogOf<FractionF64> for LogPolynomialF64 {
    fn log_of(argument: FractionF64) -> Result<Self> {
        LogPolynomialF64::log_of(argument.0)
    }

    fn n_log_n_of(argument: FractionF64) -> Result<Self> {
        LogPolynomialF64::n_log_n_of(argument.0)
    }
}

impl LogOf<&FractionF64> for LogPolynomialF64 {
    fn log_of(argument: &FractionF64) -> Result<Self> {
        LogPolynomialF64::log_of(argument.0)
    }

    fn n_log_n_of(argument: &FractionF64) -> Result<Self> {
        LogPolynomialF64::n_log_n_of(argument.0)
    }
}

// ===== enum =====

impl LogOf<FractionEnum> for LogPolynomialEnum {
    fn log_of(argument: FractionEnum) -> Result<Self>
    where
        Self: Sized,
    {
        match argument {
            FractionEnum::Exact(f) => Ok(Self::Exact(LogPolynomialExact::log_of(&f)?)),
            FractionEnum::Approx(f) => Ok(Self::Approx(LogPolynomialF64::log_of(f)?)),
            FractionEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }

    fn n_log_n_of(argument: FractionEnum) -> Result<Self>
    where
        Self: Sized,
    {
        match argument {
            FractionEnum::Exact(f) => Ok(Self::Exact(LogPolynomialExact::n_log_n_of(&f)?)),
            FractionEnum::Approx(f) => Ok(Self::Approx(LogPolynomialF64::n_log_n_of(f)?)),
            FractionEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }
}

impl LogOf<&FractionEnum> for LogPolynomialEnum {
    fn log_of(argument: &FractionEnum) -> Result<Self> {
        match argument {
            FractionEnum::Exact(f) => Ok(Self::Exact(LogPolynomialExact::log_of(f)?)),
            FractionEnum::Approx(f) => Ok(Self::Approx(LogPolynomialF64::log_of(*f)?)),
            FractionEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }

    fn n_log_n_of(argument: &FractionEnum) -> Result<Self> {
        match argument {
            FractionEnum::Exact(f) => Ok(Self::Exact(LogPolynomialExact::n_log_n_of(f)?)),
            FractionEnum::Approx(f) => Ok(Self::Approx(LogPolynomialF64::n_log_n_of(*f)?)),
            FractionEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }
}

macro_rules! log_primitive {
    ($t:ty) => {
        impl LogOf<$t> for LogPolynomialExact {
            fn log_of(argument: $t) -> Result<Self> {
                if argument.is_not_positive() {
                    return Err(anyhow!("Can only take a logarithm of a positive number."));
                }
                LogPolynomialExact::log_of(argument as u128)
            }

            fn n_log_n_of(argument: $t) -> Result<Self> {
                if argument.is_not_positive() {
                    return Err(anyhow!("Can only take a logarithm of a positive number."));
                }
                LogPolynomialExact::n_log_n_of(argument as u128)
            }
        }

        impl LogOf<$t> for LogPolynomialF64 {
            fn log_of(argument: $t) -> Result<Self> {
                if argument.is_not_positive() {
                    return Err(anyhow!("Can only take a logarithm of a positive number."));
                }
                LogPolynomialF64::log_of(argument as f64)
            }

            fn n_log_n_of(argument: $t) -> Result<Self> {
                if argument.is_not_positive() {
                    return Err(anyhow!("Can only take a logarithm of a positive number."));
                }
                LogPolynomialF64::n_log_n_of(argument as f64)
            }
        }

        impl LogOf<$t> for LogPolynomialEnum {
            fn log_of(argument: $t) -> Result<Self> {
                if argument.is_not_positive() {
                    return Err(anyhow!("Can only take a logarithm of a positive number."));
                }
                let f = FractionEnum::from(argument);
                LogPolynomialEnum::log_of(f)
            }

            fn n_log_n_of(argument: $t) -> Result<Self> {
                if argument.is_not_positive() {
                    return Err(anyhow!("Can only take a logarithm of a positive number."));
                }
                let f = FractionEnum::from(argument);
                LogPolynomialEnum::n_log_n_of(f)
            }
        }
    };
}
log_primitive!(usize);
log_primitive!(u64);
log_primitive!(u32);
log_primitive!(u16);
log_primitive!(u8);
log_primitive!(i128);
log_primitive!(i64);
log_primitive!(i32);
log_primitive!(i16);
log_primitive!(i8);

#[cfg(test)]
mod tests {
    use crate::{
        fraction::{approximate::Approximate, fraction_exact::FractionExact},
        log::LogOf,
        log_polynomial::log_polynomial_exact::LogPolynomialExact,
    };
    use malachite::Natural;

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
