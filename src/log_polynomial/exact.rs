use crate::{
    MaybeExact,
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use anyhow::{Result, anyhow};

impl MaybeExact for LogPolynomialExact {
    type Approximate = LogPolynomialF64;
    type Exact = LogPolynomialExact;

    fn is_exact(&self) -> bool {
        true
    }

    fn approx_ref(&self) -> Result<&Self::Approximate> {
        Err(anyhow!("cannot extract a float from a fraction"))
    }

    fn exact_ref(&self) -> Result<&Self::Exact> {
        Ok(self)
    }

    fn approx(self) -> Result<Self::Approximate> {
        Err(anyhow!("cannot extract a float from a fraction"))
    }

    fn exact(self) -> Result<Self::Exact> {
        Ok(self)
    }
}

impl MaybeExact for LogPolynomialF64 {
    type Approximate = LogPolynomialF64;
    type Exact = LogPolynomialExact;

    fn is_exact(&self) -> bool {
        false
    }

    fn approx_ref(&self) -> Result<&Self::Approximate> {
        Ok(self)
    }

    fn exact_ref(&self) -> Result<&Self::Exact> {
        Err(anyhow!("cannot extract a fraction from a float"))
    }

    fn approx(self) -> Result<Self::Approximate> {
        Ok(self)
    }

    fn exact(self) -> Result<Self::Exact> {
        Err(anyhow!("cannot extract a fraction from a float"))
    }
}

impl MaybeExact for LogPolynomialEnum {
    type Approximate = LogPolynomialF64;
    type Exact = LogPolynomialExact;

    fn is_exact(&self) -> bool {
        match self {
            LogPolynomialEnum::Approx(_) => false,
            LogPolynomialEnum::Exact(_) => true,
            LogPolynomialEnum::CannotCombineExactAndApprox => false,
        }
    }

    fn approx_ref(&self) -> Result<&Self::Approximate> {
        match self {
            LogPolynomialEnum::Approx(f) => Ok(f),
            LogPolynomialEnum::Exact(_) => Err(anyhow!("cannot extract a fraction from a float")),
            LogPolynomialEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }

    fn exact_ref(&self) -> Result<&<LogPolynomialEnum as MaybeExact>::Exact> {
        match self {
            LogPolynomialEnum::Approx(_) => Err(anyhow!("cannot extract a float from a fraction")),
            LogPolynomialEnum::Exact(f) => Ok(f),
            LogPolynomialEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }

    fn approx(self) -> Result<Self::Approximate> {
        match self {
            LogPolynomialEnum::Approx(f) => Ok(f),
            LogPolynomialEnum::Exact(_) => Err(anyhow!("cannot extract a fraction from a float")),
            LogPolynomialEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }

    fn exact(self) -> Result<<LogPolynomialEnum as MaybeExact>::Exact> {
        match self {
            LogPolynomialEnum::Approx(_) => Err(anyhow!("cannot extract a float from a fraction")),
            LogPolynomialEnum::Exact(f) => Ok(f),
            LogPolynomialEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine exact and approximate arithmetic."))
            }
        }
    }
}
