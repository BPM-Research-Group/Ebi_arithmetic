use anyhow::{Result, anyhow};
use malachite::rational::Rational;

use crate::{
    exact::MaybeExact,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
};

impl MaybeExact for FractionF64 {
    type Approximate = f64;
    type Exact = Rational;

    fn is_exact(&self) -> bool {
        false
    }

    fn extract_approx(&self) -> Result<&f64> {
        Ok(&self.0)
    }

    fn extract_exact(&self) -> Result<&Self::Exact> {
        Err(anyhow!("cannot extract a fraction from a float"))
    }
}

impl MaybeExact for FractionExact {
    type Approximate = f64;
    type Exact = Rational;

    fn is_exact(&self) -> bool {
        true
    }

    fn extract_approx(&self) -> Result<&f64> {
        Err(anyhow!("cannot extract a float from a fraction"))
    }

    /**
     * This is a low-level function to extract an f64. Only use if you are sure that the fraction is exact.
     * May not be available in all compilation modes.
     */
    fn extract_exact(&self) -> Result<&Rational> {
        Ok(&self.0)
    }
}

impl MaybeExact for FractionEnum {
    type Approximate = f64;
    type Exact = Rational;

    fn is_exact(&self) -> bool {
        match self {
            FractionEnum::Exact(_) => true,
            FractionEnum::Approx(_) => false,
            FractionEnum::CannotCombineExactAndApprox => false,
        }
    }

    fn extract_approx(&self) -> Result<&f64> {
        match self {
            FractionEnum::Exact(_) => Err(anyhow!("cannot extract a float from a fraction")),
            FractionEnum::Approx(f) => Ok(f),
            FractionEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine exact and approximate arithmetic"))
            }
        }
    }

    fn extract_exact(&self) -> Result<&Rational> {
        match self {
            FractionEnum::Exact(f) => Ok(f),
            FractionEnum::Approx(_) => Err(anyhow!("cannot extract a fraction from a float")),
            FractionEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine exact and approximate arithmetic"))
            }
        }
    }
}

macro_rules! approx {
    ($t: ident) => {
        impl MaybeExact for $t {
            type Approximate = $t;
            type Exact = ();

            fn is_exact(&self) -> bool {
                false
            }

            fn extract_approx(&self) -> Result<&Self::Approximate> {
                Ok(self)
            }

            fn extract_exact(&self) -> Result<&Self::Exact> {
                Err(anyhow!("cannot extract exact value"))
            }
        }
    };
}

macro_rules! exact {
    ($t: ident) => {
        impl MaybeExact for $t {
            type Approximate = ();
            type Exact = $t;

            fn is_exact(&self) -> bool {
                true
            }

            fn extract_approx(&self) -> Result<&Self::Approximate> {
                Err(anyhow!("cannot extract approximate value"))
            }

            fn extract_exact(&self) -> Result<&Self::Exact> {
                Ok(self)
            }
        }
    };
}

exact!(i128);
exact!(i64);
exact!(i32);
exact!(i16);
exact!(u128);
exact!(u64);
exact!(u32);
exact!(u16);
exact!(u8);
approx!(f64);
approx!(f32);
