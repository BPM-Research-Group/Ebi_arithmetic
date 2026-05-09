use crate::{
    exact::MaybeExact,
    is_exact_globally,
    matrix::{
        fraction_matrix_enum::FractionMatrixEnum, fraction_matrix_exact::FractionMatrixExact,
        fraction_matrix_f64::FractionMatrixF64,
    },
};
use anyhow::{Result, anyhow};

impl MaybeExact for FractionMatrixF64 {
    type Approximate = FractionMatrixF64;
    type Exact = FractionMatrixExact;

    fn is_exact(&self) -> bool {
        false
    }

    fn approx_ref(&self) -> anyhow::Result<&Self::Approximate> {
        Ok(self)
    }

    fn exact_ref(&self) -> anyhow::Result<&Self::Exact> {
        Err(anyhow!("cannot extract a fraction from a float"))
    }

    fn approx(self) -> anyhow::Result<Self::Approximate> {
        Ok(self)
    }

    fn exact(self) -> anyhow::Result<Self::Exact> {
        Err(anyhow!("cannot extract a fraction from a float"))
    }

    fn try_to_exact(_: Self::Exact) -> Result<Self> {
        Err(anyhow!("cannot put fraction in a float"))
    }

    fn try_to_approx(approx: Self::Approximate) -> Result<Self> {
        Ok(approx)
    }
}

impl MaybeExact for FractionMatrixExact {
    type Approximate = FractionMatrixF64;

    type Exact = FractionMatrixExact;

    fn is_exact(&self) -> bool {
        true
    }

    fn approx_ref(&self) -> anyhow::Result<&Self::Approximate> {
        Err(anyhow!("cannot extract a float from a fraction"))
    }

    fn exact_ref(&self) -> anyhow::Result<&Self::Exact> {
        Ok(self)
    }

    fn approx(self) -> anyhow::Result<Self::Approximate> {
        Err(anyhow!("cannot extract a float from a fraction"))
    }

    fn exact(self) -> anyhow::Result<Self::Exact> {
        Ok(self)
    }

    fn try_to_exact(exact: Self::Exact) -> Result<Self> {
        Ok(exact)
    }

    fn try_to_approx(_: Self::Approximate) -> Result<Self> {
        Err(anyhow!("cannot put float in a fraction"))
    }
}

impl MaybeExact for FractionMatrixEnum {
    type Approximate = FractionMatrixF64;

    type Exact = FractionMatrixExact;

    fn is_exact(&self) -> bool {
        true
    }

    fn approx_ref(&self) -> anyhow::Result<&Self::Approximate> {
        match self {
            FractionMatrixEnum::Approx(f) => Ok(f),
            FractionMatrixEnum::Exact(_) => Err(anyhow!("cannot extract a fraction from a float")),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine exact and approximate arithmetic"))
            }
        }
    }

    fn exact_ref(&self) -> anyhow::Result<&FractionMatrixExact> {
        match self {
            FractionMatrixEnum::Approx(_) => Err(anyhow!("cannot extract a float from a fraction")),
            FractionMatrixEnum::Exact(f) => Ok(f),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine exact and approximate arithmetic"))
            }
        }
    }

    fn approx(self) -> anyhow::Result<Self::Approximate> {
        match self {
            FractionMatrixEnum::Approx(f) => Ok(f),
            FractionMatrixEnum::Exact(_) => Err(anyhow!("cannot extract a fraction from a float")),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine exact and approximate arithmetic"))
            }
        }
    }

    fn exact(self) -> anyhow::Result<FractionMatrixExact> {
        match self {
            FractionMatrixEnum::Approx(_) => Err(anyhow!("cannot extract a float from a fraction")),
            FractionMatrixEnum::Exact(f) => Ok(f),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine exact and approximate arithmetic"))
            }
        }
    }

    fn try_to_exact(exact: <FractionMatrixEnum as MaybeExact>::Exact) -> Result<Self> {
        if is_exact_globally() {
            Ok(FractionMatrixEnum::Exact(exact))
        } else {
            Err(anyhow!("cannot put float in a fraction"))
        }
    }

    fn try_to_approx(approx: Self::Approximate) -> Result<Self> {
        if !is_exact_globally() {
            Ok(FractionMatrixEnum::Approx(approx))
        } else {
            Err(anyhow!("cannot put fraction in a float"))
        }
    }
}
