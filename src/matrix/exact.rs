use crate::{
    exact::MaybeExact,
    matrix::{fraction_matrix_exact::FractionMatrixExact, fraction_matrix_f64::FractionMatrixF64},
};
use anyhow::anyhow;

impl MaybeExact for FractionMatrixF64 {
    type Approximate = FractionMatrixF64;
    type Exact = FractionMatrixExact;

    fn is_exact(&self) -> bool {
        false
    }

    fn extract_approx(&self) -> anyhow::Result<&Self::Approximate> {
        Ok(self)
    }

    fn extract_exact(&self) -> anyhow::Result<&Self::Exact> {
        Err(anyhow!("cannot extract a fraction from a float"))
    }
}

impl MaybeExact for FractionMatrixExact {
    type Approximate = FractionMatrixF64;

    type Exact = FractionMatrixExact;

    fn is_exact(&self) -> bool {
        true
    }

    fn extract_approx(&self) -> anyhow::Result<&Self::Approximate> {
        Err(anyhow!("cannot extract a float from a fraction"))
    }

    fn extract_exact(&self) -> anyhow::Result<&Self::Exact> {
        Ok(self)
    }
}
