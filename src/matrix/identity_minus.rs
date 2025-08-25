use anyhow::{Result, anyhow};
use malachite::{base::num::basic::traits::One, rational::Rational};

use crate::{
    ebi_matrix::EbiMatrix,
    matrix::{
        fraction_matrix_enum::FractionMatrixEnum, fraction_matrix_exact::FractionMatrixExact,
        fraction_matrix_f64::FractionMatrixF64,
    }, IdentityMinus,
};

impl IdentityMinus for FractionMatrixF64 {
    fn identity_minus(&mut self) -> Result<()> {
        for i in 0..self.number_of_rows() {
            for j in 0..self.number_of_columns() {
                if i == j {
                    let idx = self.index(i, i);
                    self.values[idx] = 1f64 - self.values[idx];
                } else {
                    let idx = self.index(i, j);
                    self.values[idx] *= -1f64;
                }
            }
        }
        Ok(())
    }
}

impl IdentityMinus for FractionMatrixExact {
    fn identity_minus(&mut self) -> Result<()> {
        for i in 0..self.number_of_rows() {
            for j in 0..self.number_of_columns() {
                if i == j {
                    let idx = self.index(i, i);
                    self.values[idx] = &Rational::ONE - &self.values[idx];
                } else {
                    let idx = self.index(i, j);
                    self.values[idx] *= -Rational::ONE;
                }
            }
        }
        Ok(())
    }
}

impl IdentityMinus for FractionMatrixEnum {
    fn identity_minus(&mut self) -> Result<()> {
        match self {
            FractionMatrixEnum::Approx(m) => m.identity_minus(),
            FractionMatrixEnum::Exact(m) => m.identity_minus(),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine approximate and exact arithmetic"))
            }
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        f_en,
        fraction::fraction_enum::FractionEnum,
        matrix::{fraction_matrix_enum::FractionMatrixEnum, identity_minus::IdentityMinus},
    };

    #[test]
    fn fraction_matrix_abnormal() {
        let mut m1: FractionMatrixEnum = vec![vec![f_en!(8, 3), f_en!(3, 8)]].try_into().unwrap();

        m1.identity_minus().unwrap();

        let m3: FractionMatrixEnum = vec![vec![-f_en!(5, 3), -f_en!(3, 8)]].try_into().unwrap();

        assert_eq!(m1, m3)
    }
}
