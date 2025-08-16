use anyhow::{Result, anyhow};

use crate::matrix::{
    ebi_matrix::EbiMatrix, fraction_matrix_enum::FractionMatrixEnum,
    fraction_matrix_exact::FractionMatrixExact, fraction_matrix_f64::FractionMatrixF64,
};

pub trait IdentityMinus {
    /// For a given matrix M, computes I-M.
    ///
    /// For exact matrices, keeps the current precision level.
    fn identity_minus(&mut self) -> Result<()>;
}

impl IdentityMinus for FractionMatrixF64 {
    fn identity_minus(&mut self) -> Result<()> {
        if self.number_of_rows() != self.number_of_columns() {
            return Err(anyhow!(
                "cannot take identity-minus if rows and columns are not equal"
            ));
        }

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

macro_rules! im {
    ($types:ident, $numerators:ident, $denominators:ident) => {
        for i in 0..$numerators.len() {
            for j in 0..$numerators.len() {
                if i == j {
                    //on the diagonal

                    if $types[i][i].is_plusminus() {
                        //a normal number
                        if $numerators[i][i] >= $denominators[i][i] {
                            //1 - A if A >= 1 equals (A - 1)
                            $types[i][i] = -$types[i][i];
                            $numerators[i][i] -= &$denominators[i][i];
                        } else {
                            //1 - A if A < 1 equals 1 - A
                            $numerators[i][i] = &$denominators[i][i] - &$numerators[i][i];
                        }
                    } else {
                        //weird number; take its inverse
                        $types[i][i] = -$types[i][i];
                    }
                } else {
                    //off the diagonal
                    $types[i][j] = -$types[i][j];
                }
            }
        }
    };
}

impl IdentityMinus for FractionMatrixExact {
    fn identity_minus(&mut self) -> Result<()> {
        if self.number_of_rows() != self.number_of_columns() {
            return Err(anyhow!(
                "cannot take identity-minus if rows and columns are not equal"
            ));
        }
        match self {
            FractionMatrixExact::U64 {
                types,
                numerators,
                denominators,
                ..
            } => im!(types, numerators, denominators),
            FractionMatrixExact::BigInt {
                types,
                numerators,
                denominators,
                ..
            } => {
                im!(types, numerators, denominators)
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
        fraction_enum::FractionEnum,
        matrix::{
            ebi_matrix::EbiMatrix, fraction_matrix_enum::FractionMatrixEnum,
            identity_minus::IdentityMinus,
        },
    };

    #[test]
    fn fraction_matrix_abnormal() {
        let m1: FractionMatrixEnum = vec![
            vec![FractionEnum::infinity(), FractionEnum::neg_infinity()],
            vec![f_en!(8, 3), f_en!(3, 8)],
        ]
        .try_into()
        .unwrap();

        let mut m2 = m1.clone().reduce();

        m2.identity_minus().unwrap();

        let m3: FractionMatrixEnum = vec![
            vec![FractionEnum::neg_infinity(), FractionEnum::infinity()],
            vec![-f_en!(8, 3), f_en!(5, 8)],
        ]
        .try_into()
        .unwrap();

        assert_eq!(m2.to_vec().unwrap(), m3.to_vec().unwrap())
    }
}
