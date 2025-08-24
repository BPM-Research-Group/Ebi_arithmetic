use std::{
    fmt::{Debug, Display},
    mem,
};

use anyhow::{Error, Result, anyhow};

use crate::{
    exact::{self, MaybeExact, is_exact_globally},
    matrix::{
        ebi_matrix::EbiMatrix, fraction_matrix_exact::FractionMatrixExact,
        fraction_matrix_f64::FractionMatrixF64,
    },
};

#[derive(Clone, Debug)]
pub enum FractionMatrixEnum {
    Approx(FractionMatrixF64),
    Exact(FractionMatrixExact),
    CannotCombineExactAndApprox,
}

impl FractionMatrixEnum {
    pub fn to_vec(self) -> Result<Vec<Vec<FractionEnum>>> {
        match self {
            FractionMatrixEnum::Approx(m) => Ok(m
                .to_vec()?
                .into_iter()
                .map(|row| row.into_iter().map(|f| FractionEnum::Approx(f.0)).collect())
                .collect()),
            FractionMatrixEnum::Exact(m) => Ok(m
                .to_vec()?
                .into_iter()
                .map(|row| row.into_iter().map(|f| FractionEnum::Exact(f.0)).collect())
                .collect()),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine exact and approximate arithmetic"))
            }
        }
    }
}

impl EbiMatrix<FractionEnum> for FractionMatrixEnum {
    fn new(number_of_columns: usize) -> Self {
        if exact::is_exact_globally() {
            Self::Exact(FractionMatrixExact::new(number_of_columns))
        } else {
            Self::Approx(FractionMatrixF64::new(number_of_columns))
        }
    }

    fn number_of_rows(&self) -> usize {
        match self {
            FractionMatrixEnum::Approx(m) => m.number_of_rows(),
            FractionMatrixEnum::Exact(m) => m.number_of_rows(),
            FractionMatrixEnum::CannotCombineExactAndApprox => 0,
        }
    }

    fn number_of_columns(&self) -> usize {
        match self {
            FractionMatrixEnum::Approx(m) => m.number_of_columns(),
            FractionMatrixEnum::Exact(m) => m.number_of_columns(),
            FractionMatrixEnum::CannotCombineExactAndApprox => 0,
        }
    }

    fn reduce(self) -> Self {
        match self {
            FractionMatrixEnum::Approx(m) => FractionMatrixEnum::Approx(m.reduce()),
            FractionMatrixEnum::Exact(m) => FractionMatrixEnum::Exact(m.reduce()),
            FractionMatrixEnum::CannotCombineExactAndApprox => self,
        }
    }

    fn eq(&mut self, other: &mut Self) -> bool {
        match (self, other) {
            (FractionMatrixEnum::Approx(m1), FractionMatrixEnum::Approx(m2)) => m1.eq(m2),
            (FractionMatrixEnum::Exact(m1), FractionMatrixEnum::Exact(m2)) => m1.eq(m2),
            _ => false,
        }
    }

    fn inner_eq(&self, other: &Self) -> bool {
        match (self, other) {
            (FractionMatrixEnum::Approx(m1), FractionMatrixEnum::Approx(m2)) => m1.inner_eq(m2),
            (FractionMatrixEnum::Exact(m1), FractionMatrixEnum::Exact(m2)) => m1.inner_eq(m2),
            _ => false,
        }
    }

    fn push_columns(&mut self, number_of_columns_to_add: usize) {
        match self {
            FractionMatrixEnum::Approx(m) => m.push_columns(number_of_columns_to_add),
            FractionMatrixEnum::Exact(m) => m.push_columns(number_of_columns_to_add),
            FractionMatrixEnum::CannotCombineExactAndApprox => {}
        }
    }

    fn pop_front_columns(&mut self, number_of_columns_to_remove: usize) {
        match self {
            FractionMatrixEnum::Approx(m) => m.pop_front_columns(number_of_columns_to_remove),
            FractionMatrixEnum::Exact(m) => m.pop_front_columns(number_of_columns_to_remove),
            FractionMatrixEnum::CannotCombineExactAndApprox => {}
        }
    }

    fn get(&self, row: usize, column: usize) -> Option<FractionEnum> {
        Some(match self {
            FractionMatrixEnum::Approx(m) => FractionEnum::Approx(m.get(row, column)?.0),
            FractionMatrixEnum::Exact(m) => FractionEnum::Exact(m.get(row, column)?.0),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                FractionEnum::CannotCombineExactAndApprox
            }
        })
    }

    fn set(&mut self, row: usize, column: usize, value: FractionEnum) {
        if let FractionMatrixEnum::Approx(m) = self {
            if let FractionEnum::Approx(f) = value {
                m.set(row, column, FractionF64(f));
                return;
            }
        } else if let FractionMatrixEnum::Exact(m) = self {
            if let FractionEnum::Exact(f) = value {
                m.set(row, column, FractionExact(f));
                return;
            }
        }
        mem::swap(self, &mut FractionMatrixEnum::CannotCombineExactAndApprox);
    }

    fn set_one(&mut self, row: usize, column: usize) {
        match self {
            FractionMatrixEnum::Approx(m) => m.set_one(row, column),
            FractionMatrixEnum::Exact(m) => m.set_one(row, column),
            FractionMatrixEnum::CannotCombineExactAndApprox => {}
        }
    }
}

impl MaybeExact for FractionMatrixEnum {
    type Approximate = ();

    type Exact = FractionMatrixEnum;

    fn is_exact(&self) -> bool {
        true
    }

    fn extract_approx(&self) -> anyhow::Result<&()> {
        Err(anyhow!("cannot extract a float from a fraction"))
    }

    fn extract_exact(&self) -> anyhow::Result<&FractionMatrixEnum> {
        Ok(self)
    }
}

impl TryFrom<(usize, Vec<FractionEnum>)> for FractionMatrixEnum {
    type Error = Error;

    fn try_from(value: (usize, Vec<FractionEnum>)) -> Result<Self> {
        let (number_of_columns, values) = value;
        let number_of_rows = values.len() / number_of_columns;

        if number_of_rows * number_of_columns != values.len() {
            return Err(anyhow!("some cells of the matrix are not provided"));
        }

        if number_of_rows > 0 {
            if let Some(x) = values.iter().next() {
                //proper matrix
                if x.is_exact() {
                    //exact mode
                    let mut new_values = Vec::with_capacity(values.len());
                    for f in values {
                        match f {
                            FractionEnum::Exact(f) => new_values.push(FractionExact(f)),
                            FractionEnum::Approx(_) => {
                                return Err(anyhow!(
                                    "cannot combine approximate and exact arithmetic"
                                ));
                            }
                            FractionEnum::CannotCombineExactAndApprox => {
                                return Err(anyhow!(
                                    "cannot combine approximate and exact arithmetic"
                                ));
                            }
                        }
                    }

                    let m: FractionMatrixExact = (number_of_columns, new_values).try_into()?;
                    Ok(Self::Exact(m))
                } else {
                    //approximate mode
                    let mut new_values = Vec::with_capacity(values.len());
                    for f in values {
                        match f {
                            FractionEnum::Exact(_) => {
                                return Err(anyhow!(
                                    "cannot combine approximate and exact arithmetic"
                                ));
                            }
                            FractionEnum::Approx(f) => new_values.push(FractionF64(f)),
                            FractionEnum::CannotCombineExactAndApprox => {
                                return Err(anyhow!(
                                    "cannot combine approximate and exact arithmetic"
                                ));
                            }
                        }
                    }

                    let m: FractionMatrixF64 = (number_of_columns, new_values).try_into()?;
                    Ok(Self::Approx(m))
                }
            } else {
                //no columns
                if is_exact_globally() {
                    let new_rows = vec![vec![]; values.len()];
                    let m: FractionMatrixExact = new_rows.try_into()?;
                    Ok(Self::Exact(m))
                } else {
                    let new_rows = vec![vec![]; values.len()];
                    let m: FractionMatrixF64 = new_rows.try_into()?;
                    Ok(Self::Approx(m))
                }
            }
        } else {
            //no rows
            if is_exact_globally() {
                Ok(Self::Exact(FractionMatrixExact::new(0)))
            } else {
                Ok(Self::Approx(FractionMatrixF64::new(0)))
            }
        }
    }
}

impl TryFrom<Vec<Vec<FractionEnum>>> for FractionMatrixEnum {
    type Error = Error;

    fn try_from(value: Vec<Vec<FractionEnum>>) -> Result<Self> {
        if let Some(x) = value.iter().next() {
            if let Some(y) = x.iter().next() {
                let number_of_columns = x.len();
                //proper matrix
                if y.is_exact() {
                    //exact mode
                    let mut new_rows = Vec::with_capacity(value.len());
                    for row in value {
                        if row.len() != number_of_columns {
                            return Err(anyhow!("number of columns is not consistent"));
                        }

                        let mut new_row = Vec::with_capacity(number_of_columns);
                        for f in row {
                            match f {
                                FractionEnum::Exact(f) => new_row.push(FractionExact(f)),
                                FractionEnum::Approx(_) => {
                                    return Err(anyhow!(
                                        "cannot combine approximate and exact arithmetic"
                                    ));
                                }
                                FractionEnum::CannotCombineExactAndApprox => {
                                    return Err(anyhow!(
                                        "cannot combine approximate and exact arithmetic"
                                    ));
                                }
                            }
                        }
                        new_rows.push(new_row);
                    }

                    let m: FractionMatrixExact = new_rows.try_into()?;
                    Ok(Self::Exact(m))
                } else {
                    //approximate mode
                    let mut new_rows = Vec::with_capacity(value.len());
                    for row in value {
                        if row.len() != number_of_columns {
                            return Err(anyhow!("number of columns is not consistent"));
                        }

                        let mut new_row = Vec::with_capacity(number_of_columns);
                        for f in row {
                            match f {
                                FractionEnum::Exact(_) => {
                                    return Err(anyhow!(
                                        "cannot combine approximate and exact arithmetic"
                                    ));
                                }
                                FractionEnum::Approx(f) => new_row.push(FractionF64(f)),
                                FractionEnum::CannotCombineExactAndApprox => {
                                    return Err(anyhow!(
                                        "cannot combine approximate and exact arithmetic"
                                    ));
                                }
                            }
                        }
                        new_rows.push(new_row);
                    }

                    let m: FractionMatrixF64 = new_rows.try_into()?;
                    Ok(Self::Approx(m))
                }
            } else {
                //rows, no columns
                if is_exact_globally() {
                    let new_rows = vec![vec![]; value.len()];
                    let m: FractionMatrixExact = new_rows.try_into()?;
                    Ok(Self::Exact(m))
                } else {
                    let new_rows = vec![vec![]; value.len()];
                    let m: FractionMatrixF64 = new_rows.try_into()?;
                    Ok(Self::Approx(m))
                }
            }
        } else {
            //no rows
            if is_exact_globally() {
                Ok(Self::Exact(FractionMatrixExact::new(0)))
            } else {
                Ok(Self::Approx(FractionMatrixF64::new(0)))
            }
        }
    }
}

impl Display for FractionMatrixEnum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Approx(m) => Display::fmt(m, f),
            Self::Exact(m) => Display::fmt(m, f),
            Self::CannotCombineExactAndApprox => write!(f, "CannotCombineExactAndApprox"),
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        f_en,
        fraction_enum::FractionEnum,
        matrix::{ebi_matrix::EbiMatrix, fraction_matrix_enum::FractionMatrixEnum},
    };

    #[test]
    fn fraction_matrix_abnormal() {
        let mut m1: FractionMatrixEnum = vec![vec![
            FractionEnum::infinity(),
            FractionEnum::neg_infinity(),
            f_en!(8, 3),
        ]]
        .try_into()
        .unwrap();

        let mut m2 = m1.clone().reduce();

        assert!(m1.eq(&mut m2))
    }

    #[test]
    fn fraction_matrix_reversible() {
        let m1 = vec![vec![
            FractionEnum::infinity(),
            FractionEnum::neg_infinity(),
            f_en!(8, 3),
        ]];

        let m2: FractionMatrixEnum = m1.clone().try_into().unwrap();
        let m2 = m2.reduce();

        let m3 = m2.to_vec().unwrap();

        assert_eq!(m1, m3);
    }
}
