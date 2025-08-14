use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Result, anyhow};
use fraction::{Sign, ToPrimitive};
use num::{BigInt, BigUint, One, Zero};
use num_bigint::ToBigInt;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    ebi_matrix::EbiMatrix,
    exact::{MaybeExact, is_exact_globally},
    fraction::ToExact,
    fraction_enum::FractionEnum,
    fraction_exact::FractionExact,
};

#[derive(Clone, PartialEq, Debug)]
pub enum FractionMatrixEnum {
    F64 {
        number_of_columns: usize,
        values: Vec<Vec<f64>>,
    },
    Fractions {
        number_of_columns: usize,
        values: Vec<Vec<fraction::BigFraction>>,
    },
    I64 {
        number_of_columns: usize,
        values: Vec<Vec<i64>>,
        denominator: BigUint,
    },
    BigInt {
        number_of_columns: usize,
        values: Vec<Vec<BigInt>>,
        denominator: BigUint,
    },
    CannotCombineExactAndApprox,
}

impl FractionMatrixEnum {
    /// Obtains an element from the matrix.
    /// This may be an expensive operation.
    pub fn get(&self, row: usize, column: usize) -> FractionEnum {
        match self {
            FractionMatrixEnum::F64 { values, .. } => FractionEnum::Approx(values[row][column]),
            FractionMatrixEnum::Fractions { values, .. } => {
                FractionEnum::Exact(values[row][column].clone())
            }
            FractionMatrixEnum::I64 {
                values,
                denominator,
                ..
            } => FractionEnum::Exact(FractionExact::to_exact((
                values[row][column],
                denominator.clone(),
            ))),
            FractionMatrixEnum::BigInt {
                values,
                denominator,
                ..
            } => FractionEnum::Exact(FractionExact::to_exact((
                values[row][column].clone(),
                denominator.clone(),
            ))),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                FractionEnum::CannotCombineExactAndApprox
            }
        }
    }

    pub fn to_vec(self) -> Result<Vec<Vec<FractionEnum>>> {
        match self {
            FractionMatrixEnum::F64 { values, .. } => Ok(values
                .into_iter()
                .map(|row| row.into_iter().map(|f| FractionEnum::Approx(f)).collect())
                .collect()),
            FractionMatrixEnum::Fractions { values, .. } => Ok(values
                .into_iter()
                .map(|row| row.into_iter().map(|f| FractionEnum::Exact(f)).collect())
                .collect()),
            FractionMatrixEnum::I64 {
                values,
                denominator,
                ..
            } => Ok(values
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|f| {
                            FractionEnum::Exact(FractionExact::to_exact((f, denominator.clone())))
                        })
                        .collect()
                })
                .collect()),
            FractionMatrixEnum::BigInt {
                values,
                denominator,
                ..
            } => Ok(values
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|f| {
                            FractionEnum::Exact(FractionExact::to_exact((f, denominator.clone())))
                        })
                        .collect()
                })
                .collect()),
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("cannot combine exact and approximate arithmetic"))
            }
        }
    }

    fn lowest_common_multiple_of_denominators(
        values: &Vec<Vec<fraction::BigFraction>>,
    ) -> Option<BigUint> {
        let normal = AtomicBool::new(true);
        let denoms = values
            .par_iter()
            .flat_map(|row| {
                row.par_iter().map(|value| {
                    let x = value.denom(); //we are sure that the number is exact here, so the unwrap() should be safe.
                    if let Some(x) = x {
                        x.clone()
                    } else {
                        normal.store(false, Ordering::Relaxed);
                        BigUint::zero()
                    }
                })
            })
            //denominators
            .reduce(|| BigUint::one(), |a, b| num::integer::lcm(a, b));
        if normal.load(Ordering::Relaxed) {
            Some(denoms)
        } else {
            None
        }
    }
}

impl EbiMatrix for FractionMatrixEnum {
    fn new(number_of_columns: usize) -> Self {
        Self::Fractions {
            number_of_columns,
            values: vec![],
        }
    }

    fn number_of_rows(&self) -> usize {
        match self {
            FractionMatrixEnum::F64 { values, .. } => values.len(),
            FractionMatrixEnum::Fractions { values, .. } => values.len(),
            FractionMatrixEnum::I64 { values, .. } => values.len(),
            FractionMatrixEnum::BigInt { values, .. } => values.len(),
            FractionMatrixEnum::CannotCombineExactAndApprox => 0,
        }
    }

    fn number_of_columns(&self) -> usize {
        match self {
            FractionMatrixEnum::F64 {
                number_of_columns, ..
            }
            | FractionMatrixEnum::Fractions {
                number_of_columns, ..
            }
            | FractionMatrixEnum::I64 {
                number_of_columns, ..
            }
            | FractionMatrixEnum::BigInt {
                number_of_columns, ..
            } => *number_of_columns,
            FractionMatrixEnum::CannotCombineExactAndApprox => 0,
        }
    }

    fn reduce(self) -> Self {
        //try to transform fractions into BigInts
        let result = if let FractionMatrixEnum::Fractions {
            values,
            number_of_columns,
        } = self
        {
            //find the lowest common multiple
            if let Some(lcm) = Self::lowest_common_multiple_of_denominators(&values) {
                //there is a lowest common multiple => every number is normal

                //transform to Vec<Vec<BigInt>>
                let values = values
                    .into_par_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|cell| {
                                let num = cell.numer().unwrap() * &lcm / cell.denom().unwrap();
                                match cell.sign() {
                                    Some(Sign::Plus) => num.to_bigint().unwrap(),
                                    Some(Sign::Minus) => -num.to_bigint().unwrap(),
                                    None => num.to_bigint().unwrap(),
                                }
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();

                Self::BigInt {
                    number_of_columns,
                    values: values,
                    denominator: lcm,
                }
            } else {
                //there is no lowest common multiple; do nothing
                Self::Fractions {
                    number_of_columns,
                    values,
                }
            }
        } else {
            self
        };

        //try to transform BigInts into i64s
        let result2 = if let FractionMatrixEnum::BigInt {
            number_of_columns,
            values,
            denominator,
        } = result
        {
            let values2 = values
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|cell| {
                            if let Some(x) = cell.to_i64() {
                                Ok(x)
                            } else {
                                Err(anyhow!("bla"))
                            }
                        })
                        .collect::<Result<Vec<_>>>()
                })
                .collect::<Result<Vec<_>>>();

            if let Ok(values2) = values2 {
                Self::I64 {
                    number_of_columns,
                    values: values2,
                    denominator,
                }
            } else {
                Self::BigInt {
                    number_of_columns,
                    values,
                    denominator,
                }
            }
        } else {
            result
        };

        result2
    }
}

impl MaybeExact for FractionMatrixEnum {
    type Approximate = ();

    type Exact = FractionMatrixEnum;

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

impl From<Vec<Vec<FractionEnum>>> for FractionMatrixEnum {
    fn from(value: Vec<Vec<FractionEnum>>) -> Self {
        if let Some(x) = value.iter().next() {
            if let Some(y) = x.iter().next() {
                let number_of_columns = x.len();
                //proper matrix
                if y.is_exact() {
                    //exact mode
                    let new_values = value
                        .into_iter()
                        .map(|row| {
                            row.into_iter()
                                .map(|cell| cell.extract_exact().cloned())
                                .collect::<Result<Vec<_>>>()
                        })
                        .collect::<Result<Vec<_>>>();

                    if let Ok(new_values) = new_values {
                        Self::Fractions {
                            number_of_columns,
                            values: new_values,
                        }
                    } else {
                        Self::CannotCombineExactAndApprox
                    }
                } else {
                    //approximate mode
                    let new_values = value
                        .into_iter()
                        .map(|row| {
                            row.into_iter()
                                .map(|cell| cell.extract_approx().cloned())
                                .collect::<Result<Vec<_>>>()
                        })
                        .collect::<Result<Vec<_>>>();

                    if let Ok(new_values) = new_values {
                        Self::F64 {
                            number_of_columns,
                            values: new_values,
                        }
                    } else {
                        Self::CannotCombineExactAndApprox
                    }
                }
            } else {
                //rows, no columns
                if is_exact_globally() {
                    Self::Fractions {
                        number_of_columns: 0,
                        values: vec![vec![]; value.len()],
                    }
                } else {
                    Self::F64 {
                        number_of_columns: 0,
                        values: vec![vec![]; value.len()],
                    }
                }
            }
        } else {
            //no rows
            if is_exact_globally() {
                Self::Fractions {
                    number_of_columns: 0,
                    values: vec![],
                }
            } else {
                Self::F64 {
                    number_of_columns: 0,
                    values: vec![],
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        ebi_matrix::EbiMatrix, f_en, fraction_enum::FractionEnum,
        fraction_matrix_enum::FractionMatrixEnum,
    };

    #[test]
    fn fraction_matrix_abnormal() {
        let m1: FractionMatrixEnum = vec![vec![
            FractionEnum::infinity(),
            FractionEnum::neg_infinity(),
            f_en!(8, 3),
        ]]
        .into();

        let m2 = m1.clone().reduce();

        assert_eq!(m1, m2)
    }

    #[test]
    fn fraction_matrix_reversible() {
        let m1 = vec![vec![
            FractionEnum::infinity(),
            FractionEnum::neg_infinity(),
            f_en!(8, 3),
        ]];

        let m2: FractionMatrixEnum = m1.clone().into();
        let m2 = m2.reduce();

        let m3 = m2.to_vec().unwrap();

        assert_eq!(m1, m3);
    }
}
