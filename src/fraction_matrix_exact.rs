use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Result, anyhow};
use fraction::{Sign, ToPrimitive};
use num::{BigInt, BigUint, One, Zero};
use num_bigint::ToBigInt;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{ebi_matrix::EbiMatrix, exact::MaybeExact, fraction_exact::FractionExact};

#[derive(Clone, PartialEq, Debug)]
pub enum FractionMatrixExact {
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
}

impl FractionMatrixExact {
    /// Obtains an element from the matrix.
    /// This may be an expensive operation.
    pub fn get(&self, row: usize, column: usize) -> FractionExact {
        match self {
            FractionMatrixExact::Fractions { values, .. } => {
                FractionExact(values[row][column].clone())
            }
            FractionMatrixExact::I64 {
                values,
                denominator,
                ..
            } => FractionExact::from((values[row][column], denominator.clone())),
            FractionMatrixExact::BigInt {
                values,
                denominator,
                ..
            } => FractionExact::from((values[row][column].clone(), denominator.clone())),
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
                    let x = value.denom();
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

impl EbiMatrix for FractionMatrixExact {
    fn new(number_of_columns: usize) -> Self {
        Self::Fractions {
            number_of_columns,
            values: vec![],
        }
    }

    fn number_of_rows(&self) -> usize {
        match self {
            FractionMatrixExact::Fractions { values, .. } => values.len(),
            FractionMatrixExact::I64 { values, .. } => values.len(),
            FractionMatrixExact::BigInt { values, .. } => values.len(),
        }
    }

    fn number_of_columns(&self) -> usize {
        match self {
            FractionMatrixExact::Fractions {
                number_of_columns, ..
            }
            | FractionMatrixExact::I64 {
                number_of_columns, ..
            }
            | FractionMatrixExact::BigInt {
                number_of_columns, ..
            } => *number_of_columns,
        }
    }

    fn optimise(self) -> Self {
        //try to transform fractions into BigInts
        let result = if let FractionMatrixExact::Fractions {
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
        let result2 = if let FractionMatrixExact::BigInt {
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

impl From<Vec<Vec<FractionExact>>> for FractionMatrixExact {
    fn from(value: Vec<Vec<FractionExact>>) -> Self {
        if let Some(x) = value.iter().next() {
            let number_of_columns = x.len();
            //proper matrix

            //exact mode
            let new_values = value
                .into_iter()
                .map(|row| row.into_iter().map(|cell| cell.0).collect::<Vec<_>>())
                .collect::<Vec<_>>();

            Self::Fractions {
                number_of_columns,
                values: new_values,
            }
        } else {
            //no rows
            Self::Fractions {
                number_of_columns: 0,
                values: vec![],
            }
        }
    }
}

impl MaybeExact for FractionMatrixExact {
    type Approximate = ();

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

#[cfg(test)]
mod tests {
    use num_bigint::ToBigUint;

    use crate::{
        ebi_matrix::EbiMatrix, f_e, fraction_exact::FractionExact, fraction_matrix_exact::FractionMatrixExact,
    };

    #[test]
    fn fraction_matrix() {
        let m1: FractionMatrixExact = vec![vec![f_e!(1, 4), f_e!(2, 5), f_e!(8, 3)]].into();
        let m2 = m1.clone().optimise();
        assert_ne!(m1, m2);

        let m3 = FractionMatrixExact::I64 {
            number_of_columns: 3,
            values: vec![vec![15, 24, 160]],
            denominator: 60.to_biguint().unwrap(),
        };
        assert_eq!(m2, m3);
    }

    #[test]
    fn fraction_matrix_abnormal() {
        let m1: FractionMatrixExact = vec![vec![
            FractionExact::infinity(),
            FractionExact::neg_infinity(),
            f_e!(8, 3),
        ]]
        .into();

        let m2 = m1.clone().optimise();

        assert_eq!(m1, m2)
    }
}
