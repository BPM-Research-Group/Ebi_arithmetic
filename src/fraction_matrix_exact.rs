use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Result, anyhow};
use fraction::{Sign, ToPrimitive};
use num::{BigUint, Zero, integer::gcd};
use num_bigint::ToBigUint;

use crate::{
    ebi_matrix::EbiMatrix,
    exact::MaybeExact,
    fraction_exact::{FractionExact, Type},
};

#[derive(Clone, Debug)]
pub enum FractionMatrixExact {
    U64 {
        number_of_columns: usize,
        types: Vec<Vec<Type>>,
        numerators: Vec<Vec<u64>>,
        denominators: Vec<Vec<u64>>,
    },
    BigInt {
        number_of_columns: usize,
        types: Vec<Vec<Type>>,
        numerators: Vec<Vec<BigUint>>,
        denominators: Vec<Vec<BigUint>>,
    },
}

impl FractionMatrixExact {
    /// Obtains an element from the matrix.
    /// This may be an expensive operation; consider using to_vec() to avoid some cloning.
    pub fn get(&self, row: usize, column: usize) -> FractionExact {
        match self {
            FractionMatrixExact::U64 {
                numerators,
                denominators,
                types,
                ..
            } => FractionExact::from((
                types[row][column],
                numerators[row][column].clone(),
                denominators[row][column].clone(),
            )),
            FractionMatrixExact::BigInt {
                numerators,
                denominators,
                types,
                ..
            } => FractionExact::from((
                types[row][column],
                numerators[row][column].clone(),
                denominators[row][column].clone(),
            )),
        }
    }

    /// Obtains all elements from the matrix.
    /// This may be an expensive operation.
    pub fn to_vec(self) -> Vec<Vec<FractionExact>> {
        match self {
            FractionMatrixExact::U64 {
                types,
                numerators,
                denominators,
                ..
            } => numerators
                .into_iter()
                .zip(denominators)
                .zip(types)
                .map(|((row_num, row_den), row_types)| {
                    row_num
                        .into_iter()
                        .zip(row_den)
                        .zip(row_types)
                        .map(|((num, den), typee)| FractionExact::from((typee, num, den)))
                        .collect()
                })
                .collect(),
            FractionMatrixExact::BigInt {
                types,
                numerators,
                denominators,
                ..
            } => numerators
                .into_iter()
                .zip(denominators)
                .zip(types)
                .map(|((row_num, row_den), row_types)| {
                    row_num
                        .into_iter()
                        .zip(row_den)
                        .zip(row_types)
                        .map(|((num, den), typee)| FractionExact::from((typee, num, den)))
                        .collect()
                })
                .collect(),
        }
    }

    #[cfg(test)]
    fn inner_eq(&self, rhs: &Self) -> bool {
        match (self, rhs) {
            (
                FractionMatrixExact::U64 {
                    number_of_columns,
                    types,
                    numerators,
                    denominators,
                },
                FractionMatrixExact::U64 {
                    number_of_columns: number_of_columns2,
                    types: types2,
                    numerators: numerators2,
                    denominators: denominators2,
                },
            ) => {
                number_of_columns == number_of_columns2
                    && types == types2
                    && numerators == numerators2
                    && denominators == denominators2
            }
            (
                FractionMatrixExact::BigInt {
                    number_of_columns,
                    types,
                    numerators,
                    denominators,
                },
                FractionMatrixExact::BigInt {
                    number_of_columns: number_of_columns2,
                    types: types2,
                    numerators: numerators2,
                    denominators: denominators2,
                },
            ) => {
                number_of_columns == number_of_columns2
                    && types == types2
                    && numerators == numerators2
                    && denominators == denominators2
            }
            _ => false,
        }
    }
}

impl EbiMatrix for FractionMatrixExact {
    fn new(number_of_columns: usize) -> Self {
        Self::U64 {
            number_of_columns,
            types: vec![],
            numerators: vec![],
            denominators: vec![],
        }
    }

    fn number_of_rows(&self) -> usize {
        match self {
            FractionMatrixExact::U64 { numerators, .. } => numerators.len(),
            FractionMatrixExact::BigInt { numerators, .. } => numerators.len(),
        }
    }

    fn number_of_columns(&self) -> usize {
        match self {
            FractionMatrixExact::U64 {
                number_of_columns, ..
            }
            | FractionMatrixExact::BigInt {
                number_of_columns, ..
            } => *number_of_columns,
        }
    }

    fn reduce(mut self) -> Self {
        match self {
            FractionMatrixExact::U64 {
                ref types,
                ref mut numerators,
                ref mut denominators,
                ..
            } => {
                numerators
                    .iter_mut()
                    .zip(denominators.iter_mut())
                    .zip(types)
                    .for_each(|((row_num, row_den), row_type)| {
                        row_num
                            .iter_mut()
                            .zip(row_den.iter_mut())
                            .zip(row_type)
                            .for_each(|((num, den), typee)| match typee {
                                Type::Plus | Type::Minus => {
                                    let gcd = gcd(*num, *den);
                                    *num /= gcd;
                                    *den /= gcd;
                                }
                                _ => {}
                            })
                    });
                self
            }
            FractionMatrixExact::BigInt {
                number_of_columns,
                types,
                mut numerators,
                mut denominators,
                ..
            } => {
                let fits_u64 = AtomicBool::new(true);
                let max_u64 = u64::MAX.to_biguint().unwrap();
                numerators
                    .iter_mut()
                    .zip(denominators.iter_mut())
                    .zip(types.iter())
                    .for_each(|((row_num, row_den), row_type)| {
                        row_num
                            .iter_mut()
                            .zip(row_den.iter_mut())
                            .zip(row_type)
                            .for_each(|((num, den), typee)| match typee {
                                Type::Plus | Type::Minus => {
                                    let gcd = gcd(num.clone(), den.clone());
                                    *num /= &gcd;
                                    *den /= gcd;

                                    if fits_u64.load(Ordering::Relaxed) {
                                        //check whether the values would fit in u64
                                        if &*num > &max_u64 || &*den > &max_u64 {
                                            println!("set false {} {}", num, den);
                                            fits_u64.store(false, Ordering::Release);
                                        }
                                    }
                                }
                                _ => {}
                            })
                    });
                if fits_u64.load(Ordering::Acquire) {
                    //all values fit in u64; return a u64-matrix
                    Self::U64 {
                        number_of_columns: number_of_columns,
                        types: types,
                        numerators: numerators
                            .into_iter()
                            .map(|row| row.into_iter().map(|num| num.to_u64().unwrap()).collect())
                            .collect(),
                        denominators: denominators
                            .into_iter()
                            .map(|row| row.into_iter().map(|num| num.to_u64().unwrap()).collect())
                            .collect(),
                    }
                } else {
                    FractionMatrixExact::BigInt {
                        number_of_columns,
                        types,
                        numerators,
                        denominators,
                    }
                }
            }
        }
    }

    fn invert(&mut self) -> Result<()> {
        todo!()
    }
}

impl From<Vec<Vec<FractionExact>>> for FractionMatrixExact {
    fn from(value: Vec<Vec<FractionExact>>) -> Self {
        if let Some(x) = value.iter().next() {
            let number_of_columns = x.len();
            //has rows

            let types = value
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|v| match &v.0 {
                            fraction::GenericFraction::Rational(Sign::Plus, _) => Type::Plus,
                            fraction::GenericFraction::Rational(Sign::Minus, _) => Type::Minus,
                            fraction::GenericFraction::Infinity(Sign::Plus) => Type::Infinite,
                            fraction::GenericFraction::Infinity(Sign::Minus) => Type::NegInfinite,
                            fraction::GenericFraction::NaN => Type::NaN,
                        })
                        .collect()
                })
                .collect();

            let numerators = value
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|cell| match cell.0.numer() {
                            Some(x) => x.clone(),
                            None => BigUint::zero(),
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();

            let denominators = value
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|cell| match cell.0.denom() {
                            Some(x) => x.clone(),
                            None => BigUint::zero(),
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();

            Self::BigInt {
                number_of_columns,
                types,
                numerators,
                denominators,
            }
        } else {
            //no rows
            Self::BigInt {
                number_of_columns: 0,
                types: vec![],
                numerators: vec![],
                denominators: vec![],
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

    use crate::{
        ebi_matrix::EbiMatrix,
        f_e,
        fraction_exact::{FractionExact, Type},
        fraction_matrix_exact::FractionMatrixExact,
    };

    #[test]
    fn fraction_matrix() {
        let m1: FractionMatrixExact = vec![vec![f_e!(1, 4), f_e!(2, 5), f_e!(8, 3)]].into();
        let m2 = m1.clone().reduce();

        println!("{:?}", m2);

        if let FractionMatrixExact::U64 { .. } = m2 {
        } else {
            panic!()
        }

        let m3 = FractionMatrixExact::U64 {
            number_of_columns: 3,
            types: vec![vec![Type::Plus, Type::Plus, Type::Plus]],
            numerators: vec![vec![1, 2, 5, 8]],
            denominators: vec![vec![4, 5, 3]],
        };
        assert!(m2.inner_eq(&m3));
    }

    #[test]
    fn fraction_matrix_abnormal() {
        let m1: FractionMatrixExact = vec![vec![
            FractionExact::infinity(),
            FractionExact::neg_infinity(),
            f_e!(8, 3),
        ]]
        .into();

        let m2 = m1.clone().reduce();

        assert!(m1.inner_eq(&m2))
    }

    #[test]
    fn fraction_matrix_reversible() {
        let m1 = vec![vec![
            FractionExact::infinity(),
            FractionExact::neg_infinity(),
            f_e!(8, 3),
        ]];

        let m2: FractionMatrixExact = m1.clone().into();
        let m2 = m2.reduce();

        let m3 = m2.to_vec();

        assert_eq!(m1, m3);
    }
}
