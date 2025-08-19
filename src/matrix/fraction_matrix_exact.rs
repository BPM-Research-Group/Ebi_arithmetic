use std::{
    mem,
    sync::atomic::{AtomicBool, Ordering},
};

use anyhow::{Error, Result, anyhow};
use fraction::ToPrimitive;
use itertools::Itertools;
use num::{BigUint, Zero, integer::gcd};
use num_bigint::ToBigUint;

use crate::{
    ebi_number::One,
    exact::MaybeExact,
    fraction_exact::FractionExact,
    fraction_raw::{fraction_raw::FractionRaw, getters::FractionRawGetter},
    matrix::{ebi_matrix::EbiMatrix, loose_fraction::Type},
    pop_front_columns, push_columns,
};

#[derive(Clone, Debug)]
pub enum FractionMatrixExact {
    U64 {
        number_of_columns: usize,
        number_of_rows: usize,
        types: Vec<Type>,
        numerators: Vec<u64>,
        denominators: Vec<u64>,
    },
    BigInt {
        number_of_columns: usize,
        number_of_rows: usize,
        types: Vec<Type>,
        numerators: Vec<BigUint>,
        denominators: Vec<BigUint>,
    },
}

impl FractionMatrixExact {
    /// Obtains all elements from the matrix.
    /// This may be an expensive operation.
    pub fn to_vec(self) -> Result<Vec<Vec<FractionExact>>> {
        Ok(match self {
            FractionMatrixExact::U64 {
                types,
                numerators,
                denominators,
                number_of_columns,
                ..
            } => numerators
                .into_iter()
                .zip(denominators)
                .zip(types)
                .chunks(number_of_columns)
                .into_iter()
                .map(|row| {
                    row.map(|((num, den), typee)| FractionExact::from((typee, num, den)))
                        .collect()
                })
                .collect(),
            FractionMatrixExact::BigInt {
                types,
                numerators,
                denominators,
                number_of_columns,
                ..
            } => numerators
                .into_iter()
                .zip(denominators)
                .zip(types)
                .chunks(number_of_columns)
                .into_iter()
                .map(|row| {
                    row.map(|((num, den), typee)| FractionExact::from((typee, num, den)))
                        .collect()
                })
                .collect(),
        })
    }

    pub(crate) fn index(&self, row: usize, column: usize) -> usize {
        row * self.number_of_columns() + column
    }

    pub(crate) fn clone_to_biguint(&self) -> FractionMatrixExact {
        match self {
            FractionMatrixExact::U64 {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => FractionMatrixExact::BigInt {
                number_of_columns: *number_of_columns,
                number_of_rows: *number_of_rows,
                types: types.clone(),
                numerators: numerators.iter().map(|f| f.to_biguint().unwrap()).collect(),
                denominators: denominators
                    .iter()
                    .map(|f| f.to_biguint().unwrap())
                    .collect(),
            },
            FractionMatrixExact::BigInt { .. } => self.clone(),
        }
    }
}

impl EbiMatrix<FractionExact> for FractionMatrixExact {
    fn new(number_of_columns: usize) -> Self {
        Self::U64 {
            number_of_columns,
            number_of_rows: 0,
            types: vec![],
            numerators: vec![],
            denominators: vec![],
        }
    }

    fn number_of_rows(&self) -> usize {
        match self {
            FractionMatrixExact::U64 { number_of_rows, .. }
            | FractionMatrixExact::BigInt { number_of_rows, .. } => *number_of_rows,
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
                ref mut types,
                ref mut numerators,
                ref mut denominators,
                ..
            } => {
                numerators
                    .iter_mut()
                    .zip(denominators.iter_mut())
                    .zip(types.iter_mut())
                    .for_each(|((num, den), typee)| match typee {
                        Type::Plus | Type::Minus => {
                            if typee.is_plusminus() {
                                if den.is_zero() {
                                    *typee = Type::NaN;
                                    num.set_zero();
                                    return;
                                }
                                if num.is_zero() {
                                    den.set_one();
                                    return;
                                }
                                if num == den {
                                    num.set_one();
                                    den.set_one();
                                    return;
                                }

                                let gcd = gcd(*num, *den);
                                *num /= gcd;
                                *den /= gcd;
                            } else {
                                num.set_zero();
                                den.set_zero();
                            }
                        }
                        _ => {}
                    });
                self
            }
            FractionMatrixExact::BigInt {
                number_of_columns,
                number_of_rows,
                mut types,
                mut numerators,
                mut denominators,
                ..
            } => {
                let fits_u64 = AtomicBool::new(true);
                let max_u64 = u64::MAX.to_biguint().unwrap();
                numerators
                    .iter_mut()
                    .zip(denominators.iter_mut())
                    .zip(types.iter_mut())
                    .for_each(|((num, den), typee)| match typee {
                        Type::Plus | Type::Minus => {
                            if typee.is_plusminus() {
                                if den.is_zero() {
                                    *typee = Type::NaN;
                                    num.set_zero();
                                    return;
                                }
                                if num.is_zero() {
                                    den.set_one();
                                    return;
                                }
                                if num == den {
                                    num.set_one();
                                    den.set_one();
                                    return;
                                }
                                let gcd = gcd(num.clone(), den.clone());
                                *num /= &gcd;
                                *den /= gcd;

                                if fits_u64.load(Ordering::Relaxed) {
                                    //check whether the values would fit in u64
                                    if &*num > &max_u64 || &*den > &max_u64 {
                                        fits_u64.store(false, Ordering::Release);
                                    }
                                }
                            }
                        }
                        _ => {}
                    });
                if fits_u64.load(Ordering::Acquire) {
                    //all values fit in u64; return a u64-matrix
                    Self::U64 {
                        number_of_columns,
                        number_of_rows,
                        types: types,
                        numerators: numerators
                            .into_iter()
                            .map(|num| num.to_u64().unwrap())
                            .collect(),
                        denominators: denominators
                            .into_iter()
                            .map(|num| num.to_u64().unwrap())
                            .collect(),
                    }
                } else {
                    FractionMatrixExact::BigInt {
                        number_of_columns,
                        number_of_rows,
                        types,
                        numerators,
                        denominators,
                    }
                }
            }
        }
    }

    fn eq(&mut self, other: &mut Self) -> bool {
        let mut s = self.clone().reduce();

        let mut o = other.clone().reduce();

        mem::swap(&mut o, other);
        mem::swap(&mut s, self);

        self.inner_eq(&other)
    }

    fn inner_eq(&self, rhs: &Self) -> bool {
        match (self, rhs) {
            (
                FractionMatrixExact::U64 {
                    number_of_columns,
                    number_of_rows,
                    types,
                    numerators,
                    denominators,
                },
                FractionMatrixExact::U64 {
                    number_of_columns: number_of_columns2,
                    number_of_rows: number_of_rows2,
                    types: types2,
                    numerators: numerators2,
                    denominators: denominators2,
                },
            ) => {
                number_of_columns == number_of_columns2
                    && number_of_rows == number_of_rows2
                    && types == types2
                    && numerators == numerators2
                    && denominators == denominators2
            }
            (
                FractionMatrixExact::BigInt {
                    number_of_columns,
                    number_of_rows,
                    types,
                    numerators,
                    denominators,
                },
                FractionMatrixExact::BigInt {
                    number_of_columns: number_of_columns2,
                    number_of_rows: number_of_rows2,
                    types: types2,
                    numerators: numerators2,
                    denominators: denominators2,
                },
            ) => {
                number_of_columns == number_of_columns2
                    && number_of_rows == number_of_rows2
                    && types == types2
                    && numerators == numerators2
                    && denominators == denominators2
            }
            _ => false,
        }
    }

    fn push_columns(&mut self, number_of_columns_to_add: usize) {
        match self {
            FractionMatrixExact::U64 {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => {
                push_columns!(
                    Type::Plus,
                    number_of_columns_to_add,
                    types,
                    *number_of_rows,
                    *number_of_columns
                );
                push_columns!(
                    u64::zero(),
                    number_of_columns_to_add,
                    numerators,
                    *number_of_rows,
                    *number_of_columns
                );
                push_columns!(
                    u64::one(),
                    number_of_columns_to_add,
                    denominators,
                    *number_of_rows,
                    *number_of_columns
                );
                *number_of_columns += number_of_columns_to_add;
            }
            FractionMatrixExact::BigInt {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => {
                push_columns!(
                    Type::Plus,
                    number_of_columns_to_add,
                    types,
                    *number_of_rows,
                    *number_of_columns
                );
                push_columns!(
                    BigUint::zero(),
                    number_of_columns_to_add,
                    numerators,
                    *number_of_rows,
                    *number_of_columns
                );
                push_columns!(
                    BigUint::one(),
                    number_of_columns_to_add,
                    denominators,
                    *number_of_rows,
                    *number_of_columns
                );
                *number_of_columns += number_of_columns_to_add;
            }
        };
    }

    fn pop_front_columns(&mut self, number_of_columns_to_remove: usize) {
        match self {
            FractionMatrixExact::U64 {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => {
                pop_front_columns!(
                    number_of_columns_to_remove,
                    types,
                    *number_of_rows,
                    *number_of_columns
                );
                pop_front_columns!(
                    number_of_columns_to_remove,
                    numerators,
                    *number_of_rows,
                    *number_of_columns
                );
                pop_front_columns!(
                    number_of_columns_to_remove,
                    denominators,
                    *number_of_rows,
                    *number_of_columns
                );
                *number_of_columns -= number_of_columns_to_remove;
            }
            FractionMatrixExact::BigInt {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => {
                pop_front_columns!(
                    number_of_columns_to_remove,
                    types,
                    *number_of_rows,
                    *number_of_columns
                );
                pop_front_columns!(
                    number_of_columns_to_remove,
                    numerators,
                    *number_of_rows,
                    *number_of_columns
                );
                pop_front_columns!(
                    number_of_columns_to_remove,
                    denominators,
                    *number_of_rows,
                    *number_of_columns
                );
                *number_of_columns -= number_of_columns_to_remove;
            }
        }
    }

    fn get(&self, row: usize, column: usize) -> Option<FractionExact> {
        let idx = self.index(row, column);
        Some(match self {
            FractionMatrixExact::U64 {
                numerators,
                denominators,
                types,
                ..
            } => FractionExact::from((
                *types.get(idx)?,
                numerators.get(idx)?.clone(),
                denominators.get(idx)?.clone(),
            )),
            FractionMatrixExact::BigInt {
                numerators,
                denominators,
                types,
                ..
            } => FractionExact::from((
                *types.get(idx)?,
                numerators.get(idx)?.clone(),
                denominators.get(idx)?.clone(),
            )),
        })
    }

    fn set(&mut self, row: usize, column: usize, value: FractionExact) {
        match self {
            FractionMatrixExact::U64 {
                numerators,
                denominators,
                types,
                number_of_columns,
                ..
            } => {
                let idx = row * *number_of_columns + column;
                //check whether the new value fits in u64
                if let Some(FractionRaw(typee, nom, den)) = FractionRaw::try_u64(&value) {
                    //value fits in u64; set
                    types[idx] = typee;
                    numerators[idx] = nom;
                    denominators[idx] = den;
                } else {
                    //value does not fit in u64; tranform the entire matrix
                    let mut new_m = self.clone_to_biguint();
                    new_m.set(row, column, value);
                    std::mem::swap(self, &mut new_m);
                }
            }
            FractionMatrixExact::BigInt {
                numerators,
                denominators,
                types,
                number_of_columns,
                ..
            } => {
                let idx = row * *number_of_columns + column;
                types[idx] = (&value.0).into();
                numerators[idx] = match value.0.numer() {
                    Some(x) => x.clone(),
                    None => BigUint::zero(),
                };
                denominators[idx] = match value.0.denom() {
                    Some(x) => x.clone(),
                    None => BigUint::zero(),
                };
            }
        }
    }

    fn set_one(&mut self, row: usize, column: usize) {
        match self {
            FractionMatrixExact::U64 {
                types,
                numerators,
                denominators,
                number_of_columns,
                ..
            } => {
                let idx = row * *number_of_columns + column;
                types[idx] = Type::Plus;
                numerators[idx] = 1;
                denominators[idx] = 1;
            }
            FractionMatrixExact::BigInt {
                types,
                numerators,
                denominators,
                number_of_columns,
                ..
            } => {
                let idx = row * *number_of_columns + column;
                types[idx] = Type::Plus;
                numerators[idx] = BigUint::one();
                denominators[idx] = BigUint::one();
            }
        }
    }
}

impl TryFrom<(usize, Vec<FractionExact>)> for FractionMatrixExact {
    type Error = Error;

    fn try_from(value: (usize, Vec<FractionExact>)) -> Result<Self> {
        let (number_of_columns, values) = value;
        let number_of_rows = values.len() / number_of_columns;

        if number_of_rows * number_of_columns != values.len() {
            return Err(anyhow!("some cells of the matrix are not provided"));
        }

        if number_of_rows != 0 {
            //has rows

            let mut types = Vec::with_capacity(number_of_rows * number_of_columns);
            for v in values.iter() {
                types.push((&v.0).into());
            }

            let numerators = values
                .iter()
                .map(|cell| match cell.0.numer() {
                    Some(x) => x.clone(),
                    None => BigUint::zero(),
                })
                .collect::<Vec<_>>();

            let denominators = values
                .into_iter()
                .map(|cell| match cell.0.denom() {
                    Some(x) => x.clone(),
                    None => BigUint::zero(),
                })
                .collect::<Vec<_>>();

            Ok(Self::BigInt {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            })
        } else {
            //no rows
            Ok(Self::BigInt {
                number_of_columns: 0,
                number_of_rows: 0,
                types: vec![],
                numerators: vec![],
                denominators: vec![],
            })
        }
    }
}

impl TryFrom<Vec<Vec<FractionExact>>> for FractionMatrixExact {
    type Error = Error;

    fn try_from(value: Vec<Vec<FractionExact>>) -> Result<Self> {
        let number_of_rows = value.len();
        if let Some(x) = value.iter().next() {
            let number_of_columns = x.len();
            //has rows

            let mut types = Vec::with_capacity(number_of_rows * number_of_columns);
            for row in value.iter() {
                if row.len() != number_of_columns {
                    return Err(anyhow!("number of columns is not consistent"));
                }

                for v in row {
                    types.push((&v.0).into());
                }
            }

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
                .flatten()
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
                .flatten()
                .collect::<Vec<_>>();

            Ok(Self::BigInt {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            })
        } else {
            //no rows
            Ok(Self::BigInt {
                number_of_columns: 0,
                number_of_rows: 0,
                types: vec![],
                numerators: vec![],
                denominators: vec![],
            })
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

impl std::fmt::Display for FractionMatrixExact {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{{{")?;

        match self {
            FractionMatrixExact::U64 {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => todo!(),
            FractionMatrixExact::BigInt {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => {
                for i in 0..types.len() {
                    let frac = BigUint::get_ref(i, types, numerators, denominators);
                    write!(f, "{}", frac)?;

                    if i + 1 % number_of_columns == 0 {
                        //comma or something
                        if i + 1 == types.len() {
                            //end of matrix
                            write!(f, "}}")?;
                        } else {
                            //end of row
                            write!(f, "}},\n {{")?;
                        }
                    } else {
                        write!(f, ", ")?;
                    }
                }
            }
        }

        write!(f, "}}}}")
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        f_e,
        fraction_exact::FractionExact,
        matrix::{
            ebi_matrix::EbiMatrix, fraction_matrix_exact::FractionMatrixExact, loose_fraction::Type,
        },
    };

    #[test]
    fn fraction_matrix() {
        let m1: FractionMatrixExact = vec![vec![f_e!(1, 4), f_e!(2, 5), f_e!(8, 3)]]
            .try_into()
            .unwrap();
        let m2 = m1.clone().reduce();

        if let FractionMatrixExact::U64 { .. } = m2 {
        } else {
            panic!()
        }

        let m3 = FractionMatrixExact::U64 {
            number_of_rows: 1,
            number_of_columns: 3,
            types: vec![Type::Plus, Type::Plus, Type::Plus],
            numerators: vec![1, 2, 8],
            denominators: vec![4, 5, 3],
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
        .try_into()
        .unwrap();
        let m2 = m1.clone().reduce();

        let m3 = FractionMatrixExact::U64 {
            number_of_rows: 1,
            number_of_columns: 3,
            types: vec![Type::Infinite, Type::NegInfinite, Type::Plus],
            numerators: vec![0, 0, 8],
            denominators: vec![0, 0, 3],
        };

        assert!(m2.inner_eq(&m3));
    }

    #[test]
    fn fraction_matrix_reversible() {
        let m1 = vec![vec![
            FractionExact::infinity(),
            FractionExact::neg_infinity(),
            f_e!(8, 3),
        ]];

        let m2: FractionMatrixExact = m1.clone().try_into().unwrap();
        let m2 = m2.reduce();

        let m3 = m2.to_vec().unwrap();

        assert_eq!(m1, m3);
    }
}
