use std::{
    ops::Mul,
    sync::atomic::{AtomicBool, Ordering},
};

use anyhow::{Result, anyhow};
use fraction::{Sign, ToPrimitive};
use itertools::iproduct;
use num::{BigUint, Zero, integer::gcd};
use num_bigint::ToBigUint;

use crate::{
    ebi_matrix::EbiMatrix,
    ebi_number::One,
    exact::MaybeExact,
    fraction_exact::FractionExact,
    loose_fraction::{self, LooseFraction, Type},
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
    pub fn to_vec(self) -> Result<Vec<Vec<FractionExact>>> {
        Ok(match self {
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
        })
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
                ref mut types,
                ref mut numerators,
                ref mut denominators,
                ..
            } => {
                numerators
                    .iter_mut()
                    .zip(denominators.iter_mut())
                    .zip(types.iter_mut())
                    .for_each(|((row_num, row_den), row_type)| {
                        row_num
                            .iter_mut()
                            .zip(row_den.iter_mut())
                            .zip(row_type.iter_mut())
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
                            })
                    });
                self
            }
            FractionMatrixExact::BigInt {
                number_of_columns,
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
                    .for_each(|((row_num, row_den), row_type)| {
                        row_num
                            .iter_mut()
                            .zip(row_den.iter_mut())
                            .zip(row_type)
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

impl Mul for &FractionMatrixExact {
    type Output = Result<FractionMatrixExact>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.number_of_columns() != rhs.number_of_rows() {
            return Err(anyhow!(
                "cannot multiply matrix of size {}x{} with a matrix of size {}x{}",
                self.number_of_rows(),
                self.number_of_columns(),
                rhs.number_of_rows(),
                rhs.number_of_columns()
            ));
        }

        match (self, rhs) {
            (
                FractionMatrixExact::U64 {
                    types,
                    numerators,
                    denominators,
                    ..
                },
                FractionMatrixExact::U64 {
                    types: types2,
                    numerators: numerators2,
                    denominators: denominators2,
                    ..
                },
            ) => {
                let n = self.number_of_rows();
                let m = self.number_of_columns();
                let p = rhs.number_of_columns();
                let mut new_types = vec![vec![Type::Plus; p]; n];
                let mut new_num = vec![vec![0; p]; n];
                let mut new_den = vec![vec![1; p]; n];

                let mut last_attempted = None;
                'outer: for i in 0..n {
                    for j in 0..p {
                        for k in 0..m {
                            if loose_fraction::checked_add_assign_mul(
                                &mut new_types[i][j],
                                &mut new_num[i][j],
                                &mut new_den[i][j],
                                types[i][k],
                                &numerators[i][k],
                                &denominators[i][k],
                                types2[k][j],
                                &numerators2[k][j],
                                &denominators2[k][j],
                            ) {
                                //no overlflow; continue
                            } else {
                                //overflow
                                println!("overflow detected at {} {}", i, j);

                                new_types[i][j] = Type::Plus;
                                last_attempted = Some((i, j));
                                break 'outer;
                            }
                        }
                    }
                }

                if let Some((r, c)) = last_attempted {
                    //overflow occurred, salvage results and finish it as a larger data type
                    let mut new_new_num = vec![vec![BigUint::zero(); p]; n];
                    let mut new_new_den = vec![vec![BigUint::one(); p]; n];

                    //first copy the already obtained results to the larger data type
                    iproduct!(0..n, 0..p).take(r * n + c).for_each(|(i, j)| {
                        println!(
                            "salvage {} {} being {}/{}",
                            i, j, new_num[i][j], new_den[i][j]
                        );
                        new_new_num[i][j] = new_num[i][j].to_biguint().unwrap();
                        new_new_den[i][j] = new_den[i][j].to_biguint().unwrap();
                    });

                    //second, finish the multiplication on the larger data type
                    iproduct!(0..n, 0..p).skip(r * n + c).for_each(|(i, j)| {
                        println!("compute {} {}", i, j);
                        for k in 0..m {
                            BigUint::add_assign_mul(
                                &mut new_types[i][j],
                                &mut new_new_num[i][j],
                                &mut new_new_den[i][j],
                                types[i][k],
                                &numerators[i][k].to_biguint().unwrap(),
                                &denominators[i][k].to_biguint().unwrap(),
                                types2[k][j],
                                &numerators2[k][j],
                                &denominators2[k][j],
                            )
                        }
                    });

                    Ok(FractionMatrixExact::BigInt {
                        number_of_columns: n,
                        types: new_types,
                        numerators: new_new_num,
                        denominators: new_new_den,
                    })
                } else {
                    //completed normally
                    Ok(FractionMatrixExact::U64 {
                        number_of_columns: n,
                        types: new_types,
                        numerators: new_num,
                        denominators: new_den,
                    })
                }
            }
            (
                FractionMatrixExact::U64 {
                    types,
                    numerators,
                    denominators,
                    ..
                },
                FractionMatrixExact::BigInt {
                    types: types2,
                    numerators: numerators2,
                    denominators: denominators2,
                    ..
                },
            ) => {
                let n = self.number_of_rows();
                let m = self.number_of_columns();
                let p = rhs.number_of_columns();
                let mut new_types = vec![vec![Type::Plus; p]; n];
                let mut new_num = vec![vec![BigUint::zero(); p]; n];
                let mut new_den = vec![vec![BigUint::one(); p]; n];

                iproduct!(0..n, 0..p).for_each(|(i, j)| {
                    for k in 0..m {
                        BigUint::add_assign_mul(
                            &mut new_types[i][j],
                            &mut new_num[i][j],
                            &mut new_den[i][j],
                            types[i][k],
                            &numerators[i][k],
                            &denominators[i][k],
                            types2[k][j],
                            &numerators2[k][j],
                            &denominators2[k][j],
                        );
                    }
                });

                Ok(FractionMatrixExact::BigInt {
                    number_of_columns: n,
                    types: new_types,
                    numerators: new_num,
                    denominators: new_den,
                })
            }
            (
                FractionMatrixExact::BigInt {
                    types,
                    numerators,
                    denominators,
                    ..
                },
                FractionMatrixExact::U64 {
                    types: types2,
                    numerators: numerators2,
                    denominators: denominators2,
                    ..
                },
            ) => {
                let n = self.number_of_rows();
                let m = self.number_of_columns();
                let p = rhs.number_of_columns();
                let mut new_types = vec![vec![Type::Plus; p]; n];
                let mut new_num = vec![vec![BigUint::zero(); p]; n];
                let mut new_den = vec![vec![BigUint::one(); p]; n];

                iproduct!(0..n, 0..p).for_each(|(i, j)| {
                    for k in 0..m {
                        BigUint::add_assign_mul(
                            &mut new_types[i][j],
                            &mut new_num[i][j],
                            &mut new_den[i][j],
                            types[i][k],
                            &numerators[i][k],
                            &denominators[i][k],
                            types2[k][j],
                            &numerators2[k][j],
                            &denominators2[k][j],
                        );
                    }
                });

                Ok(FractionMatrixExact::BigInt {
                    number_of_columns: n,
                    types: new_types,
                    numerators: new_num,
                    denominators: new_den,
                })
            }
            (
                FractionMatrixExact::BigInt {
                    types,
                    numerators,
                    denominators,
                    ..
                },
                FractionMatrixExact::BigInt {
                    types: types2,
                    numerators: numerators2,
                    denominators: denominators2,
                    ..
                },
            ) => {
                let n = self.number_of_rows();
                let m = self.number_of_columns();
                let p = rhs.number_of_columns();
                let mut new_types = vec![vec![Type::Plus; p]; n];
                let mut new_num = vec![vec![BigUint::zero(); p]; n];
                let mut new_den = vec![vec![BigUint::one(); p]; n];

                iproduct!(0..n, 0..p).for_each(|(i, j)| {
                    for k in 0..m {
                        println!("add_assign_mul {} {} {}", i, j, k);

                        BigUint::add_assign_mul(
                            &mut new_types[i][j],
                            &mut new_num[i][j],
                            &mut new_den[i][j],
                            types[i][k],
                            &numerators[i][k],
                            &denominators[i][k],
                            types2[k][j],
                            &numerators2[k][j],
                            &denominators2[k][j],
                        );
                    }
                });

                Ok(FractionMatrixExact::BigInt {
                    number_of_columns: n,
                    types: new_types,
                    numerators: new_num,
                    denominators: new_den,
                })
            }
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        ebi_matrix::EbiMatrix,
        ebi_number::{One, Zero},
        f_e, f0_e, f1_e,
        fraction_exact::FractionExact,
        fraction_matrix_exact::FractionMatrixExact,
        loose_fraction::Type,
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
            numerators: vec![vec![1, 2, 8]],
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

        let m3 = FractionMatrixExact::U64 {
            number_of_columns: 3,
            types: vec![vec![Type::Infinite, Type::NegInfinite, Type::Plus]],
            numerators: vec![vec![0, 0, 8]],
            denominators: vec![vec![0, 0, 3]],
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

        let m2: FractionMatrixExact = m1.clone().into();
        let m2 = m2.reduce();

        let m3 = m2.to_vec().unwrap();

        assert_eq!(m1, m3);
    }

    #[test]
    fn fraction_matrix_mul() {
        let m1: FractionMatrixExact = vec![
            vec![f_e!(1), f_e!(2), f_e!(3)],
            vec![f_e!(4), f_e!(5), f_e!(6)],
        ]
        .into();

        (&m1 * &m1).unwrap_err();

        let m2: FractionMatrixExact = vec![
            vec![f_e!(7), f_e!(8)],
            vec![f_e!(9), f_e!(10)],
            vec![f_e!(11), f_e!(12)],
        ]
        .into();

        (&m2 * &m2).unwrap_err();

        let prod = (&m1 * &m2).unwrap();

        println!("{:?}", prod);

        assert_eq!(prod.number_of_columns(), 2);
        assert_eq!(prod.number_of_rows(), 2);

        let m3 = vec![vec![f_e!(58), f_e!(64)], vec![f_e!(139), f_e!(154)]];

        assert_eq!(prod.clone().to_vec().unwrap(), m3);

        let m2 = m2.reduce();

        let prod = (&m1 * &m2).unwrap();
        println!("{:?}", prod);
        assert_eq!(prod.to_vec().unwrap(), m3);

        let m1 = m1.reduce();

        let prod = (&m1 * &m2).unwrap();
        println!("{:?}", prod);
        assert_eq!(prod.to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_nan() {
        let m1 = vec![vec![
            FractionExact::infinity(),
            FractionExact::neg_infinity(),
            f_e!(-8, 3),
        ]];
        let m1: FractionMatrixExact = m1.into();

        (&m1 * &m1).unwrap_err();

        let m2: FractionMatrixExact = vec![vec![f0_e!()], vec![f1_e!()], vec![f_e!(-8, 3)]].into();

        (&m2 * &m2).unwrap_err();

        let prod = (&m1 * &m2).unwrap();

        assert_eq!(prod.number_of_columns(), 1);
        assert_eq!(prod.number_of_rows(), 1);

        println!("{:?}", prod);

        let m3 = vec![vec![FractionExact::nan()]];

        assert_eq!(prod.clone().to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_1() {
        let m1: FractionMatrixExact = vec![
            vec![f_e!(u64::MAX), f_e!(2), f_e!(3)],
            vec![f_e!(4), f_e!(5), f_e!(6)],
        ]
        .into();
        let m1 = m1.reduce();
        println!("m1 {:?}", m1);

        let m2: FractionMatrixExact = vec![
            vec![f_e!(u64::MAX), f_e!(8)],
            vec![f_e!(9), f_e!(10)],
            vec![f_e!(11), f_e!(12)],
        ]
        .into();
        let m2 = m2.reduce();
        println!("m2 {:?}", m2);

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "340282366920938463426481119284349108276".parse().unwrap(),
                "147573952589676412976".parse().unwrap(),
            ],
            ["73786976294838206571".parse().unwrap(), f_e!(154)],
        ];

        assert_eq!(prod.to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_2() {
        let m1: FractionMatrixExact = vec![
            vec![f_e!(u64::MAX), f_e!(2), f_e!(3)],
            vec![f_e!(4), f_e!(5), f_e!(6)],
        ]
        .into();
        let m1 = m1.reduce();
        println!("m1 {:?}", m1);

        let m2: FractionMatrixExact = vec![
            vec![f_e!(1), f_e!(8)],
            vec![f_e!(9), f_e!(10)],
            vec![f_e!(11), f_e!(12)],
        ]
        .into();
        let m2 = m2.reduce();
        println!("m2 {:?}", m2);

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "18446744073709551666".parse().unwrap(),
                "147573952589676412976".parse().unwrap(),
            ],
            [f_e!(115), f_e!(154)],
        ];

        assert_eq!(prod.to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_3() {
        let m1: FractionMatrixExact = vec![
            vec![-f_e!(u64::MAX), f_e!(2), f_e!(3)],
            vec![f_e!(4), f_e!(5), f_e!(6)],
        ]
        .into();
        let m1 = m1.reduce();
        println!("m1 {:?}", m1);

        let m2: FractionMatrixExact = vec![
            vec![f_e!(1), f_e!(8)],
            vec![f_e!(9), f_e!(10)],
            vec![f_e!(11), f_e!(12)],
        ]
        .into();
        let m2 = m2.reduce();
        println!("m2 {:?}", m2);

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "-18446744073709551564".parse().unwrap(),
                "-147573952589676412864".parse().unwrap(),
            ],
            [f_e!(115), f_e!(154)],
        ];

        assert_eq!(prod.to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_4() {
        let m1: FractionMatrixExact = vec![
            vec![-f_e!(u64::MAX), f_e!(2), f_e!(3)],
            vec![f_e!(4), f_e!(5), f_e!(6)],
        ]
        .into();
        // let m1 = m1.reduce();

        let m2: FractionMatrixExact = vec![
            vec![f_e!(u64::MAX), f_e!(8)],
            vec![f_e!(9), f_e!(10)],
            vec![f_e!(11), f_e!(12)],
        ]
        .into();
        let m2 = m2.reduce();

        let prod = (&m1 * &m2).unwrap();
        println!("prod {:?}", prod);

        let m3 = vec![
            vec![
                "-340282366920938463426481119284349108174".parse().unwrap(),
                "-147573952589676412864".parse().unwrap(),
            ],
            vec!["73786976294838206571".parse().unwrap(), f_e!(154)],
        ];

        assert_eq!(prod.to_vec().unwrap(), m3);
    }
}
