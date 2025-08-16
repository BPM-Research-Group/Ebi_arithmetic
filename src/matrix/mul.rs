use anyhow::{Result, anyhow};
use itertools::iproduct;
use num::BigUint;
use num_bigint::ToBigUint;
use std::ops::Mul;

use crate::{
    ebi_number::{One, Zero},
    matrix::{
        ebi_matrix::EbiMatrix,
        fraction_matrix_enum::FractionMatrixEnum,
        fraction_matrix_exact::FractionMatrixExact,
        fraction_matrix_f64::FractionMatrixF64,
        loose_fraction::{self, LooseFraction, Type},
    },
};

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
                                //overflow detected
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
                        new_new_num[i][j] = new_num[i][j].to_biguint().unwrap();
                        new_new_den[i][j] = new_den[i][j].to_biguint().unwrap();
                    });

                    //second, finish the multiplication on the larger data type
                    iproduct!(0..n, 0..p).skip(r * n + c).for_each(|(i, j)| {
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

impl Mul for &FractionMatrixF64 {
    type Output = Result<FractionMatrixF64>;

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

        let n = self.number_of_rows();
        let m = self.number_of_columns();
        let p = rhs.number_of_columns();
        let mut values = vec![0f64; p * n];

        iproduct!(0..n, 0..p).for_each(|(i, j)| {
            for k in 0..m {
                let idx_ik = self.index(i, k);
                let idx_kj = rhs.index(k, j);
                values[i * n + j] += self.values[idx_ik] * rhs.values[idx_kj];
            }
        });

        Ok(FractionMatrixF64 {
            values,
            number_of_columns: n,
            number_of_rows: p,
        })
    }
}

impl Mul for &FractionMatrixEnum {
    type Output = Result<FractionMatrixEnum>;

    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (FractionMatrixEnum::Approx(m1), FractionMatrixEnum::Approx(m2)) => {
                Ok(FractionMatrixEnum::Approx((m1 * m2)?))
            }
            (FractionMatrixEnum::Exact(m1), FractionMatrixEnum::Exact(m2)) => {
                Ok(FractionMatrixEnum::Exact((m1 * m2)?))
            }
            _ => Ok(FractionMatrixEnum::CannotCombineExactAndApprox),
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        ebi_number::{One, Zero},
        f, f0, f1,
        fraction::Fraction,
        matrix::ebi_matrix::EbiMatrix,
        matrix::fraction_matrix::FractionMatrix,
    };

    #[test]
    fn fraction_matrix_mul() {
        let m1: FractionMatrix = vec![vec![f!(1), f!(2), f!(3)], vec![f!(4), f!(5), f!(6)]]
            .try_into()
            .unwrap();

        (&m1 * &m1).unwrap_err();

        let m2: FractionMatrix = vec![
            vec![f!(7), f!(8)],
            vec![f!(9), f!(10)],
            vec![f!(11), f!(12)],
        ]
        .try_into()
        .unwrap();

        (&m2 * &m2).unwrap_err();

        let prod = (&m1 * &m2).unwrap();

        assert_eq!(prod.number_of_columns(), 2);
        assert_eq!(prod.number_of_rows(), 2);

        let m3 = vec![vec![f!(58), f!(64)], vec![f!(139), f!(154)]];

        assert_eq!(prod.clone().to_vec().unwrap(), m3);

        let m2 = m2.reduce();

        let prod = (&m1 * &m2).unwrap();
        assert_eq!(prod.to_vec().unwrap(), m3);

        let m1 = m1.reduce();

        let prod = (&m1 * &m2).unwrap();
        assert_eq!(prod.to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_nan() {
        let m1 = vec![vec![
            Fraction::infinity(),
            Fraction::neg_infinity(),
            f!(-8, 3),
        ]];
        let m1: FractionMatrix = m1.try_into().unwrap();

        (&m1 * &m1).unwrap_err();

        let m2: FractionMatrix = vec![vec![f0!()], vec![f1!()], vec![f!(-8, 3)]]
            .try_into()
            .unwrap();

        (&m2 * &m2).unwrap_err();

        let prod = (&m1 * &m2).unwrap();

        assert_eq!(prod.number_of_columns(), 1);
        assert_eq!(prod.number_of_rows(), 1);

        let m3 = vec![vec![Fraction::nan()]];

        assert_eq!(prod.clone().to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_1() {
        let m1: FractionMatrix = vec![vec![f!(u64::MAX), f!(2), f!(3)], vec![f!(4), f!(5), f!(6)]]
            .try_into()
            .unwrap();
        let m1 = m1.reduce();

        let m2: FractionMatrix = vec![
            vec![f!(u64::MAX), f!(8)],
            vec![f!(9), f!(10)],
            vec![f!(11), f!(12)],
        ]
        .try_into()
        .unwrap();
        let m2 = m2.reduce();

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "340282366920938463426481119284349108276".parse().unwrap(),
                "147573952589676412976".parse().unwrap(),
            ],
            ["73786976294838206571".parse().unwrap(), f!(154)],
        ];

        assert_eq!(prod.to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_2() {
        let m1: FractionMatrix = vec![vec![f!(u64::MAX), f!(2), f!(3)], vec![f!(4), f!(5), f!(6)]]
            .try_into()
            .unwrap();
        let m1 = m1.reduce();

        let m2: FractionMatrix = vec![
            vec![f!(1), f!(8)],
            vec![f!(9), f!(10)],
            vec![f!(11), f!(12)],
        ]
        .try_into()
        .unwrap();
        let m2 = m2.reduce();

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "18446744073709551666".parse().unwrap(),
                "147573952589676412976".parse().unwrap(),
            ],
            [f!(115), f!(154)],
        ];

        assert_eq!(prod.to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_3() {
        let m1: FractionMatrix = vec![vec![-f!(u64::MAX), f!(2), f!(3)], vec![f!(4), f!(5), f!(6)]]
            .try_into()
            .unwrap();
        let m1 = m1.reduce();

        let m2: FractionMatrix = vec![
            vec![f!(1), f!(8)],
            vec![f!(9), f!(10)],
            vec![f!(11), f!(12)],
        ]
        .try_into()
        .unwrap();
        let m2 = m2.reduce();

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "-18446744073709551564".parse().unwrap(),
                "-147573952589676412864".parse().unwrap(),
            ],
            [f!(115), f!(154)],
        ];

        assert_eq!(prod.to_vec().unwrap(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_4() {
        let m1: FractionMatrix = vec![vec![-f!(u64::MAX), f!(2), f!(3)], vec![f!(4), f!(5), f!(6)]]
            .try_into()
            .unwrap();
        // let m1 = m1.reduce();

        let m2: FractionMatrix = vec![
            vec![f!(u64::MAX), f!(8)],
            vec![f!(9), f!(10)],
            vec![f!(11), f!(12)],
        ]
        .try_into()
        .unwrap();
        let m2 = m2.reduce();

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            vec![
                "-340282366920938463426481119284349108174".parse().unwrap(),
                "-147573952589676412864".parse().unwrap(),
            ],
            vec!["73786976294838206571".parse().unwrap(), f!(154)],
        ];

        assert_eq!(prod.to_vec().unwrap(), m3);
    }
}
