use anyhow::{Result, anyhow};
use itertools::iproduct;
use malachite::{base::num::basic::traits::Zero, rational::Rational};
use std::ops::Mul;

use crate::matrix::{
    ebi_matrix::EbiMatrix, fraction_matrix_enum::FractionMatrixEnum,
    fraction_matrix_exact::FractionMatrixExact, fraction_matrix_f64::FractionMatrixF64,
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

        let n = self.number_of_rows();
        let m = self.number_of_columns();
        let p = rhs.number_of_columns();
        let mut values = vec![Rational::ZERO; p * n];

        iproduct!(0..n, 0..p).for_each(|(i, j)| {
            for k in 0..m {
                let idx_ik = self.index(i, k);
                let idx_kj = rhs.index(k, j);
                values[i * n + j] += &self.values[idx_ik] * &rhs.values[idx_kj];
            }
        });

        Ok(FractionMatrixExact {
            values,
            number_of_columns: n,
            number_of_rows: p,
        })
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

    use crate::fraction::fraction::Fraction;
    use std::time::Instant;

    use rand::Rng;

    use crate::{
        exact::MaybeExact,
        f,
        fraction::{fraction_exact::FractionExact, fraction_f64::FractionF64},
        matrix::{
            ebi_matrix::EbiMatrix, fraction_matrix::FractionMatrix,
            fraction_matrix_exact::FractionMatrixExact, fraction_matrix_f64::FractionMatrixF64,
        },
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

        assert_eq!(prod.clone().to_vec(), m3);

        let prod = (&m1 * &m2).unwrap();
        assert_eq!(prod.to_vec(), m3);

        let prod = (&m1 * &m2).unwrap();
        assert_eq!(prod.to_vec(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_1() {
        let m1: FractionMatrix = vec![vec![f!(u64::MAX), f!(2), f!(3)], vec![f!(4), f!(5), f!(6)]]
            .try_into()
            .unwrap();

        let m2: FractionMatrix = vec![
            vec![f!(u64::MAX), f!(8)],
            vec![f!(9), f!(10)],
            vec![f!(11), f!(12)],
        ]
        .try_into()
        .unwrap();

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "340282366920938463426481119284349108276".parse().unwrap(),
                "147573952589676412976".parse().unwrap(),
            ],
            ["73786976294838206571".parse().unwrap(), f!(154)],
        ];

        assert_eq!(prod.to_vec(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_2() {
        let m1: FractionMatrix = vec![vec![f!(u64::MAX), f!(2), f!(3)], vec![f!(4), f!(5), f!(6)]]
            .try_into()
            .unwrap();

        let m2: FractionMatrix = vec![
            vec![f!(1), f!(8)],
            vec![f!(9), f!(10)],
            vec![f!(11), f!(12)],
        ]
        .try_into()
        .unwrap();

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "18446744073709551666".parse().unwrap(),
                "147573952589676412976".parse().unwrap(),
            ],
            [f!(115), f!(154)],
        ];

        assert_eq!(prod.to_vec(), m3);
    }

    #[test]
    fn fraction_matrix_mul_overflow_3() {
        let m1: FractionMatrix = vec![vec![-f!(u64::MAX), f!(2), f!(3)], vec![f!(4), f!(5), f!(6)]]
            .try_into()
            .unwrap();

        let m2: FractionMatrix = vec![
            vec![f!(1), f!(8)],
            vec![f!(9), f!(10)],
            vec![f!(11), f!(12)],
        ]
        .try_into()
        .unwrap();

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            [
                "-18446744073709551564".parse().unwrap(),
                "-147573952589676412864".parse().unwrap(),
            ],
            [f!(115), f!(154)],
        ];

        assert_eq!(prod.to_vec(), m3);
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

        let prod = (&m1 * &m2).unwrap();

        let m3 = vec![
            vec![
                "-340282366920938463426481119284349108174".parse().unwrap(),
                "-147573952589676412864".parse().unwrap(),
            ],
            vec!["73786976294838206571".parse().unwrap(), f!(154)],
        ];

        assert_eq!(prod.to_vec(), m3);
    }

    #[test]
    fn bench_mul() {
        let repeat = 5;
        let size = 100_usize;

        let mut rng = rand::thread_rng();
        let sqrt = 1000_u64;
        let numerators = vec![rng.gen_range(0..sqrt); size * size];
        let denominators = vec![rng.gen_range(0..sqrt); size * size];

        let matrices_f64: Vec<FractionMatrixF64> = (0..repeat)
            .into_iter()
            .map(|i| FractionMatrixF64 {
                number_of_columns: size,
                number_of_rows: size,
                values: numerators
                    .iter()
                    .zip(denominators.iter())
                    .enumerate()
                    .map(|(x, (nom, den))| {
                        if x == i {
                            FractionF64::from((*nom as i64, den + 1)).0
                        } else {
                            FractionF64::from((*nom as i64, *den)).0
                        }
                    })
                    .collect::<Vec<_>>(),
            })
            .collect();

        let matrices_exact: Vec<FractionMatrixExact> = (0..repeat)
            .into_iter()
            .map(|i| FractionMatrixExact {
                number_of_columns: size,
                number_of_rows: size,
                values: numerators
                    .iter()
                    .zip(denominators.iter())
                    .enumerate()
                    .map(|(x, (nom, den))| {
                        if x == i {
                            FractionExact::from((*nom as i64, den + 1)).0
                        } else {
                            FractionExact::from((*nom as i64, *den)).0
                        }
                    })
                    .collect(),
            })
            .collect();

        //exact
        {
            let before = Instant::now();
            for m in matrices_exact {
                let m3 = (&m * &m).unwrap();

                if !m3.is_exact() {
                    panic!()
                }
            }

            println!("exact:         {:.2?}", before.elapsed());
        }

        //f64
        {
            let before = Instant::now();
            for m in matrices_f64 {
                let m3 = (&m * &m).unwrap();

                if m3.is_exact() {
                    panic!()
                }
            }

            println!("approx f64:    {:.2?}", before.elapsed());
        }
    }
}
