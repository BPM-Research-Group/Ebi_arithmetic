use std::mem;

use crate::{
    ebi_number::{One, Zero},
    fraction_raw::{getters::FractionRawGetter, recip::Recip},
    matrix::{
        ebi_matrix::EbiMatrix, fraction_matrix_exact::FractionMatrixExact,
        fraction_matrix_f64::FractionMatrixF64, gauss_jordan::GaussJordan,
    },
};
use anyhow::{Result, anyhow};
use num::BigUint;

pub trait Inversion {
    fn invert(self) -> Result<Self>
    where
        Self: Sized;
}

impl Inversion for FractionMatrixF64 {
    fn invert(mut self) -> Result<Self> {
        if self.number_of_columns() != self.number_of_rows() {
            return Err(anyhow!("can only take the inverse of a square matrix"));
        }

        //optimisation: size-zero matrix
        if self.number_of_rows().is_zero() {
            return Ok(self);
        }

        //optimisation: size-one matrix
        if self.number_of_rows().is_one() {
            if self.values[0].is_zero() {
                return Err(anyhow!("matrix is not invertible"));
            }

            self.values[0] = self.values[0].recip();
            return Ok(self);
        }

        //optimisation: size-two matrix
        if self.number_of_rows() == 2 {
            //compute determinant
            let mut det = self.values[0];
            det *= &self.values[3];
            let mut det2 = self.values[1];
            det2 *= &self.values[2];
            det -= det2;

            if det.is_zero() {
                return Err(anyhow!("matrix is not invertible"));
            }

            // log::debug!("determinant {}", det);

            det = det.recip();

            //perform inverse
            let (m1, m2) = self.values.split_at_mut(2);
            mem::swap(&mut m1[0], &mut m2[1]);

            self.values[0] *= &det;
            self.values[3] *= &det;

            self.values[2] *= -&det;
            self.values[1] *= -det;
            return Ok(self);
        }

        // log::debug!("\t\tcompute inverse of {}", self);

        //extend the rows with the identity matrix
        self.push_columns(self.number_of_rows);
        for i in 0..self.number_of_rows {
            let idx_ii = self.index(i, self.number_of_rows + i);
            self.values[idx_ii] = f64::one();
        }

        // println!("add identity\n{}", self);

        //solve
        self = self.gauss_jordan_reduced()?;

        // println!("solved\n{}", self);

        //remove the columns
        self.pop_front_columns(self.number_of_rows);

        // log::info!("inverse done");

        Ok(self)
    }
}

macro_rules! inverse {
    () => {};
}

impl Inversion for FractionMatrixExact {
    fn invert(mut self) -> Result<Self> {
        if self.number_of_columns() != self.number_of_rows() {
            return Err(anyhow!("can only take the inverse of a square matrix"));
        }

        //optimisation: size-zero matrix
        if self.number_of_rows().is_zero() {
            return Ok(self);
        }

        match &mut self {
            FractionMatrixExact::U64 {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => {
                todo!()
            }
            FractionMatrixExact::BigInt {
                number_of_rows,
                types,
                numerators,
                denominators,
                ..
            } => {
                //optimisation: size-one matrix
                if number_of_rows.is_one() {
                    if numerators[0].is_zero() {
                        return Err(anyhow!("matrix is not invertible: determinant is zero"));
                    }

                    BigUint::get_mut(0, types, numerators, denominators).recip();
                    return Ok(self);
                }

                //optimisation: size-two matrix
                if *number_of_rows == 2 {
                    //compute determinant
                    let mut det = BigUint::get_clone(0, &types, &numerators, &denominators);
                    det *= BigUint::get_ref(3, &types, &numerators, &denominators);
                    let mut det2 = BigUint::get_clone(1, &types, &numerators, &denominators);
                    det2 *= BigUint::get_ref(2, &types, &numerators, &denominators);
                    det -= &det2;

                    if det.is_zero() {
                        return Err(anyhow!("matrix is not invertible: determinant is zero"));
                    }

                    // log::debug!("determinant {}", det);

                    det.recip();

                    //perform inverse
                    let (typee1, typee2) = types.split_at_mut(2);
                    let (num1, num2) = numerators.split_at_mut(2);
                    let (den1, den2) = denominators.split_at_mut(2);
                    mem::swap(&mut typee1[0], &mut typee2[1]);
                    mem::swap(&mut num1[0], &mut num2[1]);
                    mem::swap(&mut den1[0], &mut den2[1]);

                    //self.values[0] *= &det;
                    let mut f = BigUint::get_mut(0, types, numerators, denominators);
                    f *= &det;

                    // self.values[3] *= &det;
                    let mut f = BigUint::get_mut(3, types, numerators, denominators);
                    f *= &det;

                    //self.values[2] *= -&det;
                    let mut f = BigUint::get_mut(2, types, numerators, denominators);
                    f *= -&det;

                    //self.values[1] *= -det;
                    let mut f = BigUint::get_mut(1, types, numerators, denominators);
                    f *= -det;
                    return Ok(self);
                }
            }
        };

        // log::debug!("\t\tcompute inverse of {}", self);

        //extend the rows with the identity matrix
        self.push_columns(self.number_of_rows());
        for i in 0..self.number_of_rows() {
            self.set_one(i, self.number_of_rows() + i);
        }

        // println!("add identity\n{}", self);

        //solve
        self = self.gauss_jordan_reduced()?;

        // println!("solved\n{}", self);

        //remove the columns
        self.pop_front_columns(self.number_of_rows());

        // log::info!("inverse done");

        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        fraction_exact::FractionExact,
        fraction_f64::FractionF64,
        matrix::{
            ebi_matrix::EbiMatrix, fraction_matrix_exact::FractionMatrixExact,
            fraction_matrix_f64::FractionMatrixF64, inversion::Inversion,
        },
    };

    #[test]
    fn inverse_f64() {
        let mut m: FractionMatrixF64 = vec![
            vec![1.into(), 0.into(), 0.into(), 0.into()],
            vec![0.into(), 1.into(), 0.into(), FractionF64::from((-3, 5))],
            vec![0.into(), FractionF64::from((-3, 4)), 1.into(), 0.into()],
            vec![0.into(), 0.into(), 0.into(), 1.into()],
        ]
        .try_into()
        .unwrap();

        let i: FractionMatrixF64 = vec![
            vec![1.into(), 0.into(), 0.into(), 0.into()],
            vec![0.into(), 1.into(), 0.into(), FractionF64::from((3, 5))],
            vec![
                0.into(),
                FractionF64::from((3, 4)),
                1.into(),
                FractionF64::from((9, 20)),
            ],
            vec![0.into(), 0.into(), 0.into(), 1.into()],
        ]
        .try_into()
        .unwrap();

        m = m.invert().unwrap();

        // println!("\t\tinverted matrix    {:?}", m);
        // println!("\t\tcorrect inverse    {:?}", i);

        assert!(m.inner_eq(&i));
    }

    #[test]
    fn inverse_biguint() {
        let mut m: FractionMatrixExact = vec![
            vec![1.into(), 0.into(), 0.into(), 0.into()],
            vec![0.into(), 1.into(), 0.into(), (-3, 5).into()],
            vec![0.into(), (-3, 4).into(), 1.into(), 0.into()],
            vec![0.into(), 0.into(), 0.into(), 1.into()],
        ]
        .try_into()
        .unwrap();

        let i: FractionMatrixExact = vec![
            vec![1.into(), 0.into(), 0.into(), 0.into()],
            vec![0.into(), 1.into(), 0.into(), (3, 5).into()],
            vec![0.into(), (3, 4).into(), 1.into(), (9, 20).into()],
            vec![0.into(), 0.into(), 0.into(), 1.into()],
        ]
        .try_into()
        .unwrap();

        m = m.invert().unwrap();

        // println!("\t\tinverted matrix    {:?}", m);
        // println!("\t\tcorrect inverse    {:?}", i);

        assert!(m.inner_eq(&i));
    }
}
