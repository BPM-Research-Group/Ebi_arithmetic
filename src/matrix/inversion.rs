use std::{
    mem,
    ops::{DivAssign, MulAssign, SubAssign},
    sync::atomic::AtomicBool,
};

use crate::{
    ebi_number::{One, Zero},
    matrix::{ebi_matrix::EbiMatrix, fraction_matrix_f64::FractionMatrixF64},
};
use anyhow::{Result, anyhow};

pub trait Inversion {
    fn invert(&mut self) -> Result<()>;
}

impl Inversion for FractionMatrixF64 {
    fn invert(&mut self) -> Result<()> {
        if self.number_of_columns() != self.number_of_rows() {
            return Err(anyhow!("can only take the inverse of a square matrix"));
        }

        //optimisation: size-zero matrix
        if self.number_of_rows().is_zero() {
            return Ok(());
        }

        //optimisation: size-one matrix
        if self.number_of_rows().is_one() {
            if self.values[0].is_zero() {
                return Err(anyhow!("matrix is not invertible"));
            }

            self.values[0] = self.values[0].recip();
            return Ok(());
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
            return Ok(());
        }

        // log::debug!("\t\tcompute inverse of {}", self);

        //extend the rows with the identity matrix
        self.push_columns(self.number_of_rows);
        for i in 0..self.number_of_rows {
            let idx_ii = self.index(i, i);
            self.values[idx_ii] = f64::one();
        }

        // log::debug!("\t\tadd identity       {}", self);

        //solve
        solve(self.number_of_columns, &mut self.values)?;

        // log::debug!("\t\tsolved             {}", self);

        //remove the columns
        self.pop_front_columns(self.number_of_rows);

        // log::info!("inverse done");

        Ok(())
    }
}

pub fn solve<T>(number_of_columns: usize, values: &mut Vec<T>) -> Result<()>
where
    T: Sized
        + Zero
        + Clone
        + SubAssign
        + for<'a> MulAssign<&'a T>
        + for<'a> DivAssign<&'a T>
        + Send
        + One,
{
    if values.len() == 0 {
        return Ok(());
    }

    for row_a in 0..values.len() - 1 {
        if values[row_a * number_of_columns + row_a].is_zero() {
            continue;
        } else {
            for row_b in row_a..values.len() - 1 {
                //optimisation: do not attempt to add a factor of 0
                if !values[(row_b + 1) * number_of_columns + row_a].is_zero() {
                    let mut factor = values[(row_b + 1) * number_of_columns + row_a].clone();
                    factor /= &values[row_a * number_of_columns + row_a];
                    // let factor = &values[row_b + 1][row_a] / &values[row_a][row_a];

                    // log::debug!(
                    //     "\t\t\t\t\tfactor row_a {}, row_b {}, {}",
                    //     row_a, row_b, factor
                    // );
                    for column in row_a..number_of_columns {
                        let mut old = values[row_a * number_of_columns + column].clone();
                        old *= &factor;
                        values[(row_b + 1) + number_of_columns + column] -= old;
                    }

                    // log::debug!("\t\t\t       now {}", self);
                }
            }
        }
    }

    // log::debug!("\t\trow-reduced echelon{}", self);

    // log::info!("number of columns {}", self.get_number_of_columns());

    // log::info!("first step done");

    for i in (0..values.len()).rev() {
        if values[i * number_of_columns + i].is_zero() {
            continue;
        } else {
            for j in (0..i).rev() {
                let mut factor = values[j * number_of_columns + i].clone();
                factor /= &values[i * number_of_columns + i];
                // let factor = &values[j][i] / &values[i][i];

                for k in i..number_of_columns {
                    let mut old = values[i * number_of_columns + k].clone();
                    old *= &factor;
                    values[j * number_of_columns + k] -= old;
                }
            }
        }
    }

    // log::debug!("\t\tsecond step        {}", self);

    // log::info!("second step done");

    let failed = AtomicBool::new(false);
    let number_of_rows = values.len() / number_of_columns;

    values
        .chunks_mut(number_of_columns)
        .enumerate()
        .for_each(|(i, row)| {
            let factor = row[i].clone();
            if factor.is_zero() {
                failed.store(true, std::sync::atomic::Ordering::Relaxed);
            } else {
                for j in number_of_rows..number_of_columns {
                    row[j] /= &factor;
                }
                row[i] = T::one();
            }
        });

    if failed.load(std::sync::atomic::Ordering::Relaxed) {
        return Err(anyhow!("matrix is not invertible"));
    }

    // log::info!("third step done");

    Ok(())
}
