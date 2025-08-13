use std::{
    mem,
    ops::{DivAssign, MulAssign, Neg, SubAssign},
    sync::atomic::AtomicBool,
};

use crate::ebi_number::{EbiNumber, Fractional, One, Zero};
use anyhow::{Result, anyhow};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

pub(crate) fn invert<T>(number_of_columns: &mut usize, values: &mut Vec<Vec<T>>) -> Result<()>
where
    T: EbiNumber
        + Fractional
        + for<'a> MulAssign<&'a T>
        + for<'a> DivAssign<&'a T>
        + SubAssign
        + Neg<Output = T>
        + Send,
    //for<'a> &'a T: Div<Output = T>
{
    let number_of_rows = values.len();
    if number_of_columns != &number_of_rows {
        return Err(anyhow!("can only take the inverse of a square matrix"));
    }

    //optimisation: size-zero matrix
    if number_of_rows.is_zero() {
        return Ok(());
    }

    //optimisation: size-one matrix
    if number_of_rows.is_one() {
        if values[0][0].is_zero() {
            return Err(anyhow!("matrix is not invertible"));
        }

        values[0][0] = values[0][0].recip();
        return Ok(());
    }

    //optimisation: size-two matrix
    if number_of_rows == 2 {
        //compute determinant
        let mut det = values[0][0].clone();
        det *= &values[1][1];
        let mut det2 = values[0][1].clone();
        det2 *= &values[1][0];
        det -= det2;

        if det.is_zero() {
            return Err(anyhow!("matrix is not invertible"));
        }

        // log::debug!("determinant {}", det);

        det = det.recip();

        //perform inverse
        let (row1, row2) = values.split_at_mut(1);
        mem::swap(&mut row1[0][0], &mut row2[0][1]);

        values[0][0] *= &det;
        values[1][1] *= &det;

        let mindet = det.clone().neg();
        values[1][0] *= &mindet;
        values[0][1] *= &mindet;

        return Ok(());
    }

    // log::debug!("\t\tcompute inverse of {}", self);

    //extend the rows with the identity matrix
    let n = number_of_rows;
    for (r, row) in values.iter_mut().enumerate() {
        row.extend(vec![T::zero(); n]);
        row[n + r] = T::one();
    }
    *number_of_columns += n;

    // log::debug!("\t\tadd identity       {}", self);

    //solve
    solve(number_of_columns, values)?;

    // log::debug!("\t\tsolved             {}", self);

    //reduce the rows
    for (_, row) in values.iter_mut().enumerate() {
        row.drain(0..n);
    }
    *number_of_columns = n;

    // log::info!("inverse done");

    Ok(())
}

pub fn solve<T>(number_of_columns: &mut usize, values: &mut Vec<Vec<T>>) -> Result<()>
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
        if values[row_a][row_a].is_zero() {
            continue;
        } else {
            for row_b in row_a..values.len() - 1 {
                //optimisation: do not attempt to add a factor of 0
                if !values[row_b + 1][row_a].is_zero() {
                    let mut factor = values[row_b + 1][row_a].clone();
                    factor /= &values[row_a][row_a];
                    // let factor = &values[row_b + 1][row_a] / &values[row_a][row_a];

                    // log::debug!(
                    //     "\t\t\t\t\tfactor row_a {}, row_b {}, {}",
                    //     row_a, row_b, factor
                    // );
                    for column in row_a..*number_of_columns {
                        let mut old = values[row_a][column].clone();
                        old *= &factor;
                        values[row_b + 1][column] -= old;
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
        if values[i][i].is_zero() {
            continue;
        } else {
            for j in (0..i).rev() {
                let mut factor = values[j][i].clone();
                factor /= &values[i][i];
                // let factor = &values[j][i] / &values[i][i];

                for k in i..values[0].len() {
                    let mut old = values[i][k].clone();
                    old *= &factor;
                    values[j][k] -= old;
                }
            }
        }
    }

    // log::debug!("\t\tsecond step        {}", self);

    // log::info!("second step done");

    let failed = AtomicBool::new(false);
    let number_of_columns = values[0].len();
    let number_of_rows = values.len();

    if number_of_rows > 100 {
        values.par_iter_mut().enumerate().for_each(|(i, row)| {
            solve_step_3(row, i, &failed, number_of_rows, number_of_columns);
        });
    } else {
        values.iter_mut().enumerate().for_each(|(i, row)| {
            solve_step_3(row, i, &failed, number_of_rows, number_of_columns);
        });
    }

    if failed.load(std::sync::atomic::Ordering::Relaxed) {
        return Err(anyhow!("matrix is not invertible"));
    }

    // log::info!("third step done");

    Ok(())
}

fn solve_step_3<T>(
    row: &mut Vec<T>,
    i: usize,
    failed: &AtomicBool,
    number_of_rows: usize,
    number_of_columns: usize,
) where
    T: One + Zero + Clone + for<'a> DivAssign<&'a T>,
{
    let factor = row[i].clone();
    if factor.is_zero() {
        failed.store(true, std::sync::atomic::Ordering::Relaxed);
    } else {
        for j in number_of_rows..number_of_columns {
            row[j] /= &factor;
        }
        row[i] = T::one();
    }
}
