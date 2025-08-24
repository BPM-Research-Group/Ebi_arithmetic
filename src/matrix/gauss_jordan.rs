use anyhow::{Result, anyhow};
use num::BigUint;
use std::sync::atomic::AtomicBool;

use crate::{
    fraction::ebi_number::{One, Zero},
    fraction_raw::{getters::FractionRawGetter, one::SetOne, zero::IsZero},
    matrix::{
        ebi_matrix::EbiMatrix, fraction_matrix_enum::FractionMatrixEnum,
        fraction_matrix_exact::FractionMatrixExact, fraction_matrix_f64::FractionMatrixF64,
    },
};

pub trait GaussJordan {
    /// Applies Gaussian elimination to obtain a matrix in row echelon form.
    fn gauss_jordan(&mut self);

    /// Applies Gaussian elimination to obtain a matrix in reduced row echelon form.
    fn gauss_jordan_reduced(self) -> Result<Self>
    where
        Self: Sized;
}

impl GaussJordan for FractionMatrixF64 {
    fn gauss_jordan(&mut self) {
        let number_of_rows = self.number_of_rows();
        let number_of_columns = self.number_of_columns();

        if number_of_rows == 0 || number_of_columns == 0 {
            return;
        }

        for row_a in 0..number_of_rows - 1 {
            if self.values[row_a * number_of_columns + row_a].is_zero() {
                continue;
            } else {
                for row_b in row_a..number_of_rows - 1 {
                    //optimisation: do not attempt to add a factor of 0
                    if !self.values[(row_b + 1) * number_of_columns + row_a].is_zero() {
                        let mut factor =
                            self.values[(row_b + 1) * number_of_columns + row_a].clone();
                        factor /= &self.values[row_a * number_of_columns + row_a];
                        // let factor = &values[row_b + 1][row_a] / &values[row_a][row_a];

                        // println!(
                        //     "\t\t\t\t\tfactor row_a {}, row_b {}, {}",
                        //     row_a, row_b, factor
                        // );
                        for column in row_a..number_of_columns {
                            let mut old = self.values[row_a * number_of_columns + column].clone();
                            old *= &factor;
                            self.values[(row_b + 1) * number_of_columns + column] -= old;
                        }

                        // log::debug!("\t\t\t       now {}", self);
                    }
                }
            }
        }

        // println!("row-reduced echelon\n{:?}", values);

        // log::info!("number of columns {}", self.get_number_of_columns());

        // log::info!("first step done");

        for i in (0..number_of_rows).rev() {
            if self.values[i * number_of_columns + i].is_zero() {
                continue;
            } else {
                for j in (0..i).rev() {
                    let mut factor = self.values[j * number_of_columns + i].clone();
                    factor /= &self.values[i * number_of_columns + i];
                    // let factor = &values[j][i] / &values[i][i];

                    for k in i..number_of_columns {
                        let mut old = self.values[i * number_of_columns + k].clone();
                        old *= &factor;
                        self.values[j * number_of_columns + k] -= old;
                    }
                }
            }
        }

        // log::debug!("\t\tsecond step        {}", self);

        // log::info!("second step done");
    }

    fn gauss_jordan_reduced(mut self) -> Result<Self> {
        self.gauss_jordan();

        let number_of_rows = self.number_of_rows();
        let number_of_columns = self.number_of_columns();

        let failed = AtomicBool::new(false);

        self.values
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
                    row[i] = f64::one();
                }
            });

        if failed.load(std::sync::atomic::Ordering::Relaxed) {
            return Err(anyhow!("matrix has no reduced row-echelon form"));
        }

        // log::info!("third step done");

        Ok(self)
    }
}

macro_rules! gauss_jordan_reduced {
    ($number_of_rows:expr, $number_of_columns:expr, $types:expr, $numerators: expr, $denominators:expr, $t:ident) => {
        for row in 0..*$number_of_rows {
            let factor = $t::get_clone(
                row * *$number_of_columns + row,
                &$types,
                &$numerators,
                &$denominators,
            );
            if factor.is_zero() {
                return Err(anyhow!("matrix has no reduced row-echelon form"));
            } else {
                for j in *$number_of_rows..*$number_of_columns {
                    let mut f = $t::get_mut(
                        row * *$number_of_columns + j,
                        $types,
                        $numerators,
                        $denominators,
                    );
                    f /= &factor;
                }

                let mut f = $t::get_mut(
                    row * *$number_of_columns + row,
                    $types,
                    $numerators,
                    $denominators,
                );
                f.set_one();
            }
        }
    };
}

impl GaussJordan for FractionMatrixExact {
    fn gauss_jordan(&mut self) {
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
                if *number_of_rows == 0 || *number_of_columns == 0 {
                    return;
                }

                for row_a in 0..*number_of_rows - 1 {
                    if BigUint::get_ref(
                        row_a * *number_of_columns + row_a,
                        types,
                        numerators,
                        denominators,
                    )
                    .is_zero()
                    {
                        continue;
                    } else {
                        for row_b in row_a..*number_of_rows - 1 {
                            //optimisation: do not attempt to add a factor of 0
                            if !BigUint::get_ref(
                                (row_b + 1) * *number_of_columns + row_a,
                                types,
                                numerators,
                                denominators,
                            )
                            .is_zero()
                            {
                                let mut factor = BigUint::get_clone(
                                    (row_b + 1) * *number_of_columns + row_a,
                                    types,
                                    numerators,
                                    denominators,
                                );
                                factor /= BigUint::get_ref(
                                    row_a * *number_of_columns + row_a,
                                    types,
                                    numerators,
                                    denominators,
                                );
                                // let factor = &values[row_b + 1][row_a] / &values[row_a][row_a];

                                // println!(
                                //     "\t\t\t\t\tfactor row_a {}, row_b {}, {}",
                                //     row_a, row_b, factor
                                // );
                                for column in row_a..*number_of_columns {
                                    let mut old = BigUint::get_clone(
                                        row_a * *number_of_columns + column,
                                        types,
                                        numerators,
                                        denominators,
                                    );
                                    old *= &factor;
                                    let mut f = BigUint::get_mut(
                                        (row_b + 1) * *number_of_columns + column,
                                        types,
                                        numerators,
                                        denominators,
                                    );
                                    f -= old;
                                }

                                // log::debug!("\t\t\t       now {}", self);
                            }
                        }
                    }
                }

                // println!("row-reduced echelon\n{:?}", values);

                // log::info!("number of columns {}", self.get_number_of_columns());

                // log::info!("first step done");

                for i in (0..*number_of_rows).rev() {
                    if BigUint::get_ref(
                        i * *number_of_columns + i,
                        types,
                        &numerators,
                        &denominators,
                    )
                    .is_zero()
                    {
                        continue;
                    } else {
                        for j in (0..i).rev() {
                            let mut factor = BigUint::get_clone(
                                j * *number_of_columns + i,
                                types,
                                &numerators,
                                &denominators,
                            );
                            factor /= BigUint::get_ref(
                                i * *number_of_columns + i,
                                types,
                                &numerators,
                                &denominators,
                            );
                            // let factor = &values[j][i] / &values[i][i];

                            for k in i..*number_of_columns {
                                let mut old = BigUint::get_clone(
                                    i * *number_of_columns + k,
                                    types,
                                    &numerators,
                                    &denominators,
                                );
                                old *= &factor;
                                let mut f = BigUint::get_mut(
                                    j * *number_of_columns + k,
                                    types,
                                    numerators,
                                    denominators,
                                );
                                f -= old;
                            }
                        }
                    }
                }
            }
        }
    }

    fn gauss_jordan_reduced(mut self) -> Result<Self> {
        self.gauss_jordan();

        match &mut self {
            FractionMatrixExact::U64 {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => {
                gauss_jordan_reduced!(
                    number_of_rows,
                    number_of_columns,
                    types,
                    numerators,
                    denominators,
                    u64
                );
            }
            FractionMatrixExact::BigInt {
                number_of_columns,
                number_of_rows,
                types,
                numerators,
                denominators,
            } => {
                gauss_jordan_reduced!(
                    number_of_rows,
                    number_of_columns,
                    types,
                    numerators,
                    denominators,
                    BigUint
                );

                // log::info!("third step done");
            }
        };
        Ok(self)
    }
}

impl GaussJordan for FractionMatrixEnum {
    fn gauss_jordan(&mut self) {
        match self {
            FractionMatrixEnum::Approx(m) => m.gauss_jordan(),
            FractionMatrixEnum::Exact(m) => m.gauss_jordan(),
            FractionMatrixEnum::CannotCombineExactAndApprox => {}
        }
    }

    fn gauss_jordan_reduced(self) -> Result<Self> {
        match self {
            FractionMatrixEnum::Approx(m) => {
                Ok(FractionMatrixEnum::Approx(m.gauss_jordan_reduced()?))
            }
            FractionMatrixEnum::Exact(m) => {
                Ok(FractionMatrixEnum::Exact(m.gauss_jordan_reduced()?))
            }
            FractionMatrixEnum::CannotCombineExactAndApprox => {
                return Err(anyhow!("cannot combine exact and approximate arithmetic"));
            }
        }
    }
}
