use std::fmt::Display;

use crate::{
    ebi_number::Zero, exact::MaybeExact, fraction::EPSILON, fraction_f64::FractionF64,
    matrix::ebi_matrix::EbiMatrix, pop_front_columns, push_columns,
};
use anyhow::{Error, Result, anyhow};
use itertools::Itertools;

#[derive(Clone, Debug)]
pub struct FractionMatrixF64 {
    pub(crate) values: Vec<f64>,
    pub(crate) number_of_rows: usize,
    pub(crate) number_of_columns: usize,
}

impl FractionMatrixF64 {
    /// Obtains an element from the matrix.
    /// This may be an expensive operation depending on the compile mode.
    pub fn get(&self, row: usize, column: usize) -> Option<FractionF64> {
        Some(FractionF64(*self.values.get(self.index(row, column))?))
    }

    /// Transforms the matrix into a multidimensional vector
    pub fn to_vec(self) -> Result<Vec<Vec<FractionF64>>> {
        Ok(self
            .values
            .into_iter()
            .chunks(self.number_of_columns)
            .into_iter()
            .map(|row| row.into_iter().map(|f| FractionF64(f)).collect())
            .collect())
    }

    pub(crate) fn index(&self, row: usize, column: usize) -> usize {
        row * self.number_of_columns + column
    }
}

impl EbiMatrix<FractionF64> for FractionMatrixF64 {
    fn new(number_of_columns: usize) -> Self {
        Self {
            values: vec![],
            number_of_columns,
            number_of_rows: 0,
        }
    }

    fn number_of_rows(&self) -> usize {
        self.number_of_rows
    }

    fn number_of_columns(&self) -> usize {
        self.number_of_columns
    }

    fn reduce(self) -> Self {
        self
    }

    fn eq(&mut self, other: &mut Self) -> bool {
        self.inner_eq(other)
    }

    fn inner_eq(&self, other: &Self) -> bool {
        self.number_of_columns == other.number_of_columns
            && self.number_of_rows == other.number_of_rows
            && self
                .values
                .iter()
                .zip(other.values.iter())
                .all(|(a, b)| (a - b).abs() < EPSILON)
    }

    fn push_columns(&mut self, number_of_columns_to_add: usize) {
        push_columns!(
            f64::zero(),
            number_of_columns_to_add,
            self.values,
            self.number_of_rows,
            self.number_of_columns
        );
        self.number_of_columns += number_of_columns_to_add;
    }

    fn pop_front_columns(&mut self, number_of_columns_to_remove: usize) {
        pop_front_columns!(
            number_of_columns_to_remove,
            self.values,
            self.number_of_rows,
            self.number_of_columns
        );
        self.number_of_columns -= number_of_columns_to_remove;
    }

    fn get(&self, row: usize, column: usize) -> Option<FractionF64> {
        let idx = self.index(row, column);
        Some(FractionF64(*self.values.get(idx)?))
    }

    fn set(&mut self, row: usize, column: usize, value: FractionF64) {
        let idx = self.index(row, column);
        self.values[idx] = value.0;
    }

    fn set_one(&mut self, row: usize, column: usize) {
        let idx = self.index(row, column);
        self.values[idx] = 1f64;
    }
}

impl TryFrom<Vec<Vec<FractionF64>>> for FractionMatrixF64 {
    type Error = Error;
    fn try_from(values: Vec<Vec<FractionF64>>) -> Result<Self> {
        let number_of_columns = if let Some(x) = values.iter().next() {
            x.len()
        } else {
            0
        };
        let number_of_rows = values.len();

        let mut new_values = Vec::with_capacity(values.len() * number_of_columns);
        for row in values {
            if row.len() != number_of_columns {
                return Err(anyhow!("number of columns is not consistent"));
            }

            for v in row {
                new_values.push(v.0);
            }
        }

        Ok(Self {
            values: new_values,
            number_of_rows,
            number_of_columns,
        })
    }
}

impl MaybeExact for FractionMatrixF64 {
    type Approximate = FractionMatrixF64;
    type Exact = ();

    fn is_exact(&self) -> bool {
        false
    }

    fn extract_approx(&self) -> anyhow::Result<&Self::Approximate> {
        Ok(self)
    }

    fn extract_exact(&self) -> anyhow::Result<&Self::Exact> {
        Err(anyhow!("cannot extract a fraction from a float"))
    }
}

impl std::fmt::Display for FractionMatrixF64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{{{")?;
        for (i, row) in self.values.chunks(self.number_of_columns).enumerate() {
            for (j, fraction) in row.iter().enumerate() {
                write!(f, "{}", fraction.to_string())?;
                if j < row.len() - 1 {
                    write!(f, ", ")?;
                }
            }
            if i < self.number_of_rows - 1 {
                write!(f, "}},\n {{")?;
            }
        }
        write!(f, "}}}}")
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        f_a,
        fraction_f64::FractionF64,
        matrix::{ebi_matrix::EbiMatrix, fraction_matrix_f64::FractionMatrixF64},
    };

    #[test]
    fn fraction_matrix_reversible() {
        let m1 = vec![vec![
            FractionF64::infinity(),
            FractionF64::neg_infinity(),
            f_a!(8, 3),
        ]];

        let m2: FractionMatrixF64 = m1.clone().try_into().unwrap();
        let m2 = m2.reduce();

        let m3 = m2.to_vec().unwrap();

        assert_eq!(m1, m3);
    }
}
