use crate::{ebi_matrix::EbiMatrix, exact::MaybeExact, fraction_f64::FractionF64};
use anyhow::anyhow;

#[derive(Clone)]
pub struct FractionMatrixF64 {
    values: Vec<Vec<FractionF64>>,
    number_of_columns: usize,
}

impl FractionMatrixF64 {
    /// Obtains an element from the matrix.
    /// This may be an expensive operation.
    pub fn get(&self, row: usize, column: usize) -> FractionF64 {
        self.values[row][column]
    }
}

impl EbiMatrix for FractionMatrixF64 {
    fn new(number_of_columns: usize) -> Self {
        Self {
            values: vec![],
            number_of_columns: number_of_columns,
        }
    }

    fn number_of_rows(&self) -> usize {
        self.values.len()
    }

    fn number_of_columns(&self) -> usize {
        self.number_of_columns
    }

    fn optimise(self) -> Self {
        self
    }
}

impl From<Vec<Vec<FractionF64>>> for FractionMatrixF64 {
    fn from(values: Vec<Vec<FractionF64>>) -> Self {
        let number_of_columns = if let Some(x) = values.iter().next() {
            x.len()
        } else {
            0
        };
        Self {
            values,
            number_of_columns,
        }
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
