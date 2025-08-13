use crate::{
    ebi_matrix::EbiMatrix, exact::MaybeExact, fraction_f64::FractionF64, inversion::invert,
};
use anyhow::{Result, anyhow};

#[derive(Clone)]
pub struct FractionMatrixF64 {
    values: Vec<Vec<f64>>,
    number_of_columns: usize,
}

impl FractionMatrixF64 {
    /// Obtains an element from the matrix.
    /// This may be an expensive operation.
    pub fn get(&self, row: usize, column: usize) -> FractionF64 {
        FractionF64(self.values[row][column])
    }

    pub fn to_vec(self) -> Result<Vec<Vec<FractionF64>>> {
        Ok(self
            .values
            .into_iter()
            .map(|row| row.into_iter().map(|f| FractionF64(f)).collect())
            .collect())
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

    fn reduce(self) -> Self {
        self
    }

    fn invert(&mut self) -> Result<()> {
        invert(&mut self.number_of_columns, &mut self.values)
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
            values: values
                .into_iter()
                .map(|row| row.into_iter().map(|f| f.0).collect())
                .collect(),
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

#[cfg(test)]
mod tests {
    use crate::{
        ebi_matrix::EbiMatrix, f_a, fraction_f64::FractionF64,
        fraction_matrix_f64::FractionMatrixF64,
    };

    #[test]
    fn fraction_matrix_reversible() {
        let m1 = vec![vec![
            FractionF64::infinity(),
            FractionF64::neg_infinity(),
            f_a!(8, 3),
        ]];

        let m2: FractionMatrixF64 = m1.clone().into();
        let m2 = m2.reduce();

        let m3 = m2.to_vec().unwrap();

        assert_eq!(m1, m3);
    }
}
