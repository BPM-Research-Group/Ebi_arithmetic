use crate::{ebi_matrix::EbiMatrix, exact::MaybeExact, fraction_f64::FractionF64};
use anyhow::{Error, Result, anyhow};

#[derive(Clone, PartialEq, Debug)]
pub struct FractionMatrixF64 {
    pub(crate) values: Vec<Vec<f64>>,
    pub(crate) number_of_columns: usize,
}

impl FractionMatrixF64 {
    /// Obtains an element from the matrix.
    /// This may be an expensive operation depending on the compile mode.
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

    fn eq(&mut self, other: &mut Self) -> bool {
        self == other
    }

    fn inner_eq(&self, other: &Self) -> bool {
        self == other
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

        let mut new_values = Vec::with_capacity(values.len());
        for row in values {
            if row.len() != number_of_columns {
                return Err(anyhow!("number of columns is not consistent"));
            }

            let mut new_row = Vec::with_capacity(number_of_columns);
            for v in row {
                new_row.push(v.0);
            }

            new_values.push(new_row);
        }

        Ok(Self {
            values: new_values,
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

        let m2: FractionMatrixF64 = m1.clone().try_into().unwrap();
        let m2 = m2.reduce();

        let m3 = m2.to_vec().unwrap();

        assert_eq!(m1, m3);
    }
}
