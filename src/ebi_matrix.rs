use crate::exact::MaybeExact;

pub trait EbiMatrix: Clone + MaybeExact {
    /// Creates a new reduced matrix with no rows and the given number of columns.
    fn new(number_of_columns: usize) -> Self;

    /// Returns the number of rows
    fn number_of_rows(&self) -> usize;

    /// Returns the number of columns
    fn number_of_columns(&self) -> usize;

    /// Reduces all fractions in the matrix, and stores them more compactly if possible.
    /// This may be a costly call, but has the potential to make subsequent operations more efficient.
    ///
    /// Has forseeably no effect on approximate arithmetic.
    fn reduce(self) -> Self;

    /// Tests for equivalence.
    /// Reduces both matrices first, so may be an expensive operation.
    fn eq(&mut self, other: &mut Self) -> bool;

    /// Tests for internal equivalence.
    /// Returns true if both matrices are internally the same.
    ///
    /// use `eq` to test for numerical equivalence.
    fn inner_eq(&self, other: &Self) -> bool;
}
