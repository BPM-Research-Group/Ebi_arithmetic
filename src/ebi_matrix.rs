use crate::exact::MaybeExact;

pub trait EbiMatrix: Clone + MaybeExact {
    /// Creates a new matrix with no rows and the given number of columns.
    fn new(number_of_columns: usize) -> Self;

    /// Returns the number of rows
    fn number_of_rows(&self) -> usize;

    /// Returns the number of columns
    fn number_of_columns(&self) -> usize;

    /// Optimises the matrix for repeated operations.
    /// This may be a costly one-time call, but has the potential to make subsequent operations much more efficient.
    ///
    /// For exact numbers, this will extract the lowest common multiple of the denominators, such that only the nominators need to be stored.
    /// Has foreseeably no effect on approximate arithmetic.
    fn optimise(self) -> Self;
}
