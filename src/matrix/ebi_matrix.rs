use crate::{exact::MaybeExact, matrix::identity_minus::IdentityMinus};

pub trait EbiMatrix: Clone + MaybeExact + IdentityMinus {
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

    /// Adds a number of columns to the right side of the matrix.
    /// Initially, these added columns will be filled with zeroes.
    fn push_columns(&mut self, number_of_columns_to_add: usize);

    /// Removes columsn from the left of the matrix.
    fn pop_front_columns(&mut self, number_of_columns_to_remove: usize);
}

// shared macros
#[macro_export]
macro_rules! push_columns {
    ($zero:expr, $number_of_columns_to_add:expr, $values:expr, $number_of_rows:expr, $number_of_columns:expr) => {
        for row in (0..$number_of_rows).rev() {
            $values.splice(
                row * $number_of_columns + $number_of_columns
                    ..row * $number_of_columns + $number_of_columns,
                vec![$zero; $number_of_columns_to_add],
            );
        }
    };
}

#[macro_export]
macro_rules! pop_front_columns {
    ($number_of_columns_to_remove:expr, $values:expr, $number_of_rows:expr, $number_of_columns:expr) => {
        for row in (0..$number_of_rows).rev() {
            $values.drain(
                row * $number_of_columns..row * $number_of_columns + $number_of_columns_to_remove,
            );
        }
    };
}
