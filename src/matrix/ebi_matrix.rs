use crate::{
    exact::MaybeExact,
    matrix::{gauss_jordan::GaussJordan, identity_minus::IdentityMinus},
};

pub trait EbiMatrix<T>: Clone + MaybeExact + IdentityMinus + GaussJordan {
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

    /// Gets a particular value of the matrix, if it exists.
    /// This may be an expensive operation.
    fn get(&self, row: usize, column: usize) -> Option<T>;

    /// Sets a particular value of the matrix, if the row and column exist.
    /// If row and column do not exist, behaviour is undefined, and may panic.
    /// Prefer set_one if one is the value you're setting.
    fn set(&mut self, row: usize, column: usize, value: T);

    /// Sets a particular value of the matrix to one, if the row and column exist.
    /// If row and column do not exist, behaviour is undefined, and may panic.
    fn set_one(&mut self, row: usize, column: usize);
}
