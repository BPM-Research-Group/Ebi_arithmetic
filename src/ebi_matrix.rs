use crate::{
    exact::MaybeExact
};
use anyhow::Result;

pub trait EbiMatrix<T>:
    Clone + MaybeExact + IdentityMinus + GaussJordan + TryFrom<Vec<Vec<T>>> + Eq
where
    T: Clone,
{
    /// Creates a new matrix with each value initialised to the given value.
    fn new(number_of_rows: usize, number_of_columns: usize, value: T) -> Result<Self>;

    /// Returns the number of rows
    fn number_of_rows(&self) -> usize;

    /// Returns the number of columns
    fn number_of_columns(&self) -> usize;

    /// Gets a particular value of the matrix, if it exists.
    /// This may be an expensive operation.
    fn get(&self, row: usize, column: usize) -> Option<T>;

    /// Sets a particular value of the matrix, if the row and column exist.
    /// If row and column do not exist, behaviour is undefined, and may panic.
    /// Prefer set_one and set_zero if possible.
    fn set(&mut self, row: usize, column: usize, value: T);

    /// Sets a particular value of the matrix to zero, if the row and column exist.
    /// If row and column do not exist, behaviour is undefined, and may panic.
    fn set_zero(&mut self, row: usize, column: usize);

    /// Sets a particular value of the matrix to one, if the row and column exist.
    /// If row and column do not exist, behaviour is undefined, and may panic.
    fn set_one(&mut self, row: usize, column: usize);

    /// Adds a number of columns to the right side of the matrix.
    /// Initially, these added columns will be filled with zeroes.
    fn push_columns(&mut self, number_of_columns_to_add: usize);

    /// Removes columsn from the left of the matrix.
    fn pop_front_columns(&mut self, number_of_columns_to_remove: usize);

    /// Returns a vector of the matrix
    fn to_vec(self) -> Vec<Vec<T>>;
}

pub trait IdentityMinus {
    /// For a given matrix M, computes I-M.
    /// The matrix does not need to be squared.
    fn identity_minus(&mut self) -> Result<()>;
}

pub trait Inversion {
    fn invert(self) -> Result<Self>
    where
        Self: Sized;
}

pub trait GaussJordan {
    /// Applies Gaussian elimination to obtain a matrix in row echelon form.
    fn gauss_jordan(&mut self);

    /// Applies Gaussian elimination to obtain a matrix in reduced row echelon form.
    fn gauss_jordan_reduced(self) -> Result<Self>
    where
        Self: Sized;
}