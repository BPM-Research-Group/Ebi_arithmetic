use crate::log_polynomial::log_polynomial::LogPolynomial;
use anyhow::Result;

pub trait Log {
    /// Returns the 2-logarithm of the given argument: log_2(argument).
    /// Returns an error if the argument is not positive.
    ///
    /// This is a potentially expensive operation, as the prime factors of the argument may be computed.
    /// May return an error if the argument is too large, that is, for now, cannot be represented by an u128.
    ///
    /// If a multiplication with argument is foreseen, then the n_log_n function is more efficient.
    fn log(&self) -> Result<LogPolynomial>
    where
        Self: Sized;

    /// Returns the value argument * log_2(argument).
    /// Returns an error if the argument is not positive.
    ///
    /// This is a potentially expensive operation, as the prime factors of the argument may be computed.
    /// May return an error if the argument is too large, that is, for now, cannot be represented by an u128.
    fn n_log_n(&self) -> Result<LogPolynomial>
    where
        Self: Sized;
}

impl<T> Log for T
where
    LogPolynomial: LogOf<T>,
{
    fn log(&self) -> Result<LogPolynomial>
    where
        Self: Sized,
    {
        LogPolynomial::log_of(self)
    }

    fn n_log_n(&self) -> Result<LogPolynomial>
    where
        Self: Sized,
    {
        LogPolynomial::n_log_n_of(self)
    }
}

pub trait LogOf<T> {
    /// Returns the 2-logarithm of the given argument: log_2(argument).
    /// Returns an error if the argument is not positive.
    ///
    /// This is a potentially expensive operation, as the prime factors of the argument may be computed.
    /// May return an error if the argument is too large, that is, for now, cannot be represented by an u128.
    ///
    /// If a multiplication with argument is foreseen, then the n_log_n function is more efficient.
    fn log_of(argument: &T) -> Result<Self>
    where
        Self: Sized;

    /// Returns the value argument * log_2(argument).
    /// Returns an error if the argument is not positive.
    ///
    /// This is a potentially expensive operation, as the prime factors of the argument may be computed.
    /// May return an error if the argument is too large, that is, for now, cannot be represented by an u128.
    fn n_log_n_of(argument: &T) -> Result<Self>
    where
        Self: Sized;
}
