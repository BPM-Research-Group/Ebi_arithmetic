use anyhow::Result;

pub trait Log<T> {
    /// Returns the logarithm of the given argument: log(argument).
    /// Returns an error if the argument is not positive.
    fn log(argument: T) -> Result<Self>
    where
        Self: Sized;
}
