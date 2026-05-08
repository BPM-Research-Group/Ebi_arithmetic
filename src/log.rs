use anyhow::Result;

pub trait Log<T> {
    /// Returns the logarithm of the given argument: log(argument).
    /// Returns an error if the argument is not positive.
    /// 
    /// This is a potentially expensive operation, as the prime factors of the argument may be computed.
    fn log(argument: T) -> Result<Self>
    where
        Self: Sized;
}
