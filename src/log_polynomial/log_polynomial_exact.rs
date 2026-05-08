use fnv::FnvBuildHasher;
use malachite::{
    Natural, Rational,
    base::num::basic::traits::{One, Two},
};
use std::{
    collections::{HashMap, hash_map::Entry},
    fmt::Display,
};

#[derive(Debug, Eq, PartialEq, Clone)]
/// A log polynomial is a sum of coefficient * log_2(argument) (R * log(N)) terms.
/// In the background, each term's N is prime (which makes the representation unique).
pub struct LogPolynomialExact {
    pub(crate) argument2coefficient: HashMap<Natural, Rational, FnvBuildHasher>,
}

impl LogPolynomialExact {
    pub fn zero() -> Self {
        Self {
            argument2coefficient: HashMap::<_, _, FnvBuildHasher>::default(),
        }
    }

    pub fn one() -> Self {
        let mut argument2coefficient = HashMap::<_, _, FnvBuildHasher>::default();
        argument2coefficient.insert(Natural::TWO, Rational::ONE);
        Self {
            argument2coefficient,
        }
    }

    /// Constructs a new LogPolynomial with representing coefficient * log_2(argument).
    /// The caller must ensure that argument is prime.
    pub fn add_raw(&mut self, coefficient: Rational, argument: Natural) {
        match self.argument2coefficient.entry(argument) {
            Entry::Occupied(mut occupied_entry) => *occupied_entry.get_mut() += coefficient,
            Entry::Vacant(vacant_entry) => {
                vacant_entry.insert(coefficient);
            }
        };
    }
}

impl Display for LogPolynomialExact {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.argument2coefficient
                .iter()
                .map(|(argument, coefficient)| format!("{coefficient} log({argument})"))
                .collect::<Vec<_>>()
                .join(" + ")
        )
    }
}
