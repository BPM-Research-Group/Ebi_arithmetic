use crate::{Signed, fraction::approximate::Approximate};
use anyhow::Result;
use fnv::FnvBuildHasher;
use malachite::{
    Natural, Rational,
    base::num::basic::traits::{NegativeOne, One, Two},
};
use std::{
    collections::{HashMap, hash_map::Entry},
    fmt::Display,
    io::Write,
};

#[derive(Debug, Eq, PartialEq, Clone)]
/// A log polynomial is a sum of coefficient * log_2(argument) (R * log(N)) terms.
///
/// Internal invariant: no coefficient is zero.
/// Internal invariant: each argument is prime (which makes the representation unique).
pub struct LogPolynomialExact {
    pub(crate) argument2coefficient: HashMap<Natural, Rational, FnvBuildHasher>,
}

impl LogPolynomialExact {
    pub fn export(&self, f: &mut dyn Write) -> Result<()> {
        writeln!(f, "{}", self)?;
        Ok(writeln!(
            f,
            "Approximately {}",
            self.clone().approximate()?
        )?)
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
        let mut list = self.argument2coefficient.iter().collect::<Vec<_>>();
        list.sort();

        if list.is_empty() {
            return write!(f, "0");
        }

        let display = |argument, coefficient| {
            if argument == &Natural::TWO {
                format!("{coefficient}")
            } else if coefficient == &Rational::ONE {
                format!("log({argument})")
            } else if coefficient == &Rational::NEGATIVE_ONE {
                format!("-log({argument})")
            } else {
                format!("{coefficient} * log({argument})")
            }
        };

        let (arg0, coe0) = list.remove(0);
        write!(f, "{}", display(arg0, coe0))?;

        for (arg, coe) in list {
            if coe.is_not_negative() {
                write!(f, " + {}", display(arg, coe))?;
            } else {
                let mut s = display(arg, coe);
                s.remove(0);
                write!(f, " - {}", s)?;
            }
        }

        write!(f, "")
    }
}
