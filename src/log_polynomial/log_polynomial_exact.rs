use malachite::{Natural, Rational};
use std::{collections::HashMap, fmt::Display};

#[derive(Debug)]
pub struct LogPolynomialExact {
    pub(crate) argument2coefficient: HashMap<Natural, Rational>,
}

impl LogPolynomialExact {
    pub fn zero() -> Self {
        Self {
            argument2coefficient: HashMap::new(),
        }
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
