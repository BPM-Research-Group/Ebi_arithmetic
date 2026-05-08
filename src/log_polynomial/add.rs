use crate::log_polynomial::log_polynomial_exact::LogPolynomialExact;
use std::{collections::hash_map::Entry, ops::AddAssign};

impl AddAssign for LogPolynomialExact {
    fn add_assign(&mut self, rhs: Self) {
        for (argument, coefficient) in rhs.argument2coefficient.into_iter() {
            match self.argument2coefficient.entry(argument) {
                Entry::Occupied(mut occupied_entry) => *occupied_entry.get_mut() += coefficient,
                Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(coefficient);
                }
            };
        }
    }
}
