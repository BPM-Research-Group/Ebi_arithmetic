use crate::{
    Zero, fraction::fraction_exact::FractionExact,
    log_polynomial::log_polynomial_exact::LogPolynomialExact,
};
use malachite::{Natural, base::num::basic::traits::Two};
use std::{collections::hash_map::Entry, ops::AddAssign};

impl AddAssign for LogPolynomialExact {
    fn add_assign(&mut self, rhs: Self) {
        for (argument, coefficient) in rhs.argument2coefficient.into_iter() {
            match self.argument2coefficient.entry(argument) {
                Entry::Occupied(mut occupied_entry) => {
                    *occupied_entry.get_mut() += coefficient;
                    if occupied_entry.get().is_zero() {
                        occupied_entry.remove_entry();
                    }
                }
                Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(coefficient);
                }
            };
        }
    }
}

impl AddAssign<&LogPolynomialExact> for LogPolynomialExact {
    fn add_assign(&mut self, rhs: &LogPolynomialExact) {
        for (argument, coefficient) in rhs.argument2coefficient.iter() {
            match self.argument2coefficient.entry(argument.clone()) {
                Entry::Occupied(mut occupied_entry) => {
                    *occupied_entry.get_mut() += coefficient;
                    if occupied_entry.get().is_zero() {
                        occupied_entry.remove_entry();
                    }
                }
                Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(coefficient.clone());
                }
            };
        }
    }
}

impl AddAssign<FractionExact> for LogPolynomialExact {
    fn add_assign(&mut self, rhs: FractionExact) {
        if !rhs.is_zero() {
            match self.argument2coefficient.entry(Natural::TWO) {
                Entry::Occupied(mut occupied_entry) => {
                    *occupied_entry.get_mut() += rhs.0;
                    if occupied_entry.get().is_zero() {
                        occupied_entry.remove_entry();
                    }
                }
                Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(rhs.0);
                }
            };
        }
    }
}

impl AddAssign<&FractionExact> for LogPolynomialExact {
    fn add_assign(&mut self, rhs: &FractionExact) {
        if !rhs.is_zero() {
            match self.argument2coefficient.entry(Natural::TWO) {
                Entry::Occupied(mut occupied_entry) => {
                    *occupied_entry.get_mut() += &rhs.0;
                    if occupied_entry.get().is_zero() {
                        occupied_entry.remove_entry();
                    }
                }
                Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(rhs.0.clone());
                }
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::log_polynomial::log_polynomial_exact::LogPolynomialExact;

    #[test]
    fn add() {
        let mut a = LogPolynomialExact::zero();
        let mut b = LogPolynomialExact::one();
        let c = LogPolynomialExact::from(2);

        a += &b;
        assert_eq!(a, b);

        b += &a;
        assert_eq!(b, c);
    }
}
