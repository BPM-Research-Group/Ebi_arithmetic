use malachite::{Natural, Rational, base::num::basic::traits::Two};

use crate::{
    Zero,
    fraction::{fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64},
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact, log_polynomial_f64::LogPolynomialF64
    },
};
use std::{collections::hash_map::Entry, ops::SubAssign};

impl SubAssign for LogPolynomialExact {
    fn sub_assign(&mut self, rhs: Self) {
        for (argument, coefficient) in rhs.argument2coefficient.into_iter() {
            match self.argument2coefficient.entry(argument) {
                Entry::Occupied(mut occupied_entry) => {
                    *occupied_entry.get_mut() -= coefficient;
                    if occupied_entry.get().is_zero() {
                        occupied_entry.remove_entry();
                    }
                }
                Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(-coefficient);
                }
            };
        }
    }
}

impl SubAssign<FractionExact> for LogPolynomialExact {
    fn sub_assign(&mut self, rhs: FractionExact) {
        match self.argument2coefficient.entry(Natural::TWO) {
            Entry::Occupied(mut occupied_entry) => {
                *occupied_entry.get_mut() -= rhs.0;
                if occupied_entry.get().is_zero() {
                    occupied_entry.remove_entry();
                }
            }
            Entry::Vacant(vacant_entry) => {
                vacant_entry.insert(-rhs.0);
            }
        };
    }
}

impl SubAssign<&FractionExact> for LogPolynomialExact {
    fn sub_assign(&mut self, rhs: &FractionExact) {
        match self.argument2coefficient.entry(Natural::TWO) {
            Entry::Occupied(mut occupied_entry) => {
                *occupied_entry.get_mut() -= &rhs.0;
                if occupied_entry.get().is_zero() {
                    occupied_entry.remove_entry();
                }
            }
            Entry::Vacant(vacant_entry) => {
                vacant_entry.insert(-&rhs.0);
            }
        };
    }
}

impl SubAssign<Rational> for LogPolynomialExact {
    fn sub_assign(&mut self, rhs: Rational) {
        match self.argument2coefficient.entry(Natural::TWO) {
            Entry::Occupied(mut occupied_entry) => {
                *occupied_entry.get_mut() -= &rhs;
                if occupied_entry.get().is_zero() {
                    occupied_entry.remove_entry();
                }
            }
            Entry::Vacant(vacant_entry) => {
                vacant_entry.insert(-rhs);
            }
        };
    }
}

impl SubAssign<&Rational> for LogPolynomialExact {
    fn sub_assign(&mut self, rhs: &Rational) {
        match self.argument2coefficient.entry(Natural::TWO) {
            Entry::Occupied(mut occupied_entry) => {
                *occupied_entry.get_mut() -= rhs;
                if occupied_entry.get().is_zero() {
                    occupied_entry.remove_entry();
                }
            }
            Entry::Vacant(vacant_entry) => {
                vacant_entry.insert(-rhs);
            }
        };
    }
}

impl SubAssign for LogPolynomialF64 {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0
    }
}

impl SubAssign<FractionF64> for LogPolynomialF64 {
    fn sub_assign(&mut self, rhs: FractionF64) {
        self.0 -= rhs.0
    }
}

impl SubAssign<&FractionF64> for LogPolynomialF64 {
    fn sub_assign(&mut self, rhs: &FractionF64) {
        self.0 -= rhs.0
    }
}

impl SubAssign<f64> for LogPolynomialF64 {
    fn sub_assign(&mut self, rhs: f64) {
        self.0 -= rhs
    }
}

impl SubAssign<&f64> for LogPolynomialF64 {
    fn sub_assign(&mut self, rhs: &f64) {
        self.0 -= rhs
    }
}

impl SubAssign for LogPolynomialEnum {
    fn sub_assign(&mut self, rhs: Self) {
        if self.matches(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), LogPolynomialEnum::Approx(f)) => lp.sub_assign(f),
                (LogPolynomialEnum::Exact(lp), LogPolynomialEnum::Exact(f)) => lp.sub_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

impl SubAssign<FractionEnum> for LogPolynomialEnum {
    fn sub_assign(&mut self, rhs: FractionEnum) {
        if self.matches_fraction(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), FractionEnum::Approx(f)) => lp.sub_assign(f),
                (LogPolynomialEnum::Exact(lp), FractionEnum::Exact(f)) => lp.sub_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

impl SubAssign<&FractionEnum> for LogPolynomialEnum {
    fn sub_assign(&mut self, rhs: &FractionEnum) {
        if self.matches_fraction(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), FractionEnum::Approx(f)) => lp.sub_assign(f),
                (LogPolynomialEnum::Exact(lp), FractionEnum::Exact(f)) => lp.sub_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}