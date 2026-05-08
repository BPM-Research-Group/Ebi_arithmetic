use crate::{
    Zero,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use malachite::{Natural, Rational, base::num::basic::traits::Two};
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

impl AddAssign<Rational> for LogPolynomialExact {
    fn add_assign(&mut self, rhs: Rational) {
        if !rhs.is_zero() {
            match self.argument2coefficient.entry(Natural::TWO) {
                Entry::Occupied(mut occupied_entry) => {
                    *occupied_entry.get_mut() += rhs;
                    if occupied_entry.get().is_zero() {
                        occupied_entry.remove_entry();
                    }
                }
                Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(rhs);
                }
            };
        }
    }
}

impl AddAssign<&Rational> for LogPolynomialExact {
    fn add_assign(&mut self, rhs: &Rational) {
        if !rhs.is_zero() {
            match self.argument2coefficient.entry(Natural::TWO) {
                Entry::Occupied(mut occupied_entry) => {
                    *occupied_entry.get_mut() += rhs;
                    if occupied_entry.get().is_zero() {
                        occupied_entry.remove_entry();
                    }
                }
                Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(rhs.clone());
                }
            };
        }
    }
}

impl AddAssign for LogPolynomialF64 {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}

impl AddAssign<&LogPolynomialF64> for LogPolynomialF64 {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += rhs.0
    }
}

impl AddAssign<FractionF64> for LogPolynomialF64 {
    fn add_assign(&mut self, rhs: FractionF64) {
        self.0 += rhs.0
    }
}

impl AddAssign<&FractionF64> for LogPolynomialF64 {
    fn add_assign(&mut self, rhs: &FractionF64) {
        self.0 += rhs.0
    }
}

impl AddAssign<f64> for LogPolynomialF64 {
    fn add_assign(&mut self, rhs: f64) {
        self.0 += rhs
    }
}

impl AddAssign<&f64> for LogPolynomialF64 {
    fn add_assign(&mut self, rhs: &f64) {
        self.0 += rhs
    }
}

impl AddAssign for LogPolynomialEnum {
    fn add_assign(&mut self, rhs: Self) {
        if self.matches(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), LogPolynomialEnum::Approx(f)) => lp.add_assign(f),
                (LogPolynomialEnum::Exact(lp), LogPolynomialEnum::Exact(f)) => lp.add_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

impl AddAssign<&LogPolynomialEnum> for LogPolynomialEnum {
    fn add_assign(&mut self, rhs: &Self) {
        if self.matches(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), LogPolynomialEnum::Approx(f)) => lp.add_assign(f),
                (LogPolynomialEnum::Exact(lp), LogPolynomialEnum::Exact(f)) => lp.add_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

impl AddAssign<FractionEnum> for LogPolynomialEnum {
    fn add_assign(&mut self, rhs: FractionEnum) {
        if self.matches_fraction(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), FractionEnum::Approx(f)) => lp.add_assign(f),
                (LogPolynomialEnum::Exact(lp), FractionEnum::Exact(f)) => lp.add_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

impl AddAssign<&FractionEnum> for LogPolynomialEnum {
    fn add_assign(&mut self, rhs: &FractionEnum) {
        if self.matches_fraction(&rhs) {
            match (self, rhs) {
                (LogPolynomialEnum::Approx(lp), FractionEnum::Approx(f)) => lp.add_assign(f),
                (LogPolynomialEnum::Exact(lp), FractionEnum::Exact(f)) => lp.add_assign(f),
                _ => {}
            }
        } else {
            *self = LogPolynomialEnum::CannotCombineExactAndApprox
        }
    }
}

macro_rules! addassign_primitive {
    ($t:ty) => {
        impl AddAssign<$t> for LogPolynomialExact {
            fn add_assign(&mut self, rhs: $t) {
                self.add_assign(Rational::from(rhs))
            }
        }

        impl AddAssign<$t> for LogPolynomialF64 {
            fn add_assign(&mut self, rhs: $t) {
                self.add_assign(rhs as f64)
            }
        }

        impl AddAssign<$t> for LogPolynomialEnum {
            fn add_assign(&mut self, rhs: $t) {
                match self {
                    LogPolynomialEnum::Exact(lp) => lp.add_assign(rhs),
                    LogPolynomialEnum::Approx(lp) => lp.add_assign(rhs),
                    LogPolynomialEnum::CannotCombineExactAndApprox => {}
                }
            }
        }
    };
}
addassign_primitive!(usize);
addassign_primitive!(u128);
addassign_primitive!(u64);
addassign_primitive!(u32);
addassign_primitive!(u16);
addassign_primitive!(u8);
addassign_primitive!(i128);
addassign_primitive!(i64);
addassign_primitive!(i32);
addassign_primitive!(i16);
addassign_primitive!(i8);

#[cfg(test)]
mod tests {
    use crate::{One, Zero, log_polynomial::log_polynomial::LogPolynomial};

    #[test]
    fn add() {
        let mut a = LogPolynomial::zero();
        let mut b = LogPolynomial::one();
        let c = LogPolynomial::from(2);

        a += &b;
        assert_eq!(a, b);

        b += &a;
        assert_eq!(b, c);
    }
}
