use crate::{
    Random,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
    is_exact_globally,
};
use malachite::{
    Rational,
    base::{
        num::basic::traits::{One, Zero},
        random::Seed,
    },
    rational::random::random_rational_range,
};
use rand::Rng;

impl Random for FractionExact {
    fn random_non_zero_probability(bit_length: u64, seed: Seed) -> Self {
        let mut range = random_rational_range(
            seed,
            Rational::ZERO,
            Rational::ONE,
            bit_length,
            1,
            bit_length,
            1,
        );
        Self(Rational::ONE - range.next().unwrap())
    }
}

impl Random for FractionF64 {
    fn random_non_zero_probability(_bit_length: u64, seed: Seed) -> Self {
        Self(1.0 - seed.get_rng().random_range(0.0..1.0))
    }
}

impl Random for FractionEnum {
    fn random_non_zero_probability(bit_length: u64, seed: Seed) -> Self {
        if is_exact_globally() {
            let mut range = random_rational_range(
                seed,
                Rational::ZERO,
                Rational::ONE,
                bit_length,
                1,
                bit_length,
                1,
            );
            Self::Exact(Rational::ONE - range.next().unwrap())
        } else {
            Self::Approx(1.0 - seed.get_rng().random_range(0.0..1.0))
        }
    }
}
