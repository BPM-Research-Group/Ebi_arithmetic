use fraction::BigFraction;

use crate::ebi_number::{EbiNumber, Fractional, Infinite, Normal, One, Round, Signed, Zero};

impl EbiNumber for BigFraction {}

impl Zero for BigFraction {
    fn zero() -> Self {
        num::Zero::zero()
    }

    fn is_zero(&self) -> bool {
        num::Zero::is_zero(self)
    }
}

impl One for BigFraction {
    fn one() -> Self {
        num::One::one()
    }

    fn is_one(&self) -> bool {
        num::One::is_one(self)
    }
}

impl Signed for BigFraction {
    fn abs(&self) -> Self {
        BigFraction::abs(&self)
    }

    fn is_positive(&self) -> bool {
        self.is_sign_positive() && !Zero::is_zero(self)
    }

    fn is_negative(&self) -> bool {
        self.is_sign_negative() && !Zero::is_zero(self)
    }
}

impl Infinite for BigFraction {
    fn is_infinite(&self) -> bool {
        BigFraction::is_infinite(&self)
    }
}

impl Normal for BigFraction {
    fn is_nan(&self) -> bool {
        BigFraction::is_nan(&self)
    }
}

impl Round for BigFraction {
    fn floor(self) -> Self {
        BigFraction::floor(&self)
    }

    fn ceil(self) -> Self {
        BigFraction::ceil(&self)
    }
}

impl Fractional for BigFraction {
    fn recip(&self) -> Self {
        BigFraction::recip(&self)
    }
}