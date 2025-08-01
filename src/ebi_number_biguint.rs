use num::{BigUint, One as NumOne, Zero as NumZero};

use crate::ebi_number::{EbiNumber, Infinite, Normal, One, Round, Signed, Zero};

impl EbiNumber for BigUint {}

impl Zero for BigUint {
    fn zero() -> Self {
        NumZero::zero()
    }

    fn is_zero(&self) -> bool {
        NumZero::is_zero(self)
    }
}

impl One for BigUint {
    fn one() -> Self {
        NumOne::one()
    }

    fn is_one(&self) -> bool {
        NumOne::is_one(self)
    }
}

impl Signed for BigUint {
    fn abs(&self) -> Self {
        self.clone()
    }

    fn is_positive(&self) -> bool {
        !Zero::is_zero(self)
    }

    fn is_negative(&self) -> bool {
        false
    }
}

impl Infinite for BigUint {
    fn is_infinite(&self) -> bool {
        false
    }
}

impl Normal for BigUint {
    fn is_nan(&self) -> bool {
        false
    }
}

impl Round for BigUint {
    fn floor(self) -> Self {
        self
    }

    fn ceil(self) -> Self {
        self
    }
}
