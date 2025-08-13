use num::{BigInt, One as NumOne, Signed as NumSigned, Zero as NumZero};

use crate::ebi_number::{EbiNumber, Infinite, Normal, One, Round, Signed, Zero};

impl EbiNumber for BigInt {}

impl One for BigInt {
    fn one() -> Self {
        <BigInt as NumOne>::one()
    }

    fn is_one(&self) -> bool {
        <BigInt as NumOne>::is_one(&self)
    }
}

impl Zero for BigInt {
    fn zero() -> Self {
        <BigInt as NumZero>::zero()
    }

    fn is_zero(&self) -> bool {
        <BigInt as NumZero>::is_zero(&self)
    }
}

impl Signed for BigInt {
    fn abs(&self) -> Self {
        <BigInt as NumSigned>::abs(&self)
    }

    fn is_positive(&self) -> bool {
        <BigInt as NumSigned>::is_positive(&self)
    }

    fn is_negative(&self) -> bool {
        <BigInt as NumSigned>::is_negative(&self)
    }

    fn is_not_negative(&self) -> bool {
        !Signed::is_negative(self)
    }

    fn is_not_positive(&self) -> bool {
        !Signed::is_positive(self)
    }
}

impl Infinite for BigInt {
    fn is_infinite(&self) -> bool {
        false
    }
}

impl Normal for BigInt {
    fn is_nan(&self) -> bool {
        false
    }
}

impl Round for BigInt {
    fn floor(self) -> Self {
        self
    }

    fn ceil(self) -> Self {
        self
    }
}
