// ============ implementations ============

use fraction::BigFraction;
use num::{BigInt, BigUint, Float, One as NumOne, Signed as NumSigned, Zero as NumZero};

use crate::{fraction::EPSILON, traits::{Infinite, One, Signed, Zero}};

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

impl Zero for BigUint {
    fn zero() -> Self {
        num::Zero::zero()
    }

    fn is_zero(&self) -> bool {
        num::Zero::is_zero(self)
    }
}

impl One for BigUint {
    fn one() -> Self {
        num::One::one()
    }

    fn is_one(&self) -> bool {
        num::One::is_one(self)
    }
}

impl One for f64 {
    fn one() -> Self {
        1.0
    }

    fn is_one(&self) -> bool {
        (self - 1.0).abs() - &EPSILON < 0.0
    }
}

impl Zero for f64 {
    fn zero() -> Self {
        0.0
    }

    fn is_zero(&self) -> bool {
        <f64 as Signed>::abs(&self) - &EPSILON < 0.0
    }
}

impl Signed for f64 {
    fn abs(&self) -> Self {
        <f64 as Float>::abs(*self)
    }

    fn is_positive(&self) -> bool {
        *self != 0f64 && self > &EPSILON
    }

    fn is_negative(&self) -> bool {
        *self != 0f64 && self < &-EPSILON
    }

    fn is_not_negative(&self) -> bool {
        self > &-EPSILON
    }

    fn is_not_positive(&self) -> bool {
        self < &EPSILON
    }
}

impl Infinite for f64 {
    fn is_infinite(&self) -> bool {
        Float::is_infinite(*self)
    }
}

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

macro_rules! ttype_signed {
    ($t:ident) => {
        impl Zero for $t {
            fn zero() -> Self {
                0
            }

            fn is_zero(&self) -> bool {
                num::Zero::is_zero(self)
            }
        }

        impl One for $t {
            fn one() -> Self {
                1
            }

            fn is_one(&self) -> bool {
                num::One::is_one(self)
            }
        }

        impl Signed for $t {
            fn abs(&self) -> Self {
                NumSigned::abs(&self)
            }

            fn is_positive(&self) -> bool {
                NumSigned::is_positive(self)
            }

            fn is_negative(&self) -> bool {
                NumSigned::is_negative(self)
            }

            fn is_not_negative(&self) -> bool {
                !Signed::is_negative(self)
            }

            fn is_not_positive(&self) -> bool {
                !Signed::is_positive(self)
            }
        }

        impl Infinite for $t {
            fn is_infinite(&self) -> bool {
                false
            }
        }
    };
}

macro_rules! ttype {
    ($t:ident) => {
        impl Zero for $t {
            fn zero() -> Self {
                0
            }

            fn is_zero(&self) -> bool {
                num::Zero::is_zero(self)
            }
        }

        impl One for $t {
            fn one() -> Self {
                1
            }

            fn is_one(&self) -> bool {
                num::One::is_one(self)
            }
        }

        impl Signed for $t {
            fn abs(&self) -> Self {
                *self
            }

            fn is_positive(&self) -> bool {
                *self > 0
            }

            fn is_negative(&self) -> bool {
                false
            }

            fn is_not_negative(&self) -> bool {
                !self.is_negative()
            }

            fn is_not_positive(&self) -> bool {
                !self.is_positive()
            }
        }
    };
}

ttype!(usize);
ttype!(u128);
ttype!(u16);
ttype!(u32);
ttype!(u64);
ttype!(u8);
ttype_signed!(i128);
ttype_signed!(i16);
ttype_signed!(i32);
ttype_signed!(i64);
ttype_signed!(i8);
