// ============ implementations ============
use num::Signed as NumSigned;

use crate::ebi_number::{EbiNumber, Infinite, Normal, One, Signed, Zero};

macro_rules! ttype_signed {
    ($t:ident) => {
        impl EbiNumber for $t {}

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
        }

        impl Infinite for $t {
            fn is_infinite(&self) -> bool {
                false
            }
        }

        impl Normal for $t {
            fn is_nan(&self) -> bool {
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
        }

        impl Infinite for $t {
            fn is_infinite(&self) -> bool {
                false
            }
        }

        impl Normal for $t {
            fn is_nan(&self) -> bool {
                false
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
