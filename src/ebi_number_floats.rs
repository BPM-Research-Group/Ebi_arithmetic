use num::Float;

use crate::{
    ebi_number::{EbiNumber, Fractional, Infinite, Normal, One, Round, Signed, Zero},
    fraction::EPSILON,
};

macro_rules! float {
    ($t: ident, $e: expr) => {
        impl EbiNumber for $t {}

        impl One for $t {
            fn one() -> Self {
                1.0
            }

            fn is_one(&self) -> bool {
                (self - 1.0).abs() - &$e < 0.0
            }
        }

        impl Zero for $t {
            fn zero() -> Self {
                0.0
            }

            fn is_zero(&self) -> bool {
                <$t as Signed>::abs(&self) - &$e < 0.0
            }
        }

        impl Signed for $t {
            fn abs(&self) -> Self {
                <$t as Float>::abs(*self)
            }

            fn is_positive(&self) -> bool {
                *self != 0.0 && self > &$e
            }

            fn is_negative(&self) -> bool {
                *self != 0.0 && self < &-$e
            }

            fn is_not_negative(&self) -> bool {
                self > &-$e
            }

            fn is_not_positive(&self) -> bool {
                self < &$e
            }
        }

        impl Infinite for $t {
            fn is_infinite(&self) -> bool {
                Float::is_infinite(*self)
            }
        }

        impl Normal for $t {
            fn is_nan(&self) -> bool {
                Float::is_nan(*self)
            }
        }

        impl Round for $t {
            fn floor(self) -> $t {
                $t::floor(self)
            }

            fn ceil(self) -> $t {
                $t::ceil(self)
            }
        }

        impl Fractional for $t {
            fn recip(&self) -> Self {
                $t::recip(*self)
            }
        }
    };
}

float!(f32, f32::EPSILON);
float!(f64, EPSILON);
