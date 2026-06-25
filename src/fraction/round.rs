use malachite::{
    base::num::{
        arithmetic::traits::{Ceiling, Floor},
        basic::traits::OneHalf,
    },
    rational::Rational,
};

use crate::{
    Signed,
    ebi_number::Round,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
};

impl Round for FractionF64 {
    fn floor(self) -> Self {
        FractionF64(self.0.floor())
    }

    fn ceil(self) -> Self {
        FractionF64(self.0.ceil())
    }

    fn round_half_away_from_zero(self) -> Self {
        FractionF64(self.0.round())
    }
}

impl Round for FractionExact {
    fn floor(self) -> Self {
        Self(Round::floor(self.0))
    }

    fn ceil(self) -> Self {
        Self(Round::ceil(self.0))
    }

    fn round_half_away_from_zero(self) -> Self {
        Self(Round::round_half_away_from_zero(self.0))
    }
}

impl Round for FractionEnum {
    fn floor(self) -> Self {
        match self {
            Self::Exact(f) => Self::Exact(Round::floor(f)),
            Self::Approx(f) => Self::Approx(f.floor()),
            Self::CannotCombineExactAndApprox => Self::CannotCombineExactAndApprox,
        }
    }

    fn ceil(self) -> Self {
        match self {
            Self::Exact(f) => Self::Exact(f.ceil()),
            Self::Approx(f) => Self::Approx(f.ceil()),
            Self::CannotCombineExactAndApprox => Self::CannotCombineExactAndApprox,
        }
    }

    fn round_half_away_from_zero(self) -> Self {
        match self {
            Self::Exact(f) => Self::Exact(f.round_half_away_from_zero()),
            Self::Approx(f) => Self::Approx(f.round()),
            Self::CannotCombineExactAndApprox => Self::CannotCombineExactAndApprox,
        }
    }
}

impl Round for Rational {
    fn floor(self) -> Self {
        Floor::floor(self).into()
    }

    fn ceil(self) -> Self {
        Ceiling::ceiling(self).into()
    }

    fn round_half_away_from_zero(self) -> Self {
        if self.is_positive() {
            Floor::floor(self + Rational::ONE_HALF).into()
        } else {
            Ceiling::ceiling(self - Rational::ONE_HALF).into()
        }
    }
}

macro_rules! float {
    ($t: ident, $e: expr) => {
        impl Round for $t {
            fn floor(self) -> $t {
                $t::floor(self)
            }

            fn ceil(self) -> $t {
                $t::ceil(self)
            }

            fn round_half_away_from_zero(self) -> $t {
                $t::round(self)
            }
        }
    };
}

float!(f32, f32::EPSILON);
float!(f64, EPSILON);

macro_rules! ttype {
    ($t:ident) => {
        impl Round for $t {
            fn floor(self) -> Self {
                self
            }

            fn ceil(self) -> Self {
                self
            }

            fn round_half_away_from_zero(self) -> Self {
                self
            }
        }
    };
}

ttype!(usize);
ttype!(u128);
ttype!(u64);
ttype!(u32);
ttype!(u16);
ttype!(u8);
ttype!(i128);
ttype!(i64);
ttype!(i32);
ttype!(i16);
ttype!(i8);

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use crate::{One, Round, Zero, fraction::fraction_exact::FractionExact};

    #[test]
    fn rounding() {
        let zero = FractionExact::zero();
        assert_eq!(zero, zero.clone().round_half_away_from_zero());

        assert_eq!(
            FractionExact::from_str("1/2")
                .unwrap()
                .round_half_away_from_zero(),
            FractionExact::one()
        );

        assert_eq!(
            FractionExact::from_str("-1/2")
                .unwrap()
                .round_half_away_from_zero(),
            -FractionExact::one()
        );

        assert_eq!(
            FractionExact::from_str("-2/3")
                .unwrap()
                .round_half_away_from_zero(),
            -FractionExact::one()
        );
    }
}
