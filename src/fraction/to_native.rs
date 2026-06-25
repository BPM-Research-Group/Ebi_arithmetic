use crate::{
    Signed, ToNative,
    fraction::{
        fraction_enum::FractionEnum, fraction_exact::FractionExact, fraction_f64::FractionF64,
    },
};
use malachite::{
    Rational,
    base::num::{
        arithmetic::traits::{Ceiling, Floor},
        basic::traits::OneHalf,
    },
};
use std::usize;

impl ToNative for Rational {
    fn to_usize(&self) -> usize {
        let int = if self.is_positive() {
            Floor::floor(self + Rational::ONE_HALF)
        } else {
            Ceiling::ceiling(self - Rational::ONE_HALF)
        };

        match usize::try_from(&int) {
            Ok(x) => x,
            Err(_) => {
                if int.is_negative() {
                    0
                } else {
                    usize::MAX
                }
            }
        }
    }
}

impl ToNative for FractionF64 {
    fn to_usize(&self) -> usize {
        self.0.to_usize()
    }
}

impl ToNative for FractionExact {
    fn to_usize(&self) -> usize {
        self.0.to_usize()
    }
}

impl ToNative for FractionEnum {
    fn to_usize(&self) -> usize {
        match self {
            FractionEnum::Exact(rational) => rational.to_usize(),
            FractionEnum::Approx(f) => f.to_usize(),
            FractionEnum::CannotCombineExactAndApprox => usize::MAX,
        }
    }
}

impl ToNative for usize {
    fn to_usize(&self) -> usize {
        *self
    }
}

macro_rules! floats {
    ($t:ty) => {
        impl ToNative for $t {
            fn to_usize(&self) -> usize {
                if self.is_nan() {
                    return usize::MAX;
                }
                let int = self.round();
                if Signed::is_negative(&int) {
                    return 0;
                }
                if int.is_infinite() {
                    return usize::MAX;
                }
                int as usize
            }
        }
    };
}

floats!(f64);
floats!(f32);

macro_rules! ttype {
    ($t:ty) => {
        impl ToNative for $t {
            fn to_usize(&self) -> usize {
                #[allow(unused_comparisons)]
                if *self < 0 {
                    return 0;
                }
                match usize::try_from(*self) {
                    Ok(x) => x,
                    Err(_) => usize::MAX,
                }
            }
        }
    };
}

ttype!(u8);
ttype!(u16);
ttype!(u32);
ttype!(u64);
ttype!(u128);
ttype!(i8);
ttype!(i16);
ttype!(i32);
ttype!(i64);
ttype!(i128);
