use num::BigUint;

use crate::{
    ebi_number::{One, Zero},
    fraction_raw::fraction_raw::{FractionRaw, FractionRawMut, FractionRawRef},
    matrix::loose_fraction::Type,
};

pub trait SetZero {
    fn set_zero(&mut self);
}

pub trait IsZero {
    fn is_zero(&self) -> bool;
}

macro_rules! zero {
    ($typee:expr, $num:expr) => {
        match $typee {
            Type::Plus | Type::Minus => $num.is_zero(),
            Type::NaN => false,
            Type::Infinite => false,
            Type::NegInfinite => false,
        }
    };
}

macro_rules! one_2 {
    ($t:ident) => {
        impl IsZero for FractionRaw<$t> {
            fn is_zero(&self) -> bool {
                let FractionRaw(typee, num, _) = self;
                zero!(typee, num)
            }
        }

        impl SetZero for FractionRaw<$t> {
            fn set_zero(&mut self) {
                let FractionRaw(typee, num, den) = self;
                *typee = Type::Plus;
                *num = $t::zero();
                *den = $t::one();
            }
        }

        impl<'a> IsZero for FractionRawRef<'a, $t> {
            fn is_zero(&self) -> bool {
                let FractionRawRef(typee, num, _) = self;
                zero!(typee, num)
            }
        }

        impl<'a> SetZero for FractionRawMut<'a, $t> {
            fn set_zero(&mut self) {
                let FractionRawMut(typee, num, den) = self;
                **typee = Type::Plus;
                **num = $t::zero();
                **den = $t::one();
            }
        }
    };
}

one_2!(BigUint);
one_2!(u64);
