use num::BigUint;

use crate::{
    ebi_number::One,
    fraction_raw::fraction_raw::{FractionRaw, FractionRawMut},
    matrix::loose_fraction::Type,
};

pub trait SetOne {
    fn set_one(&mut self);
}

macro_rules! one {
    ($typee:expr, $num:expr, $den:expr) => {
        match $typee {
            Type::Plus => $num == $den,
            Type::Minus => false,
            Type::NaN => false,
            Type::Infinite => false,
            Type::NegInfinite => false,
        }
    };
}

macro_rules! one_2 {
    ($t:ident) => {
        impl One for FractionRaw<$t> {
            fn one() -> Self {
                FractionRaw(Type::Plus, $t::one(), $t::one())
            }

            fn is_one(&self) -> bool {
                let FractionRaw(typee, num, den) = self;
                one!(typee, num, den)
            }
        }

        impl SetOne for FractionRaw<$t> {
            fn set_one(&mut self) {
                let FractionRaw(typee, num, den) = self;
                *typee = Type::Plus;
                *num = $t::one();
                *den = $t::one();
            }
        }

        impl<'a> SetOne for FractionRawMut<'a, $t> {
            fn set_one(&mut self) {
                let FractionRawMut(typee, num, den) = self;
                **typee = Type::Plus;
                **num = $t::one();
                **den = $t::one();
            }
        }
    };
}

one_2!(BigUint);
one_2!(u64);
