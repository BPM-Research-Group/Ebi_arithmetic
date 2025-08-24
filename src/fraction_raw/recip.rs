use num::BigUint;

use crate::{
    fraction::ebi_number::Zero,
    fraction_raw::fraction_raw::{FractionRaw, FractionRawMut},
    matrix::loose_fraction::Type,
};

pub trait Recip {
    fn recip(&mut self);
}

macro_rules! recip {
    ($typee:expr, $num: expr, $den: expr) => {
        match $typee {
            Type::Plus | Type::Minus => {
                if $num.is_zero() {
                    *$typee = Type::NaN
                } else {
                    std::mem::swap($num, $den)
                }
            }
            Type::NaN => {}
            Type::Infinite => {
                *$typee = Type::NaN;
            }
            Type::NegInfinite => {
                *$typee = Type::NaN;
            }
        };
    };
}

macro_rules! recip_2 {
    ($t:ident) => {
        impl Recip for FractionRaw<$t> {
            fn recip(&mut self) {
                let FractionRaw(typee, num, den) = self;

                recip!(typee, num, den);
            }
        }

        impl<'a> Recip for FractionRawMut<'a, $t> {
            fn recip(&mut self) {
                let FractionRawMut(typee, num, den) = self;

                recip!(*typee, *num, *den);
            }
        }
    };
}

recip_2!(BigUint);
recip_2!(u64);
