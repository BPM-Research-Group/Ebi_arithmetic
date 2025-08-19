use std::ops::Neg;

use num::BigUint;

use crate::fraction_raw::fraction_raw::FractionRaw;

macro_rules! neg {
    ($t:ident) => {
        impl Neg for FractionRaw<$t> {
            type Output = Self;

            fn neg(mut self) -> Self::Output {
                self.0 = -self.0;
                self
            }
        }

        impl Neg for &FractionRaw<$t> {
            type Output = FractionRaw<$t>;

            fn neg(self) -> Self::Output {
                let mut x = self.clone();
                x.0 = -x.0;
                x
            }
        }
    };
}

neg!(BigUint);
neg!(u64);
