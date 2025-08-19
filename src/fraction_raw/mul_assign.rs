use std::ops::MulAssign;

use num::BigUint;

use crate::fraction_raw::fraction_raw::{FractionRaw, FractionRawMut, FractionRawRef};

macro_rules! mul_assign {
    ($type1:expr, $num1:expr, $den1:expr, $type2:expr, $num2:expr, $den2:expr) => {
        *$type1 = *$type1 * $type2;

        if $type1.is_plusminus() {
            *$num1 *= $num2;
            *$den1 *= $den2;
        }
    };
}

macro_rules! mul_assign_2 {
    ($t:ident) => {
        //raw

        impl<'a> MulAssign<FractionRaw<$t>> for FractionRaw<$t> {
            fn mul_assign(&mut self, rhs: FractionRaw<$t>) {
                let FractionRaw(type1, num1, den1) = self;
                let FractionRaw(type2, num2, den2) = rhs;

                mul_assign!(type1, num1, den1, type2, num2, den2);
            }
        }

        impl<'a> MulAssign<&FractionRaw<$t>> for FractionRaw<$t> {
            fn mul_assign(&mut self, rhs: &FractionRaw<$t>) {
                let FractionRaw(type1, num1, den1) = self;
                let FractionRaw(type2, num2, den2) = rhs;

                mul_assign!(type1, num1, den1, *type2, num2, den2);
            }
        }

        impl<'a> MulAssign<FractionRawRef<'a, $t>> for FractionRaw<$t> {
            fn mul_assign(&mut self, rhs: FractionRawRef<'a, $t>) {
                let FractionRaw(type1, num1, den1) = self;
                let FractionRawRef(type2, num2, den2) = rhs;

                mul_assign!(type1, num1, den1, type2, num2, den2);
            }
        }

        impl<'a> MulAssign<FractionRawMut<'a, $t>> for FractionRaw<$t> {
            fn mul_assign(&mut self, rhs: FractionRawMut<'a, $t>) {
                let FractionRaw(type1, num1, den1) = self;
                let FractionRawMut(type2, num2, den2) = rhs;

                mul_assign!(type1, num1, den1, *type2, &*num2, &*den2);
            }
        }

        //mut

        impl<'a> MulAssign<FractionRaw<$t>> for FractionRawMut<'a, $t> {
            fn mul_assign(&mut self, rhs: FractionRaw<$t>) {
                let FractionRawMut(type1, num1, den1) = self;
                let FractionRaw(type2, num2, den2) = rhs;

                mul_assign!(*type1, *num1, *den1, type2, num2, den2);
            }
        }

        impl<'a> MulAssign<&FractionRaw<$t>> for FractionRawMut<'a, $t> {
            fn mul_assign(&mut self, rhs: &FractionRaw<$t>) {
                let FractionRawMut(type1, num1, den1) = self;
                let FractionRaw(type2, num2, den2) = rhs;

                mul_assign!(*type1, *num1, *den1, *type2, num2, den2);
            }
        }

        impl<'a> MulAssign<FractionRawRef<'a, $t>> for FractionRawMut<'a, $t> {
            fn mul_assign(&mut self, rhs: FractionRawRef<'a, $t>) {
                let FractionRawMut(type1, num1, den1) = self;
                let FractionRawRef(type2, num2, den2) = rhs;

                mul_assign!(*type1, *num1, *den1, type2, num2, den2);
            }
        }

        impl<'a> MulAssign<FractionRawMut<'a, $t>> for FractionRawMut<'a, $t> {
            fn mul_assign(&mut self, rhs: FractionRawMut<'a, $t>) {
                let FractionRawMut(type1, num1, den1) = self;
                let FractionRawMut(type2, num2, den2) = rhs;

                mul_assign!(*type1, *num1, *den1, *type2, &*num2, &*den2);
            }
        }
    };
}

mul_assign_2!(BigUint);
mul_assign_2!(u64);

#[cfg(test)]
mod tests {
    use crate::{fraction_raw::fraction_raw::FractionRaw, matrix::loose_fraction::Type};

    #[test]
    fn mul() {
        let mut a = FractionRaw(Type::Minus, 2u64, 2);
        let b = FractionRaw(Type::Minus, 3u64, 4);

        a *= b;

        let FractionRaw(type1, num1, den1) = a;

        assert_eq!(type1, Type::Plus);
        assert_eq!(num1, 6);
        assert_eq!(den1, 8);
    }
}
