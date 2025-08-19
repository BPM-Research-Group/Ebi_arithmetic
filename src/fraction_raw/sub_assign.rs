use std::ops::SubAssign;

use num::BigUint;

use crate::{
    fraction_raw::fraction_raw::{FractionRaw, FractionRawMut, FractionRawRef},
    matrix::loose_fraction::Type,
};

macro_rules! sub_assign {
    ($type1:expr, $num1:expr, $den1:expr, $type2:expr, $num2:expr, $den2:expr) => {
        if $type1.is_plusminus() && $type2.is_plusminus() {
            //first, make the fractions comparable
            //numerator
            *$num1 *= $den2;
            let mut num2 = $num2 * &*$den1;

            //denominator
            *$den1 *= $den2;

            //second, perform the subtraction
            match (*$type1, $type2) {
                (Type::Plus, Type::Plus) => {
                    if &*$num1 >= &num2 {
                        *$num1 -= num2;
                    } else {
                        std::mem::swap($num1, &mut num2);
                        *$num1 -= num2;
                        *$type1 = Type::Minus;
                    }
                }
                (Type::Plus, Type::Minus) => *$num1 += num2,
                (Type::Minus, Type::Plus) => *$num1 += num2,
                (Type::Minus, Type::Minus) => {
                    if &num2 >= &*$num1 {
                        num2 -= &*$num1;
                        std::mem::swap($num1, &mut num2);
                        *$type1 = Type::Plus;
                    } else {
                        *$num1 -= num2;
                    }
                }
                _ => unreachable!(),
            }
        } else {
            *$type1 = match (*$type1, $type2) {
                (Type::Plus, Type::Infinite) => Type::NegInfinite,
                (Type::Plus, Type::NegInfinite) => Type::Infinite,
                (Type::Minus, Type::Infinite) => Type::NegInfinite,
                (Type::Minus, Type::NegInfinite) => Type::Infinite,
                (Type::Infinite, Type::Plus) => Type::Infinite,
                (Type::Infinite, Type::Minus) => Type::Infinite,
                (Type::Infinite, Type::Infinite) => Type::NaN,
                (Type::Infinite, Type::NegInfinite) => Type::Infinite,
                (Type::NegInfinite, Type::Plus) => Type::NegInfinite,
                (Type::NegInfinite, Type::Minus) => Type::NegInfinite,
                (Type::NegInfinite, Type::Infinite) => Type::NegInfinite,
                (Type::NegInfinite, Type::NegInfinite) => Type::NaN,
                (Type::NaN, _) => Type::NaN,
                (_, Type::NaN) => Type::NaN,
                _ => unreachable!(),
            }
        }
    };
}

macro_rules! sub_assign_2 {
    ($t:ident) => {
        impl SubAssign<FractionRaw<$t>> for FractionRaw<$t> {
            fn sub_assign(&mut self, rhs: FractionRaw<$t>) {
                let FractionRaw(type1, num1, den1) = self;
                let FractionRaw(type2, num2, den2) = rhs;

                sub_assign!(type1, num1, den1, type2, num2, &den2)
            }
        }

        impl SubAssign<&FractionRaw<$t>> for FractionRaw<$t> {
            fn sub_assign(&mut self, rhs: &FractionRaw<$t>) {
                let FractionRaw(type1, num1, den1) = self;
                let FractionRaw(type2, num2, den2) = rhs;

                sub_assign!(type1, num1, den1, type2, num2, den2)
            }
        }

        impl<'a> SubAssign<FractionRawMut<'a, $t>> for FractionRaw<$t> {
            fn sub_assign(&mut self, rhs: FractionRawMut<'a, $t>) {
                let FractionRaw(type1, num1, den1) = self;
                let FractionRawMut(type2, num2, den2) = rhs;

                sub_assign!(type1, num1, den1, type2, &*num2, &*den2)
            }
        }

        impl<'a> SubAssign<FractionRawRef<'a, $t>> for FractionRaw<$t> {
            fn sub_assign(&mut self, rhs: FractionRawRef<'a, $t>) {
                let FractionRaw(type1, num1, den1) = self;
                let FractionRawRef(type2, num2, den2) = rhs;

                sub_assign!(type1, num1, den1, type2, num2, den2)
            }
        }

        impl<'a> SubAssign<FractionRaw<$t>> for FractionRawMut<'a, $t> {
            fn sub_assign(&mut self, rhs: FractionRaw<$t>) {
                let FractionRawMut(type1, num1, den1) = self;
                let FractionRaw(type2, num2, den2) = rhs;

                sub_assign!(*type1, *num1, *den1, type2, num2, &den2)
            }
        }

        impl<'a> SubAssign<&FractionRaw<$t>> for FractionRawMut<'a, $t> {
            fn sub_assign(&mut self, rhs: &FractionRaw<$t>) {
                let FractionRawMut(type1, num1, den1) = self;
                let FractionRaw(type2, num2, den2) = rhs;

                sub_assign!(*type1, *num1, *den1, type2, num2, den2)
            }
        }

        impl<'a> SubAssign<FractionRawMut<'a, $t>> for FractionRawMut<'a, $t> {
            fn sub_assign(&mut self, rhs: FractionRawMut<'a, $t>) {
                let FractionRawMut(type1, num1, den1) = self;
                let FractionRawMut(type2, num2, den2) = rhs;

                sub_assign!(*type1, *num1, *den1, type2, &*num2, &*den2)
            }
        }

        impl<'a> SubAssign<FractionRawRef<'a, $t>> for FractionRawMut<'a, $t> {
            fn sub_assign(&mut self, rhs: FractionRawRef<'a, $t>) {
                let FractionRawMut(type1, num1, den1) = self;
                let FractionRawRef(type2, num2, den2) = rhs;

                sub_assign!(*type1, *num1, *den1, type2, num2, den2)
            }
        }
    };
}

sub_assign_2!(BigUint);
sub_assign_2!(u64);
