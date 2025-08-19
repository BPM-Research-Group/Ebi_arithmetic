use std::ops::{Neg, SubAssign};

use crate::{
    ebi_number::{One, Zero},
    fraction_exact::FractionExact,
    matrix::loose_fraction::Type,
};
use fraction::ToPrimitive;
use num::BigUint;

#[derive(Clone)]
pub struct FractionRaw<T>(pub Type, pub T, pub T)
where
    T: Clone;
pub struct FractionRawRef<'a, T>(pub Type, pub &'a T, pub &'a T);
pub struct FractionRawMut<'a, T>(pub &'a mut Type, pub &'a mut T, pub &'a mut T);

impl FractionRaw<u64> {
    pub fn try_u64(value: &FractionExact) -> Option<Self> {
        let typee = (&value.0).into();
        let num = match value.0.numer() {
            Some(num) => match num.to_u64() {
                Some(x) => x,
                None => {
                    return None;
                }
            },
            None => 0,
        };
        let den = match value.0.denom() {
            Some(den) => match den.to_u64() {
                Some(x) => x,
                None => {
                    return None;
                }
            },
            None => 0,
        };

        Some(Self(typee, num, den))
    }
}

impl<'a> FractionRaw<BigUint> {
    pub fn recip(&mut self) {
        let FractionRaw(typee, num, den) = self;
        match typee {
            Type::Plus | Type::Minus => {
                if num.is_zero() {
                    *typee = Type::NaN
                } else {
                    std::mem::swap(num, den)
                }
            }
            Type::NaN => {}
            Type::Infinite => {
                *typee = Type::NaN;
            }
            Type::NegInfinite => {
                *typee = Type::NaN;
            }
        };
    }
}

impl Zero for FractionRaw<BigUint> {
    fn zero() -> Self {
        FractionRaw(Type::Plus, BigUint::zero(), BigUint::one())
    }

    fn is_zero(&self) -> bool {
        let FractionRaw(typee, num, _) = self;
        match typee {
            Type::Plus | Type::Minus => num.is_zero(),
            Type::NaN => false,
            Type::Infinite => false,
            Type::NegInfinite => false,
        }
    }
}

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

impl SubAssign<&FractionRaw<BigUint>> for FractionRaw<BigUint> {
    fn sub_assign(&mut self, rhs: &FractionRaw<BigUint>) {
        let FractionRaw(type1, num1, den1) = self;
        let FractionRaw(type2, num2, den2) = rhs;

        sub_assign!(type1, num1, den1, type2, num2, den2)
    }
}

impl<'a> SubAssign<FractionRawRef<'a, BigUint>> for FractionRaw<BigUint> {
    fn sub_assign(&mut self, rhs: FractionRawRef<'a, BigUint>) {
        let FractionRaw(type1, num1, den1) = self;
        let FractionRawRef(type2, num2, den2) = rhs;

        sub_assign!(type1, num1, den1, type2, num2, den2)
    }
}

impl Neg for FractionRaw<BigUint> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.0 = -self.0;
        self
    }
}

impl Neg for &FractionRaw<BigUint> {
    type Output = FractionRaw<BigUint>;

    fn neg(self) -> Self::Output {
        let mut x = self.clone();
        x.0 = -x.0;
        x
    }
}
