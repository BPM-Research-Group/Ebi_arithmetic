use std::fmt::Display;

use crate::{fraction_exact::FractionExact, matrix::loose_fraction::Type};
use fraction::ToPrimitive;

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

macro_rules! display {
    ($self:expr, $f:expr) => {
        match $self.0 {
            Type::Plus => write!($f, "{}/{}", $self.1, $self.2),
            Type::Minus => write!($f, "-{}/{}", $self.1, $self.2),
            Type::NaN => write!($f, "NaN"),
            Type::Infinite => write!($f, "∞"),
            Type::NegInfinite => write!($f, "-∞"),
        }
    };
}

impl<T> Display for FractionRaw<T>
where
    T: Display + Clone,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        display!(self, f)
    }
}

impl<'a, T> Display for FractionRawRef<'a, T>
where
    T: Display + Clone,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        display!(self, f)
    }
}

impl<'a, T> Display for FractionRawMut<'a, T>
where
    T: Display + Clone,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        display!(self, f)
    }
}
