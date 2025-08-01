use crate::{
    ebi_number::{One, Zero},
    fraction_enum::FractionEnum,
    fraction_exact::FractionExact,
    fraction_f64::FractionF64,
};

macro_rules! fraction {
    ($t: ident) => {
        impl fraction::Zero for $t {
            fn zero() -> Self {
                Zero::zero()
            }

            fn is_zero(&self) -> bool {
                Zero::is_zero(self)
            }
        }

        impl fraction::One for $t {
            fn one() -> Self {
                One::one()
            }

            fn is_one(&self) -> bool {
                One::is_one(self)
            }
        }
    };
}
fraction!(FractionEnum);
fraction!(FractionExact);
fraction!(FractionF64);
