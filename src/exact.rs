use anyhow::{Result, anyhow};
use num::{BigInt, BigUint};
use std::sync::atomic::AtomicBool;

//======================== exactness tools ========================//

static EXACT: AtomicBool = AtomicBool::new(true);

/// Enables or disables exact arithmetic globally.
/// Exact arithmetic cannot be combined with approximate arithmetic.
pub fn set_exact_globally(exact: bool) {
    EXACT.store(exact, std::sync::atomic::Ordering::Relaxed);
}

pub fn is_exact_globally() -> bool {
    if cfg!(any(
        all(
            feature = "exactarithmetic",
            feature = "approximatearithmetic"
        ),
        all(
            not(feature = "exactarithmetic"),
            not(feature = "approximatearithmetic")
        )
    )) {
        EXACT.load(std::sync::atomic::Ordering::Relaxed)
    } else if cfg!(feature = "exactarithmetic") {
        true
    } else {
        false
    }
}

//======================== exactness trait ========================//
pub trait MaybeExact {
    type Approximate;
    type Exact;

    fn is_exact(&self) -> bool;

    /**
     * This is a low-level function to extract an approximate value. Will only succeed if the fraction is approximate.
     */
    fn extract_approx(&self) -> Result<Self::Approximate>;

    /**
     * This is a low-level function to extract an exact value. Will only succeed if the fraction is exact.
     */
    fn extract_exact(&self) -> Result<&Self::Exact>;
}

macro_rules! approx {
    ($t: ident) => {
        impl MaybeExact for $t {
            type Approximate = $t;
            type Exact = ();

            fn is_exact(&self) -> bool {
                false
            }

            fn extract_approx(&self) -> Result<Self::Approximate> {
                Ok(*self)
            }

            fn extract_exact(&self) -> Result<&Self::Exact> {
                Err(anyhow!("cannot extract exact value"))
            }
        }
    };
}

macro_rules! exact {
    ($t: ident) => {
        impl MaybeExact for $t {
            type Approximate = ();
            type Exact = $t;

            fn is_exact(&self) -> bool {
                true
            }

            fn extract_approx(&self) -> Result<Self::Approximate> {
                Err(anyhow!("cannot extract approximate value"))
            }

            fn extract_exact(&self) -> Result<&Self::Exact> {
                Ok(self)
            }
        }
    };
}

exact!(i128);
exact!(i64);
exact!(i32);
exact!(i16);
exact!(u128);
exact!(u64);
exact!(u32);
exact!(u16);
exact!(u8);
exact!(BigUint);
exact!(BigInt);
approx!(f64);
approx!(f32);
