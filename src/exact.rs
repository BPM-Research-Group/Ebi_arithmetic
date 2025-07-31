use anyhow::{anyhow,Result};
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

macro_rules! prim {
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

prim!(i128);
prim!(i64);
prim!(i32);
prim!(i16);
prim!(u128);
prim!(u64);
prim!(u32);
prim!(u16);
prim!(u8);