//======================== set type alias based on compile flags ========================//
/// The fraction matrix is a matrix of fractions.
/// It applies two strategies to potentially save time for exact arithmetic matrices:
/// - it postpones reduction of fractions to the moment of export, or a user-chosen moment.
/// - it attempts to store values in primitives rather than BigUints at each reduction.
#[cfg(any(
    all(
        not(feature = "exactarithmetic"),
        not(feature = "approximatearithmetic")
    ),
    all(feature = "exactarithmetic", feature = "approximatearithmetic")
))]
pub type LogPolynomial = super::log_polynomial_enum::LogPolynomialEnum;

#[cfg(all(not(feature = "exactarithmetic"), feature = "approximatearithmetic"))]
pub type LogPolynomial = super::log_polynomial_f64::LogPolynomialF64;

#[cfg(all(feature = "exactarithmetic", not(feature = "approximatearithmetic")))]
pub type LogPolynomial = super::log_polynomial_exact::LogPolynomialExact;
