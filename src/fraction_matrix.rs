//======================== set type alias based on compile flags ========================//

#[cfg(any(
    all(
        not(feature = "exactarithmetic"),
        not(feature = "approximatearithmetic")
    ),
    all(feature = "exactarithmetic", feature = "approximatearithmetic")
))]
pub type FractionMatrix = super::fraction_matrix_enum::FractionMatrixEnum;

#[cfg(all(not(feature = "exactarithmetic"), feature = "approximatearithmetic"))]
pub type FractionMatrix = super::fraction_matrix_f64::FractionMatrixF64;

#[cfg(all(feature = "exactarithmetic", not(feature = "approximatearithmetic")))]
pub type FractionMatrix = super::fraction_matrix_exact::FractionMatrixExact;

//======================== tests ========================//
#[cfg(test)]
mod tests {
    use crate::{ebi_matrix::EbiMatrix, f, fraction::Fraction, fraction_matrix::FractionMatrix};

    #[test]
    fn fraction_matrix() {
        let m: FractionMatrix = vec![vec![f!(1, 4), f!(2, 5), f!(8, 3)]].into();

        let _ = m.reduce();
    }

    #[test]
    fn fraction_matrix_abnormal() {
        let m = vec![vec![
            Fraction::infinity(),
            Fraction::neg_infinity(),
            f!(8, 3),
        ]];
        let m: FractionMatrix = m.into();
        let _ = m.reduce();
    }

    #[test]
    fn fraction_matrix_empty() {
        let m = vec![vec![]];
        let m: FractionMatrix = m.into();
        assert_eq!(m.number_of_rows(), 1);
        assert_eq!(m.number_of_columns(), 0);
        let _ = m.reduce();
    }

    #[test]
    fn fraction_matrix_get() {
        let m: FractionMatrix = vec![vec![f!(1, 4), f!(2, 5), f!(8, 3)]].into();

        assert_eq!(m.get(0, 0), f!(1, 4));
    }
}
