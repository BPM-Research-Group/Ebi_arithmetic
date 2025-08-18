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
pub type FractionMatrix = super::fraction_matrix_enum::FractionMatrixEnum;

#[cfg(all(not(feature = "exactarithmetic"), feature = "approximatearithmetic"))]
pub type FractionMatrix = super::fraction_matrix_f64::FractionMatrixF64;

#[cfg(all(feature = "exactarithmetic", not(feature = "approximatearithmetic")))]
pub type FractionMatrix = super::fraction_matrix_exact::FractionMatrixExact;

//======================== tests ========================//
#[cfg(test)]
mod tests {
    use crate::{
        ebi_number::Zero, f, f0, fraction::Fraction, matrix::ebi_matrix::EbiMatrix,
        matrix::fraction_matrix::FractionMatrix,
    };

    #[test]
    fn fraction_matrix() {
        let m: FractionMatrix = vec![vec![f!(1, 4), f!(2, 5), f!(8, 3)]].try_into().unwrap();

        let _ = m.reduce();
    }

    #[test]
    fn fraction_matrix_abnormal() {
        let m = vec![vec![
            Fraction::infinity(),
            Fraction::neg_infinity(),
            f!(8, 3),
        ]];
        let m: FractionMatrix = m.try_into().unwrap();
        let _ = m.reduce();
    }

    #[test]
    fn fraction_matrix_empty() {
        let m = vec![vec![]];
        let m: FractionMatrix = m.try_into().unwrap();
        assert_eq!(m.number_of_rows(), 1);
        assert_eq!(m.number_of_columns(), 0);
        let _ = m.reduce();
    }

    #[test]
    fn fraction_matrix_get() {
        let m: FractionMatrix = vec![vec![f!(1, 4), f!(2, 5), f!(8, 3)]].try_into().unwrap();

        assert!(m.get(10, 10).is_none());

        assert_eq!(m.get(0, 0).unwrap(), f!(1, 4));
    }

    #[test]
    #[should_panic]
    fn fraction_matrix_incomplete() {
        let _: FractionMatrix = vec![vec![f!(1, 4), f!(2, 5)], vec![f!(8, 3)]]
            .try_into()
            .unwrap();
    }

    #[test]
    fn fraction_matrix_pop_front() {
        let mut m1: FractionMatrix = vec![vec![f!(1, 4), f!(2, 5), f!(8, 3)]].try_into().unwrap();

        m1.pop_front_columns(1);

        let m3: FractionMatrix = vec![vec![f!(2, 5), f!(8, 3)]].try_into().unwrap();

        println!("{:?}", m1);
        println!("{:?}", m3);

        assert!(m1.inner_eq(&m3));
    }

    #[test]
    fn fraction_matrix_push_columns() {
        let mut m1: FractionMatrix = vec![vec![f!(1, 4), f!(2, 5), f!(8, 3)]].try_into().unwrap();

        m1.push_columns(1);

        let m3: FractionMatrix = vec![vec![f!(1, 4), f!(2, 5), f!(8, 3), f0!()]]
            .try_into()
            .unwrap();

        println!("{:?}", m1);
        println!("{:?}", m3);

        assert!(m1.inner_eq(&m3));
    }
}
