use num::BigUint;

use crate::{
    fraction_raw::fraction_raw::{FractionRaw, FractionRawMut, FractionRawRef},
    matrix::loose_fraction::Type,
};

pub trait FractionRawGetter<T: Clone> {
    fn get_clone(
        index: usize,
        types: &[Type],
        numerators: &[T],
        denominators: &[T],
    ) -> FractionRaw<T>;

    fn get_ref<'a>(
        index: usize,
        types: &'a [Type],
        numerators: &'a [T],
        denominators: &'a [T],
    ) -> FractionRawRef<'a, T>;

    fn get_mut<'a>(
        index: usize,
        types: &'a mut [Type],
        numerators: &'a mut [T],
        denominators: &'a mut [T],
    ) -> FractionRawMut<'a, T>;
}

impl FractionRawGetter<BigUint> for BigUint {
    fn get_clone(
        index: usize,
        types: &[Type],
        numerators: &[BigUint],
        denominators: &[BigUint],
    ) -> FractionRaw<BigUint> {
        FractionRaw(
            types[index],
            numerators[index].clone(),
            denominators[index].clone(),
        )
    }

    fn get_ref<'a>(
        index: usize,
        types: &'a [Type],
        numerators: &'a [BigUint],
        denominators: &'a [BigUint],
    ) -> FractionRawRef<'a, BigUint> {
        FractionRawRef(types[index], &numerators[index], &denominators[index])
    }

    fn get_mut<'a>(
        index: usize,
        types: &'a mut [Type],
        numerators: &'a mut [BigUint],
        denominators: &'a mut [BigUint],
    ) -> FractionRawMut<'a, BigUint> {
        FractionRawMut(
            &mut types[index],
            &mut numerators[index],
            &mut denominators[index],
        )
    }
}

impl FractionRawGetter<u64> for u64 {
    fn get_clone(
        index: usize,
        types: &[Type],
        numerators: &[u64],
        denominators: &[u64],
    ) -> FractionRaw<u64> {
        FractionRaw(
            types[index],
            numerators[index].clone(),
            denominators[index].clone(),
        )
    }

    fn get_ref<'a>(
        index: usize,
        types: &'a [Type],
        numerators: &'a [u64],
        denominators: &'a [u64],
    ) -> FractionRawRef<'a, u64> {
        FractionRawRef(types[index], &numerators[index], &denominators[index])
    }

    fn get_mut<'a>(
        index: usize,
        types: &'a mut [Type],
        numerators: &'a mut [u64],
        denominators: &'a mut [u64],
    ) -> FractionRawMut<'a, u64> {
        FractionRawMut(
            &mut types[index],
            &mut numerators[index],
            &mut denominators[index],
        )
    }
}
