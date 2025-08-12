use anyhow::{Context, Result, anyhow};
use fraction::{GenericFraction, Ratio, Sign};
use num::BigUint;
use num_bigint::RandBigInt;
use rand::Rng;

use crate::{
    exact::{MaybeExact, is_exact_globally},
    fraction_enum::FractionEnum,
    fraction_exact::FractionExact,
    fraction_f64::FractionF64,
    ebi_number::Zero,
};

pub trait ChooseRandomly {
    type Cache;
    /**
     * Return a random index from 0 (inclusive) to the length of the list (exclusive).
     * The likelihood of each index to be returned is proportional to the value of the fraction at that index.
     *
     * The fractions do not need to sum to 1.
     */
    fn choose_randomly(fractions: &Vec<Self>) -> Result<usize>
    where
        Self: Sized;

    fn choose_randomly_create_cache<'a>(
        fractions: impl Iterator<Item = &'a Self>,
    ) -> Result<Self::Cache>
    where
        Self: Sized,
        Self: 'a;

    fn choose_randomly_cached(cache: &Self::Cache) -> usize
    where
        Self: Sized;
}

#[cfg(any(
    all(
        not(feature = "exactarithmetic"),
        not(feature = "approximatearithmetic")
    ),
    all(feature = "exactarithmetic", feature = "approximatearithmetic")
))]
pub type FractionRandomCache = FractionRandomCacheEnum;

#[cfg(all(not(feature = "exactarithmetic"), feature = "approximatearithmetic"))]
pub type FractionRandomCache = FractionRandomCacheF64;

#[cfg(all(feature = "exactarithmetic", not(feature = "approximatearithmetic")))]
pub type FractionRandomCache = FractionRandomCacheExact;

pub enum FractionRandomCacheEnum {
    Exact(Vec<fraction::BigFraction>, BigUint),
    Approx(Vec<f64>),
}

impl ChooseRandomly for FractionEnum {
    type Cache = FractionRandomCacheEnum;

    fn choose_randomly(fractions: &Vec<FractionEnum>) -> Result<usize> {
        if fractions.is_empty() {
            return Err(anyhow!("cannot take an element of an empty list"));
        }

        //normalise the probabilities
        let mut probabilities: Vec<FractionEnum> = fractions.iter().cloned().collect();
        let sum = probabilities
            .iter()
            .fold(FractionEnum::zero(), |x, y| &x + y);
        if sum == FractionEnum::CannotCombineExactAndApprox {
            return Err(anyhow!("cannot combine exact and approximate arithmetic"));
        }
        probabilities.retain_mut(|v| {
            *v /= &sum;
            true
        });

        //select a random value
        let mut rng = rand::thread_rng();
        let rand_val = if sum.is_exact() {
            //strategy: the highest denominator determines how much precision we need
            let temp_zero = BigUint::zero();
            let max_denom = probabilities
                .iter()
                .map(|f| {
                    if let FractionEnum::Exact(e) = f {
                        e.denom().unwrap()
                    } else {
                        &temp_zero
                    }
                })
                .max()
                .unwrap();
            //Generate a random value with the number of bits of the highest denominator. Repeat until this value is <= the max denominator.
            let mut rand_val = rng.gen_biguint(max_denom.bits());
            while &rand_val > max_denom {
                rand_val = rng.gen_biguint(max_denom.bits());
            }
            //create the fraction from the random nominator and the max denominator
            FractionEnum::try_from((rand_val, max_denom.clone())).unwrap()
        } else {
            FractionEnum::Approx(rng.gen_range(0.0..=1.0))
        };

        let mut cum_prob = FractionEnum::zero();
        for (index, value) in probabilities.iter().enumerate() {
            cum_prob += value;
            if rand_val < cum_prob {
                return Ok(index);
            }
        }
        Ok(probabilities.len() - 1)
    }

    fn choose_randomly_create_cache<'a>(
        mut fractions: impl Iterator<Item = &'a Self>,
    ) -> Result<FractionRandomCacheEnum>
    where
        Self: Sized,
        Self: 'a,
    {
        if is_exact_globally() {
            //exact mode
            if let Some(first) = fractions.next() {
                let mut cumulative_probabilities = vec![
                    first
                        .extract_exact()
                        .with_context(|| "cannot combine exact and approximate arithmetic")?
                        .clone(),
                ];
                let mut highest_denom = first.extract_exact()?.denom().unwrap();

                while let Some(fraction) = fractions.next() {
                    highest_denom = highest_denom.max(fraction.extract_exact()?.denom().unwrap());

                    let mut x = fraction
                        .extract_exact()
                        .with_context(|| "cannot combine exact and approximate arithmetic")?
                        .clone();
                    x += cumulative_probabilities.last().unwrap();
                    cumulative_probabilities.push(x);
                }
                let highest_denom = highest_denom.clone();

                Ok(FractionRandomCacheEnum::Exact(
                    cumulative_probabilities,
                    highest_denom,
                ))
            } else {
                Err(anyhow!("cannot take an element of an empty list"))
            }
        } else {
            //approximate mode
            if let Some(first) = fractions.next() {
                let mut cumulative_probabilities = vec![
                    *first
                        .extract_approx()
                        .with_context(|| "cannot combine exact and approximate arithmetic")?,
                ];

                while let Some(fraction) = fractions.next() {
                    cumulative_probabilities.push(
                        fraction
                            .extract_approx()
                            .with_context(|| "cannot combine exact and approximate arithmetic")?
                            + cumulative_probabilities.last().unwrap(),
                    );
                }

                Ok(FractionRandomCacheEnum::Approx(cumulative_probabilities))
            } else {
                Err(anyhow!("cannot take an element of an empty list"))
            }
        }
    }

    fn choose_randomly_cached(cache: &FractionRandomCacheEnum) -> usize
    where
        Self: Sized,
    {
        match cache {
            FractionRandomCacheEnum::Exact(cumulative_probabilities, highest_denom) => {
                //select a random value
                let mut rng = rand::thread_rng();
                let rand_val = {
                    //strategy: the highest denominator determines how much precision we need

                    //Generate a random value with the number of bits of the highest denominator. Repeat until this value is <= the max denominator.
                    let mut rand_val = rng.gen_biguint(highest_denom.bits());
                    while &rand_val > highest_denom {
                        rand_val = rng.gen_biguint(highest_denom.bits());
                    }
                    //create the fraction from the random nominator and the max denominator
                    GenericFraction::Rational(
                        Sign::Plus,
                        Ratio::new(rand_val, highest_denom.clone()),
                    )
                };

                match cumulative_probabilities.binary_search(&rand_val) {
                    Ok(index) | Err(index) => index,
                }
            }
            FractionRandomCacheEnum::Approx(cumulative_probabilities) => {
                //select a random value
                let mut rng = rand::thread_rng();
                let rand_val = rng.gen_range(0.0..=*cumulative_probabilities.last().unwrap());

                match cumulative_probabilities.binary_search_by(|probe| probe.total_cmp(&rand_val))
                {
                    Ok(index) | Err(index) => index,
                }
            }
        }
    }
}

pub struct FractionRandomCacheExact {
    cumulative_probabilities: Vec<FractionExact>,
    highest_denom: BigUint,
}

impl ChooseRandomly for FractionExact {
    type Cache = FractionRandomCacheExact;

    fn choose_randomly(fractions: &Vec<FractionExact>) -> Result<usize> {
        if fractions.is_empty() {
            return Err(anyhow!("cannot take an element of an empty list"));
        }

        //normalise the probabilities
        let mut probabilities: Vec<FractionExact> = fractions.iter().cloned().collect();
        let sum = probabilities
            .iter()
            .fold(FractionExact::zero(), |x, y| &x + y);
        probabilities.retain_mut(|v| {
            *v /= &sum;
            true
        });

        //select a random value
        let mut rng = rand::thread_rng();
        let rand_val = {
            //strategy: the highest denominator determines how much precision we need
            let max_denom = probabilities
                .iter()
                .map(|f| match f {
                    FractionExact(e) => e.denom().unwrap(),
                })
                .max()
                .unwrap();
            //Generate a random value with the number of bits of the highest denominator. Repeat until this value is <= the max denominator.
            let mut rand_val = rng.gen_biguint(max_denom.bits());
            while &rand_val > max_denom {
                rand_val = rng.gen_biguint(max_denom.bits());
            }
            //create the fraction from the random nominator and the max denominator
            FractionExact::try_from((rand_val, max_denom.clone())).unwrap()
        };

        let mut cum_prob = FractionExact::zero();
        for (index, value) in probabilities.iter().enumerate() {
            cum_prob += value;
            if rand_val < cum_prob {
                return Ok(index);
            }
        }
        Ok(probabilities.len() - 1)
    }

    fn choose_randomly_create_cache<'a>(
        mut fractions: impl Iterator<Item = &'a Self>,
    ) -> Result<FractionRandomCacheExact>
    where
        Self: Sized,
        Self: 'a,
    {
        if let Some(first) = fractions.next() {
            let mut cumulative_probabilities = vec![first.clone()];
            let mut highest_denom = first.0.denom().unwrap();

            while let Some(fraction) = fractions.next() {
                highest_denom = highest_denom.max(fraction.0.denom().unwrap());

                cumulative_probabilities.push(fraction + cumulative_probabilities.last().unwrap());
            }
            let highest_denom = highest_denom.clone();

            Ok(FractionRandomCacheExact {
                cumulative_probabilities,
                highest_denom,
            })
        } else {
            Err(anyhow!("cannot take an element of an empty list"))
        }
    }

    fn choose_randomly_cached(cache: &FractionRandomCacheExact) -> usize
    where
        Self: Sized,
    {
        //select a random value
        let mut rng = rand::thread_rng();
        let rand_val = {
            //strategy: the highest denominator determines how much precision we need

            //Generate a random value with the number of bits of the highest denominator. Repeat until this value is <= the max denominator.
            let mut rand_val = rng.gen_biguint(cache.highest_denom.bits());
            while rand_val > cache.highest_denom {
                rand_val = rng.gen_biguint(cache.highest_denom.bits());
            }
            //create the fraction from the random nominator and the max denominator
            FractionExact::try_from((rand_val, cache.highest_denom.clone())).unwrap()
        };

        match cache.cumulative_probabilities.binary_search(&rand_val) {
            Ok(index) | Err(index) => index,
        }
    }
}

pub struct FractionRandomCacheF64 {
    cumulative_probabilities: Vec<FractionF64>,
}

impl ChooseRandomly for FractionF64 {
    type Cache = FractionRandomCacheF64;

    fn choose_randomly(fractions: &Vec<FractionF64>) -> Result<usize> {
        if fractions.is_empty() {
            return Err(anyhow!("cannot take an element of an empty list"));
        }

        //normalise the probabilities
        let mut probabilities: Vec<FractionF64> = fractions.iter().cloned().collect();
        let sum = probabilities
            .iter()
            .fold(FractionF64::zero(), |x, y| &x + y);
        probabilities.retain_mut(|v| {
            *v /= &sum;
            true
        });

        //select a random value
        let mut rng = rand::thread_rng();
        let rand_val = FractionF64(rng.gen_range(0.0..=1.0));

        let mut cum_prob = FractionF64::zero();
        for (index, value) in probabilities.iter().enumerate() {
            cum_prob += value;
            if rand_val < cum_prob {
                return Ok(index);
            }
        }
        Ok(probabilities.len() - 1)
    }

    fn choose_randomly_create_cache<'a>(
        mut fractions: impl Iterator<Item = &'a Self>,
    ) -> Result<FractionRandomCacheF64>
    where
        Self: Sized,
        Self: 'a,
    {
        if let Some(first) = fractions.next() {
            let mut cumulative_probabilities = vec![*first];

            while let Some(fraction) = fractions.next() {
                cumulative_probabilities.push(fraction + cumulative_probabilities.last().unwrap());
            }

            Ok(FractionRandomCacheF64 {
                cumulative_probabilities,
            })
        } else {
            Err(anyhow!("cannot take an element of an empty list"))
        }
    }

    fn choose_randomly_cached(cache: &FractionRandomCacheF64) -> usize
    where
        Self: Sized,
    {
        //select a random value
        let mut rng = rand::thread_rng();
        let rand_val = FractionF64::from(
            rng.gen_range(
                0.0..=*cache
                    .cumulative_probabilities
                    .last()
                    .unwrap()
                    .extract_approx()
                    .unwrap(),
            ),
        );

        match cache.cumulative_probabilities.binary_search(&rand_val) {
            Ok(index) | Err(index) => index,
        }
    }
}
