# Ebi_arithmetic

This is a package with arithmetic options for the crate Ebi.

The aim of this package is to provide exact and approximate arithmetic in an as transparent as possible fashion.

This package is still subject to change and may break compatibility in minor releases.

The most common entry point is the `Fraction' type alias, which will point to a fraction type based on the given compile flags.

To conveniently declare a constant, use a macro:

```
f0!() = Fraction::zero();
f1!() = Fraction::one();
f!(0) = Fraction::from(0);
f!(1, 2) = Fraction::from((1, 2))
```