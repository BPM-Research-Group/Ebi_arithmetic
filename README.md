# Ebi_arithmetic

This is a package with arithmetic structs for the process mining crate Ebi.

The aim of this package is to provide both exact and approximate arithmetic in an as transparent as possible fashion. That is, users can interact solely with the `Fraction` and `FractionMatrix` aliases.

These aliases are changed to one of several structs, based on the compile features. If the crate is compiled without compile features, there is an atomic boolean that selects whether exact or approximate arithmetic is used, such that the mode can be selected at runtime. If the `exactarithmetic` or `approximatearithmetic` compile feature is provided, then the crate will function entirely in exact or approximate mode. A mode selected by compilation feature is generally faster than a runtime-selected mode.

To conveniently declare a constant `Fraction`, a macro can be used:

```
f0!() = Fraction::zero();
f1!() = Fraction::one();
f!(0) = Fraction::from(0);
f!(1, 2) = Fraction::from((1, 2))
```

This package is still subject to change and may break compatibility in minor releases.