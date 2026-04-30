# Concepts

ACMG uses a conductor-first model:

```text
N -> Q(zeta_N) -> exact S/T/F/R data -> Galois and finite-field reductions
```

Main modules:

- Cyclotomics: exact cyclotomic contexts and reductions.
- ModularData: S/T data, validation, and Verlinde extraction.
- FR: F/R equation infrastructure and exact data containers.
- Gauge: gauge parameters and gauge-fixing records.
- BraidRepresentations: braid matrices from F/R data.
- Search: stratum and block-U enumeration.

Representative functions include `modular_data`, `validate_modular_data`,
`fr_equation_system`, `gauge_fix`, and `braid_representation`.

Notes: finite-field and diagnostic paths are useful for experiments, but exact
verification is still required before treating a result as a mathematical
classification statement.
