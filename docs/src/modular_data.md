# Modular Data

The modular-data layer constructs and validates exact S/T data over a
`CyclotomicContext`.

Representative functions:

- `modular_data`
- `semion_modular_data`, `fibonacci_modular_data`, `ising_modular_data`
- `validate_modular_data`, `validate_exact_modular_data`
- `check_modular_relations`, `check_unitarity`
- `verlinde_coefficient`, `extract_fusion_rule_Fp`

Notes: built-in examples are stable public API.  Low-level finite-field
extraction helpers are useful for tests and examples, but users should prefer
the higher-level constructors and validation functions when possible.
