# Finite Fields

Finite-field support reduces cyclotomic and modular data modulo split primes
and provides prototype finite-field F/R experiments.

Representative functions:

- `reduce_mod_p`
- `find_zeta_in_Fp`, `cyclotomic_to_Fp`
- `frobenius_metadata`
- `solve_FR_mod_p`
- `solve_finite_field`, `cyclotomic_reconstruct`

Notes: exact modular-data reduction is part of the main workflow.  Finite-field
F/R solving and reconstruction are experimental APIs; residues are evidence
that require exact lifting and verification.
