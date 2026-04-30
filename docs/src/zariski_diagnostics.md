# Zariski Diagnostics

The Zariski diagnostics layer inspects finite-field braid images using matrix
algebra, commutant, determinant, projective-order, and sampled relation data.

Representative functions:

- `generated_matrix_algebra`
- `commutant`
- `generated_subgroup`
- `finite_group_diagnostics`
- `zariski_closure_diagnostics`

Notes: these APIs are experimental.  Their return values are computational
evidence, not mathematical theorems, and they do not compute a full
characteristic-zero Zariski closure.
