"""
Quantum dimension enumeration (STUB).

Quantum dimensions d_i must be totally positive algebraic integers
in a real cyclotomic field Q(ζ_N + ζ_N^{-1}), with d_0 = 1 and
D² = Σ d_i² bounded.

To be implemented:
- enumerate_quantum_dims(rank, D_max; field_info) -> Vector{Vector{R}}
- lift to F_p via root-of-unity substitution

For prototype, we provide known d-tuples for specific MTCs for testing.
"""

# Known quantum dimensions for test MTCs (abstract form, to be realized in F_p)
const KNOWN_QUANTUM_DIMS = Dict(
    "trivial"  => (rank=1, d_squared_sum=1, description="d = [1]"),
    "semion"   => (rank=2, d_squared_sum=2, description="d = [1, 1]"),
    "fibonacci"=> (rank=2, d_squared_sum_formula="1 + φ²", description="d = [1, φ], φ = (1+√5)/2"),
    "ising"    => (rank=3, d_squared_sum=4, description="d = [1, √2, 1]"),
    "su2_3"    => (rank=3, d_squared_sum_formula="1 + φ² + 1 = 3 + φ²", description="d = [1, φ, 1] or similar"),
    "z3"       => (rank=3, d_squared_sum=3, description="d = [1, 1, 1]"),
)
