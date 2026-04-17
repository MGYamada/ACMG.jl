"""
    ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category classification via p-fix modular data enumeration
over F_p, with subsequent F/R symbol construction.

Design philosophy:
- Primary object: modular data (S, T) over F_p
- Fusion rule N_ij^k is derived via Verlinde formula
- Classification proceeds prime-by-prime, union over admissible primes
- F, R symbols are solved from fusion rule directly (reusing v5 HC solver)

Module organization:
- FpArith:       F_p arithmetic primitives (inverse, sqrt, n-th roots)
- Dimensions:    Quantum dimension enumeration (small totally real algebraic integers)
- ModularData:   (S, T) enumeration and validation over F_p
- FusionExtract: Verlinde-based fusion rule extraction and integrality lift
- Admissibility: BNRW checks (Cauchy, FS) — stubs initially, to be filled from v5
"""
module ACMG

using LinearAlgebra
using Primes

# Core arithmetic layer
include("FpArith.jl")

# Modular data and fusion rule types
include("Types.jl")

# Quantum dimension enumeration
include("Dimensions.jl")

# (S, T) enumeration over F_p
include("ModularData.jl")

# Fusion rule extraction
include("FusionExtract.jl")

# Top-level enumeration driver
include("Enumerator.jl")

# Exports
export ModularDatumFp, FusionRule
export validate_modular_data, build_modular_datum, compute_alpha, compute_charge_conjugation
export extract_fusion_rule_Fp, lift_fusion_to_Z, extract_and_lift
export verlinde_coefficient
export is_square, sqrt_mod, primitive_root, root_of_unity, roots_of_unity
export matmul_mod, matpow_mod, diagmul_right, diagmul_left, lift_symmetric
export fusion_isomorphic, fusion_matrix, validate

end # module ACMG
