"""
    ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category classification by fixing the conductor `N`,
running an SL(2, ℤ/N) stratum + block-U enumeration, multi-prime F_p
sweep + CRT reconstruction, and pentagon/hexagon `(F, R)` solution.

Top-level pipeline (Phase 5):

    using ACMG
    classified = classify_mtcs_at_conductor(5; max_rank = 2,
                                             primes = [11, 31, 41, 61])

returns a `Vector{ClassifiedMTC}`, one per Galois sector per stratum.

Module organisation:
- FpArith:         F_p arithmetic primitives (Tonelli-Shanks, primitive
                   roots, matmul, lift_symmetric, ...)
- Types:           `ModularDatumFp`, `FusionRule` with axiom validation
- Dimensions:      quantum dimension enumeration (stub)
- ModularData:     (S, T) axiom checking over F_p
- FusionExtract:   Verlinde-based fusion rule extraction
- Enumerator:      top-level legacy driver (stub)
- SL2Reps:         SL(2, ℤ/N) irreducible representation catalog
                   (via Oscar + GAP/SL2Reps)  [Phase 0]
- StratumEnum:     combinatorial partition `Σ m_λ d_λ = r`  [Phase 1]
- BlockU:          general O(n) Cayley + single-prime MTC driver
                   (`find_mtcs_at_prime`)  [Phase 2]
- CRT:             multi-prime CRT reconstruction + Galois-aware
                   grouping  [Phase 3]
- Phase4:          submodule — pentagon/hexagon (F, R) solver
                   + ribbon verify  [Phase 4]
- Phase5:          this file's `classify_mtcs_at_conductor` + friends,
                   wiring everything into an `N → List[ClassifiedMTC]`
                   pipeline  [Phase 5]
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

# Top-level enumeration driver (legacy)
include("Enumerator.jl")

# SL(2, ℤ/N) irrep catalog (Oscar + GAP/SL2Reps) — Phase 0
include("SL2Reps.jl")

# Stratum enumeration (combinatorial partition of rank by irrep dimensions) — Phase 1
include("StratumEnum.jl")

# Block-U parametrisation and MTC reconstruction — Phase 2
include("BlockU.jl")

# Multi-prime CRT reconstruction — Phase 3
include("CRT.jl")

# Phase 4: Pentagon/Hexagon solver for (F, R) classification — submodule
include("Phase4/Phase4.jl")

# Phase 5: end-to-end pipeline driver — uses Phase 4 submodule + all of the above
include("Phase5.jl")

# ============================================================
#  Exports
# ============================================================

# Core types and Fp arithmetic
export ModularDatumFp, FusionRule
export validate_modular_data, build_modular_datum, compute_alpha, compute_charge_conjugation
export extract_fusion_rule_Fp, lift_fusion_to_Z, extract_and_lift
export verlinde_coefficient
export is_square, sqrt_mod, primitive_root, root_of_unity, roots_of_unity
export matmul_mod, matpow_mod, diagmul_right, diagmul_left, lift_symmetric
export fusion_isomorphic, fusion_matrix, validate

# Phase 0: atomic catalog
export AtomicIrrep, build_atomic_catalog, all_divisors

# Phase 1: strata
export Stratum, enumerate_strata, count_strata, describe_stratum, find_unit_indices

# Phase 2: block-U
export MTCCandidate, build_block_diagonal, reduce_matrix_to_Fp, reduce_vector_to_Fp
export find_zeta_in_Fp, cyclotomic_to_Fp
export t_eigenspace_decomposition, parameter_dim
export o2_circle_points, apply_o2_block, verlinde_find_unit
export cayley_so_n, inverse_mod_p, enumerate_so_n_Fp, enumerate_o_n_Fp, apply_block_U
export find_mtcs_at_prime, signed_Fp

# Phase 3: CRT reconstruction
export acmg_crt, crt2, rational_reconstruct, compute_sqrt_d_mod_p
export compute_sqrt3_cyclotomic_mod_p, compute_sqrt2_cyclotomic_mod_p
export fusion_signature, group_mtcs_by_fusion, group_mtcs_galois_aware
export reconstruct_rational, reconstruct_in_Z_sqrt_d
export reconstruct_matrix_in_Z_sqrt_d, reconstruct_S_matrix
export verify_reconstruction, describe_matrix

# Phase 4: (F, R) classification — re-exported at ACMG top level.
#
# Users can now write `ACMG.verify_mtc(...)`, `ACMG.get_pentagon_system(...)`
# etc. directly without `ACMG.Phase4.` prefixing. The `Phase4` submodule
# is still accessible as `ACMG.Phase4` for those who want explicit
# namespacing.
using .Phase4: get_pentagon_system,
               solve_pentagon_newton, solve_pentagon_homotopy,
               refine_solution_newton,
               hexagon_equations, get_hexagon_system,
               solve_hexagon_homotopy,
               assign_F_to_associator!, invert_associator_numeric,
               DiscreteLogTable, lift_T_Fp_to_complex,
               lift_S_sqrtd_to_complex, lift_mtc_candidate,
               pentagon_residuals, hexagon_residuals,
               extract_R_block, block_positions_R,
               ribbon_residuals, VerifyReport, verify_mtc

export get_pentagon_system
export solve_pentagon_newton, solve_pentagon_homotopy, refine_solution_newton
export hexagon_equations, get_hexagon_system
export solve_hexagon_homotopy
export assign_F_to_associator!, invert_associator_numeric
export DiscreteLogTable, lift_T_Fp_to_complex, lift_S_sqrtd_to_complex
export lift_mtc_candidate
export pentagon_residuals, hexagon_residuals
export extract_R_block, block_positions_R
export ribbon_residuals, VerifyReport, verify_mtc

# Phase 5: end-to-end pipeline
export ClassifiedMTC
export compute_FR_from_ST, classify_from_group, classify_mtcs_at_conductor

end # module ACMG
