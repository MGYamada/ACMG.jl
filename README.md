# ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category classification via p-fix modular data enumeration.

## Overview

This is the v6 rewrite of the MTC classification pipeline, built on the
**p-fix** strategy: instead of fixing a conductor N and working over Z[ζ_N],
we fix a prime p and enumerate modular data over F_p. Multiple primes are
run in parallel, and the resulting candidates are lifted to Z and classified
by fusion rule.

## Current status (v0.1.0-prototype)

**Implemented:**
- F_p arithmetic primitives (`FpArith.jl`):
  - Tonelli-Shanks square root
  - Euler criterion (QR test)
  - Primitive roots, n-th roots of unity
  - Matrix operations mod p
- Core types (`Types.jl`):
  - `ModularDatumFp`: (S, T) over F_p with metadata
  - `FusionRule`: ring-agnostic fusion data with axiom validation
- Modular data validation (`ModularData.jl`):
  - Axiom checking for (S, T)
  - α (central charge phase) extraction from (ST)³ = α · S²
  - Charge conjugation extraction from S²
- Verlinde extraction (`FusionExtract.jl`):
  - F_p fusion coefficients
  - Lift to Z with non-negativity check
  - Full pipeline: ModularDatumFp → FusionRule

**Stubbed (not yet implemented):**
- Quantum dimension enumeration (`Dimensions.jl`)
- (S, T) enumeration driver (`Enumerator.jl`)

## Tests

Tests are direct encodings of the hand calculations for Fibonacci and Ising:

- `test_fparith.jl`: F_p primitive correctness
- `test_fibonacci.jl`: Fibonacci at p=41 (N=5, rank 2)
  - ζ_5 = 10, √5 = 13, φ = 7, D² = 9, D = 3
  - S = [14 16; 16 27], T = [1, 18]
  - (ST)³ = 21·I, α = 21 = ζ_{20}^7 ↔ c = 14/5
- `test_ising.jl`: Ising at p=17 (N=16, rank 3)
  - ζ_16 = 3, √2 = 11, D² = 4, D = 2
  - S = [9 14 9; 14 0 3; 9 3 9], T = [1, 3, 16]
  - (ST)³ = 3·I, α = 3 = ζ_{16} ↔ c = 1/2

Run with:
```julia
] test ACMG
```

## Design notes

### Primary object: modular data, not fusion rule

Earlier discussion considered enumerating fusion rules directly as a
scheme F_r over Z. This was rejected because:
1. Fusion rings are a much larger class than MTCs (most are non-modular)
2. Modular condition is hard to check from fusion rule alone
3. v5's admissibility checks (BNRW) are already built around (S, T)

We keep (S, T) as primary and obtain fusion rule as a projection.

### Prime selection for rank ≤ 4

Target primes for rank ≤ 4 enumeration:
- p = 41: covers N | 40 = 2³ · 5
- p = 61: covers N | 60 = 2² · 3 · 5
- p = 97: covers N | 96 = 2⁵ · 3
- p = 113: covers N | 112 = 2⁴ · 7
- p = 181: covers N | 180
- p = 241: covers N | 240

A prime p is **admissible** for an MTC C if:
- N | p-1 (ζ_N ∈ F_p)
- D² is a square in F_p (so D ∈ F_p; else work in F_{p²})
- All N_ij^k lift to small non-negative integers

The square-of-D² condition is the restrictive one: Chebotarev gives density
≥ 1/2 among splitting primes, so several primes are needed.

### Next steps

1. Implement `Dimensions.jl`: enumerate totally real algebraic integer tuples
   (d_0=1, d_1, ..., d_{r-1}) with Σd_i² ≤ D_max² and d_i ∈ μ_∞-integral.
2. Implement S-solver: given T and {d_i}, solve (ST)² = α·T⁻¹·S using
   row-by-row pruning via SS† = D²·C.
3. Wire up `Enumerator.jl` as top-level driver.
4. Validate against NRW NsdGOL4.g for rank ≤ 4.
5. Extend to rank 5 and beyond (target NRW NsdGOL5.g with (3,2) type at N=24).
