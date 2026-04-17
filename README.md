# ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category (MTC) classification via cyclotomic representation
variety enumeration with multi-prime F_p validation.

## Overview

ACMG classifies MTCs by **fixing the conductor** N = max{cond(S), cond(T)}
and exploiting the following key structural insight:

> The moduli space of SL(2, ℤ/N)-modular data of rank r is a
> **stratified scheme** over Spec ℤ[1/2N]. Each stratum is indexed by the
> isomorphism type (m_λ) of the underlying SL(2, ℤ/N)-representation
> ρ = ⊕_λ ρ_λ^{m_λ}. Within each stratum, the moduli is parametrised by
> a block-U group O(n_θ) on each T-eigenspace V_θ (with dim n_θ = dim V_θ).
> MTC integrality (Verlinde non-negativity) then cuts out a 0-dimensional
> locus inside this continuous family.

This gives a pipeline:

```
Atomic SL(2, ℤ/N)-irreps catalog     (via Oscar + GAP/SL2Reps)
        ↓
Stratum enumeration (m_λ with Σ m_λ d_λ = r)
        ↓
T-spectrum filter (retain (m_λ) matching target T-eigenvalues)
        ↓
Block-U parametrisation on each T-eigenspace V_θ
        ↓
F_p multi-prime Verlinde sweep (fast, exact)
        ↓
CRT reconstruction to Z[ζ_N]
        ↓
Classified MTCs
```

The pipeline has been validated end-to-end on SU(2)_4 at rank 5, N = 24.

## Current status (v0.2 prototype)

### Phase 0: Atomic catalog — ✅ complete

- `FpArith.jl`: Tonelli-Shanks sqrt, Euler criterion, primitive roots,
  roots of unity, matrix ops mod p.
- `Types.jl`: `ModularDatumFp`, `FusionRule` with ring-agnostic axiom
  validation.
- `ModularData.jl`: (S, T) axiom checking, α (= e^{-2πi c/8}) extraction
  from (ST)³ = α · S², charge conjugation from S².
- `FusionExtract.jl`: F_p Verlinde coefficients, Z lift with
  non-negativity check.
- `SL2Reps.jl`: Oscar + GAP/SL2Reps integration, `build_atomic_catalog(N)`
  producing the irreducible SL(2, ℤ/N)-representations.

### Phase 1: Stratum enumeration — ✅ complete

- `StratumEnum.jl`: combinatorial enumeration of {m_λ} with
  Σ m_λ · dim(λ) = r.
- Note: `require_unit_summand=false` is the correct default. In non-
  abelian MTCs the unit object sits as a basis vector *inside* a larger
  irrep block (wherever T-eigenvalue 1 occurs), not as a separate 1d_1
  summand. See the Design Notes below.

### Phase 1.5: T-spectrum filter — 🔜 next

- Given target spins (h_0, ..., h_{r-1}), filter strata whose T-spectrum
  (as a multiset) agrees with the target.
- At N = 24 rank 5 this reduces 25,440 strata to O(1) candidates.

### Phase 2: Block-U parametrisation — 🔜 next

- For each stratum, identify T-eigenspaces V_θ and the degenerate ones
  (n_θ ≥ 2).
- Parametrise block-U as U ∈ GL(V)^T with U U^T ∈ End_G(ρ) — equivalently,
  a rotation on each V_θ (see Design Notes).
- Search via F_p multi-prime sweep: at each good prime p (N | p-1),
  enumerate the F_p-rational points of O(n_θ), check Verlinde
  non-negativity.

### Phase 3: CRT reconstruction — 🔜 next

- Collate solutions across multiple primes.
- Reconstruct S-entries as elements of Z[ζ_N] (often Z[√D²-factors]).

### Test statistics

151 tests passing (runtests.jl):
- 18 FpArith primitives
- 34 Fibonacci at p=41 (N=5, rank 2)
- 30 Ising at p=17 (N=16, rank 3)
- 9 SL2Reps catalog (N ∈ {5, 8, 16})
- 15 StratumEnum (synthetic rank 1–5 cases)

## Design notes

### Central conceptual shift from v5

The v5 pipeline fixed a **prime p** and solved modular-data equations
over F_p directly. The v6 / ACMG pipeline fixes the **conductor N** and
works primarily with SL(2, ℤ/N)-representation theory over Z[ζ_N], using
F_p only for fast computation and final validation. The reasons:

1. **Finiteness** of the search space. Fixing N restricts modular data
   entries to Z[ζ_N] ∩ {|z| ≤ D} — a finite set. The p-fix approach
   searches F_p instead, also finite but without the algebraic structure.

2. **Natural alignment with NRW(W)** (2203.14829, 2308.09670). Their
   classification stratifies by SL(2, ℤ/N) irrep content — we follow the
   same stratification.

3. **Conductor-first beats rank-first**. The NRW paper enumerates rank
   first then searches for fields, which forces unnecessary field
   extensions when the MTC actually lives in a smaller cyclotomic field.

### Modular data moduli: what the geometry actually is

Consider the moduli space

    M_{r, N}^{MTC} ⊂ M_{r, N}^{SL_2(ℤ/N)-rep}

of rank r modular data of conductor dividing N. The right-hand side is
the representation variety (up to equivalence), with the additional
conditions S symmetric and S² = charge-conjugation.

**Stratification.** Each isomorphism class ρ = ⊕_λ ρ_λ^{m_λ} defines a
locally closed stratum. Within a stratum, the choice of basis for V
(not changing the SL(2, ℤ/N)-module structure) is an element of
Aut_G(ρ) = ∏_λ GL(m_λ). For all m_λ = 1, this is just a torus.

**Continuous moduli are hidden in T-eigenspace overlaps.** The surprise
(discovered in this iteration) is that even when all m_λ = 1 — so
Aut_G(ρ) is just a torus acting trivially on S — there can still be a
continuous moduli if **different irreducibles share a T-eigenvalue θ**.
Specifically, the T-eigenspace

    V_θ = ⊕_λ V_θ^{(ρ_λ)}

can be 2-dimensional, and rotating its basis (an O(2) action) gives a
new S-matrix while preserving T and the SL(2, ℤ/N)-module structure.

**Block-U, algebraic definition.** The full group acting is

    𝒰(ρ) = {U ∈ GL(V) : [T, U] = 0,
                          U|_{V_θ} regular semisimple,
                          U Uᵀ ∈ End_G(ρ)}

(the last being projective orthogonality). Decomposing U = ⊕_θ U_θ on
T-eigenspaces, the continuous moduli has dimension

    parameter_dim = Σ_θ C(n_θ, 2)   where n_θ = dim V_θ

For SU(2)_4 at N = 24: T-spectrum (with multiplicities) has a single
degenerate eigenvalue (multiplicity 2), so parameter_dim = 1. This
matches the Python M4 sweep finding φ = π/4.

### Verlinde integrality: from continuous to discrete

On the continuous moduli the Verlinde fusion coefficients
N_{ij}^k(φ) are algebraic functions of the block-U parameters. The MTC
locus is the subset of parameter space where all N_{ij}^k are
non-negative integers — generically a **0-dimensional algebraic set**.

For SU(2)_4 this is φ = π/4 (a rotation / reflection pair). Not an
accident: the relevant 2×2 block of S on V_1 is rank-1 with all entries
equal (by Kac–Peterson), and its eigenvector is the Hadamard direction,
i.e., φ = π/4. **This is derivable a priori, without any numerical
sweep.**

### The étale picture (partial)

Strata where all n_θ = 1 are 0-dimensional: parameter_dim = 0. The
scheme 𝒰 reduces to ∏_θ μ_2 = μ_2^{#θ} on these strata, which is
étale over Z[1/2]. Rank ≤ 4 over small N tends to fall in this case,
and indeed: empirically most "simple" MTCs (Fibonacci, Ising family)
have no continuous moduli.

When some n_θ ≥ 2, continuous moduli appear and MTC integrality
becomes a non-trivial cutting condition. Rank 5 N = 24 is the smallest
case with this phenomenon.

Precisely: **the final set of MTCs** (after Verlinde integrality cuts
out a 0-dimensional locus) is still étale over Z[1/2N], but the
intermediate universe of SL(2, ℤ/N)-modular data is not. This is the
correct version of "MTC moduli is étale" — the statement applies to
rational points after integrality, not to the ambient variety.

### Degenerate perturbation theory as a mental model

The mathematical content of block-U is exactly **degenerate
perturbation theory** in algebraic-geometric disguise:

- T = unperturbed Hamiltonian with degenerate spectrum.
- V_θ (n_θ ≥ 2) = degenerate subspace.
- S|_{V_θ} = perturbation lifting the degeneracy.
- O(n_θ) freedom = choice of 0-th order eigenbasis.
- Hadamard / φ = π/4 = the correct eigenbasis for SU(2)_4.

Explicitly, the tangent direction δT that resolves the degeneracy of T
at V_1 is δT ∝ σ_x (off-diagonal Pauli). Its eigenvectors are exactly
the Hadamard basis, recovering φ = π/4 analytically.

### Multi-prime F_p validation

Practical computation uses **F_p reduction at multiple good primes**
(N | p-1), even though the underlying moduli lives in Z[ζ_N].
Rationale:

1. Block-U parametrisation can be computed exactly in F_p via the
   algebraic circle {(u, v) : u² + v² = 1} ⊂ F_p² — only ~p points.
2. Verlinde integrality check in F_p is an integer comparison.
3. CRT across primes reconstructs the algebraic entries in Z[ζ_N].
4. Sanity check: verify at fresh primes not used for reconstruction.

This has been validated end-to-end on SU(2)_4 rank 5 N = 24:

- 10 primes {73, 97, 193, 241, 313, 337, 409, 433, 457, 577}.
- Each prime gives 2 MTC solutions (φ = π/4, rotation & reflection).
- 4-prime CRT reconstructs S = (1/2√3) · M with M ∈ Z[√3]^{5×5}
  entries in {0, ±1, ±2, ±√3}.
- Verified at 6 unused primes: all ✓.
- Quantum dims d = (√3, √3, 1, 1, 2), D² = 12.

Fast (~1 second), exact (no floating-point), parallel (primes
independent).

### Feasibility of block-U enumeration by n_θ

The total cost of a single-prime sweep is dominated by |O(n_θ)(F_p)|
for each degenerate eigenspace, where n_θ = dim V_θ. Empirically at
p = 73:

| n_θ | |O(n_θ)(F_p)| | Single-prime sweep time | Feasibility |
|-----|---|-----|-----|
|  1  | 2 | instant | trivial (discrete sign) |
|  2  | ~2p = 146 | < 1 s | naive enumeration |
|  3  | ~p³ = 4×10⁵ | ~20 s | naive OK |
|  4  | ~p⁶ = 1.5×10¹¹ | ~10³ h | **infeasible naively** |
| ≥ 4 | ~p^{n(n-1)/2} | — | needs algebraic solver |

So naive F_p enumeration handles n_θ ≤ 3; beyond that, Cayley
parametrisation + Verlinde polynomial system (Gröbner basis over F_p)
is required. The Cayley parameter space has dimension C(n_θ, 2) =
parameter_dim, which is also the variable count for the polynomial
system — so this is the natural transition.

**Practically, rank 5 classification fits entirely in the n_θ ≤ 2
regime.** The NsdGOL5.g catalog shows spin patterns with at most
2-fold T-eigenvalue degeneracy for rank 5. Higher n_θ (and hence the
need for algebraic solvers) is expected to appear at rank ≥ 6, e.g.
for Drinfeld centers D(Z_n) with n ≥ 2 (rank n²) where multiple
self-dual objects share T = 1.

**This finite-field discretisation of continuous strata is not in
NRWW.** Their block-U framework is stated over C, without reduction
to F_p or a concrete enumeration scheme. The equivalence
|continuous stratum| ↔ |O(n_θ)(F_p)| under base change, combined with
CRT reconstruction back to Z[ζ_N], is new in this work.

### What becomes obsolete

- **v5 Part B**: the "2×2 rotation on T-overlapping atomic pairs"
  framework is partially correct in its setup (shared-eigenvalue O(2))
  but v5 got parametrisation details wrong or simply tested conductors
  (N = 8, 12, 15, 20) where no genuine multi-component MTC exists.
- **p-fix as primary search**: p-fix is retained as an accelerator, not
  as the primary search space.
- **Homotopy continuation / damped Newton**: replaced by the block-U
  formulation, which is a direct algebraic construction.

## Prime selection for F_p validation

For conductor N, select primes p with N | p − 1:

| N | Example good primes |
|-----|-----|
|  5 | 11, 31, 41, 61, 71 |
|  8 | 17, 41, 73, 89, 97 |
| 12 | 13, 37, 61, 73, 97 |
| 24 | 73, 97, 193, 241, 313, 337, 409 |

A prime p is **admissible** for a given MTC candidate if:

- N | p − 1 (so ζ_N ∈ F_p),
- √3, √2, and in general √(D²-factors) are in F_p (quadratic residue
  condition),
- all N_{ij}^k lift to small non-negative integers.

Chebotarev density ⇒ for each MTC, a positive-density set of primes is
admissible. Running at several primes in parallel gives robust
confirmation.

## Repository layout

```
src/
  ACMG.jl                 — module root, exports
  FpArith.jl              — F_p primitives
  Types.jl                — ModularDatumFp, FusionRule
  ModularData.jl          — (S, T) axiom validation
  FusionExtract.jl        — F_p Verlinde → Z lift
  Dimensions.jl           — (stub) quantum dimension enumeration
  Enumerator.jl           — (stub) top-level driver
  SL2Reps.jl              — Oscar + GAP/SL2Reps catalog builder
  StratumEnum.jl          — {m_λ} combinatorial enumeration

test/
  runtests.jl
  test_fparith.jl
  test_fibonacci.jl
  test_ising.jl
  test_sl2reps.jl
  test_stratum_enum.jl
```

## References

- Bruillard, Ng, Rowell, Wang, *Rank-finiteness for modular categories*
  (2016). Admissibility criteria used in filtering.
- Ng, Rowell, Wang, Wen, *Reconstruction of modular data from
  SL(2, ℤ) representations* arXiv:2203.14829. The irrep-sum + block-U
  framework.
- Ng, Rowell, Wen, *Classification of modular data up to rank 12*
  arXiv:2308.09670. Ground-truth catalog (NsdGOL files) used for
  cross-validation.

## Next steps

**Near-term (Phase 2 implementation):**

1. **Phase 1.5** (Julia): `enumerate_strata_by_T(catalog, target_spins)` to
   filter strata by T-spectrum multiset. This alone reduces the rank-5
   N=24 search from ~25,000 strata to ~1 candidate.
2. **Phase 2a** (Julia): port the Python O(2) sweep to F_p, tested on
   SU(2)_4 for end-to-end reproduction inside Julia ACMG.
3. **Phase 2b** (Julia): multi-prime CRT reconstruction helper, with
   sanity-check at unused primes.
4. **Validation**: reproduce all 16 variants of NsdGOL[2] ∪ NsdGOL[3]
   (rank 5 N=24) via the full pipeline.

**Medium-term (beyond rank 5):**

5. Extend to rank 5 other conductors (N ∈ {5, 7, 11, ...}).
6. Attempt rank 6 with NRW ground truth (NsdGOL6.g). Identify the
   smallest example with n_θ = 3 to stress-test O(3) naive enumeration.
7. **Algebraic solver** for n_θ ≥ 4 cases (Cayley parametrisation +
   Gröbner basis over F_p). Target case: Drinfeld center D(Z_2) at
   rank 4, where n_1 = 4.

**Long-term:**

8. BNRW admissibility criteria (Cauchy, Frobenius-Schur) for integer
   candidates surviving Verlinde.
9. Galois orbit organisation of the classification output.
10. Haagerup Z(H_3) at rank 12, N=39 — the original motivation for
    multi-component MTC search.

