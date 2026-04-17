# ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category (MTC) classification via block-U moduli enumeration
over cyclotomic fields and finite-field reductions.

## Overview

ACMG is a computational pipeline for classifying Modular Tensor Categories
(MTCs) of rank $r$ and conductor $N$. The core idea is to parametrize MTCs as
rational points of a stratified scheme constructed from the representation
variety of $SL_2(\mathbb{Z}/N)$, with the Verlinde integrality condition
cutting out a 0-dimensional locus.

This is the v6 rewrite, replacing the earlier "p-fix + homotopy continuation"
approach of v5 with a cleaner scheme-theoretic pipeline whose final step uses
multi-prime $\mathbb{F}_p$ arithmetic and CRT reconstruction.

## Current status (v0.2, Phase 0 + Phase 1 complete)

**Implemented:**
- `FpArith.jl`: $\mathbb{F}_p$ primitives (Tonelli-Shanks, primitive roots,
  matrix ops mod $p$)
- `Types.jl`, `ModularData.jl`, `FusionExtract.jl`: core types and validators
  for $(S, T)$ over $\mathbb{F}_p$, axiom checking, Verlinde extraction
- `SL2Reps.jl`: Oscar.jl + GAP/SL2Reps bridge; builds atomic catalog of
  $SL_2(\mathbb{Z}/N)$ irreducible representations via cyclotomic arithmetic
- `StratumEnum.jl`: combinatorial enumeration of irreducible decomposition
  patterns $\{m_\lambda\}$ with $\sum_\lambda m_\lambda \dim(\rho_\lambda) = r$

**Tests:** 151 passing (as of 2026-04-17)
- F_p arithmetic: 18 tests
- Fibonacci (N=5 at p=41): 34 tests
- Ising (N=16 at p=17): 30 tests
- SL2Reps bridge (N=5, 8, 16 catalogs): 9 tests
- StratumEnum: 15 tests across 6 testsets

**Not yet implemented:**
- Phase 1.5: T-spectrum filter (reduce ~$10^4$ strata to single target)
- Phase 2: Block-U parametrization and MTC reconstruction
- Phase 3: BNRW admissibility checks (Cauchy, Frobenius-Schur)
- Phase 4: Galois orbit and central charge classification

## Design notes

### The moduli-theoretic viewpoint

MTCs of rank $r$ and conductor $N$ correspond to rational points of a
stratified scheme
$$\mathcal{M}_{r,N}^{\text{MTC}} \subset \mathcal{M}_{r,N}^{\text{rep}}$$
where $\mathcal{M}_{r,N}^{\text{rep}}$ is the representation variety of
$SL_2(\mathbb{Z}/N)$ in dimension $r$, and the MTC locus is cut out by
Verlinde non-negative integrality of fusion coefficients.

The stratification is indexed by tuples $(m_\lambda)$ giving the multiplicity
of each atomic irrep $\rho_\lambda$ in the decomposition
$\rho = \bigoplus_\lambda \rho_\lambda^{m_\lambda}$ with $\sum m_\lambda \dim \rho_\lambda = r$.

### Block-U structure

Fix a stratum $(m_\lambda)$ and diagonalize $T = \rho(\mathfrak{t})$. Each
$T$-eigenspace $V_\theta$ has dimension $n_\theta = \sum_\lambda m_\lambda \cdot (\text{mult of } \theta \text{ in } \rho_\lambda)$.

**Block-U group.** The subspace of $GL(V)$ commuting with $T$ and preserving
the $SL_2(\mathbb{Z}/N)$-module structure up to scalar on each irrep block is
$$\mathcal{U}(\rho) = \{U \in GL(V) : TU = UT,\ U|_{V_\theta} \text{ regular semisimple},\ UU^T \in \text{End}_G(\rho)\}$$
This is the "projective orthogonal group of the representation", and acts on
modular data by $S \mapsto U^T S U$.

**Continuous moduli.** When a $T$-eigenvalue has multiplicity $n_\theta \geq 2$,
$V_\theta$ carries a continuous $O(n_\theta)$-action. The parameter dimension
of the moduli is
$$\text{parameter\_dim} = \sum_\theta \binom{n_\theta}{2}$$
Strata with all $n_\theta = 1$ give purely discrete MTC moduli (étale over
$\mathbb{Z}[\frac{1}{2N}]$); strata with $n_\theta \geq 2$ have continuous
representation-theoretic moduli, discretized by Verlinde integrality.

### Degenerate perturbation-theoretic interpretation

The block-U moduli is a geometric avatar of **degenerate perturbation theory**:
- $T$ is the unperturbed Hamiltonian with degenerate spectrum $\{V_\theta\}$
- $S$ is the perturbation, breaking the degeneracy
- $O(n_\theta)$ parametrizes the choice of 0-th order eigenbasis on $V_\theta$
- Correct eigenbasis = eigenvectors of $S|_{V_\theta}$
- Verlinde integrality selects special perturbation directions

**Example (SU(2)_4 at rank 5, N=24):** the degenerate eigenspace $V_1$ has
dimension 2, spanned by two objects with topological spin 0. The block $S|_{V_1}$
is rank 1 with uniform entries in the (1,1) direction, forcing its eigenbasis
to be the Hadamard basis — equivalently, a 45° reflection, $\phi = \pi/4$.
This is derivable *a priori* from the Kac-Peterson formula, without sweeping.

### The classification pipeline

```
Phase 0: Atomic catalog
  Use Oscar.jl + GAP/SL2Reps to enumerate all irreps of SL_2(Z/N)
  Each irrep: (dimension, level, parity, S-matrix, T-diagonal) in Z[ζ_N]

Phase 1: Stratum enumeration
  Enumerate partitions {m_λ} with Σ m_λ · dim(ρ_λ) = r
  (Combinatorial, by dimension only)

Phase 1.5: T-spectrum filter (planned)
  For a target spin profile, filter strata by predicted T-eigenvalue multiset
  Empirically: reduces ~25,000 rank-5 strata at N=24 to 1 candidate

Phase 2: Block-U parametrization & reconstruction (planned)
  For the filtered stratum, parametrize the block-U group
  Execute multi-prime F_p sweep for Verlinde-integer points
  CRT-reconstruct to Z[ζ_N]

Phase 3: BNRW admissibility (planned)
  Cauchy theorem, Frobenius-Schur indicators
  (Reuse from v5 where possible)

Phase 4: Galois orbit organization (planned)
  Collect MTCs into Galois orbits, assign central charge mod 8
```

### The multi-prime F_p + CRT strategy

The final block-U step is implemented in $\mathbb{F}_p$ for speed and
exactness:

1. **Choose good primes** $p_1, \ldots, p_k$ with $N \mid p_i - 1$
2. **Reduce atomic catalog** to each $\mathbb{F}_{p_i}$ (exact modular arithmetic)
3. **Build block-diagonal atomic $(S_{\text{atomic}}, T_{\text{atomic}})$**
4. **Enumerate $\mathbb{F}_{p_i}$-rational points of the block-U group**
   — typically $O(p_i^{\text{parameter\_dim}})$, small for each prime
5. **For each point, compute $S' = U^T S_{\text{atomic}} U$** and check
   Verlinde integrality by explicit computation of $N_{ij}^k$
6. **Multi-prime CRT** on resulting $(S, T)$ entries to reconstruct
   algebraic integers in $\mathbb{Z}[\zeta_N]$
7. **Cross-validate** at additional primes not used in reconstruction

**Empirical validation (2026-04-17).** For SU(2)_4 at rank 5, $N = 24$:
- 103 atomic irreps; 25,440 rank-5 strata; 1 candidate after T-spectrum filter
- 2 block-U solutions at each prime (reflection $\phi = \pi/4$ + rotation)
- 10 primes $p \in \{73, 97, 193, 241, 313, 337, 409, 433, 457, 577\}$
- CRT from 4 primes reconstructed $S' = \frac{1}{2\sqrt{3}} M$ with
  $M \in \mathbb{Z}[\sqrt{3}]^{5 \times 5}$ having entries in $\{0, \pm 1, \pm 2, \pm\sqrt{3}\}$
- Independently verified at remaining 6 primes
- Quantum dimensions $d = (\sqrt{3}, \sqrt{3}, 1, 1, 2)$, $D^2 = 12$

This matches the ground truth from NRW's `NsdGOL5.g` (Galois orbit GO[2]).

### Why not $\mathbb{F}_p$ as primary search space?

The earlier v5 strategy of searching directly in $\mathbb{F}_p$ without the
representation-theoretic structure hit obstructions:
- Cross-block $F_p$ search for S (with T fixed) is expensive without the
  atomic decomposition
- Homotopy continuation over $\mathbb{F}_p$ has no topology (no damping/line
  search), unlike over $\mathbb{C}$
- False positives at rank ≥ 4 require BNRW filtering

The v6 design uses $\mathbb{F}_p$ only for the final block-U sweep, where the
search space is both small and structured.

### Why not $\mathbb{Q}(\zeta_N)$ directly?

In principle the entire pipeline could run over $\mathbb{Q}(\zeta_N)$ (e.g. in
Oscar.jl), and this is a valid implementation option. But:
- $\mathbb{Q}(\zeta_N)$ arithmetic is $\sim 10\text{-}100\times$ slower than $\mathbb{F}_p$
- The block-U sweep needs to enumerate a variety of rational points, for which
  $\mathbb{F}_p$ gives a finite exhaustive search
- Multi-prime with CRT gives both speed and exactness

We retain $\mathbb{Z}[\zeta_N]$ (via Oscar) for the **atomic catalog** (Phase 0)
where arithmetic is done once and reduced to many primes downstream.

### Key theoretical results (informal)

1. **Parameter dimension formula:**
   $$\text{parameter\_dim} = \sum_\theta \binom{n_\theta}{2}$$
   where $n_\theta = \dim V_\theta$. This counts the continuous block-U moduli.

2. **Atomic decomposition:** Every MTC with conductor $N$ has the form
   $\rho = \bigoplus_\lambda \rho_\lambda^{m_\lambda}$ as $SL_2(\mathbb{Z}/N)$-representation,
   with basis change by a block-U element connecting the canonical
   representation-theoretic basis to the physical MTC basis.

3. **Étale vs. continuous strata:** Strata with all $n_\theta = 1$ give étale
   (discrete, rigid) MTC moduli over $\mathbb{Z}[\frac{1}{2N}]$. Strata with
   $n_\theta \geq 2$ have continuous $O(n_\theta)$ moduli, with the MTC locus
   being a 0-dimensional subvariety cut out by Verlinde integrality.

4. **Multi-prime / CRT soundness:** If $(S', T') \in \mathbb{F}_p$-$MTC$ for a
   sufficient set of primes, and CRT gives integer entries in $\mathbb{Z}[\zeta_N]$
   consistent across unused primes, the result is the unique lift with those
   reductions.

## Key references

- Bruillard, Ng, Rowell, Wang (2016): admissibility criteria (BNRW)
- Ng, Rowell, Wang, Wen (2022): SL_2(Z) representation-based classification,
  ancillary data NsdGOL
- Ng, Rowell, Wen (2023): rank ≤ 12 classification up to modular data
- Morrison-Snyder: Haagerup modular data

## Environment

- Julia 1.12+
- Oscar.jl (for $\mathbb{Z}[\zeta_N]$ arithmetic and GAP bridge)
- SL2Reps 1.1 GAP package (installed via `Oscar.GAP.Packages.install("SL2Reps")`)

## Running tests

```julia
cd /path/to/ACMG
julia --project=. test/runtests.jl
```

Expected output (as of 2026-04-17):
```
Test Summary:  | Pass  Total  Time
ACMG Prototype |  151    151  2.3s
```
