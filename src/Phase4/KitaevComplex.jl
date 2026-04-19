"""
    KitaevComplex

Davydov-Yetter / Kitaev tangent cohomology for multiplicity-free fusion
categories. Implements the cochain complex

    C⁰ →δ⁰ C¹ →δ¹ C² →δ² C³ →δ³ C⁴

with contracting homotopy χⁿ : Cⁿ → Cⁿ⁻¹ satisfying the Kitaev identity

    χⁿ⁺¹ δⁿ + δⁿ⁻¹ χⁿ = id_{Cⁿ}       (Kitaev Eq. 251; EGNO Thm 9.1.3)

# Basis of Cⁿ

For a multiplicity-free fusion category,

    Cⁿ = ⊕_{X₁,…,Xₙ} End(X₁ ⊗ … ⊗ Xₙ)

Using left-associative fusion, a basis element of `End(X₁ ⊗ … ⊗ Xₙ)` is
specified by the sequence of intermediate charges c₂ = X₁⊗X₂, c₃ = c₂⊗X₃,
…, cₙ = cₙ₋₁⊗Xₙ = Z. We represent a basis element as the full tuple

    (X₁, X₂, …, Xₙ, c₂, c₃, …, cₙ)

with Nᵢⱼᵏ = 1 constraints at every step. For n = 1 this collapses to
(X, X). For n = 0 there is a single basis element (the unit).

A cochain f ∈ Cⁿ is a scalar f[basis_tuple] on each basis element.

# Status

This module is a clean implementation designed to let the Kitaev identity
`verify_homotopy` be the sole criterion for ansatz correctness. No
pentagon / slice functionality yet — that comes after the identity checks
pass.
"""
module KitaevComplex

using LinearAlgebra
using Oscar
using TensorCategories
using ACMG: FusionRule

export C_basis, C_dim
export delta_matrix, chi_matrix
export verify_homotopy

# ============================================================
# Basis of Cⁿ
# ============================================================

"""
    C_basis(N::Array{Int,3}, r::Int, n::Int) -> Vector{Tuple}

Enumerate basis of Cⁿ. Each basis element is a tuple of length 2n:

    (X₁, X₂, …, Xₙ, c₂, c₃, …, cₙ)    (for n ≥ 2)
    (X₁, X₁)                           (for n = 1, charge = X₁)
    ((),)                              (for n = 0, the unit)

with Nᵢⱼᵏ = 1 on every fusion step c_{k+1} ∈ cₖ ⊗ Xₖ₊₁.
Multiplicity-free assumption enforced (N ∈ {0, 1}).
"""
function C_basis(N::Array{Int,3}, r::Int, n::Int)
    n == 0 && return [((),)]
    n == 1 && return [(X, X) for X in 1:r]

    # n ≥ 2: iterative build
    # start with n = 2: (X₁, X₂, c₂) with N^{X₁ X₂}_{c₂} = 1
    curr = Vector{Vector{Int}}()
    for X1 in 1:r, X2 in 1:r, c2 in 1:r
        N[X1, X2, c2] == 1 && push!(curr, [X1, X2, c2])
    end

    for step in 3:n
        next = Vector{Vector{Int}}()
        for seq in curr
            # seq layout: [X₁, …, X_{step-1}, c₂, …, c_{step-1}]
            # length = 2(step - 1)
            prev_charge = seq[end]  # c_{step-1}
            # Add new X_step and new c_step
            for Xnew in 1:r, c_new in 1:r
                if N[prev_charge, Xnew, c_new] == 1
                    new_seq = copy(seq)
                    # Insert Xnew in the middle (after position step-1),
                    # append c_new at the end.
                    insert!(new_seq, step, Xnew)     # X₁, …, X_{step-1}, X_step, c₂, …
                    push!(new_seq, c_new)             # append c_step
                    push!(next, new_seq)
                end
            end
        end
        curr = next
    end

    return [tuple(s...) for s in curr]
end

"Dimension of Cⁿ."
C_dim(N, r, n) = length(C_basis(N, r, n))

# ============================================================
# Index utility
# ============================================================

function _basis_index(basis, tuple_key)
    for (i, b) in enumerate(basis)
        b == tuple_key && return i
    end
    return 0
end

# ============================================================
# Hochschild differential δⁿ : Cⁿ → Cⁿ⁺¹
# ============================================================

"""
    delta_matrix(N, d, D2, n) -> Matrix{Float64}

Hochschild differential δⁿ : Cⁿ → Cⁿ⁺¹ as a (|Cⁿ⁺¹| × |Cⁿ|) matrix.

Formula (Kitaev Eq. 247 + Eq. 248): δⁿ = Σ_{k=0}^{n+1} (-1)^k fₖⁿ.

Face maps (as `Cⁿ → Cⁿ⁺¹`):
- f₀ⁿ: prepend a new leg X₀, charge extended by X₀ ⊗ (original Z) → Z_out.
       Contribution weight: 1 (identity extension).
- fₙ₊₁ⁿ: append a new leg Xₙ₊₁, charge extended by Z ⊗ Xₙ₊₁ → Z_out.
         Contribution weight: 1.
- fₖⁿ for k = 1..n: fuse legs Xₖ and Xₖ₊₁ via an intermediate c,
  unfolding  (…, Xₖ, Xₖ₊₁, …) ↔ (…, Xₖ', Xₖ₊₁', …) with a c-channel.
  This is the subtle case.

WARNING: only f₀ and fₙ₊₁ are implemented. Interior faces (fₖ, k = 1..n)
require careful handling of intermediate-charge insertion and are
stubbed. For identity verification with n = 0 and n = 1, only endpoints
matter, so the identity check at those levels is meaningful.
"""
function delta_matrix(N::Array{Int,3}, d::AbstractVector, D2::Real, n::Int)
    r = size(N, 1)
    C_in  = C_basis(N, r, n)
    C_out = C_basis(N, r, n + 1)
    mat = zeros(Float64, length(C_out), length(C_in))

    for (col_idx, in_tuple) in enumerate(C_in)
        # --- f₀ⁿ: prepend X₀ ---
        # For n = 0: input = ((),), extend to (X₀, X₀) — any X₀
        # For n ≥ 1: input tuple layout depends on n.
        if n == 0
            # Cᵒ has one element ((),). f₀ gives id_{X₀} for any X₀.
            for X0 in 1:r
                out_idx = _basis_index(C_out, (X0, X0))
                out_idx == 0 && continue
                mat[out_idx, col_idx] += 1.0  # (-1)^0 = +1
            end
            # fₙ₊₁ = f₁: append X₁ on the right — same effect as f₀ for n=0
            # f₁: X ↦ id_{X}, sign = (-1)^1 = -1
            for X1 in 1:r
                out_idx = _basis_index(C_out, (X1, X1))
                out_idx == 0 && continue
                mat[out_idx, col_idx] += -1.0
            end

        elseif n == 1
            # input: (X, X); output basis tuple has length 3: (X₁, X₂, c₂)
            X = in_tuple[1]
            # f₀: prepend X₀, output = (X₀, X, c₂) with c₂ ∈ X₀ ⊗ X
            for X0 in 1:r, c2 in 1:r
                N[X0, X, c2] == 1 || continue
                out_idx = _basis_index(C_out, (X0, X, c2))
                out_idx == 0 && continue
                mat[out_idx, col_idx] += 1.0  # sign +1
            end
            # f₂: append X₂, output = (X, X₂, c₂) with c₂ ∈ X ⊗ X₂
            for X2 in 1:r, c2 in 1:r
                N[X, X2, c2] == 1 || continue
                out_idx = _basis_index(C_out, (X, X2, c2))
                out_idx == 0 && continue
                mat[out_idx, col_idx] += 1.0  # sign = (-1)^(n+1) = (-1)^2 = +1 for n=1
            end
            # f₁ (interior): fuse X and ... — but input has only one X.
            # For n = 1 there's only f₀ and f₂, no interior. (n+2 = 3 faces total.)
            # Actually for Cⁿ → Cⁿ⁺¹ with n = 1, we have k = 0, 1, 2 (=n+1).
            # k = 1 is interior: but our input has one variable X, interior face
            # would "fuse X with itself" which doesn't make categorical sense
            # — in fact for n=1 the interior face is the multiplication map:
            #   f₁ f (x₀, x₁) = f(x₀ · x₁)
            # In our basis: input (X, X) represents an endo on X; f₁ gives
            # output basis tuple at (X₀, X₁, c₂) with contribution
            # depending on how the endo on X extends to endo on X₀ ⊗ X₁
            # when X₀ ⊗ X₁ → c₂ ∋ X. This requires a SUM over c₂ with
            # weights d_c₂ / (d_{X₀} d_{X₁}) (Eq. 248).
            # Sign of f₁ is (-1)^1 = -1.
            # Implementation:
            for X0 in 1:r, X1 in 1:r
                # Only contribute when X₀ ⊗ X₁ ∋ X (so N[X0, X1, X] = 1)
                N[X0, X1, X] == 1 || continue
                # Eq. 248 weight: d_X / (d_{X0} d_{X1}) ... but we need
                # the OUTPUT basis c₂ to be specified. Actually from the
                # structure of Eq. 248, (fₖ f)_{a₁...aₙ₊₁}  sums over c
                # at position k, contributing scalar f_{a₁...c...aₙ₊₁} ×
                # (d_c / (d_{aₖ} d_{aₖ₊₁})). The OUTPUT charge c₂ varies,
                # but the contribution to OUR chosen output basis requires
                # c₂ = X (the input charge being fused).
                # So for each output (X0, X1, c₂ = X):
                out_idx = _basis_index(C_out, (X0, X1, X))
                out_idx == 0 && continue
                mat[out_idx, col_idx] += -1.0 * (d[X] / (d[X0] * d[X1]))
            end

        else
            # n ≥ 2: general interior + endpoints (deferred for now)
            error("delta_matrix for n ≥ 2 not yet implemented")
        end
    end

    return mat
end

# ============================================================
# Contracting homotopy χⁿ : Cⁿ → Cⁿ⁻¹
# ============================================================

"""
    chi_matrix(N, d, D2, n; ansatz=:kitaev_first_leg) -> Matrix{Float64}

Kitaev contracting homotopy χⁿ : Cⁿ → Cⁿ⁻¹ as a (|Cⁿ⁻¹| × |Cⁿ|) matrix.

# Ansatz :kitaev_first_leg  (Kitaev Eq. 250 reading)

    (χⁿ f)_{(X₁, …, Xₙ₋₁; charges)} = (1/D²) Σ_c d_c · f_{(c, X₁, …, Xₙ₋₁; …)}

where the input-basis charges are constrained so that closing the c-leg
yields the same output-basis charges as the input.

For n = 1: output is C⁰ (1D, unit). χ¹ sends any f ∈ C¹ to a scalar:
    χ¹ f = (1/D²) Σ_X d_X · f_{(X, X)}

For n = 2: output basis (X₁, X₁). Input basis (c, X₁, ch=X₁) where
closing c gives back X₁.
"""
function chi_matrix(N::Array{Int,3}, d::AbstractVector, D2::Real, n::Int;
                    ansatz::Symbol = :kitaev_first_leg)
    n ≥ 1 || error("chi_matrix requires n ≥ 1")
    r = size(N, 1)
    C_in  = C_basis(N, r, n)
    C_out = C_basis(N, r, n - 1)
    mat = zeros(Float64, length(C_out), length(C_in))

    if ansatz === :kitaev_first_leg
        if n == 1
            # input (X, X); output ((),)
            # (χ¹ f)_unit = (1/D²) Σ_X d_X · f_{(X,X)}
            out_idx = _basis_index(C_out, ((),))
            for (col_idx, in_tuple) in enumerate(C_in)
                X = in_tuple[1]
                mat[out_idx, col_idx] += d[X] / D2
            end

        elseif n == 2
            # input: (X₁, X₂, c) — which we read as (c=X₁, X₁'=X₂, ch=c) by
            # interpreting the "first leg" to trace = X₁ — but then output
            # basis (X, X) must have X = X₂ and charge X₂; and closing X₁
            # requires the original charge c to reduce to X₂, i.e. the
            # c-loop closes to give a c ⊗ X₂ → X₂ path, which needs N[c, X₂, X₂]=1.
            #
            # Hmm — this requires deciding which leg to call "first". Two
            # natural conventions exist; pick the simplest.
            #
            # CONVENTION: "first leg" = position 1 of the tuple (X₁).
            # So input basis (X₁, X₂, c) has X₁ as the first leg to trace.
            # Trace over X₁ means: sum over X₁ with weight d_{X₁}, with the
            # closing constraint that the charge c must reduce to X₂ (the
            # remaining charge in C¹).
            #
            # Closing requires: N[X₁, X₂, c] (which is already true by C² basis)
            # AND the output basis is (X₂, X₂), so Z_out = X₂.
            # The c-loop forces Z_in = c to collapse to Z_out = X₂.
            # This is only valid if c = X₂ (closing the X₁ loop returns to X₂).
            #
            # More precisely: close X₁ means we multiply the endomorphism on
            # X₁ ⊗ X₂ by the projector onto the identity X₁-loop times d_{X₁}.
            # The result lives in End(X₂). Only the c = X₂ component of the
            # original endomorphism contributes after the loop closure.
            # (Because fusing X₁ with itself back to vacuum, composed with
            # X₁ ⊗ X₂ → c, gives a non-zero contribution only when c = X₂.)
            #
            # So: (χ² f)_{X₂, X₂} = (1/D²) Σ_{X₁ : N[X₁, X₂, X₂] = 1} d_{X₁} · f_{(X₁, X₂, X₂)}
            for (col_idx, in_tuple) in enumerate(C_in)
                X1, X2, c = in_tuple
                c == X2 || continue
                N[X1, X2, X2] == 1 || continue
                out_idx = _basis_index(C_out, (X2, X2))
                out_idx == 0 && continue
                mat[out_idx, col_idx] += d[X1] / D2
            end

        else
            error("chi_matrix ansatz $ansatz not implemented for n ≥ 3 yet")
        end

    else
        error("Unknown ansatz: $ansatz")
    end

    return mat
end

# ============================================================
# Identity check: χⁿ⁺¹ δⁿ + δⁿ⁻¹ χⁿ = id_Cⁿ
# ============================================================

"""
    verify_homotopy(N, d, D2, n; ansatz=:kitaev_first_leg, atol=1e-10)
        -> (err::Float64, M::Matrix{Float64})

Compute M = χⁿ⁺¹ δⁿ + δⁿ⁻¹ χⁿ  and report max|M - I|.

For n = 0: only χ¹ δ⁰ contributes (χⁿ for n < 1 undefined). Expect
`M = id_{C⁰}` (a 1×1 identity).

For n = 1: both terms contribute. Expect `M = id_{C¹}`.
"""
function verify_homotopy(N::Array{Int,3}, d::AbstractVector, D2::Real, n::Int;
                         ansatz::Symbol = :kitaev_first_leg, atol::Real = 1e-10)
    r = size(N, 1)
    M = zeros(Float64, C_dim(N, r, n), C_dim(N, r, n))

    if n ≥ 0
        δn  = delta_matrix(N, d, D2, n)
        χn1 = chi_matrix(N, d, D2, n + 1; ansatz = ansatz)
        M  += χn1 * δn
    end

    if n ≥ 1
        χn  = chi_matrix(N, d, D2, n; ansatz = ansatz)
        δn1 = delta_matrix(N, d, D2, n - 1)
        M  += δn1 * χn
    end

    I_n = Matrix{Float64}(I, size(M)...)
    err = maximum(abs.(M - I_n))
    return err, M
end

end # module KitaevComplex
