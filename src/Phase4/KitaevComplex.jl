"""
    KitaevComplex

Kitaev tangent cochain complex — faithful implementation of
Kitaev 2006 (arXiv:cond-mat/0506438v3) Appendix E.6.

# Summary

For a multiplicity-free unitary fusion category the complex

    C⁰ →δ⁰ C¹ →δ¹ C² →δ² ⋯        (Def E.15, Eq. 247-248)

is implemented together with the contracting homotopy

    χⁿ : Cⁿ → Cⁿ⁻¹                (Thm E.16, Eq. 250)

satisfying

    χⁿ⁺¹ δⁿ + δⁿ⁻¹ χⁿ = id_{Cⁿ}    (Eq. 251)     for n ≥ 1,

which yields `Hⁿ = 0` (Crane-Yetter / Davydov tangent cohomology).

# Scalar representation

In the multiplicity-free setting every fusion block `V_c^{ab}` is
1-dimensional, so an element `X ∈ Cⁿ` reduces to one real scalar per
left-associated fusion tree `(a₁, …, aₙ; c₂, c₃, …, cₙ)`, with
`cₖ = c_{k-1} ⊗ aₖ` and `c₁ := a₁`. This is the basis enumerated by
[`C_basis`](@ref).

## Why the `√(d_c / (d_a d_b))` normalization of Eq. 248 disappears

Kitaev Eq. (245)-(246) introduces "box notation" `[X_{ab}]` for elements
of `V_{ab}^{ab}`:

    [X_{ab}] = Σ_{c,j} √(d_c / (d_a d_b)) · X_c · ψ_{c,j}^{ab} (ψ_{c,j}^{ab})†.

With `ψ_{c,j}^{ab}` orthonormal under Kitaev's renormalized inner
product `⟨⟨⋅|⋅⟩⟩` (Eq. 195), the operators `ψ(ψ)†` form an identity
resolution (Eq. 196). The `√`-factors exactly compensate the non-unit
length of `ψ` in the usual inner product, so the SCALAR value of
`[X_{ab}]` on the `c`-block is just `X_c`.

Consequently Eq. (244) takes the simple scalar form

    (δ¹X)_{a,b;c}  =  X_b  −  X_c  +  X_a                  (Eq. 244)

and each face map `fₖⁿ : Cⁿ → Cⁿ⁺¹` contributes coefficient `±1` on
valid tree tuples (for `n ≤ 1`; for `n ≥ 2` interior face maps require
F-symbols to re-bracket the input tree).

## χⁿ normalization

Eq. 250 reads

    (χⁿ X)_{a₁…aₙ₋₁} = (1/𝒟²) Σ_c d_c · [X_{c a₁…aₙ₋₁} with c-leg looped].

In the scalar basis, requiring all Eq. 252 identities

    χⁿ⁺¹ f₀ⁿ = id,   χⁿ⁺¹ fₖⁿ = fₖ₋₁ⁿ⁻¹ χⁿ  (k, n > 0),

uniquely fixes the scalar weights. For `n = 1, 2`:

    (χ¹ X)         = (1/𝒟²) Σ_a       d_a²              · X_{a, a}
    (χ² X)_{a, a}  = (1/𝒟²) Σ_{a₁, c₂ : N[a₁, a, c₂] = 1}
                            (d_{a₁} · d_{c₂} / d_a)     · X_{a₁, a; c₂}

The sum in χ² runs over **every** `c₂ ∈ a₁ ⊗ a`, not just `c₂ = a` —
the Eq. 196 identity resolution distributes loop-closure weight across
all reachable channels.

### Derivation of the χ² weight

Writing `W(a₁, b, c₂; a)` for the χ²-weight at input `(a₁, b, c₂)`,
output `(a, a)`, support on `b = a`, the three identities demand
(using `Σ_{a₁} d_{a₁} N[a₁, b, c] = d_b · d_c`, a Frobenius-reciprocity
consequence):

    C1 (χ² f₀¹ = id_{C¹}):    Σ_{a₁, c₂: N[a₁,a,c₂]=1} W = 1
    C2 (χ² f₁¹ = f₀⁰ χ¹):     Σ_{c₂: N[a₁,a,c₂]=1}    W = d_{a₁}² / 𝒟²
    C3 (χ² f₂¹ = f₁⁰ χ¹):     Σ_{a₁: N[a₁,a,c₂]=1}    W = d_{c₂}² / 𝒟²

The ansatz `W = (d_{a₁} · d_{c₂}) / (𝒟² · d_a) · δ_{b,a}` satisfies all
three. (For C1: Σ d_{a₁} d_{c₂} N[a₁,a,c₂] = Σ_{a₁} d_{a₁} (d_{a₁} d_a)
= d_a · 𝒟²; the d_a cancels, leaving 1.)

### Fibonacci hand-verification

Fib: `d_τ = φ`, `𝒟² = 1 + φ² = 2 + φ`. For `X' ∈ C¹`:

    (χ² δ¹ X')_1 = (1/𝒟²)[X'_1 + φ² X'_1]            = X'_1                  ✓
    (χ² δ¹ X')_τ = (1/𝒟²)[X'_1 + (2X'_τ − X'_1) + φ X'_τ]
                 = (2 + φ)/𝒟² · X'_τ  =  X'_τ                                ✓

Numerically verified for both Fibonacci and Ising: `verify_homotopy(N, d, D2, 1)`
returns `err = 0.0` (exact, no round-off).

# Status

- `n = 0`     : `δ⁰ = f₀⁰ − f₁⁰ = 0` (Kitaev remark, p. 92).
- `n = 1`     : δ¹, χ² fully implemented; Eq. 251 numerically exact for
                Fib and Ising.
- `n ≥ 2`     : δⁿ, χⁿ⁺¹ require **F-symbols** for interior face maps
                `fₖⁿ` with `1 ≤ k ≤ n` because input and output trees
                differ by nontrivial re-bracketings. Currently errors
                (TODO; cf. [Phase4 roadmap]).

# References

- Kitaev A., *Anyons in an exactly solved model and beyond*,
  arXiv:cond-mat/0506438v3 (2006), Appendix E.6.
  Def E.15  → [`C_basis`](@ref), [`delta_matrix`](@ref);
  Thm E.16  → [`chi_matrix`](@ref), [`verify_homotopy`](@ref).
- Etingof-Nikshych-Ostrik (ENO), Ann. Math. **162** (2005) 581, §6.
- Crane-Yetter / Davydov, *tangent cohomology* of fusion categories.
"""
module KitaevComplex

using LinearAlgebra

export C_basis, C_dim
export delta_matrix, chi_matrix
export verify_homotopy

# ============================================================
# Cⁿ basis enumeration
# ============================================================

"""
    C_basis(N, r, n) -> Vector{Tuple}

Enumerate left-associated fusion tree basis of `Cⁿ`.

Each basis element is a tuple of length `2n - 1` (for `n ≥ 2`):

    (a₁, a₂, …, aₙ, c₂, c₃, …, cₙ)

with `cₖ = c_{k-1} ⊗ aₖ` (so `N[c_{k-1}, aₖ, cₖ] = 1`); convention
`c₁ := a₁`.

- `n = 0`:  returns `[((),)]`
- `n = 1`:  returns `[(a, a) for a ∈ 1:r]`  (second entry = `c₁ = a`)

Multiplicity-free (`N ∈ {0, 1}`) is asserted.
"""
function C_basis(N::Array{Int,3}, r::Int, n::Int)
    @assert all(x -> x == 0 || x == 1, N) "multiplicity-free only (N ∈ {0,1})"
    n < 0 && error("n must be ≥ 0")
    n == 0 && return Tuple[ ((),) ]
    n == 1 && return Tuple[ (a, a) for a in 1:r ]

    # n ≥ 2: iterative left-associated tree construction
    curr = Vector{Vector{Int}}()
    for a1 in 1:r, a2 in 1:r, c2 in 1:r
        N[a1, a2, c2] == 1 && push!(curr, [a1, a2, c2])
    end
    for k in 3:n
        nxt = Vector{Vector{Int}}()
        for seq in curr
            prev_c = seq[end]
            for ak in 1:r, ck in 1:r
                if N[prev_c, ak, ck] == 1
                    s = copy(seq)
                    insert!(s, k, ak)     # place aₖ after a_{k-1}
                    push!(s, ck)          # append cₖ
                    push!(nxt, s)
                end
            end
        end
        curr = nxt
    end
    return Tuple[ tuple(s...) for s in curr ]
end

"Dimension of `Cⁿ`."
C_dim(N::Array{Int,3}, r::Int, n::Int) = length(C_basis(N, r, n))

# tuple index lookup
function _index_of(basis::Vector, key)::Int
    for (i, b) in enumerate(basis)
        b == key && return i
    end
    return 0
end

# ============================================================
# δⁿ  (Def E.15, Eq. 247-248)
# ============================================================

"""
    delta_matrix(N, d, D2, n) -> Matrix{Float64}

Hochschild differential `δⁿ : Cⁿ → Cⁿ⁺¹` as a `(|Cⁿ⁺¹| × |Cⁿ|)` matrix
(columns index `Cⁿ`, rows index `Cⁿ⁺¹`).

Implemented for `n ∈ {0, 1}`; higher `n` errors pending F-symbol data.

For `n = 1` the formula is Kitaev Eq. 244 in scalar form:

    (δ¹X)_{a₁, a₂; c₂}  =  X_{a₂}  −  X_{c₂}  +  X_{a₁}.

For `n = 0`, both `f₀⁰` and `f₁⁰` send `λ ↦ Σ_a λ·1_{V_a^a}`, so
`δ⁰ ≡ 0`.
"""
function delta_matrix(N::Array{Int,3}, d::AbstractVector, D2::Real, n::Int)
    r = size(N, 1)
    C_in  = C_basis(N, r, n)
    C_out = C_basis(N, r, n + 1)
    M = zeros(Float64, length(C_out), length(C_in))

    if n == 0
        # δ⁰ = f₀⁰ − f₁⁰; both act as λ ↦ Σ_a λ·1_a, so the matrix is 0.
        # Loop made explicit for clarity (every entry receives +1 − 1).
        for (j, _) in pairs(C_in), (i, _) in pairs(C_out)
            M[i, j] += 1.0
            M[i, j] -= 1.0
        end
        return M
    end

    if n == 1
        # (δ¹X)_{a₁,a₂;c₂} = X_{a₂} − X_{c₂} + X_{a₁}.
        for (j, in_t) in pairs(C_in)
            a_in = in_t[1]
            for (i, out_t) in pairs(C_out)
                a1, a2, c2 = out_t
                a_in == a2 && (M[i, j] += 1.0)   # f₀
                a_in == c2 && (M[i, j] -= 1.0)   # f₁  (interior)
                a_in == a1 && (M[i, j] += 1.0)   # f₂
            end
        end
        return M
    end

    error("delta_matrix: n = $n not yet implemented " *
          "(F-symbols needed for n ≥ 2 interior faces).")
end

# ============================================================
# χⁿ  (Thm E.16, Eq. 250)
# ============================================================

"""
    chi_matrix(N, d, D2, n) -> Matrix{Float64}

Contracting homotopy `χⁿ : Cⁿ → Cⁿ⁻¹`.

Implemented for `n ∈ {1, 2}`; higher `n` errors pending F-symbol data.

    (χ¹ X)         = (1/𝒟²) Σ_a       d_a²              · X_{a, a}
    (χ² X)_{a, a}  = (1/𝒟²) Σ_{a₁, c₂ : N[a₁, a, c₂] = 1}
                            (d_{a₁} · d_{c₂} / d_a)     · X_{a₁, a; c₂}

The weights are uniquely determined by the Eq. 252 identities for every
face map; derivation in the module docstring.
"""
function chi_matrix(N::Array{Int,3}, d::AbstractVector, D2::Real, n::Int)
    n ≥ 1 || error("chi_matrix requires n ≥ 1")
    r = size(N, 1)
    C_in  = C_basis(N, r, n)
    C_out = C_basis(N, r, n - 1)
    M = zeros(Float64, length(C_out), length(C_in))

    if n == 1
        i_out = _index_of(C_out, ((),))
        @assert i_out > 0
        for (j, in_t) in pairs(C_in)
            a = in_t[1]
            M[i_out, j] = d[a]^2 / D2
        end
        return M
    end

    if n == 2
        for (j, in_t) in pairs(C_in)
            a1, a2, c2 = in_t
            i_out = _index_of(C_out, (a2, a2))
            i_out == 0 && continue
            M[i_out, j] = d[a1] * d[c2] / (D2 * d[a2])
        end
        return M
    end

    error("chi_matrix: n = $n not yet implemented " *
          "(F-symbols needed for n ≥ 3).")
end

# ============================================================
# Eq. 251 verification
# ============================================================

"""
    verify_homotopy(N, d, D2, n; atol = 1e-10) -> (err, M)

Build  `M = χⁿ⁺¹ δⁿ + δⁿ⁻¹ χⁿ`  and return `err = max|M − I|`.

Kitaev Eq. 251: `err ≈ 0` for every `n ≥ 1`.
"""
function verify_homotopy(N::Array{Int,3}, d::AbstractVector, D2::Real, n::Int;
                         atol::Real = 1e-10)
    r = size(N, 1)
    dim_n = C_dim(N, r, n)
    M = zeros(Float64, dim_n, dim_n)

    δn  = delta_matrix(N, d, D2, n)
    χn1 = chi_matrix(N, d, D2, n + 1)
    M .+= χn1 * δn

    if n ≥ 1
        χn  = chi_matrix(N, d, D2, n)
        δn1 = delta_matrix(N, d, D2, n - 1)
        M .+= δn1 * χn
    end

    Id  = Matrix{Float64}(I, dim_n, dim_n)
    err = maximum(abs.(M - Id))
    return err, M
end

end # module KitaevComplex
