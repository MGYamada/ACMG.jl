"""
    KitaevComplex

Kitaev 2006 App E.6 tangent cohomology, implemented in the
**pentagon-variable coordinate system** of TensorCategories.jl.

# Design

The previous version of this module attempted an F-coordinate
(FKey-indexed) gauge analysis and produced spurious gauge directions that
violated pentagon (see diagnostic scripts). That approach has been
entirely discarded.

This version follows Kitaev's definition directly: the contracting
homotopy `χ^n : C^n → C^{n-1}` is a canonical linear operator defined
purely from fusion data (d_a, D², N^{ab}_c). The slice constraint for
pentagon gauge fixing is simply

    χ³ F = 0

expressed as a set of **linear equations in the pentagon variables**,
living alongside the pentagon polynomials themselves.

# Status

First iteration. Current responsibilities:
  - `C2_basis(Nijk, r)`          : basis of C² as (a, b, y) triples
  - `fkey_to_xvar_map(...)`      : TensorCategories-compatible F-key
                                   → pentagon variable index
  - `chi3_matrix(...)`           : χ³ as a numerical (|C²| × nvars) matrix
  - `chi3_constraints(...)`      : χ³ F = 0 as Oscar linear polynomials

The `χ³` convention used here is an ANSATZ that must be verified on a
known solution (e.g. Fibonacci's exact F-symbol) before use. The
verification is:
  1.  `chi3_matrix * F_base_vector` should be zero (Kitaev slice
      contains the canonical F as a point).
  2.  Long-term: χδ + δχ = id as an identity on C^n, which confirms
      the chain-level Hodge decomposition and hence Ocneanu rigidity.

If (1) fails for Fibonacci, the convention in `chi3_matrix` is wrong and
should be refined (candidate fixes: add d_f / d_y factors, sum over e
more broadly, include multiplicity weights, etc.).
"""
module KitaevComplex

using LinearAlgebra
using Oscar
using TensorCategories
using ACMG: FusionRule

export C2_basis, fkey_to_xvar_map
export chi3_matrix, chi3_constraints

# ============================================================================
# Basis of C² (multiplicity-free)
# ============================================================================

"""
    C2_basis(Nijk::Array{Int,3}, r::Int) -> Vector{NTuple{3,Int}}

Basis of C² = ⊕_{a₁, a₂} V^{a₁ a₂}_{a₁ a₂} as a list of (a₁, a₂, y)
triples with N^{a₁ a₂}_y = 1.  Assumes multiplicity-free fusion at the
two-leg level.
"""
function C2_basis(Nijk::Array{Int,3}, r::Int)
    basis = NTuple{3,Int}[]
    for a in 1:r, b in 1:r, y in 1:r
        if Nijk[a, b, y] == 1
            push!(basis, (a, b, y))
        end
    end
    return basis
end

# ============================================================================
# F-key ↔ pentagon variable map (mirrors TensorCategories)
# ============================================================================

"""
    fkey_to_xvar_map(Nijk, r, one_vec) -> Dict{NTuple{6,Int}, Int}

Map (i, j, k, o, e, f) → pentagon variable index `p` such that x_p is the
associator entry F^{ijk}_{o; e, f}. Matches the traversal and pop! order
of `assign_F_to_associator!` in `HexagonEquations.jl`.
"""
function fkey_to_xvar_map(Nijk::Array{Int,3}, r::Int, one_vec::Vector{Int})
    m = r
    flat = NTuple{6,Int}[]
    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(one_vec[[i, j, k]]) > 0 && continue
        rows_e = Int[]
        for e in 1:m
            if Nijk[i, j, e] ≥ 1 && Nijk[e, k, o] ≥ 1
                push!(rows_e, e)
            end
        end
        cols_f = Int[]
        for f in 1:m
            if Nijk[j, k, f] ≥ 1 && Nijk[i, f, o] ≥ 1
                push!(cols_f, f)
            end
        end
        (isempty(rows_e) || isempty(cols_f)) && continue
        for e in rows_e, f in cols_f
            push!(flat, (i, j, k, o, e, f))
        end
    end
    n = length(flat)
    return Dict{NTuple{6,Int}, Int}(f => n - p + 1 for (p, f) in enumerate(flat))
end

# ============================================================================
# χ³ : C³ → C²  (first ansatz — to be verified on Fibonacci)
# ============================================================================
#
# Kitaev Eq. 250 for n = 3:
#
#     (χ³ X)^{a b} = (1/D²) Σ_c d_c · X^{a b c}_{[c-trace]}
#
# In the multiplicity-free associator basis, X^{abc}_{z; e, f} with rows
# e ∈ a⊗b and cols f ∈ b⊗c (total charge z).
#
# Ansatz for "c-trace": evaluate at e = y (where y is the target C² charge),
# z = y (c-leg closes the overall charge back to y), and sum over all legal
# f's with weight d_c / D².
#
#     (χ³ X)^{ab}_y = (1/D²) Σ_{c, f : legal} d_c · X^{abc}_{y; y, f}
#
# legal = N^{ab}_y · N^{yc}_y · N^{bc}_f · N^{af}_y ≠ 0.
# ============================================================================

"""
    chi3_matrix(Nijk, r, d, D2, fkey_map)
        -> (Matrix{Float64}, Vector{NTuple{3,Int}})

Build χ³ as `(|C²|_effective × nvars)` matrix, plus the list of C² basis
elements (a, b, y) giving the row ordering. Only rows with at least one
non-zero entry are kept.

`d` is the quantum-dimension vector (length `r`); `D2 = Σ d_i²`.
`fkey_map` is from `fkey_to_xvar_map`.
"""
function chi3_matrix(Nijk::Array{Int,3}, r::Int,
                     d::AbstractVector, D2::Real,
                     fkey_map::Dict{NTuple{6,Int}, Int})
    basis_full = C2_basis(Nijk, r)
    nvars = length(fkey_map)

    rows = Vector{Vector{Float64}}()
    kept_basis = NTuple{3,Int}[]

    for (a, b, y) in basis_full
        row = zeros(Float64, nvars)
        for c in 1:r
            # Ansatz: F^{abc}_{y; y, f} with c-trace at (e=y, o=y)
            Nijk[a, b, y] ≥ 1 || continue
            Nijk[y, c, y] ≥ 1 || continue
            for f in 1:r
                Nijk[b, c, f] ≥ 1 || continue
                Nijk[a, f, y] ≥ 1 || continue
                key = (a, b, c, y, y, f)
                haskey(fkey_map, key) || continue
                p = fkey_map[key]
                row[p] += d[c] / D2
            end
        end
        if any(abs.(row) .> 1e-14)
            push!(rows, row)
            push!(kept_basis, (a, b, y))
        end
    end

    M = isempty(rows) ? zeros(Float64, 0, nvars) :
                        reduce(vcat, (reshape(r_, 1, :) for r_ in rows))
    return M, kept_basis
end

# ============================================================================
# χ³ as Oscar polynomial constraints
# ============================================================================

"""
    chi3_constraints(Nijk, r, d, D2, R, fkey_map) -> Vector{QQMPolyRingElem}

χ³ F = 0 as linear Oscar polynomials in the pentagon-variable ring `R`.
Coefficients are rationalised from floating-point values of `d_c / D²`.

If `d_c` or `D²` are irrational (e.g. Fibonacci's φ, 1 + φ²), the
rationalisation introduces a rational approximation. For algebraic
precision one should run the same construction in `Q(ζ_N)` or similar,
but for numerical Newton / HC solving at `ComplexF64` precision the
approximation is harmless.
"""
function chi3_constraints(Nijk::Array{Int,3}, r::Int,
                          d::AbstractVector, D2::Real,
                          R::MPolyRing,
                          fkey_map::Dict{NTuple{6,Int}, Int})
    M, _ = chi3_matrix(Nijk, r, d, D2, fkey_map)
    xs = gens(R)

    polys = QQMPolyRingElem[]
    for row_idx in axes(M, 1)
        row = M[row_idx, :]
        any(abs.(row) .> 1e-14) || continue
        poly = zero(R)
        for (p, c) in enumerate(row)
            abs(c) < 1e-14 && continue
            q = rationalize(BigInt, c; tol = 1e-12)
            poly += q * xs[p]
        end
        push!(polys, poly)
    end
    return polys
end

end # module KitaevComplex
