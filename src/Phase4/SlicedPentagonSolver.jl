"""
    SlicedPentagonSolver

Combine pentagon equations with Kitaev `χ³ F = 0` gauge-slice constraints.

Both the pentagon equations (from TensorCategories via
`get_pentagon_system`) and the slice constraints (from
`KitaevComplex.chi3_constraints`) live in the same Oscar polynomial ring
and are appended as a single augmented system. The solver uses
damped-Newton with KrylovKit.

# Why this redesign

The previous iteration (F-coordinate `Δ_gauge` / `GaugeAnalysis`)
produced 'gauge' directions that violated pentagon (direct pentagon
evaluation gave O(0.4) residuals at 'converged' solutions). That
approach has been deleted. The Kitaev slice `χ³ F = 0` is the canonical
linear constraint in pentagon-variable space, built from fusion data
alone.

# Main API

- `dims_D2(Nijk, r)`                → (d, D²)
- `get_sliced_pentagon_system(Nijk, r)`
      → (R, pent_eqs, slice_polys, n)
- `solve_pentagon_newton_with_slice(Nijk, r, base_F_func; …)`
      → Vector{Vector{ComplexF64}}
"""
module SlicedPentagonSolver

using LinearAlgebra
using SparseArrays
using KrylovKit
using Oscar
using ACMG: FusionRule
using ..KitaevComplex
using ..PentagonEquations
using ..PentagonSolver

export dims_D2
export get_sliced_pentagon_system
export solve_pentagon_newton_with_slice

# ============================================================
# Quantum dimensions
# ============================================================

"""
    dims_D2(Nijk, r) -> (d::Vector{Float64}, D2::Float64)

Compute quantum dimensions by power-iterating the fusion matrix of
object 2 (or object 1 if the category has rank 1, trivially). Normalised
so `d[1] = 1`. Returns also `D² = Σ d_a²`.
"""
function dims_D2(Nijk::Array{Int,3}, r::Int)
    a = r ≥ 2 ? 2 : 1
    Na = zeros(Float64, r, r)
    for b in 1:r, c in 1:r
        Na[b, c] = Nijk[a, b, c]
    end

    v = ones(Float64, r)
    for _ in 1:500
        v = Na * v
        nv = norm(v)
        nv > 0 && (v ./= nv)
    end

    d1 = v[1]
    abs(d1) < 1e-14 && error("d[1] is zero; fusion rule may be degenerate")
    d = v ./ d1
    d[1] = 1.0
    D2 = sum(d .^ 2)
    return d, D2
end

# ============================================================
# Sliced system
# ============================================================

"""
    get_sliced_pentagon_system(Nijk, r)
        -> (R, pent_eqs, slice_polys, n)

Pentagon + Kitaev χ³ slice in the same Oscar ring.
"""
function get_sliced_pentagon_system(Nijk::Array{Int,3}, r::Int)
    R, pent_eqs, n = PentagonEquations.get_pentagon_system(Nijk, r)

    one_vec = zeros(Int, r)
    one_vec[1] = 1
    fkey_map = KitaevComplex.fkey_to_xvar_map(Nijk, r, one_vec)
    @assert length(fkey_map) == n "FKey map size $(length(fkey_map)) ≠ pentagon var count $n"

    d, D2 = dims_D2(Nijk, r)
    slice_polys = KitaevComplex.chi3_constraints(Nijk, r, d, D2, R, fkey_map)

    return R, pent_eqs, slice_polys, n
end

# ============================================================
# Newton with slice
# ============================================================

"""
    solve_pentagon_newton_with_slice(Nijk, r, base_F_func;
                                     initial_points=nothing,
                                     max_trials=5, max_iter=200,
                                     tol=1e-12, perturb_scale=0.05,
                                     verbose=false)
        -> Vector{Vector{ComplexF64}}
"""
function solve_pentagon_newton_with_slice(Nijk::Array{Int,3}, r::Int,
        base_F_func::Function;
        initial_points::Union{Nothing, Vector{Vector{ComplexF64}}} = nothing,
        max_trials::Int = 5,
        max_iter::Int = 200,
        tol::Real = 1e-12,
        perturb_scale::Real = 0.05,
        verbose::Bool = false)

    R, pent_eqs, slice_polys, n = get_sliced_pentagon_system(Nijk, r)
    all_eqs = vcat(pent_eqs, slice_polys)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in all_eqs]

    one_vec = zeros(Int, r)
    one_vec[1] = 1
    fkey_map = KitaevComplex.fkey_to_xvar_map(Nijk, r, one_vec)

    F_base_vec = zeros(ComplexF64, n)
    for (key, pidx) in fkey_map
        F_base_vec[pidx] = ComplexF64(base_F_func(key...))
    end

    if initial_points === nothing
        initial_points = Vector{Vector{ComplexF64}}()
        push!(initial_points, copy(F_base_vec))
        for _ in 2:max_trials
            push!(initial_points,
                  F_base_vec + perturb_scale * randn(ComplexF64, n))
        end
    end

    solutions = Vector{Vector{ComplexF64}}()
    for (trial, x0) in enumerate(initial_points)
        x = copy(x0)
        for iter in 1:max_iter
            F_val = ComplexF64[PentagonSolver.eval_poly_complex(eq, x) for eq in all_eqs]
            res = maximum(abs.(F_val))
            if res < tol
                is_dup = any(norm(x - s) < 1e-8 for s in solutions)
                is_dup || push!(solutions, copy(x))
                verbose && println("  trial $trial: converged at iter $iter, res=$res")
                break
            end

            J = PentagonSolver.sparse_jacobian(all_eqs, derivs, x, n)
            delta, _ = linsolve(v -> J' * (J * v), J' * F_val;
                                ishermitian = true, isposdef = true, verbosity = 0)

            alpha = 1.0
            for _ in 1:20
                x_new = x - alpha * delta
                F_new = ComplexF64[PentagonSolver.eval_poly_complex(eq, x_new) for eq in all_eqs]
                if maximum(abs.(F_new)) < res
                    x = x_new
                    break
                end
                alpha *= 0.5
            end
        end
    end

    verbose && println("  Newton+slice: $(length(solutions)) unique solutions across $(length(initial_points)) trials")
    return solutions
end

end # module SlicedPentagonSolver
