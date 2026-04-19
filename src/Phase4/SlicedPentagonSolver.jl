"""
    SlicedPentagonSolver

Wire Kitaev/Hodge gauge slice into the pentagon HomotopyContinuation solver.

Given a fusion rule and a prior F-symbol solution (used as a base point
for linearization), builds:
  1. The KitaevComplex F-coordinate gauge analysis,
  2. A mapping from KitaevComplex FKeys to TensorCategories pentagon
     variables x_1, ..., x_n,
  3. Linear slice constraints in the pentagon variables,
  4. An augmented HomotopyContinuation `System` with Oscar pentagon
     equations + the slice constraints.

The slice kills the continuous gauge orbit, reducing the HC mixed
volume correspondingly.

# Note on coefficient rings

The pentagon equations from `get_pentagon_system` live in an Oscar
`QQMPolyRing` (rational coefficients). The Kitaev slice constraints have
coefficients like 1/√2 or 1/√φ — irrational in general. We therefore
bypass the Oscar polynomial ring for the slice and construct the
augmented system directly at the HomotopyContinuation level, where
arbitrary real/complex coefficients are natively supported.
"""
module SlicedPentagonSolver

using LinearAlgebra
using SparseArrays
using KrylovKit
using Oscar
using TensorCategories
using ACMG: FusionRule
using ..KitaevComplex
using ..PentagonEquations
using ..PentagonSolver
import HomotopyContinuation
const HC = HomotopyContinuation

export build_fkey_to_xvar_map
export slice_constraints_as_hc
export build_sliced_hc_system
export get_sliced_pentagon_system_hc
export build_slice_linear_system
export solve_pentagon_newton_with_slice
export solve_pentagon_newton_with_slice

# ============================================================
# F-key ↔ pentagon variable map
# ============================================================

"""
    build_fkey_to_xvar_map(Nijk, r, one_vec) -> Dict{NTuple{6,Int}, Int}

Build a Dict mapping F-symbol key (i,j,k,o,e,f) (1-based) to pentagon
variable index p such that x_p corresponds to F^{ijk}_{o; e, f}.

Mirrors the traversal of `assign_F_to_associator!` in HexagonEquations.jl:
  - Nested loops (i,j,k,o) ∈ 1:m^4, skipping blocks with any of i,j,k = unit.
  - Each block has row basis e (with N^{ij}_e ≥ 1, N^{ek}_o ≥ 1)
    and column basis f (with N^{jk}_f ≥ 1, N^{if}_o ≥ 1).
  - Entries iterated in row-major order, pushed onto a stack.
  - Pentagon variables are assigned via pop! from the stack's end.

Hence the first (i,j,k,o,e,f) encountered (row-major, first non-trivial
block) gets pentagon variable index n, the second gets n-1, and so on.
"""
function build_fkey_to_xvar_map(Nijk::Array{Int,3}, r::Int, one_vec::Vector{Int})
    m = r
    flat_fkeys = NTuple{6,Int}[]
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
        isempty(rows_e) && continue
        isempty(cols_f) && continue
        for e in rows_e, f in cols_f
            push!(flat_fkeys, (i, j, k, o, e, f))
        end
    end

    n = length(flat_fkeys)
    fkey_to_xvar = Dict{NTuple{6,Int}, Int}()
    for (pos, fkey) in enumerate(flat_fkeys)
        fkey_to_xvar[fkey] = n - pos + 1
    end
    return fkey_to_xvar
end

# ============================================================
# Slice constraints as HC expressions
# ============================================================

"""
    slice_constraints_as_hc(ga, fkey_to_xvar, hc_vars) -> Vector{HC.Expression}

Translate each gauge direction (column of the unit-preserving effective
gauge image) into a linear HC.Expression in pentagon variables.

Returns `gauge_orbit_dim(ga)` linear expressions. Each expression
vanishes at the base F-solution by construction, so the base F remains
a solution of the augmented system.

The expression for gauge direction `v ∈ R^nF` is:

    sum_k v[k] * (x_{p_k} - F_base_at_key_k)

restricted to F-entries that correspond to pentagon variables (the
unit-axiom-fixed entries don't appear as variables in TensorCategories'
pentagon system).
"""
function slice_constraints_as_hc(ga::KitaevComplex.GaugeAnalysis,
                                 fkey_to_xvar::Dict,
                                 hc_vars::Vector{HC.Variable})
    fcs = ga.fcs
    # Orthonormal basis of im Δ_gauge_eff
    Dg_eff = ga.Delta_gauge_eff
    F = svd(Dg_eff; full = true)
    rank_g = ga.gauge_orbit_dim
    gauge_basis = F.U[:, 1:rank_g]   # nF × rank_g

    exprs = HC.Expression[]
    for col in 1:rank_g
        v = gauge_basis[:, col]
        expr = HC.Expression(0)
        for (k, key) in enumerate(fcs.vars)
            coeff = v[k]
            abs(coeff) < 1e-12 && continue
            if haskey(fkey_to_xvar, key)
                pidx = fkey_to_xvar[key]
                Fval = KitaevComplex.F_value(fcs, key)
                # Linear constraint vanishing at base F:
                expr += coeff * (hc_vars[pidx] - Fval)
            end
            # Unit-axiom-fixed entries: no variable, but v[k] contributes
            # v[k] * (F_const - F_const) = 0, so they can be skipped safely.
        end
        push!(exprs, expr)
    end
    return exprs
end

# ============================================================
# Augmented HC system
# ============================================================

"""
    build_sliced_hc_system(Nijk, r, base_F_func)
        -> (sys::HC.System, hc_vars::Vector{HC.Variable},
            n::Int, n_slice::Int, ga::GaugeAnalysis)

Build the augmented HomotopyContinuation System:
  - Pentagon equations from TensorCategories, converted via
    `oscar_poly_to_hc`,
  - Kitaev gauge slice constraints appended as HC expressions.

Returns the combined system, the HC variables, the number of pentagon
variables `n`, the number of slice constraints added `n_slice`, and
the GaugeAnalysis object (for later inspection).

Note: the resulting system may have `n + n_slice > n_eqs_original`.
HomotopyContinuation handles over-determined systems, but the useful
regime is when slice constraints are exactly gauge-orbit-dimension many
(one for Fibonacci, one for Ising) so the augmented system is
well-balanced with respect to the reduced mixed volume.
"""
function build_sliced_hc_system(Nijk::Array{Int,3}, r::Int,
                                base_F_func::Function)
    # 1. Oscar pentagon system
    R, eqs, n = PentagonEquations.get_pentagon_system(Nijk, r)

    # 2. KitaevComplex analysis
    fr = FusionRule(Nijk)
    fcs = KitaevComplex.build_F_coord_space(fr, base_F_func)
    ga  = KitaevComplex.analyze_gauge(fcs)

    # 3. Build HC variables and convert Oscar eqs
    hc_vars = [HC.Variable(Symbol("x", i)) for i in 1:n]
    hc_exprs = HC.Expression[]
    for eq in eqs
        (isa(eq, Integer) || iszero(eq)) && continue
        push!(hc_exprs, PentagonSolver.oscar_poly_to_hc(eq, hc_vars))
    end

    # 4. FKey ↔ pentagon variable mapping
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    fkey_to_xvar = build_fkey_to_xvar_map(Nijk, r, one_vec)

    # 5. Slice constraints as HC expressions
    slice_exprs = slice_constraints_as_hc(ga, fkey_to_xvar, hc_vars)
    for e in slice_exprs
        push!(hc_exprs, e)
    end

    sys = HC.System(hc_exprs; variables = hc_vars)
    return sys, hc_vars, n, length(slice_exprs), ga
end

"""
    get_sliced_pentagon_system_hc(Nijk, r, base_F_func) -> NamedTuple

Convenience wrapper returning a NamedTuple with fields
  (sys, hc_vars, n, n_slice, ga, fkey_to_xvar)
for downstream HC solving and inspection.
"""
function get_sliced_pentagon_system_hc(Nijk::Array{Int,3}, r::Int,
                                        base_F_func::Function)
    sys, hc_vars, n, n_slice, ga = build_sliced_hc_system(Nijk, r, base_F_func)
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    fkey_to_xvar = build_fkey_to_xvar_map(Nijk, r, one_vec)
    return (sys = sys, hc_vars = hc_vars, n = n, n_slice = n_slice,
            ga = ga, fkey_to_xvar = fkey_to_xvar)
end

# ============================================================
# Slice linear system as a plain numerical matrix
# ============================================================

"""
    build_slice_linear_system(ga, fkey_to_xvar, n) -> (A::Matrix{ComplexF64}, b::Vector{ComplexF64})

Produce the linear slice constraints as a numerical system

    A * x = b

where `x ∈ ℂ^n` is the pentagon variable vector (ordered by pentagon-variable
index) and `A` has `gauge_orbit_dim(ga)` rows.

Each row `k` corresponds to a gauge orbit vector `v_k` (orthonormal basis
column of `im Δ_gauge_eff`). The raw constraint in F-coord space is:

    Σ_i v_k[i] · (F_i - F_base_i) = 0.

For F-entries with a pentagon variable, `F_i = x_p`. For unit-axiom-fixed
entries, `F_i` is a constant equal to `F_base_i`, so those terms vanish.
What remains is:

    Σ_{i with var} v_k[i] · x_{p(i)}  =  Σ_{i with var} v_k[i] · F_base_i.

So `A[k, p(i)] += v_k[i]` and `b[k] += v_k[i] * F_base_i`, iterated over
F-entries that have a pentagon variable.

By construction `A · F_base_vec = b`, so the base F solution sits on the
slice.
"""
function build_slice_linear_system(ga::KitaevComplex.GaugeAnalysis,
                                   fkey_to_xvar::Dict,
                                   n::Int)
    fcs = ga.fcs
    F_svd = svd(ga.Delta_gauge_eff; full = true)
    rank_g = ga.gauge_orbit_dim
    gauge_basis = F_svd.U[:, 1:rank_g]   # nF × rank_g, columns = orthonormal gauge dirs

    A = zeros(ComplexF64, rank_g, n)
    b = zeros(ComplexF64, rank_g)

    for col in 1:rank_g
        v = gauge_basis[:, col]
        for (i, key) in enumerate(fcs.vars)
            coeff = v[i]
            abs(coeff) < 1e-12 && continue
            haskey(fkey_to_xvar, key) || continue     # skip unit-fixed entries
            pidx = fkey_to_xvar[key]
            Fval = ComplexF64(KitaevComplex.F_value(fcs, key))
            A[col, pidx] += coeff
            b[col]       += coeff * Fval
        end
    end
    return A, b
end

# ============================================================
# Newton with slice
# ============================================================

"""
    solve_pentagon_newton_with_slice(Nijk, r, base_F_func;
                                     initial_points = nothing,
                                     max_trials = 20, max_iter = 200,
                                     tol = 1e-12,
                                     perturb_scale = 0.1,
                                     slice_weight = 1.0,
                                     verbose = false)
        -> Vector{Vector{ComplexF64}}

Damped Newton solver for pentagon equations augmented with Kitaev gauge
slice constraints. Slice rows are appended to the residual vector and
Jacobian; normal equations `(J'J) δ = J'F` are solved via KrylovKit.

# Arguments

- `Nijk`, `r`, `base_F_func`: fusion data + base F-symbol (the base F
  must be a pentagon solution; it is used as the linearization point for
  Kitaev slice construction, and also as the default initial point for
  the first Newton trial).
- `initial_points`: if provided, a `Vector{Vector{ComplexF64}}` of starting
  points. If `nothing`, the base F (ordered by pentagon variable index)
  is used as the first trial, followed by `max_trials-1` random
  perturbations of scale `perturb_scale`.
- `slice_weight`: multiplier on slice residual / Jacobian rows (unused in
  normal equations under unit weight; included for experimentation).

# Returns

Vector of `ComplexF64` pentagon-variable vectors that satisfy both
pentagon and slice to tolerance `tol`.

# Rationale

Pentagon alone has continuous gauge orbit → Jacobian is rank-deficient
along gauge directions → Newton can drift along the orbit without
converging. Adding the Kitaev slice fills the rank deficiency canonically
(Hodge orthogonal to gauge orbit), allowing Newton to converge locally
near the base F and, with perturbations, to reach other gauge classes of
the same fusion ring (e.g., the 8 Ising gauge classes distinguished by
central charge).
"""
function solve_pentagon_newton_with_slice(Nijk::Array{Int,3}, r::Int,
                                          base_F_func::Function;
                                          initial_points::Union{Nothing,
                                              Vector{Vector{ComplexF64}}} = nothing,
                                          max_trials::Int = 20,
                                          max_iter::Int = 200,
                                          tol::Real = 1e-12,
                                          perturb_scale::Real = 0.1,
                                          slice_weight::Real = 1.0,
                                          verbose::Bool = false)
    # 1. Pentagon system (Oscar)
    R, eqs, n = PentagonEquations.get_pentagon_system(Nijk, r)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in eqs]

    # 2. Kitaev gauge analysis
    fr = FusionRule(Nijk)
    fcs = KitaevComplex.build_F_coord_space(fr, base_F_func)
    ga  = KitaevComplex.analyze_gauge(fcs)

    # 3. FKey ↔ pentagon var mapping
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    fkey_to_xvar = build_fkey_to_xvar_map(Nijk, r, one_vec)

    # 4. Slice as linear system A·x = b
    A_slice, b_slice = build_slice_linear_system(ga, fkey_to_xvar, n)
    n_slice = size(A_slice, 1)

    # 5. Base F vector (pentagon-variable-ordered)
    F_base_vec = zeros(ComplexF64, n)
    for (key, pidx) in fkey_to_xvar
        F_base_vec[pidx] = ComplexF64(base_F_func(key...))
    end

    # 6. Initial points
    if initial_points === nothing
        initial_points = Vector{Vector{ComplexF64}}()
        push!(initial_points, copy(F_base_vec))
        for _ in 2:max_trials
            push!(initial_points,
                  F_base_vec + perturb_scale * randn(ComplexF64, n))
        end
    end

    # 7. Newton loop
    solutions = Vector{Vector{ComplexF64}}()
    for (trial, x0) in enumerate(initial_points)
        x = copy(x0)

        for iter in 1:max_iter
            # Pentagon residual + slice residual
            F_pent = ComplexF64[PentagonSolver.eval_poly_complex(eq, x) for eq in eqs]
            F_sl   = A_slice * x - b_slice
            F_val  = vcat(F_pent, slice_weight * F_sl)
            res = maximum(abs.(F_val))

            if res < tol
                # Check solution not duplicate
                is_dup = any(norm(x - s) < 1e-8 for s in solutions)
                is_dup || push!(solutions, copy(x))
                verbose && println("  trial $trial: converged at iter $iter, res=$res")
                break
            end

            # Pentagon Jacobian + slice Jacobian (= A_slice, constant)
            J_pent = PentagonSolver.sparse_jacobian(eqs, derivs, x, n)
            J_sl   = sparse(slice_weight * A_slice)
            J = vcat(J_pent, J_sl)

            # Normal equations J' J δ = J' F (via KrylovKit)
            delta, _ = linsolve(v -> J' * (J * v), J' * F_val;
                                ishermitian = true, isposdef = true, verbosity = 0)

            # Damped line search
            alpha = 1.0
            for _ in 1:20
                x_new = x - alpha * delta
                F_pent_new = ComplexF64[PentagonSolver.eval_poly_complex(eq, x_new) for eq in eqs]
                F_sl_new   = A_slice * x_new - b_slice
                F_new      = vcat(F_pent_new, slice_weight * F_sl_new)
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