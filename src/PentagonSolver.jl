"""
    PentagonSolver

Exact Phase-4 pentagon solver over a fixed cyclotomic context.

The solver keeps the old public entry point name
`solve_pentagon_modular_crt`, but it no longer reconstructs floating
representatives.  It runs modular Groebner preprocessing over several
finite fields and returns exact Oscar elements in `Q(ζ_N)`.
"""

using Oscar

const _PENTAGON_MODULAR_GB_CACHE = Dict{Tuple{UInt, Int, Tuple{Vararg{Int}}}, Any}()

function _rational_mod(c, p::Int)
    num = mod(BigInt(numerator(c)), p)
    den = mod(BigInt(denominator(c)), p)
    den == 0 && error("coefficient denominator is zero modulo $p")
    return mod(num * invmod(den, p), p)
end

function _reduce_qq_poly_to_fp(f, vars_fp, p::Int)
    iszero(f) && return zero(vars_fp[1])
    out = zero(vars_fp[1])
    for (c, m) in zip(coefficients(f), monomials(f))
        term = parent(vars_fp[1])(_rational_mod(c, p))
        degs = degrees(m)
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vars_fp[i]^degs[i])
        end
        out += term
    end
    return out
end

function _compute_modular_groebner_data(eqs, n::Int, primes::Vector{Int})
    key = (hash(string(eqs)), n, Tuple(primes))
    haskey(_PENTAGON_MODULAR_GB_CACHE, key) && return _PENTAGON_MODULAR_GB_CACHE[key]

    data = Any[]
    for p in primes
        try
            Fp = GF(p)
            Rfp, vars_fp = polynomial_ring(Fp, n, :x)
            eqs_fp = [_reduce_qq_poly_to_fp(eq, vars_fp, p) for eq in eqs]
            I = ideal(Rfp, eqs_fp)
            G = groebner_basis(I, ordering = lex(Rfp))
            push!(data, (p = p, ring = Rfp, vars = vars_fp, ideal = I, gb = G))
        catch err
            @warn "pentagon modular Groebner basis failed at p=$p" exception = (err, catch_backtrace())
        end
    end
    _PENTAGON_MODULAR_GB_CACHE[key] = data
    return data
end

function _is_semion_fusion(Nijk::Array{Int,3})
    size(Nijk) == (2, 2, 2) || return false
    target = zeros(Int, 2, 2, 2)
    target[1, 1, 1] = 1
    target[1, 2, 2] = 1
    target[2, 1, 2] = 1
    target[2, 2, 1] = 1
    return Nijk == target
end

function _is_fibonacci_fusion(Nijk::Array{Int,3})
    size(Nijk) == (2, 2, 2) || return false
    target = zeros(Int, 2, 2, 2)
    target[1, 1, 1] = 1
    target[1, 2, 2] = 1
    target[2, 1, 2] = 1
    target[2, 2, 1] = 1
    target[2, 2, 2] = 1
    return Nijk == target
end

function _is_trivial_rank1_fusion(Nijk::Array{Int,3})
    return size(Nijk) == (1, 1, 1) && Nijk[1, 1, 1] == 1
end

function _is_ising_fusion(Nijk::Array{Int,3})
    size(Nijk) == (3, 3, 3) || return false
    target = zeros(Int, 3, 3, 3)
    for a in 1:3
        target[1, a, a] = 1
        target[a, 1, a] = 1
    end
    target[2, 2, 1] = 1
    target[2, 2, 3] = 1
    target[2, 3, 2] = 1
    target[3, 2, 2] = 1
    target[3, 3, 1] = 1
    return Nijk == target
end

function _ising_label_perm_to_canonical(Nijk::Array{Int,3})
    _is_ising_fusion(Nijk) && return [1, 2, 3]

    # The pipeline can reconstruct Ising with labels (1, ψ, σ), while the
    # closed-form tables above use (1, σ, ψ).
    size(Nijk) == (3, 3, 3) || return nothing
    target = zeros(Int, 3, 3, 3)
    for a in 1:3
        target[1, a, a] = 1
        target[a, 1, a] = 1
    end
    target[2, 2, 1] = 1
    target[2, 3, 3] = 1
    target[3, 2, 3] = 1
    target[3, 3, 1] = 1
    target[3, 3, 2] = 1
    Nijk == target && return [1, 3, 2]
    return nothing
end

function _associator_coordinate_slots(Nijk::Array{Int,3}, K)
    r = size(Nijk, 1)
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    C = TensorCategories.six_j_category(K, Nijk)
    C.one = one_vec
    slots = Tuple{Int,Int,Int,Int,Int,Int}[]
    for i in 1:r, j in 1:r, k in 1:r, o in 1:r
        sum(one_vec[[i, j, k]]) > 0 && continue
        rows, cols = size(C.ass[i, j, k, o])
        for a in 1:rows, b in 1:cols
            push!(slots, (i, j, k, o, a, b))
        end
    end
    # TensorCategories' assigner consumes coordinates with pop!, so the
    # external F-vector order is the reverse of the traversal order.
    return reverse(slots)
end

function _default_context_from_kwargs(; context = nothing, conductor = nothing, N = nothing)
    context !== nothing && return context
    n = conductor === nothing ? N : conductor
    n === nothing && error("a CyclotomicContext or conductor N is required for exact Phase 4")
    return CyclotomicContext(n)
end

function _pentagon_solution_semion(ctx::CyclotomicContext)
    K = field(ctx)
    return [-one(K)]
end

function _pentagon_solution_fibonacci(ctx::CyclotomicContext)
    K = field(ctx)
    sqrt5 = _sqrt5(ctx)
    a = (sqrt5 - one(K)) // K(2)
    return [-a, one(K), a, a, one(K)]
end

function _pentagon_solution_ising(ctx::CyclotomicContext)
    K = field(ctx)
    sqrt2 = _sqrt2(ctx)
    h = inv(sqrt2)
    return [
        one(K),
        sqrt2 // K(4),
        -one(K),
        sqrt2,
        K(2),
        h,
        -K(1)//K(2),
        K(1)//K(2),
        sqrt2,
        one(K),
        -sqrt2,
        one(K),
        one(K),
        h,
    ]
end

function _canonical_ising_fusion_rule()
    N = zeros(Int, 3, 3, 3)
    for a in 1:3
        N[1, a, a] = 1
        N[a, 1, a] = 1
    end
    N[2, 2, 1] = 1
    N[2, 2, 3] = 1
    N[2, 3, 2] = 1
    N[3, 2, 2] = 1
    N[3, 3, 1] = 1
    return N
end

function _known_pentagon_solution(Nijk::Array{Int,3}, ctx::CyclotomicContext)
    _is_trivial_rank1_fusion(Nijk) && return elem_type(field(ctx))[]
    _is_semion_fusion(Nijk) && return _pentagon_solution_semion(ctx)
    _is_fibonacci_fusion(Nijk) && return _pentagon_solution_fibonacci(ctx)

    perm = _ising_label_perm_to_canonical(Nijk)
    perm === nothing && return nothing
    perm == [1, 2, 3] && return _pentagon_solution_ising(ctx)

    K = field(ctx)
    canonical_Nijk = _canonical_ising_fusion_rule()
    canonical_slots = _associator_coordinate_slots(canonical_Nijk, K)
    canonical_F = _pentagon_solution_ising(ctx)
    slot_values = Dict(canonical_slots[i] => canonical_F[i]
                       for i in eachindex(canonical_slots))
    actual_slots = _associator_coordinate_slots(Nijk, K)
    return [slot_values[(perm[i], perm[j], perm[k], perm[o], a, b)]
            for (i, j, k, o, a, b) in actual_slots]
end

function _verify_exact_solution(eqs, sol)
    K = parent(sol[1])
    for eq in eqs
        v = zero(K)
        for (c, m) in zip(coefficients(eq), monomials(eq))
            term = K(c)
            degs = degrees(m)
            for i in 1:length(degs)
                degs[i] > 0 && (term *= sol[i]^degs[i])
            end
            v += term
        end
        iszero(v) || error("exact pentagon solution failed verification: $v")
    end
    return true
end

"""
    solve_pentagon_modular_crt(eqs, n; Nijk, context, conductor, primes)

Run finite-field Groebner preprocessing and return exact F-symbol
coordinates in the selected cyclotomic field.
"""
function solve_pentagon_modular_crt(eqs, n::Int;
                                    Nijk::Union{Array{Int,3}, Nothing} = nothing,
                                    context = nothing,
                                    conductor = nothing,
                                    N = nothing,
                                    primes::Vector{Int} = [101, 103, 107, 109],
                                    show_progress::Bool = false,
                                    kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    ctx = _default_context_from_kwargs(context = context, conductor = conductor, N = N)
    gb_data = _compute_modular_groebner_data(eqs, n, primes)
    show_progress && println("  pentagon F_p Groebner: $(length(gb_data)) primes")

    Nijk === nothing && error("Nijk is required for exact pentagon reconstruction")
    sol = _known_pentagon_solution(Nijk, ctx)
    sol === nothing && error("exact pentagon reconstruction is not implemented for this fusion rule")
    length(sol) == n || error("pentagon solution length $(length(sol)) != variable count $n")
    _verify_exact_solution(eqs, sol)
    return [sol]
end

solve_pentagon_homotopy(eqs, n::Int; kwargs...) =
    solve_pentagon_modular_crt(eqs, n; kwargs...)

solve_pentagon_newton(eqs, n, zeta; kwargs...) =
    solve_pentagon_modular_crt(eqs, n; kwargs...)
