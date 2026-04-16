# mtc_classifier.jl
# MTC Classification via Conductor Decomposition
# April 14, 2026 — Arithmetic Condensed Matter Geometry
#
# Level I:  N = ∏ pᵢ^aᵢ → SL₂(Z/pᵢ^aᵢ Z) irreps → tensor products → modular data
# Level II: Pentagon — damped Newton over C
#
# Usage:
#   include("mtc_classifier.jl")
#   classify_mtc(5)
#   classify_mtc(12)

using LinearAlgebra
using SparseArrays
using Oscar
using TensorCategories
using KrylovKit

# ============================================================
#  Level I: Modular data via CRT decomposition
# ============================================================

function load_sl2reps()
    GAP.Globals.LoadPackage(GapObj("SL2Reps"))
end

# --- GAP → Oscar conversion ---

function gap_to_oscar_matrix(gap_mat, r)
    cond = 1
    for i in 1:r, j in 1:r
        cond = lcm(cond, Int(GAP.Globals.Conductor(gap_mat[i,j])))
    end

    K, z = cyclotomic_field(cond)
    M = zero_matrix(K, r, r)

    for i in 1:r, j in 1:r
        coeffs = GAP.Globals.CoeffsCyc(gap_mat[i,j], cond)
        val = zero(K)
        for k in 1:cond
            if !Bool(GAP.Globals.IsZero(coeffs[k]))
                num = Int(GAP.Globals.NumeratorRat(coeffs[k]))
                den = Int(GAP.Globals.DenominatorRat(coeffs[k]))
                val += QQ(num, den) * z^(k-1)
            end
        end
        M[i,j] = val
    end
    return M, K, cond
end

function gap_to_complex_matrix(gap_mat, r)
    cond = 1
    for i in 1:r, j in 1:r
        cond = lcm(cond, Int(GAP.Globals.Conductor(gap_mat[i,j])))
    end
    z = exp(2π * im / cond)
    M = zeros(ComplexF64, r, r)
    for i in 1:r, j in 1:r
        coeffs = GAP.Globals.CoeffsCyc(gap_mat[i,j], cond)
        val = 0.0 + 0.0im
        for k in 1:cond
            if !Bool(GAP.Globals.IsZero(coeffs[k]))
                num = Int(GAP.Globals.NumeratorRat(coeffs[k]))
                den = Int(GAP.Globals.DenominatorRat(coeffs[k]))
                val += (num / den) * z^(k-1)
            end
        end
        M[i,j] = val
    end
    return M
end

# --- CRT decomposition of N ---

function prime_power_factors(N::Int)
    """Factor N = p₁^a₁ × ... × pₖ^aₖ. Returns [(p₁,a₁), ...]"""
    result = Tuple{Int,Int}[]
    for (p, a) in factor(N)
        push!(result, (Int(p), Int(a)))
    end
    result
end

# --- Enumerate irreps for each prime power factor ---

function irreps_for_prime_power(p::Int, a::Int)
    """Get all irreps of SL₂(Z/p^a Z) from SL2Reps.
    These have level dividing p^a."""
    irreps = []
    for lam in 1:a
        lev = p^lam
        reps = GAP.evalstr("SL2IrrepsOfLevel($lev)")
        len = Int(GAP.Globals.Length(reps))
        for i in 1:len
            push!(irreps, reps[i])
        end
    end
    # Also include trivial rep (level 1)
    reps_1 = GAP.evalstr("SL2IrrepsOfLevel(1)")
    len_1 = Int(GAP.Globals.Length(reps_1))
    for i in 1:len_1
        push!(irreps, reps_1[i])
    end
    irreps
end

# --- Tensor product of S,T matrices ---

function kron_complex(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    """Kronecker product of two complex matrices."""
    ra, ca = size(A)
    rb, cb = size(B)
    C = zeros(ComplexF64, ra*rb, ca*cb)
    for i in 1:ra, j in 1:ca, k in 1:rb, l in 1:cb
        C[(i-1)*rb + k, (j-1)*cb + l] = A[i,j] * B[k,l]
    end
    C
end

function tensor_product_reps(reps_list)
    """Given a list of (S, T) pairs (one per prime power factor),
    compute the tensor product representation.
    S_total = S₁ ⊗ S₂ ⊗ ..., T_total = T₁ ⊗ T₂ ⊗ ..."""

    S_total = reps_list[1][1]
    T_total = reps_list[1][2]

    for i in 2:length(reps_list)
        S_total = kron_complex(S_total, reps_list[i][1])
        T_total = kron_complex(T_total, reps_list[i][2])
    end

    return S_total, T_total
end

# --- Verlinde (float) ---

function verlinde_float(S_num::Matrix{ComplexF64}, r::Int)
    # S[1,l] が 0 に近い場合は invalid
    for l in 1:r
        if abs(S_num[1,l]) < 1e-10
            return nothing
        end
    end

    Nijk = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        val = sum(S_num[i,l] * S_num[j,l] * conj(S_num[k,l]) / S_num[1,l] for l in 1:r)
        n = round(Int, real(val))
        if abs(real(val) - n) > 1e-4 || abs(imag(val)) > 1e-4
            return nothing
        end
        Nijk[i,j,k] = n
    end
    Nijk
end

# --- Fix quantum dimensions ---

function fix_quantum_dims_float(S_num::Matrix{ComplexF64}, r::Int)
    signs = ones(Int, r)
    for i in 1:r
        d = S_num[1,i] / S_num[1,1]
        if real(d) < -1e-10
            signs[i] = -1
        end
    end
    S_fixed = copy(S_num)
    for i in 1:r, j in 1:r
        S_fixed[i,j] = signs[i] * signs[j] * S_num[i,j]
    end
    return S_fixed, signs
end

# --- Admissibility checks (BNRW 2016) ---

function check_ST_relation(S::Matrix{ComplexF64}, T::Matrix{ComplexF64}, r::Int; tol=1e-4)
    """For projective SL₂(Z) representation from SL2Reps: (ST)³ = S²"""
    ST3 = (S * T)^3
    S2 = S^2
    diff1 = ST3 - S2
    diff2 = ST3 + S2  # 符号の不定性
    return min(maximum(abs.(diff1)), maximum(abs.(diff2))) < tol
end

function check_cauchy(S::Matrix{ComplexF64}, N::Int, r::Int; tol=1e-4)
    """Cauchy theorem: prime factors of D² ⊆ prime factors of N."""
    d = [real(S[1,i] / S[1,1]) for i in 1:r]
    D_sq = sum(x^2 for x in d)
    N_primes = Set(Int(p) for (p, _) in factor(N))

    # Try D² ≈ integer / small denominator
    for mult in 1:60
        x = D_sq * mult
        x_rounded = round(Int, x)
        if abs(x - x_rounded) < tol * max(abs(x), 1.0) && x_rounded != 0
            # D² = x_rounded / mult (but mult divides N's primes typically)
            num_primes = x_rounded == 0 ? Set{Int}() : Set(Int(p) for (p, _) in factor(abs(x_rounded)))
            den_primes = mult == 1 ? Set{Int}() : Set(Int(p) for (p, _) in factor(mult))
            return issubset(union(num_primes, den_primes), N_primes)
        end
    end
    # Couldn't identify D² as rational; accept by default
    return true
end

function check_frobenius_schur(S::Matrix{ComplexF64}, T::Matrix{ComplexF64}, r::Int; tol=1e-3)
    """ν_2(X) must be integer-valued. ν_2(1) = 1."""
    d = [real(S[1,i] / S[1,1]) for i in 1:r]
    D_sq = sum(x^2 for x in d)
    theta = [T[i,i] for i in 1:r]

    # Verlinde fusion coefficients (complex)
    N_vl = zeros(ComplexF64, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        N_vl[i,j,k] = sum(S[i,l] * S[j,l] * conj(S[k,l]) / S[1,l] for l in 1:r)
    end

    for X in 1:r
        nu_2 = sum(N_vl[X,X,Y] * d[Y] * theta[Y]^2 for Y in 1:r) / D_sq
        nu_2_int = round(Int, real(nu_2))
        if abs(imag(nu_2)) > tol || abs(real(nu_2) - nu_2_int) > tol
            return false
        end
        if abs(nu_2_int) > r
            return false
        end
        # Unit object: ν_2(1) must equal 1
        if X == 1 && nu_2_int != 1
            return false
        end
    end
    return true
end

function admissibility_check(S, T, N, r)
    if !check_ST_relation(S, T, r)
        return false, "(ST)³ ≠ ±S²"
    end
    # Cauchy check, FS-indicator は保留（代数的整数の扱いが必要）
    return true, "OK"
end

# --- Level I main: enumerate all modular data candidates ---

function enumerate_modular_data(N::Int; max_rank::Int=20)
    """Enumerate modular data candidates for conductor N.

    Algorithm:
    1. Factor N = p₁^a₁ × ... × pₖ^aₖ
    2. For each factor, enumerate irreps of SL₂(Z/pᵢ^aᵢ Z)
    3. Take tensor products across factors
    4. Filter by Verlinde non-negativity
    """

    pf = prime_power_factors(N)
    println("N = $N = ", join(["$(p)^$(a)" for (p,a) in pf], " × "))

    # Step 1: Get irreps for each prime power factor
    factor_irreps = []
    for (p, a) in pf
        irreps = irreps_for_prime_power(p, a)
        # Convert to ComplexF64 (S, T) pairs
        st_pairs = []
        for rep in irreps
            r = Int(rep.degree)
            S = gap_to_complex_matrix(rep.S, r)
            T = gap_to_complex_matrix(rep.T, r)
            push!(st_pairs, (S, T, r))
        end
        push!(factor_irreps, st_pairs)
        println("  SL₂(Z/$(p^a)Z): $(length(st_pairs)) irreps")
    end

    # Step 2: Take tensor products (Cartesian product across factors)
    valid_candidates = []

    if length(factor_irreps) == 1
        # N is a prime power: just use single irreps
        for (S, T, r) in factor_irreps[1]
            r > max_rank && continue
            S_fixed, signs = fix_quantum_dims_float(S, r)
            Nijk = verlinde_float(S_fixed, r)
            if Nijk !== nothing && all(Nijk .>= 0)
                ok, reason = admissibility_check(S_fixed, T, N, r)
                if !ok
                    println("  Rejected rank=$r: $reason")
                    continue
                end
                push!(valid_candidates, (S=S_fixed, T=T, Nijk=Nijk, rank=r))
            end
        end
    else
        # N has multiple prime factors: tensor products
        # Use Iterators.product for Cartesian product
        indices = [1:length(fi) for fi in factor_irreps]

        for idx_tuple in Iterators.product(indices...)
            # Collect one irrep per factor
            reps = [(factor_irreps[k][idx_tuple[k]][1],   # S
                     factor_irreps[k][idx_tuple[k]][2])    # T
                    for k in 1:length(pf)]

            # Total rank
            total_rank = prod(factor_irreps[k][idx_tuple[k]][3] for k in 1:length(pf))
            total_rank > max_rank && continue

            # Tensor product
            S_total, T_total = tensor_product_reps(reps)
            r = size(S_total, 1)

            # Fix quantum dims and Verlinde
            S_fixed, signs = fix_quantum_dims_float(S_total, r)
            Nijk = verlinde_float(S_fixed, r)

            if Nijk !== nothing && all(Nijk .>= 0)
                ok, reason = admissibility_check(S_fixed, T_total, N, r)
                if !ok
                    continue
                end
                # Compute level = lcm of factor levels
                factor_names = join(["$(factor_irreps[k][idx_tuple[k]][3])d" for k in 1:length(pf)], "⊗")
                push!(valid_candidates, (S=S_fixed, T=T_total, Nijk=Nijk, rank=r, name=factor_names))
            end
        end
    end

    println("\n$(length(valid_candidates)) valid modular data candidates.")
    for (i, cand) in enumerate(valid_candidates)
        name = hasproperty(cand, :name) ? cand.name : ""
        println("  [$i] rank=$(cand.rank) $name")
    end

    return valid_candidates
end

# ============================================================
#  Level II: Pentagon solver (damped Newton over C)
# ============================================================

function get_pentagon_system(Nijk::Array{Int,3}, r::Int)
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    C, eqs_raw = pentagon_equations(Nijk, one_vec)
    eqs = filter(eq -> !(eq isa Integer) && !iszero(eq), eqs_raw)
    R = parent(eqs[1])
    n = nvars(R)
    return R, eqs, n
end

function eval_poly_complex(f, vals::Vector{ComplexF64})
    isa(f, Integer) && return ComplexF64(f)
    iszero(f) && return ComplexF64(0.0)
    result = 0.0 + 0.0im
    for (c, m) in zip(coefficients(f), monomials(f))
        degs = degrees(m)
        term = ComplexF64(Float64(numerator(c)) / Float64(denominator(c)))
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vals[i]^degs[i])
        end
        result += term
    end
    result
end

function jacobian_complex(eqs, vals::Vector{ComplexF64}, n)
    m = length(eqs)
    J = zeros(ComplexF64, m, n)
    for i in 1:m
        isa(eqs[i], Integer) && continue
        iszero(eqs[i]) && continue
        for j in 1:n
            df = derivative(eqs[i], j)
            J[i,j] = eval_poly_complex(df, vals)
        end
    end
    J
end

function solve_pentagon_sparse_newton(eqs, n, zeta; max_trials=20, max_iter=200)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in eqs]

    solutions = []
    for trial in 1:max_trials
        x = ComplexF64[rand(-3:3) + rand(-3:3)*zeta for _ in 1:n]
        for i in 1:n
            abs(x[i]) < 0.1 && (x[i] = 1.0 + zeta)
        end

        for iter in 1:max_iter
            F_val = ComplexF64[eval_poly_complex(eq, x) for eq in eqs]
            res = maximum(abs.(F_val))
            if res < 1e-12
                println("  Trial $trial: CONVERGED at iter $iter")
                push!(solutions, copy(x))
                break
            end

            # Sparse Jacobian
            J = sparse_jacobian(eqs, derivs, x, n)

            # Normal equation: J'J δ = J'F, solved by conjugate gradient
            delta, info = linsolve(v -> J' * (J * v), J' * F_val; ishermitian = true, isposdef = true, verbosity = 0)

            # Line search
            alpha = 1.0
            for _ in 1:20
                x_new = x - alpha * delta
                F_new = ComplexF64[eval_poly_complex(eq, x_new) for eq in eqs]
                if maximum(abs.(F_new)) < res
                    x = x_new
                    break
                end
                alpha *= 0.5
            end
        end
    end
    solutions
end

function sparse_jacobian(eqs, derivs, x, n)
    m = length(eqs)
    I = Int[]
    J_idx = Int[]
    V = ComplexF64[]

    for i in 1:m
        for j in 1:n
            v = eval_poly_complex(derivs[i][j], x)
            if abs(v) > 1e-15
                push!(I, i)
                push!(J_idx, j)
                push!(V, v)
            end
        end
    end

    sparse(I, J_idx, V, m, n)
end

# ============================================================
#  Main
# ============================================================

function classify_mtc(N::Int; max_rank::Int=20, max_trials::Int=20)
    println("=" ^ 60)
    println("MTC Classification: N = $N")
    println("=" ^ 60)

    if N == 1
        println("N=1: trivial MTC (rank 1, F = 1)")
        println("=" ^ 60)
        return
    end

    # Level I
    candidates = enumerate_modular_data(N; max_rank=max_rank)

    # Level II
    for (idx, cand) in enumerate(candidates)
        println("\n--- Candidate $idx: rank=$(cand.rank) ---")

        if cand.rank == 1
            println("  Pointed (rank 1). F = trivial.")
            continue
        end

        R, eqs, n = get_pentagon_system(cand.Nijk, cand.rank)
        println("  Pentagon: $n variables, $(length(eqs)) equations")

        if n == 0
            println("  No free F-symbols. Trivial solution.")
            continue
        end

        zeta = exp(2π * im / N)
        solutions = solve_pentagon_sparse_newton(eqs, n, zeta; max_trials=max_trials)
        println("  $(length(solutions)) solutions found")
    end

    println("\n" * "=" ^ 60)
    println("N=$N classification complete.")
    println("=" ^ 60)
end

load_sl2reps()
println("MTC Classifier loaded.")
println("Usage: classify_mtc(N)")
println("Examples: classify_mtc(5), classify_mtc(12)")
