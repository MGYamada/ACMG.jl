# mtc_classifier.jl
# MTC Classification via Conductor Decomposition
# April 14, 2026 — Arithmetic Condensed Matter Geometry
#
# Algorithm:
#   Level I:   SL₂(Z/NZ) reps → modular data (S,T)
#   Level II: Pentagon/hexagon over F_p — linear propagation — CRT
#
# Key insight: pentagon with 4 known F-symbols is LINEAR in the 5th.
#   Non-homogeneous → solve uniquely (ax = b)
#   Homogeneous     → gauge freedom, fix x = 1
#   No backtracking needed.
#
# Usage:
#   using Oscar
#   include("mtc_classifier.jl")
#   classify_mtc(5)

using Oscar

# ============================================================
#  Level I: Modular data
# ============================================================

function load_sl2reps()
    GAP.Globals.LoadPackage(GapObj("SL2Reps"))
end

function sl2_irreps(N::Int)
    GAP.evalstr("SL2IrrepsOfLevel($N)")
end

function gap_to_oscar_matrix(gap_mat, r)
    # Step 1: 全成分の conductor の lcm
    cond = 1
    for i in 1:r, j in 1:r
        cond = lcm(cond, Int(GAP.Globals.Conductor(gap_mat[i,j])))
    end
    
    K, z = cyclotomic_field(cond)
    M = zero_matrix(K, r, r)
    
    for i in 1:r, j in 1:r
        # GAP 側で直接 cond に統一して係数を取る
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

function fix_quantum_dims(S, K, r, cond)
    z = exp(2π * im / cond)
    
    signs = ones(Int, r)
    for i in 1:r
        d = S[1,i] * inv(S[1,1])
        # AbsSimpleNumFieldElem の多項式係数を取って数値評価
        deg = degree(K)
        val = sum(Float64(coeff(d, k)) * z^k for k in 0:deg-1)
        if real(val) < -1e-10
            signs[i] = -1
        end
    end
    
    S_fixed = copy(S)
    for i in 1:r, j in 1:r
        S_fixed[i,j] = signs[i] * signs[j] * S[i,j]
    end
    return S_fixed, signs
end

function verlinde(S, r::Int, K)
    aut = complex_conjugation(K)
    Nijk = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        val = sum(S[i,l] * S[j,l] * aut(S[k,l]) * inv(S[1,l]) for l in 1:r)
        # AbsSimpleNumFieldElem → Rational → Int
        Nijk[i,j,k] = Int(coeff(val, 0))
    end
    Nijk
end

function is_valid_fusion(Nijk, r)
    all(Nijk[i,j,k] >= 0 for i in 1:r, j in 1:r, k in 1:r)
end

# ============================================================
#  Level II: Pentagon over F_p
# ============================================================

function make_system(Nijk::Array{Int,3}, r::Int)
    values = Dict{NTuple{6,Int}, Int}()
    free = NTuple{6,Int}[]
    gauge_used = Set{NTuple{3,Int}}()

    for a in 1:r, b in 1:r, c in 1:r, d in 1:r, e in 1:r, f in 1:r
        (Nijk[a,b,e]*Nijk[e,c,d]*Nijk[b,c,f]*Nijk[a,f,d] == 0) && continue
        sext = (a,b,c,d,e,f)
        if a == 1 || b == 1 || c == 1
            values[sext] = 1
        elseif (a,b,e) ∉ gauge_used
            values[sext] = 1
            push!(gauge_used, (a,b,e))
        else
            push!(free, sext)
        end
    end
    (values=values, free=free, fusion=Nijk, rank=r)
end

# --- TODO: Pentagon evaluator (your convention) ---

function F_get(values, sext, p)
    # Look up F-symbol value, return mod p
    mod(get(values, sext, 0), p)
end

function pentagon_eval(a, b, c, d, g, values, Nijk, r, p)
    # Evaluate pentagon equation for (a,b,c,d,g).
    # Returns LHS - RHS mod p.
    #
    # Pentagon (Turaev convention):
    #   Σ_h F[a,b,c; g; e,h] * F[a,h,c,d; g; ?,?] * F[b,c,d; ?; h,?]
    #     = F[b,c,d; g; ?,?] * F[a,?,d; g; ?,?]
    #
    # TODO: fill in exact index contraction for your convention
    0  # placeholder
end

# --- Core: Linear propagation ---

function solve_over_Fp(system, p::Int)
    values = Dict{NTuple{6,Int},Int}(k => mod(v,p) for (k,v) in system.values)
    unsolved = Set(system.free)
    r = system.rank
    Nijk = system.fusion

    changed = true
    while changed
        changed = false
        for a in 1:r, b in 1:r, c in 1:r, d in 1:r, g in 1:r
            result = try_solve_linear(a, b, c, d, g, values, unsolved, Nijk, r, p)
            if result !== nothing
                var, val = result
                values[var] = val
                delete!(unsolved, var)
                changed = true
            end
        end
    end

    # Remaining unsolved = gauge freedom → fix to 1
    for var in unsolved
        values[var] = 1
    end

    # Verify
    if verify_all(values, Nijk, r, p)
        return values
    end
    nothing
end

function try_solve_linear(a, b, c, d, g, values, unsolved, Nijk, r, p)
    # Evaluate pentagon (a,b,c,d,g). 
    # If exactly one unknown remains and equation is linear:
    #   non-homogeneous (coeff ≠ 0) → return (var, solution)
    #   homogeneous (coeff = 0)     → return nothing (skip)
    # Otherwise → return nothing
    #
    # TODO: implement with your pentagon convention.
    # This is the ONLY function you need to write carefully.
    nothing
end

function verify_all(values, Nijk, r, p)
    for a in 1:r, b in 1:r, c in 1:r, d in 1:r, g in 1:r
        if pentagon_eval(a, b, c, d, g, values, Nijk, r, p) % p != 0
            return false
        end
    end
    true
end

# --- Hexagon (linear in R given F) ---

function solve_hexagon_Fp(F_values, Nijk, r, p)
    # With F known, hexagon is linear in R. Solve directly.
    # TODO: implement
    Dict{NTuple{3,Int}, Int}()
end

# ============================================================
#  CRT Reconstruction
# ============================================================

function good_primes(N::Int, count::Int)
    result = Int[]
    p = N + 1
    while length(result) < count
        is_prime(p) && p % N == 1 && push!(result, p)
        p += 1
    end
    result
end

function crt_reconstruct(residues, N)
    K, zeta = cyclotomic_field(N)
    # TODO: recover Z[ζ_N] element from F_p images
    zero(K)
end

# ============================================================
#  Main
# ============================================================

function classify_mtc(N::Int; num_primes::Int=5)
    println("MTC Classification: N = $N")
    primes = good_primes(N, num_primes)
    println("Good primes (p ≡ 1 mod $N): $primes")

    # SL2Reps
    candidates = sl2_irreps(N)

    # Verlinde
    for rep in candidates
        S_gap, T_gap = rep.S, rep.T
        r = rep.degree
        S, K, cond = gap_to_oscar_matrix(S_gap, r)
        T, = gap_to_oscar_matrix(T_gap, r)
        S_fixed, _ = fix_quantum_dims(S, K, r, cond)
        Nijk = verlinde(S_fixed, r, K)
        is_valid_fusion(Nijk, r) || continue
        println("Valid fusion rules found for rep of degree $(r).")
    end

    # solve_over_Fp

    # CRT
end

load_sl2reps()
println("MTC Classifier loaded. Usage: classify_mtc(N)")
