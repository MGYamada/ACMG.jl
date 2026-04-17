"""
Phase 2: Block-U parametrisation and MTC reconstruction.

This module implements the core of the Phase 2 pipeline:

1. `build_block_diagonal`: given a stratum {m_λ} and atomic catalog,
   build the block-diagonal atomic (S, T) on V = ⊕_λ V_λ^{m_λ}.

2. `reduce_to_Fp`: reduce an Oscar Q(ζ_N)-matrix to F_p (requires N | p-1
   and a chosen primitive N-th root of unity in F_p).

3. `t_eigenspace_decomposition`: given T (diagonal over F_p), group
   indices by eigenvalue. Returns a Dict{Int, Vector{Int}} mapping
   each T-eigenvalue to the list of basis indices.

4. `sweep_O2`: for a 2-dimensional degenerate eigenspace with indices
   (i, j), enumerate the O(2)(F_p)-rational circle points (u, v) with
   u² + v² = 1, and for each compute S' = U^T S U.

5. `verlinde_check`: compute fusion coefficients N_{ij}^k from an F_p
   modular datum (S', T'), check all are small non-negative integers.

The top-level driver is `find_mtcs_at_prime(catalog, stratum, p)`,
which runs the whole pipeline at a single prime.

Follow-up functions (in a separate CRT module, TODO) will combine
results across multiple primes.
"""

using Oscar

# ============================================================
#  Step 1: Build block-diagonal atomic data from a stratum
# ============================================================

"""
    build_block_diagonal(catalog::Vector{AtomicIrrep}, stratum::Stratum)
        -> (S_atomic, T_atomic, K, blocks)

Given a stratum (catalog index → multiplicity) and an atomic catalog,
construct the block-diagonal atomic modular data as matrices over a
common cyclotomic field K = Q(ζ_N).

Returns:
- S_atomic: block-diagonal S matrix over K (size r × r where r = stratum.total_dim)
- T_atomic: diagonal T as a Vector (length r) over K
- K: the shared cyclotomic field
- blocks: Vector of (start_idx, end_idx, catalog_idx, copy_num) tuples
  describing the block structure

Currently assumes all atomic irreps in the catalog use the SAME cyclotomic
field K. This is true when `build_atomic_catalog(N)` is used consistently.
"""
function build_block_diagonal(catalog::Vector{AtomicIrrep}, stratum::Stratum)
    # Find the common cyclotomic field
    K = nothing
    for (idx, _) in stratum.multiplicities
        if K === nothing
            K = catalog[idx].K
        else
            K === catalog[idx].K || error("Atomic irreps use different cyclotomic fields; not supported yet")
        end
    end
    K === nothing && error("Empty stratum")

    r = stratum.total_dim

    # Build block-diagonal S as r × r zero matrix, then fill blocks
    S_atomic = zero_matrix(K, r, r)
    T_atomic = Vector{elem_type(K)}(undef, r)
    blocks = Tuple{Int, Int, Int, Int}[]

    pos = 1
    # Sort for determinism
    for (idx, m) in sort(collect(stratum.multiplicities); by = first)
        atom = catalog[idx]
        d = atom.dim
        for copy_num in 1:m
            # Place atom.S at positions (pos:pos+d-1, pos:pos+d-1)
            for i in 1:d
                for j in 1:d
                    S_atomic[pos + i - 1, pos + j - 1] = atom.S[i, j]
                end
            end
            for i in 1:d
                T_atomic[pos + i - 1] = atom.T[i, i]
            end
            push!(blocks, (pos, pos + d - 1, idx, copy_num))
            pos += d
        end
    end

    return S_atomic, T_atomic, K, blocks
end

# ============================================================
#  Step 2: F_p reduction of Q(ζ_N) elements
# ============================================================

"""
    find_zeta_in_Fp(N::Int, p::Int) -> Int

Find a primitive N-th root of unity in F_p. Requires N | p-1.
Returns the F_p representative as an Int.

Uses: pick a primitive root g of F_p, then ζ_N = g^{(p-1)/N}.
"""
function find_zeta_in_Fp(N::Int, p::Int)
    (p - 1) % N == 0 || error("N = $N does not divide p - 1 = $(p-1)")
    g = primitive_root(p)
    exponent = div(p - 1, N)
    return powermod(g, exponent, p)
end

"""
    cyclotomic_to_Fp(x, zeta_N_Fp::Int, p::Int) -> Int

Reduce an Oscar cyclotomic field element x ∈ Q(ζ_N) to F_p, given
the value of ζ_N in F_p.

Works by:
1. Extracting coefficients of x in the basis {1, ζ, ζ², ..., ζ^{φ(N)-1}}
2. Reducing each rational coefficient mod p (requires denominator coprime to p)
3. Summing as linear combination of powers of zeta_N_Fp
"""
function cyclotomic_to_Fp(x, zeta_N_Fp::Int, p::Int)
    # Get the coefficient vector of x in the power basis of Q(ζ_N)
    coeffs = Oscar.coefficients(x)
    result = 0
    for (k, c) in enumerate(coeffs)
        # c is a rational number, k-1 is the power of ζ
        num = Int(numerator(c))
        den = Int(denominator(c))
        den % p != 0 || error("Denominator $den divisible by p = $p")
        coeff_Fp = mod(num * invmod(den, p), p)
        zeta_pow = powermod(zeta_N_Fp, k - 1, p)
        result = mod(result + coeff_Fp * zeta_pow, p)
    end
    return result
end

"""
    reduce_matrix_to_Fp(M, zeta_N_Fp::Int, p::Int) -> Matrix{Int}

Reduce an Oscar matrix over Q(ζ_N) to an Int matrix representing values mod p.
"""
function reduce_matrix_to_Fp(M, zeta_N_Fp::Int, p::Int)
    n, m = size(M)
    result = zeros(Int, n, m)
    for i in 1:n
        for j in 1:m
            result[i, j] = cyclotomic_to_Fp(M[i, j], zeta_N_Fp, p)
        end
    end
    return result
end

"""
    reduce_vector_to_Fp(v::Vector, zeta_N_Fp::Int, p::Int) -> Vector{Int}
"""
function reduce_vector_to_Fp(v::Vector, zeta_N_Fp::Int, p::Int)
    return [cyclotomic_to_Fp(x, zeta_N_Fp, p) for x in v]
end

# ============================================================
#  Step 3: T-eigenspace decomposition
# ============================================================

"""
    t_eigenspace_decomposition(T_Fp::Vector{Int}, p::Int)
        -> Dict{Int, Vector{Int}}

Given T as a vector of F_p eigenvalues (on the basis), return a Dict
mapping each distinct eigenvalue to the list of indices (1-based) where
it appears.
"""
function t_eigenspace_decomposition(T_Fp::Vector{Int}, p::Int)
    groups = Dict{Int, Vector{Int}}()
    for (i, val) in enumerate(T_Fp)
        if haskey(groups, val)
            push!(groups[val], i)
        else
            groups[val] = [i]
        end
    end
    return groups
end

"""
    parameter_dim(eigenspaces::Dict{Int, Vector{Int}}) -> Int

Compute parameter_dim = Σ_θ C(n_θ, 2) where n_θ = |eigenspaces[θ]|.
This is the continuous dimension of the block-U moduli.
"""
function parameter_dim(eigenspaces::Dict{Int, Vector{Int}})
    total = 0
    for (_, indices) in eigenspaces
        n = length(indices)
        total += div(n * (n - 1), 2)
    end
    return total
end

# ============================================================
#  Step 4: O(2) sweep for 2-dimensional degenerate eigenspace
# ============================================================

"""
    o2_circle_points(p::Int) -> Vector{Tuple{Int, Int}}

Enumerate all (u, v) ∈ F_p × F_p with u² + v² ≡ 1 (mod p).
For p ≡ 1 (mod 4), there are 2(p-1) such pairs; else 2(p+1).

For small p (p < 500 or so) a brute enumeration is instant.
"""
function o2_circle_points(p::Int)
    pts = Tuple{Int, Int}[]
    for u in 0:(p-1)
        u2 = mod(u * u, p)
        target = mod(1 - u2, p)
        # Is `target` a QR? Find sqrt
        for v in 0:(p-1)
            if mod(v * v, p) == target
                push!(pts, (u, v))
            end
        end
    end
    return pts
end

"""
    apply_o2_block(S::Matrix{Int}, idx_pair::Tuple{Int, Int},
                   u::Int, v::Int, det_sign::Int, p::Int) -> Matrix{Int}

Apply a 2×2 O(2) transformation to S in the 2-dimensional subspace
indexed by (i, j) = idx_pair, leaving all other rows/columns unchanged.

The transformation is U^T · S · U where U is the 2×2 block:
- det_sign = +1: U = [[u, -v], [v, u]]   (rotation)
- det_sign = -1: U = [[u,  v], [v, -u]]  (reflection)

In both cases (u, v) satisfies u² + v² ≡ 1 (mod p).
"""
function apply_o2_block(S::Matrix{Int}, idx_pair::Tuple{Int, Int},
                        u::Int, v::Int, det_sign::Int, p::Int)
    n = size(S, 1)
    U = Matrix{Int}(I_n(n))
    (i, j) = idx_pair
    if det_sign == +1
        U[i, i] = u
        U[i, j] = mod(-v, p)
        U[j, i] = v
        U[j, j] = u
    elseif det_sign == -1
        U[i, i] = u
        U[i, j] = v
        U[j, i] = v
        U[j, j] = mod(-u, p)
    else
        error("det_sign must be ±1")
    end
    UT = transpose_mod(U, p)
    return matmul_mod(matmul_mod(UT, S, p), U, p)
end

# Local helpers (small utilities, will move to FpArith if reused broadly)
function I_n(n::Int)
    M = zeros(Int, n, n)
    for i in 1:n
        M[i, i] = 1
    end
    return M
end

function transpose_mod(A::Matrix{Int}, p::Int)
    n, m = size(A)
    B = zeros(Int, m, n)
    for i in 1:m
        for j in 1:n
            B[i, j] = A[j, i]
        end
    end
    return B
end

# ============================================================
#  Step 5: Verlinde integrality check
# ============================================================

"""
    signed_Fp(x::Int, p::Int) -> Int

Lift x ∈ [0, p) to the symmetric interval [-⌊p/2⌋, ⌊p/2⌋].
"""
function signed_Fp(x::Int, p::Int)
    return x <= p ÷ 2 ? x : x - p
end

"""
    verlinde_find_unit(S_Fp::Matrix{Int}, p::Int; threshold::Int = 5)
        -> Union{Nothing, Tuple{Int, Array{Int,3}}}

Given an F_p modular datum S, try each basis index u as the candidate
"unit object". For unit u to be valid:
1. S[u, m] ≠ 0 for all m (so 1/S[u,m] is defined).
2. N_{u, i}^{k} = δ_{i,k} (unit fusion is identity).
3. All N_{i,j}^k small non-negative integers (|n| ≤ threshold).

If found, returns (u, N) where N is the r×r×r tensor. Otherwise nothing.

The Verlinde formula (assuming S^2 = I, self-dual MTC):
    N_{ij}^k = Σ_m (S[i,m] · S[j,m] · S[k,m]) / S[u,m]
"""
function verlinde_find_unit(S_Fp::Matrix{Int}, p::Int; threshold::Int = 5)
    r = size(S_Fp, 1)

    for u in 1:r
        S_u_row = [S_Fp[u, m] for m in 1:r]
        any(x -> x == 0, S_u_row) && continue

        try
            S_u_inv = [invmod(x, p) for x in S_u_row]

            # First check unit axiom: N[u, i, k] = δ_{i,k}
            is_unit = true
            for i in 1:r
                for k in 1:r
                    val = 0
                    for m in 1:r
                        term = mod(S_Fp[u, m] * S_Fp[i, m], p)
                        term = mod(term * S_Fp[k, m], p)
                        term = mod(term * S_u_inv[m], p)
                        val = mod(val + term, p)
                    end
                    expected = (i == k) ? 1 : 0
                    if val != expected
                        is_unit = false
                        break
                    end
                end
                is_unit || break
            end
            is_unit || continue

            # Compute full tensor, check non-negative integer
            N = zeros(Int, r, r, r)
            ok = true
            for i in 1:r
                for j in 1:r
                    for k in 1:r
                        val = 0
                        for m in 1:r
                            term = mod(S_Fp[i, m] * S_Fp[j, m], p)
                            term = mod(term * S_Fp[k, m], p)
                            term = mod(term * S_u_inv[m], p)
                            val = mod(val + term, p)
                        end
                        s = signed_Fp(val, p)
                        if s < 0 || abs(s) > threshold
                            ok = false
                            break
                        end
                        N[i, j, k] = s
                    end
                    ok || break
                end
                ok || break
            end
            ok && return (u, N)
        catch
            continue
        end
    end
    return nothing
end

# ============================================================
#  Step 6: Top-level single-prime driver
# ============================================================

"""
    MTCCandidate

Result of a successful block-U sweep at a single prime.

Fields:
- p: the prime used
- U_params: parameters of the block-U used (e.g. (u, v, det) for O(2))
- S_Fp: the S matrix mod p after block-U transformation
- T_Fp: the T eigenvalues mod p (unchanged by block-U)
- unit_index: which basis index is the unit object
- N: fusion tensor N[i][j][k] (signed integers)
- d: quantum dimensions d_i = S[unit, i] / S[unit, unit]
- D2: total quantum dimension squared
"""
struct MTCCandidate
    p::Int
    U_params::Any
    S_Fp::Matrix{Int}
    T_Fp::Vector{Int}
    unit_index::Int
    N::Array{Int, 3}
    d::Vector{Int}
    D2::Int
end

function Base.show(io::IO, c::MTCCandidate)
    d_signed = [signed_Fp(x, c.p) for x in c.d]
    print(io, "MTCCandidate(p=$(c.p), unit=$(c.unit_index), ",
          "d=$d_signed, D²=$(signed_Fp(c.D2, c.p)), ",
          "params=$(c.U_params))")
end

"""
    find_mtcs_at_prime(catalog::Vector{AtomicIrrep}, stratum::Stratum,
                       p::Int; verlinde_threshold::Int = 3)
        -> Vector{MTCCandidate}

Top-level Phase 2 driver at a single prime.

Steps:
1. Build block-diagonal atomic (S, T) from stratum + catalog.
2. Reduce to F_p (requires N | p-1 where N = catalog[].N).
3. Decompose T into eigenspaces; compute parameter_dim.
4. For each degenerate eigenspace of dimension 2, sweep O(2)(F_p).
   (Only one degenerate 2-dim eigenspace is currently supported;
   higher dimensions or multiple degenerate spaces TODO.)
5. For each (u, v, det) combination, apply block-U to S, then check
   Verlinde integrality.
6. Return all MTC candidates found.
"""
function find_mtcs_at_prime(catalog::Vector{AtomicIrrep}, stratum::Stratum,
                            p::Int; verlinde_threshold::Int = 3)
    # Common N
    N = catalog[first(keys(stratum.multiplicities))].N

    # Build block-diagonal atomic data
    S_K, T_K, K, blocks = build_block_diagonal(catalog, stratum)

    # Reduce to F_p
    zeta_N_Fp = find_zeta_in_Fp(N, p)
    S_atomic = [cyclotomic_to_Fp(S_K[i, j], zeta_N_Fp, p) for i in 1:size(S_K, 1), j in 1:size(S_K, 2)]
    T_atomic = [cyclotomic_to_Fp(t, zeta_N_Fp, p) for t in T_K]

    # T-eigenspace decomposition
    eigenspaces = t_eigenspace_decomposition(T_atomic, p)
    p_dim = parameter_dim(eigenspaces)

    # Classify strata by degeneracy structure
    # Currently supports:
    # - All n_θ = 1 (no sweep needed, just check atomic directly)
    # - Exactly one n_θ = 2 (one O(2) sweep)
    # Higher cases TODO

    degenerate = [(theta, indices) for (theta, indices) in eigenspaces if length(indices) >= 2]

    if isempty(degenerate)
        # No degeneracy - atomic S should already be (close to) MTC
        # (but signs and identifications may differ)
        r = size(S_atomic, 1)
        result = verlinde_find_unit(S_atomic, p; threshold = verlinde_threshold)
        if result !== nothing
            (u, N) = result
            Suu_inv = invmod(S_atomic[u, u], p)
            d = [mod(S_atomic[u, i] * Suu_inv, p) for i in 1:r]
            D2 = sum(mod(d[i] * d[i], p) for i in 1:r) % p
            return [MTCCandidate(p, :atomic, copy(S_atomic), copy(T_atomic),
                                 u, N, d, D2)]
        else
            return MTCCandidate[]
        end
    end

    length(degenerate) == 1 || error(
        "Only one degenerate eigenspace supported for now; got $(length(degenerate))")
    (theta_deg, indices_deg) = degenerate[1]

    if length(indices_deg) != 2
        error("Only n_θ = 2 supported in this prototype; got n_θ = $(length(indices_deg))")
    end

    i, j = indices_deg[1], indices_deg[2]
    idx_pair = (i, j)

    candidates = MTCCandidate[]
    for (u, v) in o2_circle_points(p)
        for det_sign in [+1, -1]
            S_prime = apply_o2_block(S_atomic, idx_pair, u, v, det_sign, p)
            result = verlinde_find_unit(S_prime, p; threshold = verlinde_threshold)
            if result !== nothing
                (u_idx, N) = result
                r = size(S_prime, 1)
                Suu_inv = invmod(S_prime[u_idx, u_idx], p)
                d = [mod(S_prime[u_idx, i] * Suu_inv, p) for i in 1:r]
                D2 = sum(mod(d[m] * d[m], p) for m in 1:r) % p
                push!(candidates, MTCCandidate(
                    p, (u, v, det_sign),
                    copy(S_prime), copy(T_atomic),
                    u_idx, N, d, D2))
            end
        end
    end
    return candidates
end
