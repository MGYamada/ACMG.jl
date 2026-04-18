"""
Diagnostic: inspect S_Fp values across primes to understand the
Galois/permutation structure.
"""

using Oscar
using ACMG

catalog = ACMG.build_atomic_catalog(24; max_rank = 5, verbose = false)

d3_idx = 81
d2_idx = 49
stratum = ACMG.Stratum(Dict(d3_idx => 1, d2_idx => 1), 5)

test_primes = [73, 97, 193, 241, 313, 337, 409]

# Collect candidates at each prime, extract (S_Fp, T_Fp, d, params)
println("=== Per-prime candidate comparison ===\n")

all_candidates = Dict{Int, Vector{ACMG.MTCCandidate}}()
for p in test_primes
    cands = ACMG.find_mtcs_at_prime(catalog, stratum, p; verlinde_threshold = 3)
    all_candidates[p] = cands
end

# For each prime, print candidates
for p in test_primes
    cands = all_candidates[p]
    println("p = $p: $(length(cands)) candidates")
    for (i, c) in enumerate(cands)
        d_signed = [ACMG.signed_Fp(x, c.p) for x in c.d]
        T_signed = [ACMG.signed_Fp(x, c.p) for x in c.T_Fp]
        params_signed = (ACMG.signed_Fp(c.U_params[1], c.p),
                         ACMG.signed_Fp(c.U_params[2], c.p),
                         c.U_params[3])
        println("  [$i] unit=$(c.unit_index), params=$params_signed, d=$d_signed")
        println("       T=$T_signed")
    end
    println()
end

# Now: for two primes, compute a test entry M[i,j] = 2·√3·S[i,j]
# using cyclotomic √3, and compare across primes

println("=== Check one entry across primes (cyclotomic √3) ===\n")
println("Entry M[0,0] (Python-style 0-indexed = Julia [1,1]) = 2·√3·S[1,1]")

# Collect group 1 at each prime using fusion-tensor matching
# Group 1 rep (from our earlier run): unit=1, some specific N
rep = all_candidates[73][1]

for p in test_primes
    # Find matching candidate
    matching = nothing
    for c in all_candidates[p]
        if c.N == rep.N && c.unit_index == rep.unit_index
            matching = c
            break
        end
    end

    if matching === nothing
        println("  p=$p: no matching fusion tensor")
        continue
    end

    s3 = ACMG.compute_sqrt3_cyclotomic_mod_p(p)
    two_s3 = mod(2 * s3, p)
    entry = mod(two_s3 * matching.S_Fp[1, 1], p)
    entry_signed = ACMG.signed_Fp(entry, p)
    params_signed = (ACMG.signed_Fp(matching.U_params[1], p),
                     ACMG.signed_Fp(matching.U_params[2], p),
                     matching.U_params[3])
    println("  p=$p: 2√3·S[1,1] = $entry (signed: $entry_signed), params=$params_signed")
end

# Also: compare full S_Fp · (2√3) (signed) across primes
println("\n=== Full 2√3·S matrix, signed, across primes ===")
for p in test_primes
    matching = nothing
    for c in all_candidates[p]
        if c.N == rep.N && c.unit_index == rep.unit_index
            matching = c
            break
        end
    end
    matching === nothing && continue
    s3 = ACMG.compute_sqrt3_cyclotomic_mod_p(p)
    println("\np=$p (√3 cyc = $(ACMG.signed_Fp(s3, p))):")
    M = [mod(2 * s3 * matching.S_Fp[i, j], p) for i in 1:5, j in 1:5]
    for i in 1:5
        row = [ACMG.signed_Fp(M[i, j], p) for j in 1:5]
        println("  ", row)
    end
end
