"""
Phase 5 end-to-end pipeline demonstration.

Shows the `N → List[ClassifiedMTC]` driver `classify_mtcs_at_conductor`
in two modes:

1. `compute_FR_from_ST` on Fibonacci (rank 2, ~5 F-vars):
   Phase 4 alone — demonstrates the `(S, T) → (F, R)` wrapper.

2. `classify_mtcs_at_conductor` on SU(2)_4 (N=24, rank 5):
   Full Phase 0 → 3 + ℂ-lift, with `skip_FR=true` because rank-5
   pentagon HC (238 F-vars) is infeasible.

By default the SU(2)_4 demo passes the stratum explicitly to keep
runtime short (~1 minute). Pass `--full-enumerate` to enumerate all
strata at rank ≤ 5 (this takes much longer; ~25k strata × 7 primes).

Run:
    julia --project=. scripts/phase5_demo.jl
    julia --project=. scripts/phase5_demo.jl --full-enumerate
"""

using Oscar
using ACMG

full_enumerate = any(a -> a == "--full-enumerate", ARGS)

# ============================================================
#  Demo 1: Fibonacci (F, R) from (S, T)
# ============================================================

println("=" ^ 72)
println("  Phase 5 demo — Fibonacci (compute_FR_from_ST only)")
println("=" ^ 72)

Nijk_fib = zeros(Int, 2, 2, 2)
Nijk_fib[1, 1, 1] = 1
Nijk_fib[1, 2, 2] = 1
Nijk_fib[2, 1, 2] = 1
Nijk_fib[2, 2, 1] = 1
Nijk_fib[2, 2, 2] = 1
T_fib = ComplexF64[1.0, exp(4π * im / 5)]

println("\nFusion ring (Fibonacci): τ ⊗ τ = 1 ⊕ τ")
println("T = diag(1, exp(4πi/5))\n")

res_fib = ACMG.compute_FR_from_ST(Nijk_fib, T_fib;
                                   ribbon_atol = 1e-8,
                                   verbose = true)

println("\nResult:")
println("  pentagon HC: $(res_fib.n_pentagon) solutions")
println("  ribbon matches: $(res_fib.n_matches) of $(res_fib.n_tried) (F,R) pairs")
if res_fib.F !== nothing
    println("  chosen pair: F[$(res_fib.f_idx)] × R[$(res_fib.r_idx)]")
    println("  $(res_fib.report)")
else
    println("  ✗ no match found")
end

# ============================================================
#  Demo 2: SU(2)_4 full pipeline
# ============================================================

println("\n")
println("=" ^ 72)
println("  Phase 5 demo — SU(2)_4 full pipeline (N=24, rank 5, skip_FR)")
println("=" ^ 72)

primes_su24 = [73, 97, 193, 241, 313, 337, 409]

# Build catalog for pretty-printing (and, in quick mode, for stratum lookup).
catalog_su24 = ACMG.build_atomic_catalog(24; max_rank = 5, verbose = false)

strata_arg = if full_enumerate
    println("\n[Full enumeration of all strata at rank ≤ 5. This may take a long while.]")
    nothing  # `classify_mtcs_at_conductor` will enumerate
else
    println("\n[Quick mode: using known SU(2)_4 stratum (ρ_3, ρ_2).")
    println(" Pass --full-enumerate for the full sweep.]")
    d3_idx = 81
    d2_idx = 49
    [ACMG.Stratum(Dict(d3_idx => 1, d2_idx => 1), 5)]
end

classified = ACMG.classify_mtcs_at_conductor(24;
                                              max_rank = 5,
                                              primes = primes_su24,
                                              strata = strata_arg,
                                              scale_d = 3,
                                              scale_factor = 2,
                                              verlinde_threshold = 3,
                                              skip_FR = true,
                                              verbose = true)

# ============================================================
#  Summary
# ============================================================

println("\n" * "=" ^ 72)
println("  Summary")
println("=" ^ 72)
println("Total ClassifiedMTCs at N=24 (max_rank=5): $(length(classified))")

for (i, c) in enumerate(classified)
    println("\n[$i] $c")
    println("    stratum: $(ACMG.describe_stratum(c.stratum, catalog_su24))")
    println("    S (in ℤ[√3]), divide by $(c.scale_factor)·√3 for the true S:")
    println(ACMG.describe_matrix(c.S_Zsqrtd, 3))
end

println("\nExpected at N=24, rank 5: 2 Galois-sector SU(2)_4 MTCs.")
if full_enumerate
    println("Full enumeration may additionally surface lower-rank MTCs at N | 24.")
end
