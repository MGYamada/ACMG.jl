"""
scripts/verify_ising_pentagon_independently.jl

Independently verify pentagon equations at each Newton-reported solution,
using a hand-written evaluation that does NOT rely on Oscar's
`get_pentagon_system`. If sol 2 / sol 3 are genuine pentagon solutions,
every pentagon equation should hold.

Pentagon equation used (standard):
    [F^{fcd}_e]_{g,l} * [F^{abl}_e]_{f,k}
      = Σ_h [F^{abc}_g]_{f,h} * [F^{ahd}_e]_{g,k} * [F^{bcd}_k]_{h,l}
"""

using LinearAlgebra
using Random
using ACMG
using ACMG.Phase4

const KC  = Phase4.KitaevComplex
const SPS = Phase4.SlicedPentagonSolver
const PS  = Phase4.PentagonSolver

function ising_Nijk()
    r = 3
    N = zeros(Int, r, r, r)
    for j in 1:r
        N[1, j, j] = 1
        N[j, 1, j] = 1
    end
    N[2, 2, 1] = 1
    N[2, 3, 3] = 1
    N[3, 2, 3] = 1
    N[3, 3, 1] = 1
    N[3, 3, 2] = 1
    return N
end

function ising_F_func()
    F_sss_s = (1.0 / sqrt(2.0)) * [1.0  1.0;
                                    1.0 -1.0]
    fr = FusionRule(ising_Nijk())
    function F(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
        (fr.N[a, b, e] == 1 && fr.N[e, c, d] == 1 &&
         fr.N[b, c, f] == 1 && fr.N[a, f, d] == 1) || return 0.0
        1 in (a, b, c, d) && return 1.0
        nsig = count(==(3), (a, b, c, d))
        nsig == 0 && return 1.0
        if nsig == 2
            (a, b, c, d) == (2, 3, 2, 3) && return -1.0
            (a, b, c, d) == (3, 2, 3, 1) && return  1.0
            (a, b, c, d) == (3, 2, 3, 2) && return -1.0
            return 1.0
        end
        if nsig == 4
            d != 3 && return 0.0
            return F_sss_s[e, f]
        end
        return 1.0
    end
    return F
end

"""Get F^{abc}_{d;e,f} from the pentagon-variable vector `Fvec`, using the
fkey_map. Falls back to 1.0 for unit-axiom entries not in the variable set,
and 0.0 for structurally disallowed entries."""
function F_of(Fvec::Vector{ComplexF64}, fr::FusionRule,
              fkey_map::Dict, a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
    # Allowed check
    (fr.N[a, b, e] ≥ 1 && fr.N[e, c, d] ≥ 1 &&
     fr.N[b, c, f] ≥ 1 && fr.N[a, f, d] ≥ 1) || return ComplexF64(0)
    # Unit-axiom: if any of (a,b,c,d) is unit (1), value is 1
    1 in (a, b, c, d) && return ComplexF64(1)
    key = (a, b, c, d, e, f)
    haskey(fkey_map, key) || return ComplexF64(0)
    return Fvec[fkey_map[key]]
end

"""Evaluate all pentagon equations at Fvec. Returns vector of residuals."""
function pentagon_residuals_direct(Fvec::Vector{ComplexF64}, fr::FusionRule,
                                    fkey_map::Dict, r::Int)
    residuals = ComplexF64[]
    for a in 1:r, b in 1:r, c in 1:r, d in 1:r, e in 1:r,
        f in 1:r, g in 1:r, k in 1:r, l in 1:r

        # LHS: F^{fcd}_{e;g,l} · F^{abl}_{e;f,k}
        lhs1 = F_of(Fvec, fr, fkey_map, f, c, d, e, g, l)
        lhs1 == 0 && continue
        lhs2 = F_of(Fvec, fr, fkey_map, a, b, l, e, f, k)
        lhs2 == 0 && continue
        lhs = lhs1 * lhs2

        # RHS: Σ_h F^{abc}_{g;f,h} · F^{ahd}_{e;g,k} · F^{bcd}_{k;h,l}
        rhs = ComplexF64(0)
        for h in 1:r
            x1 = F_of(Fvec, fr, fkey_map, a, b, c, g, f, h)
            x1 == 0 && continue
            x2 = F_of(Fvec, fr, fkey_map, a, h, d, e, g, k)
            x2 == 0 && continue
            x3 = F_of(Fvec, fr, fkey_map, b, c, d, k, h, l)
            x3 == 0 && continue
            rhs += x1 * x2 * x3
        end

        push!(residuals, lhs - rhs)
    end
    return residuals
end

function main()
    Nijk = ising_Nijk()
    F_fn = ising_F_func()
    r = 3
    fr = FusionRule(Nijk)

    # Reproduce Test 4
    Random.seed!(20260418)
    sols = SPS.solve_pentagon_newton_with_slice(Nijk, r, F_fn;
        max_trials    = 10,
        max_iter      = 2000,
        perturb_scale = 0.3,
        tol           = 1e-14,
        verbose       = false)
    println("Found $(length(sols)) solutions via Newton+slice.")

    # fkey map for direct evaluation
    one_vec = [1, 0, 0]
    fkey_map = SPS.build_fkey_to_xvar_map(Nijk, r, one_vec)
    n = length(fkey_map)

    # Also get Oscar's pentagon eqs for comparison
    R, oscar_eqs, n_vars = Phase4.get_pentagon_system(Nijk, r)
    println("Oscar pentagon system: $(length(oscar_eqs)) equations, $n_vars variables.")

    println("\n" * "="^68)
    println("Direct pentagon verification (hand-written, not via Oscar)")
    println("="^68)

    for (i, s) in enumerate(sols)
        println("\n--- Solution $i ---")

        # Direct pentagon residual
        direct_res = pentagon_residuals_direct(s, fr, fkey_map, r)
        println("  Direct pentagon:  #eqs = $(length(direct_res))")
        println("                    max |residual| = $(maximum(abs.(direct_res)))")

        # Oscar pentagon residual
        oscar_res = ComplexF64[PS.eval_poly_complex(eq, s) for eq in oscar_eqs]
        println("  Oscar pentagon:   #eqs = $(length(oscar_res))")
        println("                    max |residual| = $(maximum(abs.(oscar_res)))")

        # How many direct eqs significantly violate
        n_violate_direct = count(r -> abs(r) > 1e-10, direct_res)
        n_violate_oscar  = count(r -> abs(r) > 1e-10, oscar_res)
        println("  violations (|res| > 1e-10): direct=$n_violate_direct, oscar=$n_violate_oscar")

        # Show first few direct violations with their (a,b,c,d,e,f,g,k,l)
        if n_violate_direct > 0
            println("  First few direct violations:")
            count_shown = 0
            idx = 0
            for a in 1:r, b in 1:r, c in 1:r, d in 1:r, e in 1:r,
                f in 1:r, g in 1:r, k in 1:r, l in 1:r
                idx += 1
                idx > length(direct_res) && break
                res = direct_res[idx]
                if abs(res) > 1e-10
                    count_shown += 1
                    count_shown > 5 && break
                    # Note: direct_res doesn't carry these labels; this is approximate
                    println("     residual: $res")
                end
            end
        end
    end

    # Now: how many F entries are unit-axiom-fixed (1.0) in Ising?
    # That is, how many F-variables does Oscar skip that our direct code uses?
    println("\n" * "="^68)
    println("Ising F-entry accounting")
    println("="^68)
    n_allowed = 0
    n_unit_fixed = 0
    for a in 1:r, b in 1:r, c in 1:r, d in 1:r, e in 1:r, f in 1:r
        if (fr.N[a, b, e] ≥ 1 && fr.N[e, c, d] ≥ 1 &&
            fr.N[b, c, f] ≥ 1 && fr.N[a, f, d] ≥ 1)
            n_allowed += 1
            if 1 in (a, b, c, d)
                n_unit_fixed += 1
            end
        end
    end
    println("  Total allowed F-entries:      $n_allowed")
    println("  Unit-axiom-fixed:             $n_unit_fixed")
    println("  Variable F-entries:           $(n_allowed - n_unit_fixed)")
    println("  Pentagon variables (Oscar):   $n_vars")
    println("  Pentagon vars (fkey_map):     $n")
end

main()
