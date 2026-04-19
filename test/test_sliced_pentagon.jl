"""
Tests for `ACMG.Phase4.SlicedPentagonSolver`.

Strategy:
  1. Fibonacci FKey ↔ pentagon variable mapping sanity.
  2. Fibonacci: build_sliced_hc_system returns a valid HC system,
     base F satisfies every equation in it.
  3. Ising: same checks. Reports the system size for info.

Actual HC solving (timing / mixed volume) is deferred to a separate
script to avoid inflating test runtime.
"""

using Test
using LinearAlgebra
using ACMG
using ACMG.Phase4
import HomotopyContinuation as HC

const KC  = Phase4.KitaevComplex
const SPS = Phase4.SlicedPentagonSolver

# ---------------- Fibonacci helpers ----------------

function fib_Nijk()
    r = 2
    N = zeros(Int, r, r, r)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    N[2, 2, 2] = 1
    return N
end

function fib_F_func()
    φ = (1 + sqrt(5)) / 2
    sφ = sqrt(φ)
    F_τττ_τ = [1/φ    1/sφ;
               1/sφ  -1/φ]
    fr = FusionRule(fib_Nijk())
    function F(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
        (fr.N[a, b, e] == 1 && fr.N[e, c, d] == 1 &&
         fr.N[b, c, f] == 1 && fr.N[a, f, d] == 1) || return 0.0
        1 in (a, b, c, d) && return 1.0
        return F_τττ_τ[e, f]
    end
    return F
end

# ---------------- Ising helpers ----------------

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

@testset "SlicedPentagonSolver" begin

    # --------------------------------------------------------------------
    @testset "Fibonacci FKey ↔ pentagon-var map" begin
        Nijk = fib_Nijk()
        one_vec = [1, 0]
        fkey_map = SPS.build_fkey_to_xvar_map(Nijk, 2, one_vec)

        @test length(fkey_map) > 0
        vals = collect(values(fkey_map))
        @test all(1 .≤ vals .≤ 5)            # Fibonacci has 5 pentagon vars
        @test length(unique(vals)) == length(vals)

        for key in keys(fkey_map)
            i, j, k, _o, _e, _f = key
            @test !(1 in (i, j, k))
        end
    end

    # --------------------------------------------------------------------
    @testset "Fibonacci sliced HC system" begin
        Nijk = fib_Nijk()
        F_fn = fib_F_func()

        result = SPS.get_sliced_pentagon_system_hc(Nijk, 2, F_fn)
        @test result.n == 5
        @test result.n_slice == 1

        # Base F as a ComplexF64 vector, ordered by pentagon variable index
        F_vec = zeros(ComplexF64, result.n)
        for (key, pidx) in result.fkey_to_xvar
            F_vec[pidx] = F_fn(key...)
        end

        # Evaluate HC system at base F — should all be ~0
        residuals = result.sys(F_vec)
        max_resid = maximum(abs.(residuals))
        @test max_resid < 1e-10

        @info "Fibonacci sliced HC system" n=result.n n_slice=result.n_slice n_equations=length(residuals) max_residual_at_base=max_resid
    end

    # --------------------------------------------------------------------
    @testset "Ising sliced HC system" begin
        Nijk = ising_Nijk()
        F_fn = ising_F_func()

        result = SPS.get_sliced_pentagon_system_hc(Nijk, 3, F_fn)
        @test result.n > 0
        @test result.n_slice == 1             # Ising effective gauge = 1

        F_vec = zeros(ComplexF64, result.n)
        for (key, pidx) in result.fkey_to_xvar
            F_vec[pidx] = F_fn(key...)
        end

        residuals = result.sys(F_vec)
        max_resid = maximum(abs.(residuals))
        @test max_resid < 1e-10

        @info "Ising sliced HC system" n=result.n n_slice=result.n_slice n_equations=length(residuals) max_residual_at_base=max_resid
    end

    # --------------------------------------------------------------------
    @testset "Newton with slice: base F converges instantly" begin
        # Fibonacci base F is a pentagon solution; starting Newton from it
        # with slice should converge immediately.
        sols = SPS.solve_pentagon_newton_with_slice(fib_Nijk(), 2, fib_F_func();
                                                     max_trials = 1,
                                                     max_iter = 50,
                                                     perturb_scale = 0.0,
                                                     tol = 1e-12)
        @test length(sols) == 1

        # Ising: same
        sols = SPS.solve_pentagon_newton_with_slice(ising_Nijk(), 3, ising_F_func();
                                                     max_trials = 1,
                                                     max_iter = 50,
                                                     perturb_scale = 0.0,
                                                     tol = 1e-12)
        @test length(sols) == 1
    end

    # --------------------------------------------------------------------
    @testset "Newton with slice: small perturbation converges back" begin
        # Fibonacci: starting from base F + small noise, slice-augmented
        # Newton should recover the base F (unique within the slice class).
        sols = SPS.solve_pentagon_newton_with_slice(fib_Nijk(), 2, fib_F_func();
                                                     max_trials = 5,
                                                     max_iter = 200,
                                                     perturb_scale = 0.05,
                                                     tol = 1e-10)
        @test length(sols) ≥ 1

        # Ising
        sols = SPS.solve_pentagon_newton_with_slice(ising_Nijk(), 3, ising_F_func();
                                                     max_trials = 5,
                                                     max_iter = 200,
                                                     perturb_scale = 0.05,
                                                     tol = 1e-10)
        @test length(sols) ≥ 1
        @info "Ising Newton+slice small-perturb" n_sols=length(sols)
    end
end
