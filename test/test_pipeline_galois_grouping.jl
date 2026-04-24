using Test
using ACMG

@testset "Pipeline galois grouping" begin
    @testset "classify_mtcs_at_conductor on SU(2)_4 (N=24, skip_FR)" begin
        test_primes = [73, 97, 193, 241, 313, 337, 409]

        catalog = ACMG.build_atomic_catalog(24; max_rank = 5, verbose = false)
        d3_idx = 81
        d2_idx = 49
        stratum = ACMG.Stratum(Dict(d3_idx => 1, d2_idx => 1), 5)

        classified = ACMG.classify_mtcs_at_conductor(24;
                                                     max_rank = 5,
                                                     primes = test_primes,
                                                     strata = [stratum],
                                                     scale_d = 3,
                                                     scale_factor = 2,
                                                     verlinde_threshold = 3,
                                                     skip_FR = true,
                                                     verbose = false)

        @test length(classified) == 2
        for c in classified
            @test c.rank == 5
            @test c.verify_fresh
            @test c.scale_d == 3
            @test c.scale_factor == 2
            @test c.F_values === nothing
            @test c.R_values === nothing
            @test c.verify_report === nothing
            @test size(c.S_complex) == (5, 5)
            @test length(c.T_complex) == 5
            for t in c.T_complex
                @test abs(abs(t) - 1.0) < 1e-10
            end
        end

        sectors = sort([c.galois_sector for c in classified])
        @test sectors == [1, 2]
        for c in classified
            @test c.N_input == 24
            @test c.N == 24
        end
    end

    @testset "_branch_consistency_precheck resolves anchored sign contradiction" begin
        d = 6
        p_anchor = 29
        p_other = 53
        selector = ACMG.build_sqrtd_selector(d, [p_anchor, p_other], p_anchor; verbose = false)

        N1 = zeros(Int, 1, 1, 1)
        N1[1, 1, 1] = 1
        x = (a = 3, b = 1)
        s_anchor = ACMG.compute_sqrt_d_mod_p(d, p_anchor)
        s_other_raw = ACMG.compute_sqrt_d_mod_p(d, p_other)
        s_other_true = mod(-s_other_raw, p_other)

        S_anchor = mod((x.a + x.b * s_anchor) * invmod(2 * s_anchor, p_anchor), p_anchor)
        S_other = mod((x.a + x.b * s_other_true) * invmod(2 * s_other_true, p_other), p_other)

        c_anchor = ACMG.MTCCandidate(p_anchor, :dummy, reshape([S_anchor], 1, 1),
                                     [1], 1, N1, [1], 1)
        c_other = ACMG.MTCCandidate(p_other, :dummy, reshape([S_other], 1, 1),
                                    [1], 1, N1, [1], 1)
        results = Dict(p_anchor => [c_anchor], p_other => [c_other])

        contradictions = ACMG._branch_consistency_precheck(results, p_anchor, d, selector.sqrtd_fn;
                                                           branch_sign_getter = selector.branch_sign_getter,
                                                           branch_sign_setter = selector.branch_sign_setter,
                                                           verbose = false)
        @test isempty(contradictions)
        @test selector.branch_sign_getter(p_other) in (-1, 1)
    end
end
