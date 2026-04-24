using Test
using ACMG

@testset "Pipeline galois grouping" begin
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
