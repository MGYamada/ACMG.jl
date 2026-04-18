using Test
using ACMG
using ACMG.Phase4

"""
Phase 4 Verify tests.

Strategy: solve Fibonacci pentagon/hexagon via Phase 4, then check
`verify_mtc` gives near-zero residuals. This re-validates the solver
output independently.

Fibonacci data:
- rank 2, N=5
- τ ⊗ τ = 1 ⊕ τ
- all self-dual: dual = [1, 2]
- twists: θ_1 = 1, θ_τ = exp(4πi/5)

Ribbon check uses the double-braiding form
    θ_i² · R^{i*,i}_1 · R^{i,i*}_1 = 1
which is independent of the sign of θ_i, so both of Fibonacci's two
braidings must pass.
"""

@testset "Phase 4 Verify" begin
    # ----- Fibonacci setup -----
    Nijk = zeros(Int, 2, 2, 2)
    Nijk[1, 1, 1] = 1
    Nijk[1, 2, 2] = 1
    Nijk[2, 1, 2] = 1
    Nijk[2, 2, 1] = 1
    Nijk[2, 2, 2] = 1

    dual = [1, 2]
    T_expected = ComplexF64[1.0, exp(4π * im / 5)]

    # Solve pentagon and hexagon via Phase 4
    _R, eqs, n = Phase4.get_pentagon_system(Nijk, 2)
    F_sols = Phase4.solve_pentagon_homotopy(eqs, n;
                                            slice = 1, show_progress = false)
    @test !isempty(F_sols)

    F = Phase4.refine_solution_newton(eqs, F_sols[1]; tol = 1e-14)

    _R_ring, hex_eqs, n_r = Phase4.get_hexagon_system(Nijk, 2, F)
    R_sols = Phase4.solve_hexagon_homotopy(hex_eqs, n_r; show_progress = false)
    @test !isempty(R_sols)

    R = R_sols[1]

    @testset "Pentagon residuals" begin
        pent_res = Phase4.pentagon_residuals(F, Nijk)
        @test !isempty(pent_res)
        max_pent = maximum(pent_res)
        println("  Pentagon max residual: $max_pent")
        @test max_pent < 1e-8
    end

    @testset "Hexagon residuals" begin
        hex_res = Phase4.hexagon_residuals(F, R, Nijk)
        @test !isempty(hex_res)
        max_hex = maximum(hex_res)
        println("  Hexagon max residual: $max_hex")
        @test max_hex < 1e-8
    end

    @testset "R-block extraction" begin
        # For Fibonacci, R^{22}_1 and R^{22}_2 are scalars (multiplicity-free).
        # Also R^{11}_1 = 1, R^{12}_2 = 1, R^{21}_2 = 1 (trivial unit ones).
        R11_1 = Phase4.extract_R_block(R, Nijk, 1, 1, 1)
        @test size(R11_1) == (1, 1)
        # Unit-involved R should be identity (up to gauge choice).
        # We don't test the exact value since TensorCategories' gauge is
        # not pinned; we just test the extraction returns a 1x1 block.

        R22_1 = Phase4.extract_R_block(R, Nijk, 2, 2, 1)
        @test size(R22_1) == (1, 1)
        R22_2 = Phase4.extract_R_block(R, Nijk, 2, 2, 2)
        @test size(R22_2) == (1, 1)

        # No such block:
        # Nijk[1,1,2] = 0, extract should return empty
        empty_block = Phase4.extract_R_block(R, Nijk, 1, 1, 2)
        @test size(empty_block) == (0, 0)
    end

    @testset "Ribbon consistency (double braiding)" begin
        # New formula: θ_i² · R^{i*,i}_1 · R^{i,i*}_1 = 1
        # This is independent of the sign of θ_i, so BOTH Fibonacci braidings
        # (R_sols[1] and R_sols[2]) must pass.
        rib_res = Phase4.ribbon_residuals(R, T_expected, Nijk, dual)
        @test length(rib_res) == 2
        max_rib = maximum(rib_res)
        println("  Ribbon max residual (R_sols[1]): $max_rib")
        @test max_rib < 1e-8

        # Test on the OTHER braiding too — should also pass
        R_alt = R_sols[2]
        rib_alt = Phase4.ribbon_residuals(R_alt, T_expected, Nijk, dual)
        max_rib_alt = maximum(rib_alt)
        println("  Ribbon max residual (R_sols[2]): $max_rib_alt")
        @test max_rib_alt < 1e-8
    end

    @testset "verify_mtc roundtrip" begin
        report = Phase4.verify_mtc(F, R, Nijk;
                                    T = T_expected, dual = dual)
        @test report.rank == 2
        @test report.pentagon_max < 1e-8
        @test report.hexagon_max < 1e-8
        @test report.n_pentagon_eqs > 0
        @test report.n_hexagon_eqs > 0
        @test report.ribbon_max !== nothing
        @test report.ribbon_max < 1e-8
        println("  $report")
    end
end
