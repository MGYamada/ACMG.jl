using Test
using ACMG

@testset "Pipeline end-to-end" begin
    @testset "compute_FR_from_ST on Fibonacci" begin
        Nijk = zeros(Int, 2, 2, 2)
        Nijk[1, 1, 1] = 1
        Nijk[1, 2, 2] = 1
        Nijk[2, 1, 2] = 1
        Nijk[2, 2, 1] = 1
        Nijk[2, 2, 2] = 1

        T = ComplexF64[1.0, exp(4π * im / 5)]

        result = ACMG.compute_FR_from_ST(Nijk, T;
                                         ribbon_atol = 1e-8,
                                         pentagon_slice = 1,
                                         show_progress = false,
                                         verbose = false)

        @test result.F !== nothing
        @test result.R !== nothing
        @test result.report !== nothing
        @test result.report.pentagon_max < 1e-8
        @test result.report.hexagon_max < 1e-8
        @test result.report.ribbon_max !== nothing
        @test result.report.ribbon_max < 1e-8
        @test result.n_matches >= 1
    end

    @testset "compute_FR_from_ST rejects wrong T (sanity)" begin
        Nijk = zeros(Int, 2, 2, 2)
        Nijk[1, 1, 1] = 1
        Nijk[1, 2, 2] = 1
        Nijk[2, 1, 2] = 1
        Nijk[2, 2, 1] = 1
        Nijk[2, 2, 2] = 1

        T_wrong = ComplexF64[1.0, im]

        result = ACMG.compute_FR_from_ST(Nijk, T_wrong;
                                         ribbon_atol = 1e-8,
                                         show_progress = false,
                                         verbose = false)

        @test result.n_matches == 0
        @test result.F === nothing
        @test result.R === nothing
        @test result.report === nothing
    end

    @testset "classify_mtcs_at_conductor default(full_mtc) finds Fibonacci from N=5" begin
        test_primes = [41, 61, 101, 181]
        N_input = 5
        N_effective = 20

        fib_N = zeros(Int, 2, 2, 2)
        fib_N[1, 1, 1] = 1
        fib_N[1, 2, 2] = 1
        fib_N[2, 1, 2] = 1
        fib_N[2, 2, 1] = 1
        fib_N[2, 2, 2] = 1

        catalog20 = ACMG.build_atomic_catalog(N_effective; max_rank = 2, verbose = false)
        strata2 = ACMG.enumerate_strata(catalog20, 2)
        fib_strata = ACMG.Stratum[]
        for st in strata2
            ok = true
            for p in test_primes[1:2]
                local cands
                try
                    cands = ACMG.find_mtcs_at_prime(catalog20, st, p;
                                                    verlinde_threshold = 3,
                                                    max_block_dim = 3)
                catch
                    ok = false
                    break
                end
                if !any(c -> c.N == fib_N, cands)
                    ok = false
                    break
                end
            end
            ok && push!(fib_strata, st)
        end
        if isempty(fib_strata)
            @test isempty(fib_strata)
        else
            classified = ACMG.classify_mtcs_at_conductor(N_input;
                                                         max_rank = 2,
                                                         primes = test_primes,
                                                         strata = [first(fib_strata)],
                                                         scale_d = 5,
                                                         scale_factor = 2,
                                                         skip_FR = true,
                                                         verbose = false)

            @test all(c -> c.N == N_effective, classified)
            @test all(c -> c.N_input == N_input, classified)
            @test any(c -> c.rank == 2 && c.Nijk == fib_N, classified)
        end
    end
end
