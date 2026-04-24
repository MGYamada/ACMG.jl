using Test

@testset "ACMG Prototype" begin
    include("test_fparith.jl")
    include("test_fibonacci.jl")
    include("test_ising.jl")
    include("test_sl2reps.jl")
    include("test_blocku.jl")
    include("test_stratum_enum.jl")
    include("test_block_u_general.jl")
    include("test_crt.jl")
    include("test_fr_pentagon_hexagon_fibonacci.jl")
    include("test_modular_data_lift.jl")
    include("test_verify_residuals.jl")
    include("test_pipeline_primes.jl")
    include("test_pipeline_galois_grouping.jl")
    include("test_pipeline_end_to_end.jl")
end
