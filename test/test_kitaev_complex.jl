"""
Tests for `ACMG.Phase4.KitaevComplex`.

Scope (first pass):
  - C2_basis produces expected sizes
  - fkey_to_xvar_map matches pentagon variable count
  - chi3_matrix has correct shape
  - chi3_constraints produces polynomials in the pentagon ring

Semantic verification (χ³ F_base = 0 on known solutions) is done in
`test_sliced_pentagon.jl` where it's naturally co-located with the
slice usage.
"""

using Test
using LinearAlgebra
using Oscar
using ACMG
using ACMG.Phase4

const KC = Phase4.KitaevComplex

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

@testset "KitaevComplex" begin

    @testset "C2_basis" begin
        # Fibonacci: N^{11}_1, N^{12}_2, N^{21}_2, N^{22}_1, N^{22}_2 → 5 basis elements
        basis_fib = KC.C2_basis(fib_Nijk(), 2)
        @test length(basis_fib) == 5
        @test (1, 1, 1) in basis_fib
        @test (2, 2, 1) in basis_fib
        @test (2, 2, 2) in basis_fib

        # Ising: count triples (a, b, y) with N^{ab}_y = 1
        basis_ising = KC.C2_basis(ising_Nijk(), 3)
        # 1⊗1=1, 1⊗ψ=ψ, 1⊗σ=σ, ψ⊗1=ψ, ψ⊗ψ=1, ψ⊗σ=σ, σ⊗1=σ, σ⊗ψ=σ, σ⊗σ=1+ψ (2)
        # = 1+1+1+1+1+1+1+1+2 = 10
        @test length(basis_ising) == 10
    end

    @testset "fkey_to_xvar_map: Fibonacci" begin
        Nijk = fib_Nijk()
        one_vec = [1, 0]
        fkey_map = KC.fkey_to_xvar_map(Nijk, 2, one_vec)
        # Fibonacci has 5 pentagon variables
        @test length(fkey_map) == 5
        vals = collect(values(fkey_map))
        @test sort(vals) == collect(1:5)
    end

    @testset "fkey_to_xvar_map: Ising" begin
        Nijk = ising_Nijk()
        one_vec = [1, 0, 0]
        fkey_map = KC.fkey_to_xvar_map(Nijk, 3, one_vec)
        @test length(fkey_map) == 14
        vals = collect(values(fkey_map))
        @test sort(vals) == collect(1:14)
    end

    @testset "chi3_matrix shape" begin
        Nijk = fib_Nijk()
        one_vec = [1, 0]
        fkey_map = KC.fkey_to_xvar_map(Nijk, 2, one_vec)
        # Use dummy quantum dims (will test real ones in sliced_pentagon)
        d = [1.0, 1.618]
        D2 = 1.0 + 1.618^2
        M, basis_kept = KC.chi3_matrix(Nijk, 2, d, D2, fkey_map)
        @test size(M, 2) == 5           # pentagon vars
        @test size(M, 1) == length(basis_kept)
        @test size(M, 1) ≤ 5             # ≤ |C²| = 5
    end

    @testset "chi3_constraints produces Oscar polynomials" begin
        Nijk = ising_Nijk()
        r = 3
        one_vec = [1, 0, 0]
        fkey_map = KC.fkey_to_xvar_map(Nijk, r, one_vec)
        d = [1.0, 1.0, sqrt(2.0)]
        D2 = 4.0
        # Need the pentagon ring. Use the one from get_pentagon_system.
        R, _eqs, n = Phase4.get_pentagon_system(Nijk, r)
        @test n == 14
        slice_polys = KC.chi3_constraints(Nijk, r, d, D2, R, fkey_map)
        @test all(p -> parent(p) === R, slice_polys)
        @test !isempty(slice_polys)
    end
end
