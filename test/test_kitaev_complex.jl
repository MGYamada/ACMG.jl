"""
Tests for `ACMG.Phase4.KitaevComplex`.

Core principle: the Kitaev identity χⁿ⁺¹ δⁿ + δⁿ⁻¹ χⁿ = id_Cⁿ is the
SOLE correctness criterion for the ansatz. If it passes at every n for
which the matrices are implemented, the ansatz is proven correct.

Currently implemented for n = 0 (C⁰ = ℂ, trivial) and n = 1. Later:
n = 2 once interior face maps are written.
"""

using Test
using LinearAlgebra
using ACMG
using ACMG.Phase4

const KC = Phase4.KitaevComplex

# ---------------- Fusion rules and dims ----------------

function fib_Nijk()
    N = zeros(Int, 2, 2, 2)
    N[1,1,1] = 1; N[1,2,2] = 1; N[2,1,2] = 1; N[2,2,1] = 1; N[2,2,2] = 1
    return N
end

function ising_Nijk()
    N = zeros(Int, 3, 3, 3)
    for j in 1:3
        N[1,j,j] = 1; N[j,1,j] = 1
    end
    N[2,2,1] = 1; N[2,3,3] = 1; N[3,2,3] = 1; N[3,3,1] = 1; N[3,3,2] = 1
    return N
end

function fib_dims()
    φ = (1 + sqrt(5)) / 2
    d = [1.0, φ]
    D2 = 1 + φ^2
    return d, D2
end

function ising_dims()
    d = [1.0, 1.0, sqrt(2.0)]
    D2 = 4.0
    return d, D2
end

@testset "KitaevComplex — Cⁿ basis" begin
    @testset "Fibonacci" begin
        N = fib_Nijk()
        @test KC.C_dim(N, 2, 0) == 1
        @test KC.C_dim(N, 2, 1) == 2     # two simples
        @test KC.C_dim(N, 2, 2) == 5     # fusion paths: 1⊗1, 1⊗τ, τ⊗1, τ⊗τ→1, τ⊗τ→τ
    end

    @testset "Ising" begin
        N = ising_Nijk()
        @test KC.C_dim(N, 3, 0) == 1
        @test KC.C_dim(N, 3, 1) == 3
        # n = 2: basis = (X,Y,Z) with N^{XY}_Z = 1.
        # Counts per (X,Y): (1,1)→1; (1,ψ)→1; (1,σ)→1; (ψ,1)→1; (ψ,ψ)→1;
        # (ψ,σ)→1; (σ,1)→1; (σ,ψ)→1; (σ,σ)→2 (charges 1 and ψ).
        # total = 10
        @test KC.C_dim(N, 3, 2) == 10
    end
end

@testset "KitaevComplex — identity χⁿ⁺¹δⁿ + δⁿ⁻¹χⁿ = id" begin
    @testset "Fibonacci n = 0" begin
        N = fib_Nijk()
        d, D2 = fib_dims()
        err, M = KC.verify_homotopy(N, d, D2, 0)
        @info "Fibonacci n=0" err dim=size(M) M=M
        @test err < 1e-10
    end

    @testset "Fibonacci n = 1" begin
        N = fib_Nijk()
        d, D2 = fib_dims()
        err, M = KC.verify_homotopy(N, d, D2, 1)
        @info "Fibonacci n=1" err dim=size(M) M=M
        @test err < 1e-10
    end

    @testset "Ising n = 0" begin
        N = ising_Nijk()
        d, D2 = ising_dims()
        err, M = KC.verify_homotopy(N, d, D2, 0)
        @info "Ising n=0" err dim=size(M) M=M
        @test err < 1e-10
    end

    @testset "Ising n = 1" begin
        N = ising_Nijk()
        d, D2 = ising_dims()
        err, M = KC.verify_homotopy(N, d, D2, 1)
        @info "Ising n=1" err dim=size(M) M=M
        @test err < 1e-10
    end
end
