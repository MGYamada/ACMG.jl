"""
Tests for `ACMG.Phase4.KitaevComplex`.

Covers:
  * `C_basis` / `C_dim` sizing for Fibonacci and Ising
  * `δ⁰ = 0`                                             (Kitaev p. 92 remark)
  * `δ¹` matrix entries against the explicit Eq. 244
    `(δ¹X)_{a₁,a₂;c₂} = X_{a₂} − X_{c₂} + X_{a₁}`
  * `χ¹` formula  `(χ¹X) = (1/𝒟²) Σ_a d_a² X_{a,a}`
  * `χ²` formula  `(χ²X)_a = (1/𝒟²) Σ (d_{a₁}d_{c₂}/d_a) X_{a₁,a;c₂}`
  * Homotopy identity `χⁿ⁺¹ δⁿ + δⁿ⁻¹ χⁿ = id_Cⁿ` at `n = 1`
    (Kitaev Eq. 251)

Implementation is scalar-level (multiplicity-free); `n ≥ 2` operations
error out pending F-symbol data — that's asserted too.
"""

using Test
using LinearAlgebra
using ACMG
using ACMG.Phase4

const KC = Phase4.KitaevComplex

# ----- fusion data --------------------------------------------------------

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
    return ([1.0, φ], 1 + φ^2)
end

ising_dims() = ([1.0, 1.0, sqrt(2.0)], 4.0)

# ----- basis sizing -------------------------------------------------------

@testset "KitaevComplex — Cⁿ basis sizing" begin
    N = fib_Nijk()
    @test KC.C_dim(N, 2, 0) == 1
    @test KC.C_dim(N, 2, 1) == 2
    @test KC.C_dim(N, 2, 2) == 5

    N = ising_Nijk()
    @test KC.C_dim(N, 3, 0) == 1
    @test KC.C_dim(N, 3, 1) == 3
    @test KC.C_dim(N, 3, 2) == 10
end

@testset "KitaevComplex — Cⁿ basis content (Fibonacci)" begin
    N = fib_Nijk()
    @test KC.C_basis(N, 2, 0) == [((),)]
    @test Set(KC.C_basis(N, 2, 1)) == Set([(1,1), (2,2)])
    @test Set(KC.C_basis(N, 2, 2)) == Set([
        (1,1,1), (1,2,2), (2,1,2), (2,2,1), (2,2,2)
    ])
end

# ----- δ⁰ = 0 -------------------------------------------------------------

@testset "KitaevComplex — δ⁰ = 0" begin
    for (label, N, dD) in [("Fib", fib_Nijk(), fib_dims()),
                           ("Ising", ising_Nijk(), ising_dims())]
        d, D2 = dD
        δ0 = KC.delta_matrix(N, d, D2, 0)
        @test maximum(abs.(δ0)) < 1e-12
    end
end

# ----- δ¹ explicit entry check (Eq. 244) ---------------------------------

@testset "KitaevComplex — δ¹ matches Eq. 244 on arbitrary X" begin
    for (label, N, dD) in [("Fib", fib_Nijk(), fib_dims()),
                           ("Ising", ising_Nijk(), ising_dims())]
        d, D2 = dD
        r   = size(N, 1)
        δ1  = KC.delta_matrix(N, d, D2, 1)
        C1  = KC.C_basis(N, r, 1)
        C2  = KC.C_basis(N, r, 2)

        X_coef_full = Float64[0.7, -0.3, 1.9, -2.1]
        X_coef      = X_coef_full[1:length(C1)]
        y = δ1 * X_coef

        for (i, out_t) in enumerate(C2)
            a1, a2, c2 = out_t
            idx(a) = findfirst(==((a, a)), C1)
            expected = X_coef[idx(a2)] - X_coef[idx(c2)] + X_coef[idx(a1)]
            @test isapprox(y[i], expected; atol = 1e-12)
        end
    end
end

# ----- χ¹ and χ² concrete formulas ---------------------------------------

@testset "KitaevComplex — χ¹ formula" begin
    for (label, N, dD) in [("Fib", fib_Nijk(), fib_dims()),
                           ("Ising", ising_Nijk(), ising_dims())]
        d, D2 = dD
        r     = size(N, 1)
        χ1    = KC.chi_matrix(N, d, D2, 1)
        C1    = KC.C_basis(N, r, 1)

        # χ¹ is a 1 × |C¹| row matrix  [d_a² / 𝒟²] indexed by C¹
        @test size(χ1) == (1, length(C1))
        for (j, (a, _)) in enumerate(C1)
            @test isapprox(χ1[1, j], d[a]^2 / D2; atol = 1e-12)
        end

        # Sanity: χ¹ · (all-ones C¹ vector) = Σ d_a²/𝒟² = 1
        ones_C1 = ones(length(C1))
        @test isapprox((χ1 * ones_C1)[1], 1.0; atol = 1e-12)
    end
end

@testset "KitaevComplex — χ² formula" begin
    for (label, N, dD) in [("Fib", fib_Nijk(), fib_dims()),
                           ("Ising", ising_Nijk(), ising_dims())]
        d, D2 = dD
        r     = size(N, 1)
        χ2    = KC.chi_matrix(N, d, D2, 2)
        C1    = KC.C_basis(N, r, 1)
        C2    = KC.C_basis(N, r, 2)

        @test size(χ2) == (length(C1), length(C2))

        # Entry-level check against (d_{a₁} d_{c₂}) / (𝒟² d_{a₂}) · δ_{output, (a₂,a₂)}
        for (i, out_t) in enumerate(C1), (j, in_t) in enumerate(C2)
            a1, a2, c2 = in_t
            ao = out_t[1]
            expected = (ao == a2) ? d[a1] * d[c2] / (D2 * d[a2]) : 0.0
            @test isapprox(χ2[i, j], expected; atol = 1e-12)
        end
    end
end

# ----- Eq. 251 homotopy identity at n = 1 --------------------------------

@testset "KitaevComplex — Eq. 251 at n = 1" begin
    for (label, N, dD) in [("Fib", fib_Nijk(), fib_dims()),
                           ("Ising", ising_Nijk(), ising_dims())]
        d, D2 = dD
        err, M = KC.verify_homotopy(N, d, D2, 1)
        @test err < 1e-12

        # Since δ⁰ = 0, the identity reduces to χ² δ¹ = id_{C¹}.
        δ1 = KC.delta_matrix(N, d, D2, 1)
        χ2 = KC.chi_matrix(N, d, D2, 2)
        r  = size(N, 1)
        Id = Matrix{Float64}(I, KC.C_dim(N, r, 1), KC.C_dim(N, r, 1))
        @test maximum(abs.(χ2 * δ1 - Id)) < 1e-12
    end
end

# ----- Fibonacci — hand-calculated values --------------------------------

@testset "KitaevComplex — Fibonacci hand calculation" begin
    N = fib_Nijk()
    d, D2 = fib_dims()
    φ = d[2]

    # 𝒟² = 1 + φ² = 2 + φ  (since φ² = φ + 1)
    @test isapprox(D2, 2 + φ; atol = 1e-12)

    # (χ¹)_{(1,1)} = 1/𝒟²,   (χ¹)_{(τ,τ)} = φ²/𝒟²
    χ1 = KC.chi_matrix(N, d, D2, 1)
    C1 = KC.C_basis(N, 2, 1)
    i11 = findfirst(==((1,1)), C1); iττ = findfirst(==((2,2)), C1)
    @test isapprox(χ1[1, i11], 1.0 / D2; atol = 1e-12)
    @test isapprox(χ1[1, iττ], φ^2 / D2; atol = 1e-12)

    # (χ² X)_τ Eq. 250 derivation: summed over (a₁, a, c₂) with a = τ,
    # i.e. (1,τ,τ), (τ,τ,1), (τ,τ,τ).  Weights  (d_{a₁} d_{c₂}) / (𝒟² d_τ):
    #   (1·φ)/(𝒟²·φ) = 1/𝒟²    for (1,τ,τ)
    #   (φ·1)/(𝒟²·φ) = 1/𝒟²    for (τ,τ,1)
    #   (φ·φ)/(𝒟²·φ) = φ/𝒟²    for (τ,τ,τ)
    χ2 = KC.chi_matrix(N, d, D2, 2)
    C2 = KC.C_basis(N, 2, 2)
    i_out_τ = findfirst(==((2,2)), C1)
    for (in_t, expected_w) in [((1,2,2), 1.0 / D2),
                                ((2,2,1), 1.0 / D2),
                                ((2,2,2), φ   / D2)]
        j = findfirst(==(in_t), C2)
        @test isapprox(χ2[i_out_τ, j], expected_w; atol = 1e-12)
    end
end

# ----- n ≥ 2 currently unimplemented -------------------------------------

@testset "KitaevComplex — higher-n operations guarded" begin
    N = fib_Nijk()
    d, D2 = fib_dims()
    @test_throws ErrorException KC.delta_matrix(N, d, D2, 2)
    @test_throws ErrorException KC.chi_matrix(N, d, D2, 3)
end
