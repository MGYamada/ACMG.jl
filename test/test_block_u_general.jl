using Test
using ACMG

"""
Tests for general O(n) Cayley parametrisation in BlockU.jl.
"""

@testset "General O(n) Cayley" begin

    @testset "inverse_mod_p" begin
        # 2x2 example
        p = 7
        M = [3 1; 2 4]
        # det = 12 - 2 = 10 = 3 mod 7, det_inv = 5
        # adj = [[4, -1], [-2, 3]] = [[4, 6], [5, 3]] mod 7
        # M^{-1} = 5 * adj = [[20, 30], [25, 15]] = [[6, 2], [4, 1]] mod 7
        Minv = inverse_mod_p(M, p)
        @test Minv !== nothing
        # Verify: M * Minv ≡ I
        prod_M = [mod(sum(M[i, k] * Minv[k, j] for k in 1:2), p) for i in 1:2, j in 1:2]
        @test prod_M == [1 0; 0 1]

        # 3x3 identity
        I3 = [1 0 0; 0 1 0; 0 0 1]
        @test inverse_mod_p(I3, 11) == I3

        # Singular matrix
        sing = [1 2; 2 4]
        @test inverse_mod_p(sing, 7) === nothing
    end

    @testset "cayley_so_n at n=2 reproduces O(2) circle points" begin
        # Cayley with single param a → SO(2) rotation:
        # A = [[0, a], [-a, 0]], (I-A)(I+A)^{-1} = (1/(1+a²)) [[1-a², -2a], [2a, 1-a²]]
        # cos = (1-a²)/(1+a²), sin = 2a/(1+a²) (Weierstrass)
        p = 13
        for a in 0:(p-1)
            U = cayley_so_n([a], 2, p)
            U === nothing && continue  # 1+a² ≡ 0 mod p (singular)
            # U should be a 2x2 SO(2) matrix
            @test size(U) == (2, 2)
            @test mod(det_2x2(U), p) == 1  # det = +1
            # U U^T = I
            UUT = [mod(sum(U[i, k] * U[j, k] for k in 1:2), p) for i in 1:2, j in 1:2]
            @test UUT == [1 0; 0 1]
        end
    end

    @testset "enumerate_so_n_Fp counts" begin
        # SO(1) = {1}, |SO(1)(F_p)| = 1
        @test length(enumerate_so_n_Fp(1, 7)) == 1

        # |SO(2)(F_p)|: rotation matrices form circle, |circle| = p ± 1
        # Cayley misses the singular point(s), so we get p - (#singular)
        # for p=13: 1 + a² = 0 has 2 solutions if -1 is QR (p ≡ 1 mod 4) → 13-2 = 11
        # for p=11: 1 + a² = 0 has 0 solutions (p ≡ 3 mod 4) → 11
        so2_p13 = enumerate_so_n_Fp(2, 13)
        # 13 ≡ 1 (mod 4), so -1 is QR; we lose 2 Cayley-singular params
        @test length(so2_p13) == 11

        so2_p11 = enumerate_so_n_Fp(2, 11)
        # 11 ≡ 3 (mod 4), -1 is not QR, no Cayley-singular params
        @test length(so2_p11) == 11

        # All elements should be valid SO(2): det = 1, U U^T = I
        for U in so2_p13
            @test mod(det_2x2(U), 13) == 1
            UUT = [mod(sum(U[i, k] * U[j, k] for k in 1:2), 13) for i in 1:2, j in 1:2]
            @test UUT == [1 0; 0 1]
        end
    end

    @testset "enumerate_o_n_Fp size" begin
        p = 13
        o2 = enumerate_o_n_Fp(2, p)
        so2 = enumerate_so_n_Fp(2, p)
        @test length(o2) == 2 * length(so2)
        # Half should have det = +1, other half det = -1
        det_pos = count(U -> mod(det_2x2(U), p) == 1, o2)
        det_neg = count(U -> mod(det_2x2(U), p) == p - 1, o2)
        @test det_pos == length(so2)
        @test det_neg == length(so2)
    end

    @testset "apply_block_U sanity" begin
        # 4x4 S, apply 2x2 block at (1, 2)
        p = 13
        S = [1 2 3 4; 2 5 6 7; 3 6 8 9; 4 7 9 10]
        # Identity block doesn't change S
        U_id = [1 0; 0 1]
        @test apply_block_U(S, [1, 2], U_id, p) == S

        # Test on 3x3 block at (1, 3, 4)
        S5 = [10 11 12 13 14; 11 15 16 17 18; 12 16 19 20 21;
              13 17 20 22 23; 14 18 21 23 24]
        U_id3 = [1 0 0; 0 1 0; 0 0 1]
        @test apply_block_U(S5, [1, 3, 4], U_id3, p) == S5
    end
end

# Helper: 2x2 determinant
det_2x2(M) = M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]
