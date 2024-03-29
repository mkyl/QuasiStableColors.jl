using QuasiStableColors: Optimize

using Test
using Tulip

"""
Return A, b, c in the Stigler diet problem

data from https://developers.google.com/optimization/lp/stigler_diet
"""
function stigler_diet()
    # 1939 price (cents), Calories (kcal), Protein (g),
    # Calcium (g), Iron (mg), Vitamin A (KIU), Vitamin B1 (mg), Vitamin B2 (mg),
    # Niacin (mg), Vitamin C (mg)
    data = Float64[36 44.7 1411 2 365 0 55.4 33.3 441 0
        14.1 11.6 418 0.7 54 0 3.2 1.9 68 0
        24.2 11.8 377 14.4 175 0 14.4 8.8 114 0
        7.1 11.4 252 0.1 56 0 13.5 2.3 68 0
        4.6 36.0 897 1.7 99 30.9 17.4 7.9 106 0
        8.5 28.6 680 0.8 80 0 10.6 1.6 110 0
        7.5 21.2 460 0.6 41 0 2 4.8 60 0
        7.1 25.3 907 5.1 341 0 37.1 8.9 64 0
        7.9 15.0 488 2.5 115 0 13.8 8.5 126 0
        9.1 12.2 484 2.7 125 0 13.9 6.4 160 0
        9.1 12.4 439 1.1 82 0 9.9 3 66 0
        24.8 8.0 130 0.4 31 18.9 2.8 3 17 0
        15.1 12.5 288 0.5 50 0 0 0 0 0
        11 6.1 310 10.5 18 16.8 4 16 7 177
        6.7 8.4 422 15.1 9 26 3 23.5 11 60
        30.8 10.8 9 0.2 3 44.2 0 0.2 2 0
        16.1 20.6 17 0.6 6 55.8 0.2 0 0 0
        32.6 2.9 238 1.0 52 18.6 2.8 6.5 1 0
        24.2 7.4 448 16.4 19 28.1 0.8 10.3 4 0
        14.1 3.5 49 1.7 3 16.9 0.6 2.5 0 17
        17.9 15.7 661 1.0 48 0 9.6 8.1 471 0
        16.7 8.6 18 0.2 8 2.7 0.4 0.5 0 0
        20.3 20.1 0 0 0 0 0 0 0 0
        9.8 41.7 0 0 0 0.2 0 0.5 5 0
        39.6 2.9 166 0.1 34 0.2 2.1 2.9 69 0
        36.4 2.2 214 0.1 32 0.4 2.5 2.4 87 0
        29.2 3.4 213 0.1 33 0 0 2 0 0
        22.6 3.6 309 0.2 46 0.4 1 4 120 0
        14.6 8.5 404 0.2 62 0 0.9 0 0 0
        26.8 2.2 333 0.2 139 169.2 6.4 50.8 316 525
        27.6 3.1 245 0.1 20 0 2.8 3.9 86 0
        36.6 3.3 140 0.1 15 0 1.7 2.7 54 0
        30.7 3.5 196 0.2 30 0 17.4 2.7 60 0
        24.2 4.4 249 0.3 37 0 18.2 3.6 79 0
        25.6 10.4 152 0.2 23 0 1.8 1.8 71 0
        27.4 6.7 212 0.2 31 0 9.9 3.3 50 0
        16 18.8 164 0.1 26 0 1.4 1.8 0 0
        30.3 1.8 184 0.1 30 0.1 0.9 1.8 68 46
        42.3 1.7 156 0.1 24 0 1.4 2.4 57 0
        13 5.8 705 6.8 45 3.5 1 4.9 209 0
        4.4 5.8 27 0.5 36 7.3 3.6 2.7 5 544
        6.1 4.9 60 0.4 30 17.4 2.5 3.5 28 498
        26 1.0 21 0.5 14 0 0.5 0 4 952
        30.9 2.2 40 1.1 18 11.1 3.6 1.3 10 1998
        7.1 2.4 138 3.7 80 69 4.3 5.8 37 862
        3.7 2.6 125 4.0 36 7.2 9 4.5 26 5369
        4.7 2.7 73 2.8 43 188.5 6.1 4.3 89 608
        7.3 0.9 51 3.0 23 0.9 1.4 1.4 9 313
        8.2 0.4 27 1.1 22 112.4 1.8 3.4 11 449
        3.6 5.8 166 3.8 59 16.6 4.7 5.9 21 1184
        34 14.3 336 1.8 118 6.7 29.4 7.1 198 2522
        8.1 1.1 106 0 138 918.4 5.7 13.8 33 2755
        5.1 9.6 138 2.7 54 290.7 8.4 5.4 83 1912
        16.8 3.7 20 0.4 10 21.5 0.5 1 31 196
        20.4 3.0 8 0.3 8 0.8 0.8 0.8 5 81
        21.3 2.4 16 0.4 8 2 2.8 0.8 7 399
        27.7 0.4 33 0.3 12 16.3 1.4 2.1 17 272
        10 1.0 54 2 65 53.9 1.6 4.3 32 431
        7.1 7.5 364 4 134 3.5 8.3 7.7 56 0
        10.4 5.2 136 0.2 16 12 1.6 2.7 42 218
        13.8 2.3 136 0.6 45 34.9 4.9 2.5 37 370
        8.6 1.3 63 0.7 38 53.2 3.4 2.5 36 1253
        7.6 1.6 71 0.6 43 57.9 3.5 2.4 67 862
        15.7 8.5 87 1.7 173 86.8 1.2 4.3 55 57
        9 12.8 99 2.5 154 85.7 3.9 4.3 65 257
        9.4 13.5 104 2.5 136 4.5 6.3 1.4 24 136
        7.9 20.0 1367 4.2 345 2.9 28.7 18.4 162 0
        8.9 17.4 1055 3.7 459 5.1 26.9 38.2 93 0
        5.9 26.9 1691 11.4 792 0 38.4 24.6 217 0
        22.4 0 0 0 0 0 4 5.1 50 0
        17.4 0 0 0 0 0 0 2.3 42 0
        8.6 8.7 237 3 72 0 2 11.9 40 0
        16.2 8.0 77 1.3 39 0 0.9 3.4 14 0
        51.7 34.9 0 0 0 0 0 0 0 0
        13.7 14.7 0 0.5 74 0 0 0 5 0
        13.6 9.0 0 10.3 244 0 1.9 7.5 146 0
        20.5 6.4 11 0.4 7 0.2 0.2 0.4 3 0]

    # Calories (kcal), Protein (g),
    # Calcium (g), Iron (mg), Vitamin A (KIU), Vitamin B1 (mg), Vitamin B2 (mg),
    # Niacin (mg), Vitamin C (mg)
    requirements = Float64[3, 70, 0.8, 12, 5, 1.8, 2.7, 18, 75]

    A = data[:, 2:10]
    b = requirements
    c = data[:, 1]
    transpose(A), b, c
end

@testset "lifted linear optimization" begin
    @testset "stable coloring, stigler diet" begin
        model = Tulip.Optimizer()

        A, b, c = stigler_diet()
        V1 = Optimize.minimize(model, A, b, c)
        V2 = Optimize.lifted_minimize(model, A, b, c; q=0.0)
        @test c' * V1 ≈ V2
    end

    @testset "example from paper" begin
        model = Tulip.Optimizer()

        A::Matrix{Float64} = [4 8 2; 6 5 1; 7 4 2; 3 1 22; 2 3 21]
        b::Vector{Float64} = [20, 20, 21, 50, 51]
        c::Vector{Float64} = [9, 10, 50]

        z = c' * Optimize.maximize(model, A, b, c)
        @test z ≈ 128.157 atol = 0.001
        z₀ = Optimize.lifted_maximize(model, A, b, c; q=1.0)
        # meets or beats reported error
        @test abs(z₀ - z) <= abs(z - 130.199)
    end
end
