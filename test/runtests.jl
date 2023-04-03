using Test

using FBCModelTests

using Downloads, SHA
using COBREXA

# test tolerance
const TEST_TOL = 1e-6

# this loads some data
include("data.jl")

# this loads CountTests "mock" tester for testing the tests
include("testcounter.jl")

@testset "FBCModelTests test suite" begin
    @testset "FROG" begin
        include("frog.jl")
    end
    @testset "MEMOTE" begin
        include("memote.jl")
    end
end
