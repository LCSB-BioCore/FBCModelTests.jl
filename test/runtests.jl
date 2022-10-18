using Test

using FBCModelTests

using Downloads, SHA
using COBREXA, Tulip

# this loads some data
include("data.jl")

# this loads CountTests "mock" tester for testing the tests
include("testcounter.jl")

@testset "FBCModelTests test suite" begin
    include("frog.jl")
    include("memote.jl")
end
