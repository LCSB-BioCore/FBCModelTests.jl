using Test

using FBCModelTests

using Downloads, SHA
using COBREXA, Tulip

# this loads some data
include("data.jl")

@testset "FBCModelTests test suite" begin
    include("memote.jl")
end
