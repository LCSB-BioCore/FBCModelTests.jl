using Test

using FBCModelTests
using FBCModelTests.Memote
using FBCModelTests.Memote.Annotation
using FBCModelTests.Memote.Basic
using FBCModelTests.Memote.Biomass
using FBCModelTests.Memote.Consistency
using FBCModelTests.Memote.Energy
using FBCModelTests.Memote.GPRAssociation
using FBCModelTests.Memote.Network

memote_config = FBCModelTests.Memote.Config.memote_config

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
