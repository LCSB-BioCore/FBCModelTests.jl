using Test

using FBCModelTests
using FBCModelTests.Memote
using FBCModelTests.Memote.Annotation
using FBCModelTests.Memote.Basic
using FBCModelTests.Memote.Biomass
using FBCModelTests.Memote.Consistency
using FBCModelTests.Memote.Energy
using FBCModelTests.Memote.GPRAssociation
using FBCModelTests.Memote.Metabolite
using FBCModelTests.Memote.Network
using FBCModelTests.Memote.Reaction

memote_config = FBCModelTests.Config.memote_config

using Downloads
using COBREXA, Tulip

# tried and trusted E. coli core
isfile("e_coli_core.json") || Downloads.download(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
)
model = load_model("e_coli_core.json")

# this model has an energy generating cycle
isfile("iJN746.json") ||
    Downloads.download("http://bigg.ucsd.edu/static/models/iJN746.json", "iJN746.json")
iJN746 = load_model("iJN746.json")

@testset "Test tests" begin
    include("memote.jl")
end
