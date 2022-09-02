using Test

using FBCModelTests

using Downloads
using COBREXA

isfile("e_coli_core.json") ||
    Downloads.download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")
model = load_model("e_coli_core.json")

@testset "Test tests" begin
    include("memote.jl")
end
