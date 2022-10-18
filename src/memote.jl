"""
module Memote

This package contains a collection of tests based on Memote. See Lieven, C.,
Beber, M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic
model testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
"""
module Memote

using DocStringExtensions
using COBREXA
using Test
using JuMP

include(joinpath("memote", "Utils.jl")) # memote utils
include(joinpath("memote", "Config.jl")) # memote test parameters
include.(joinpath.(Ref("memote/checks"), ["Annotation.jl", "Basic.jl", "Biomass.jl", "Consistency.jl", "Energy.jl", "GPRAssociation.jl", "Metabolite.jl", "Network.jl", "Reaction.jl"]))

include(joinpath("memote", "frontend.jl"))

end # module