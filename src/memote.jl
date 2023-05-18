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
using Distributed

include(joinpath("common.jl")) # quiet tests
include(joinpath("memote", "utils.jl")) # memote utils
include(joinpath("memote", "config.jl")) # memote test parameters
include.(
    joinpath.(
        Ref("memote/checks"),
        [
            "Annotation.jl",
            "Basic.jl",
            "Biomass.jl",
            "Consistency.jl",
            "Energy.jl",
            "GPRAssociation.jl",
            "Metabolites.jl",
            "Network.jl",
            "Reactions.jl",
        ],
    )
)

include(joinpath("memote", "frontend.jl"))

end # module
