module FBCModelTests

import Distributed
import InteractiveUtils
import Pkg

using COBREXA
using JuMP
using DelimitedFiles
using DocStringExtensions
using JSON
using MD5
using SBML
using SHA
using Test
using PeriodicTable
using SparseArrays

include("version.jl")
include("common.jl")
include("structs.jl")
include("frog.jl")
include(joinpath("memote", "metabolites.jl"))
include(joinpath("memote", "reactions.jl"))
include(joinpath("memote", "basic.jl"))
include(joinpath("memote", "consistency.jl"))
include(joinpath("memote", "gpr_associations.jl"))
include(joinpath("memote", "biomass.jl"))
include(joinpath("memote", "genes.jl"))
include(joinpath("memote", "network.jl"))
include(joinpath("memote", "memote.jl")) #  the test harness

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end
