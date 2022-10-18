module FBCModelTests

import Pkg

using COBREXA
using JuMP
using DocStringExtensions
using Test
using PeriodicTable
using SparseArrays

include("version.jl")
include("common.jl")
include("frog.jl")
include("memote.jl")

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end
