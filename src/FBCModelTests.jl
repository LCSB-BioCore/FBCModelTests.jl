"""
    module FBCModelTests

A collection of tests for constraint-based flux metabolic models.
"""
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

end
