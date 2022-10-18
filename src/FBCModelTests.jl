#=
This package contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#
"""
# FBCModelTests v$(FBCModelTests.FBCMT_VERSION)

This package allows you to generate a FROG report, as well as run MeMoTe-style
tests on a constraint-based metabolic model.
"""
module FBCModelTests

import Pkg
using COBREXA
using JuMP
using DocStringExtensions
using Test
using SparseArrays

include("version.jl")
include("Utils.jl")
include("FROG.jl")
include("Memote.jl")

export FROG, Memote

end
