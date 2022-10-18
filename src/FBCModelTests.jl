"""
# FBCModelTests v$(FBCModelTests.FBCMT_VERSION)

This package allows you to generate a FROG report, as well as run MeMoTe-style
tests on a constraint-based metabolic model.
"""
module FBCModelTests

import Pkg
include("version.jl")

include("Utils.jl")
import .Utils

# include(joinpath("frog", "FROG.jl"))

include(joinpath("memote", "Memote.jl"))
import .Memote

end
