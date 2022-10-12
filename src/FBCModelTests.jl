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

include("version.jl")
include("common.jl")
include("structs.jl")
include("frog.jl")
include("memote.jl")

export frog_generate_report,
    is_model_charge_balanced,
    is_model_mass_balanced,
    has_erroneous_energy_generating_cycles,
    is_consistent

end
