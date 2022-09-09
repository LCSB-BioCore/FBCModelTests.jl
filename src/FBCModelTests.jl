module FBCModelTests

using COBREXA
using DelimitedFiles
using DocStringExtensions
using JSON
using MD5
import Distributed
import InteractiveUtils
import Pkg
using JuMP

include("version.jl")
include("frog.jl")
include("memote.jl")

export generate_frog_report,
    is_model_charge_balanced,
    is_model_mass_balanced,
    has_erroneous_energy_generating_cycles,
    is_consistent

end
