
"""
    module FROG

A collection of reproducibility checks for constraint-based metabolic models
together with a report generator and tester.

See https://www.ebi.ac.uk/biomodels/curation/fbc for details.
"""
module FROG

using COBREXA
using DelimitedFiles
using Distributed
using DocStringExtensions
using JSON
using SBML
using Test

using ..FBCModelTests: Maybe

include("frog/structs.jl")
include("frog/report.jl")
include("frog/io.jl")
include("frog/frontend.jl")

end
