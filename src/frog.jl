
"""
    module FROG

A collection of reproducibility checks for constraint-based metabolic models
together with a report generator and tester.

See https://www.ebi.ac.uk/biomodels/curation/fbc for details.
"""
module FROG

using DocStringExtensions
using SHA, MD5, JSON, SBML, DelimitedFiles, Test, COBREXA, Distributed

include("frog/structs.jl")
include("frog/report.jl")
include("frog/io.jl")
include("frog/frontend.jl")

end
