"""
$(TYPEDEF)
"""
const ObjectiveValue = Maybe{Float64}

"""
$(TYPEDEF)
"""
const FROGMetadata = Dict{String,Any}

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct FROGReactionReport
    objective_flux::ObjectiveValue
    fraction_optimum::Float64
    variability_min::ObjectiveValue
    variability_max::ObjectiveValue
    deletion::ObjectiveValue
end

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct FROGObjectiveReport
    optimum::ObjectiveValue
    reactions::Dict{String,FROGReactionReport}
    gene_deletions::Dict{String,ObjectiveValue}
end

"""
$(TYPEDEF)
"""
const FROGReportData = Dict{String,FROGObjectiveReport}
