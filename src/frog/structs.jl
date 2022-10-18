"""
$(TYPEDEF)
"""
const ObjectiveValue = Maybe{Float64}

"""
$(TYPEDEF)
"""
const FROGMetadata = Dict{String,String}

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct FROGReactionReport
    flux::ObjectiveValue
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
