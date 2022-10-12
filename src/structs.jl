
"""
$(TYPEDEF)
"""
const ObjectiveValue = Maybe{Float64}

"""
$(TYPEDEF)
"""
const ObjectiveValues = Vector{ObjectiveValue}

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct FROGObjectiveReport
    optimum::ObjectiveValue
    flux::ObjectiveValues = []
    variabilities_min::ObjectiveValues = []
    variabilities_max::ObjectiveValues = []
    reaction_knockouts::ObjectiveValues = []
    gene_knockouts::ObjectiveValues = []
end


"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct FROGReportData
    reactions::Vector{String}
    genes::Vector{String}
    objectives::Dict{String,FROGObjectiveReport}
end
