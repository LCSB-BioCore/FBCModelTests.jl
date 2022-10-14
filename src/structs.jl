
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

#=
Use these structs to set the default parameters used in the metabolic tests.
Each subtest has its own struct, which is grouped together by `MemoteConfig`.
=#

"""
$(TYPEDEF)

Parameters used by the metabolite tests.

# Fields
$(TYPEDFIELDS)
"""
mutable struct MetaboliteConfig
    formula_corner_cases::Vector{String}
    charge_corner_cases::Vector{Int64}
    medium_only_imported::Bool
    test_annotation::String
end

metabolite_config = MetaboliteConfig(["X", "x", ""], Int64[], true, "inchi_key")

"""
$(TYPEDEF)

Parameters used by the consistency tests.

# Fields
$(TYPEDFIELDS)
"""
mutable struct ConsistencyConfig
    mass_ignored_reactions::Vector{String}
    charge_ignored_reactions::Vector{String}
    consistency_ignored_reactions::Vector{String}
    energy_dissipating_metabolites::Dict{String,String}
    additional_energy_generating_reactions::Vector{Reaction}
    ignored_energy_reactions::Vector{String}
    optimizer_modifications::Vector{Function}
end

consistency_config = ConsistencyConfig(
    String[],
    String[],
    String[],
    Dict(
        "ATP" => "atp_c",
        "CTP" => "ctp_c",
        "GTP" => "gtp_c",
        "UTP" => "utp_c",
        "ITP" => "itp_c",
        "ADP" => "adp_c",
        "CDP" => "cdp_c",
        "GDP" => "gdp_c",
        "UDP" => "udp_c",
        "IDP" => "idp_c",
        "NADH" => "nadh_c",
        "NAD" => "nad_c",
        "NADPH" => "nadph_c",
        "NADP" => "nadp_c",
        "FADH2" => "fadh2_c",
        "FAD" => "fad_c",
        "FMNH2" => "fmn_c",
        "FMN" => "fmnh2_c",
        "Ubiquinol-8" => "q8h2_c",
        "Ubiquinone-8" => "q8_c",
        "Menaquinol-8" => "mql8_c",
        "Menaquinone-8" => "mqn8_c",
        "2-Demethylmenaquinol-8" => "2dmmql8_c",
        "2-Demethylmenaquinone-8" => "2dmmq8_c",
        "ACCOA" => "accoa_c",
        "COA" => "coa_c",
        "L-Glutamate" => "glu__L_c",
        "2-Oxoglutarate" => "akg_c",
        "Ammonium" => "nh4_c",
        "H" => "h_c",
        "H[external]" => "h_e",
        "H2O" => "h2o_c",
        "Phosphate" => "pi_c",
        "Acetate" => "ac_c",
    ),
    Reaction[],
    ["ATPM"],
    Function[],
)

"""
$(TYPEDEF)

Parameters used by the biomass tests.

# Fields
$(TYPEDFIELDS)
"""
mutable struct BiomassConfig
    biomass_terms :: Vector{String}
    atpm_strings :: Vector{String}
    growth_metabolites :: Dict{String, String}
end

biomass_config = BiomassConfig(
    ["BIOMASS", "biomass", "Biomass"],
    ["ATPM", "Maintenance", "maintenance"],
    Dict("atp" => "atp_c", "adp" => "adp_c", "h2o" => "h2o_c", "pi" => "pi_c")
)

"""
$(TYPEDEF)

A grouping of parameters used by the metabolic testing infrastructure.

# Fields
$(TYPEDFIELDS)
"""
mutable struct MemoteConfig
    metabolite::MetaboliteConfig
    consistency::ConsistencyConfig
    biomass::BiomassConfig
end

memote_config = MemoteConfig(metabolite_config, consistency_config, biomass_config)
