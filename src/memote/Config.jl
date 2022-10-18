"""
$(TYPEDEF)

Module housing the configuration parameters for the memote-style tests.
"""
module Config

using ..DocStringExtensions
using ..COBREXA


"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
mutable struct AnnotationConfig
    gene_annotation_keywords::Vector{String}
    gene_annotation_regexes::Dict{String,Regex}
    metabolite_annotation_keywords::Vector{String}
    metabolite_annotation_regexes::Dict{String,Regex}
    reaction_annotation_keywords::Vector{String}
    reaction_annotation_regexes::Dict{String,Regex}
    minimum_fraction_database_conformity::Float64
    minimum_fraction_database_annotations::Float64
end

annotation_config = AnnotationConfig(
    [
        "kegg.genes",
        "refseq",
        "uniprot",
        "ecogene",
        "ncbigi",
        "ncbigene",
        "ncbiprotein",
        "ccds",
        "hprd",
        "asap",
    ],
    Dict(
        "refseq" =>
            r"^((AC|AP|NC|NG|NM|NP|NR|NT|NW|XM|XP|XR|YP|ZP)_\d+|(NZ\_[A-Z]{4}\d+))(\.\d+)?$",
        "uniprot" =>
            r"^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$",
        "ecogene" => r"^EG\d+$",
        "kegg.genes" => r"^\w+:[\w\d\.-]*$",
        "ncbigi" => r"^(GI|gi)\:\d+$",
        "ncbigene" => r"^\d+$",
        "ncbiprotein" => r"^(\w+\d+(\.\d+)?)|(NP_\d+)$",
        "ccds" => r"^CCDS\d+\.\d+$",
        "hprd" => r"^\d+$",
        "asap" => r"^[A-Za-z0-9-]+$",
    ),
    [
        "pubchem.compound",
        "kegg.compound",
        "seed.compound",
        "inchi_key",
        "inchi",
        "chebi",
        "hmdb",
        "reactome.compound",
        "metanetx.chemical",
        "bigg.metabolite",
        "biocyc",
    ],
    Dict(
        "pubchem.compound" => r"^\d+$",
        "kegg.compound" => r"^C\d+$",
        "seed.compound" => r"^cpd\d+$",
        "inchi_key" => r"^[A-Z]{14}\-[A-Z]{10}(\-[A-Z])?",
        "inchi" =>
            r"^InChI\=1S?\/[A-Za-z0-9\.]+(\+[0-9]+)?(\/[cnpqbtmsih][A-Za-z0-9\-\+\(\)\,\/\?\;\.]+)*$",
        "chebi" => r"^CHEBI:\d+$",
        "hmdb" => r"^HMDB\d{5}$",
        "reactome.compound" => r"(^R-[A-Z]{3}-[0-9]+(-[0-9]+)?$)|(^REACT_\d+(\.\d+)?$)",
        "metanetx.chemical" => r"^MNXM\d+$",
        "bigg.metabolite" => r"^[a-z_A-Z0-9]+$",
        "biocyc" => r"^[A-Z-0-9]+(?<!CHEBI)(\:)?[A-Za-z0-9+_.%-]+$",
    ),
    [
        "rhea",
        "kegg.reaction",
        "seed.reaction",
        "metanetx.reaction",
        "bigg.reaction",
        "reactome",
        "ec-code",
        "brenda",
        "biocyc",
    ],
    Dict(
        "rhea" => r"^\d{5}$",
        "kegg.reaction" => r"^R\d+$",
        "seed.reaction" => r"^rxn\d+$",
        "metanetx.reaction" => r"^MNXR\d+$",
        "bigg.reaction" => r"^[a-z_A-Z0-9]+$",
        "reactome" => r"(^R-[A-Z]{3}-[0-9]+(-[0-9]+)?$)|(^REACT_\d+(\.\d+)?$)",
        "ec-code" =>
            r"^\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$",
        "brenda" =>
            r"^\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$",
        "biocyc" => r"^[A-Z-0-9]+(?<!CHEBI)(\:)?[A-Za-z0-9+_.%-]+$",
    ),
    0.9,
    0.9,
)

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
mutable struct BasicConfig
    minimum_metabolic_coverage::Float64
    minimum_growth_rate::Float64
    maximum_growth_rate::Float64
    optimizer_modifications::Vector{Function}
end

basic_config = BasicConfig(
    0.1,
    0.01,
    5.0,
    Function[],
)

"""
$(TYPEDEF)

Parameters used by the biomass tests.

# Fields
$(TYPEDFIELDS)
"""
mutable struct BiomassConfig
    biomass_strings::Vector{String}
    growth_metabolites::Dict{String,String}
    ignored_precursors::Vector{String}
    essential_precursors::Dict{String,String}
    optimizer_modifications::Vector{Function}
end

biomass_config = BiomassConfig(
    ["BIOMASS", "biomass", "Biomass"],
    Dict("atp" => "atp_c", "adp" => "adp_c", "h2o" => "h2o_c", "pi" => "pi_c"),
    ["atp_c", "h2o_c"],
    Dict(
        "trp__L" => "trp__L_c",
        "cys__L" => "cys__L_c",
        "his__L" => "his__L_c",
        "tyr__L" => "tyr__L_c",
        "met__L" => "met__L_c",
        "phe__L" => "phe__L_c",
        "ser__L" => "ser__L_c",
        "pro__L" => "pro__L_c",
        "asp__L" => "asp__L_c",
        "thr__L" => "thr__L_c",
        "gln__L" => "gln__L_c",
        "glu__L" => "glu__L_c",
        "ile__L" => "ile__L_c",
        "arg__L" => "arg__L_c",
        "lys__L" => "lys__L_c",
        "val__L" => "val__L_c",
        "leu__L" => "leu__L_c",
        "ala__L" => "ala__L_c",
        "gly" => "gly_c",
        "asn__L" => "asn__L_c",
        "datp" => "datp_c",
        "dctp" => "dctp_c",
        "dttp" => "dttp_c",
        "dgtp" => "dgtp_c",
        "atp" => "atp_c",
        "ctp" => "ctp_c",
        "utp" => "utp_c",
        "gtp" => "gtp_c",
        "nad" => "nad_c",
        "nadp" => "nadp_c",
        "amet" => "amet_c",
        "fad" => "fad_c",
        "pydx5p" => "pydx5p_c",
        "coa" => "coa_c",
        "thmpp" => "thmpp_c",
        "fmn" => "fmn_c",
        "h2o" => "h2o_c",
    ),
    Function[],
)

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
end

consistency_config = ConsistencyConfig(
    String[],
    String[],
    String[],
)

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
mutable struct EnergyConfig
    atpm_strings::Vector{String}
    energy_dissipating_metabolites::Dict{String,String}
    additional_energy_generating_reactions::Vector{Reaction}
    ignored_energy_reactions::Vector{String}
    optimizer_modifications::Vector{Function}
end

energy_config = EnergyConfig(
    ["ATPM", "Maintenance", "maintenance"],
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

# Fields
$(TYPEDFIELDS)
"""
mutable struct GPRAssociationConfig
    test_annotation :: String
end

gpra_config = GPRAssociationConfig("inchi_key")

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

# Fields
$(TYPEDFIELDS)
"""
mutable struct NetworkConfig
    condition_number::Float64
    fva_bound::Float64
    cycle_tol::Float64
    minimum_metabolite_flux::Float64
    optimizer_modifications::Vector{Function}
end

network_config = NetworkConfig(1e9, 0.01, 1e-3, 1e-3, Function[])

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
mutable struct ReactionConfig
    test_annotation::String
    ignore_annotations::Vector{String}
    bound_default::Float64
end

reaction_config = ReactionConfig("inchi_key", ["sbo", "ec-code"], 1000.0)

"""
$(TYPEDEF)

A grouping of parameters used by the metabolic testing infrastructure.

# Fields
$(TYPEDFIELDS)
"""
mutable struct MemoteConfig
    annotation::AnnotationConfig
    basic::BasicConfig
    biomass::BiomassConfig
    consistency::ConsistencyConfig
    energy::EnergyConfig
    gpra::GPRAssociationConfig
    metabolite::MetaboliteConfig
    network::NetworkConfig
    reaction::ReactionConfig
end

memote_config = MemoteConfig(
    annotation_config,
    basic_config,
    biomass_config,
    consistency_config,
    energy_config,
    gpra_config,
    metabolite_config,
    network_config,
    reaction_config,
)


end # module