"""
    module Config

Module housing the configuration parameters for the memote-style tests.
"""
module Config

using COBREXA
using DocStringExtensions

# Add all metabolite identifier constants. Use BiGG namespace.
const trp__L = "trp__L"
const cys__L = "cys__L"
const his__L = "his__L"
const tyr__L = "tyr__L"
const met__L = "met__L"
const phe__L = "phe__L"
const ser__L = "ser__L"
const pro__L = "pro__L"
const asp__L = "asp__L"
const thr__L = "thr__L"
const gln__L = "gln__L"
const glu__L = "glu__L"
const ile__L = "ile__L"
const arg__L = "arg__L"
const lys__L = "lys__L"
const val__L = "val__L"
const leu__L = "leu__L"
const ala__L = "ala__L"
const gly = "gly"
const asn__L = "asn__L"


const datp = "datp"
const dctp = "dctp"
const dttp = "dttp"
const dgtp = "dgtp"

const atp = "atp"
const ctp = "ctp"
const gtp = "gtp"
const utp = "utp"
const itp = "itp"

const udp = "udp"
const idp = "idp"
const adp = "adp"
const cdp = "cdp"
const gdp = "gdp"

const fmn = "fmn"
const fmnh2 = "fmnh2"
const nad = "nad"
const nadp = "nadp"
const nadh = "nadh"
const nadph = "nadph"
const fad = "fad"
const fadh2 = "fadh2"
const q8h2 = "q8h2"
const q8 = "q8"
const mql8 = "mql8"
const mqn8 = "mqn8"
const dmmql8 = "2dmmql8"
const dmmq8 = "2dmmq8"

const h2o = "h2o"
const phosphate = "pi"
const amet = "amet"
const pydx5p = "pydx5p"
const thmpp = "thmpp"
const accoa = "accoa"
const coa = "coa"
const akg = "akg"
const nh4 = "nh4"
const H = "h"
const ac = "ac"

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct BasicConfig
    minimum_metabolic_coverage::Float64
    minimum_growth_rate::Float64
    maximum_growth_rate::Float64
    optimizer_modifications::Vector{Function}
end

basic_config = BasicConfig(
    minimum_metabolic_coverage = 0.1,
    minimum_growth_rate = 0.01,
    maximum_growth_rate = 5.0,
    optimizer_modifications = [silence],
)

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct AnnotationConfig
    gene_annotation_keywords::Vector{String}
    gene_annotation_regexes::Dict{String,Regex}
    metabolite_annotation_keywords::Vector{String}
    metabolite_annotation_regexes::Dict{String,Regex}
    reaction_annotation_keywords::Vector{String}
    reaction_annotation_regexes::Dict{String,Regex}
    minimum_conformal_crossreferences::Int64
    minimum_crossreferences::Int64
end

annotation_config = AnnotationConfig(
    gene_annotation_keywords = [
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
    gene_annotation_regexes = Dict(
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
    metabolite_annotation_keywords = [
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
    metabolite_annotation_regexes = Dict(
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
    reaction_annotation_keywords = [
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
    reaction_annotation_regexes = Dict(
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
    minimum_conformal_crossreferences = 3,
    minimum_crossreferences = 4,
)

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct BiomassConfig
    essential_precursors::Vector{String}
end

biomass_config = BiomassConfig(
    essential_precursors = [
        trp__L,
        cys__L,
        his__L,
        tyr__L,
        met__L,
        phe__L,
        ser__L,
        pro__L,
        asp__L,
        thr__L,
        gln__L,
        glu__L,
        ile__L,
        arg__L,
        lys__L,
        val__L,
        leu__L,
        ala__L,
        gly,
        asn__L,
        datp,
        dctp,
        dttp,
        dgtp,
        atp,
        ctp,
        utp,
        gtp,
        nad,
        nadp,
        # amet,
        # fad,
        # pydx5p,
        coa,
        # thmpp,
        # fmn,
        h2o,
    ],
)

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct ConsistencyConfig
    consistency_ignored_reactions::Vector{String}
    tolerance_threshold::Float64
end

consistency_config =
    ConsistencyConfig(consistency_ignored_reactions = String[], tolerance_threshold = 1e-07)

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct EnergyConfig
    energy_dissipating_metabolites::Vector{String}
    additional_energy_generating_reactions::Vector{Reaction}
    ignored_energy_reactions::Vector{String}
    optimizer_modifications::Vector{Function}
end

energy_config = EnergyConfig(
    energy_dissipating_metabolites = [
        atp,
        ctp,
        gtp,
        utp,
        itp,
        adp,
        cdp,
        gdp,
        udp,
        idp,
        nadh,
        nad,
        nadph,
        nadp,
        fadh2,
        fad,
        fmnh2,
        fmn,
        q8h2,
        q8,
        mql8,
        mqn8,
        dmmql8,
        dmmq8,
        h2o,
        phosphate,
        accoa,
        coa,
        akg,
        H,
        ac,
    ],
    additional_energy_generating_reactions = Reaction[],
    ignored_energy_reactions = String[],
    optimizer_modifications = [silence],
)

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct MetaboliteConfig
    test_annotations::Vector{String}
end

metabolite_config =
    MetaboliteConfig(test_annotations = ["inchi_key", "chebi", "metanetx.chemical"])

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct NetworkConfig
    condition_number::Float64
    cycle_tol::Float64
    blocked_tol::Float64
    optimizer_modifications::Vector{Function}
end

network_config = NetworkConfig(
    condition_number = 1e9,
    cycle_tol = 1e-3,
    blocked_tol = 1e-3,
    optimizer_modifications = [silence],
)

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct ReactionConfig
    mass_ignored_reactions::Vector{String}
    charge_ignored_reactions::Vector{String}
end

reaction_config =
    ReactionConfig(mass_ignored_reactions = String[], charge_ignored_reactions = String[])

"""
$(TYPEDEF)

A grouping of parameters used by the metabolic testing infrastructure.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct MemoteConfig
    annotation::AnnotationConfig
    basic::BasicConfig
    biomass::BiomassConfig
    consistency::ConsistencyConfig
    energy::EnergyConfig
    metabolite::MetaboliteConfig
    network::NetworkConfig
    reaction::ReactionConfig
end

memote_config = MemoteConfig(
    annotation = annotation_config,
    basic = basic_config,
    biomass = biomass_config,
    consistency = consistency_config,
    energy = energy_config,
    metabolite = metabolite_config,
    network = network_config,
    reaction = reaction_config,
)

end # module
