#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Return a list of all boundary reactions that allow flux into the model and
create a metabolite. Assume that boundary reactions only create a single
metabolite. Use the testing config, `default.metabolite.only_imported =
false`, to also return metabolites that can be produced by the model under
default conditions.
"""
function metabolites_medium_components(model; config = memote_config)
    mets = String[]
    for (rid, lb, ub) in zip(reactions(model), bounds(model)...)
        if is_boundary(model, rid)
            mid = first(keys(reaction_stoichiometry(model, rid)))
            config.metabolite.medium_only_imported && lb < 0 && ub <= 0 && push!(mets, mid)
            !config.metabolite.medium_only_imported && lb < 0 && push!(mets, mid)
        end
    end
    mets
end

"""
$(TYPEDSIGNATURES)

List all metabolites without a formula. Use
`config.metabolite.formula_corner_cases` to specify an extra case to check for
formula's that are not properly assigned.
"""
metabolites_no_formula(model; config = memote_config) = [
    mid for mid in metabolites(model) if isnothing(metabolite_formula(model, mid)) ||
    isempty(metabolite_formula(model, mid)) ||
    any(
        in.(
            keys(metabolite_formula(model, mid)),
            Ref(config.metabolite.formula_corner_cases),
        ),
    )
]

"""
$(TYPEDSIGNATURES)

List all metabolites without a charge. Use
`config.metabolite.charge_corner_cases` to specify an extra case to check for
charge's that are not properly assigned.
"""
metabolites_no_charge(model; config = memote_config) = [
    mid for mid in metabolites(model) if isnothing(metabolite_charge(model, mid)) ||
    metabolite_charge(model, mid) in config.metabolite.charge_corner_cases
]

"""
$(TYPEDSIGNATURES)

Test if metabolites `m1` and `m2` are different by comparing their
`config.metabolite.test_annotation` field in the annotations of each
metabolite. Note, if no annotations are present for one or both of the
metabolites, then return `true`.
"""
function metabolites_are_duplicated(model, m1, m2; config = memote_config)
    k1s = get(metabolite_annotations(model, m1), config.metabolite.test_annotation, nothing)
    isnothing(k1s) && return true
    k2s = get(metabolite_annotations(model, m2), config.metabolite.test_annotation, nothing)
    isnothing(k2s) && return true
    any(in.(k1s, Ref(k2s)))
end

"""
$(TYPEDSIGNATURES)

Return a list of unique metabolites in model. Uses
[`metabolites_are_duplicated`](@ref) internally and forwards `test_annotation`
to it. The latter argument is used to determine if two metabolites are the same
by checking for any correspondence.
"""
function metabolites_unique(model; config = memote_config)
    unique_metabolites = Set{String}()
    for m1 in metabolites(model)
        duplicate = false
        for m2 in unique_metabolites
            duplicate = metabolites_are_duplicated(model, m1, m2; config)
            duplicate && break
        end
        !duplicate && push!(unique_metabolites, m1)
    end
    return unique_metabolites
end

"""
$(TYPEDSIGNATURES)

Return a dictionary of metabolites that are duplicated in their compartment.
"""
function metabolites_duplicated_in_compartment(model; config = memote_config)
    unique_metabolites = Dict{String,Set{String}}()
    for m1 in metabolites(model)
        c1 = metabolite_compartment(model, m1)
        for m2 in metabolites(model)
            c2 = metabolite_compartment(model, m2)
            if c1 == c2 && m1 != m2 && metabolites_are_duplicated(model, m1, m2; config)
                if haskey(unique_metabolites, c1)
                    push!(unique_metabolites[c1], m1)
                else
                    unique_metabolites[c1] = Set([m1])
                end
            end
        end
    end
    return unique_metabolites
end

"""
$(TYPEDSIGNATURES)

Test if the metabolites contained in the `model`:
1. are not duplicated, tested with [`metabolites_duplicated_in_compartment`](@ref)
2. all have a formula and charge associated with them, tested with [`metabolites_no_formula`](ref) and [`metabolites_no_charge`](@ref)
3. the default medium of the cell is not empty, tested with [`metabolites_medium_components`](ref).

Each of the basic functions can be run independently. THe kwargs are forwarded as indicated by the prefix.
"""
function test_metabolites(model; config = memote_config)
    @testset "Metabolite Information" begin
        @test isempty(metabolites_duplicated_in_compartment(model; config))
        @test isempty(metabolites_no_formula(model; config))
        @test isempty(metabolites_no_charge(model; config))
        @test !isempty(metabolites_medium_components(model; config))
    end
end


#TODO: change all instances of ::MetabolicModel to ::AbstractMetabolicModel once COBREXA 2.0 releases

"""
$(TYPEDSIGNATURES)

Checks for the presence of metabolite annotations.
Returns a vector of all metabolites without any annotations.
"""

function all_unannotated_metabolites(model::MetabolicModel)
    return [metabolite_id for metabolite_id in metabolites(model) if isempty(metabolite_annotations(model, metabolite_id))]
end

"""

$(TYPEDSIGNATURES)
Iterates through a model's metabolites to check if any common biochemical databases appear in the annotations field.
Returns a dictionary with the biochemical database names being checked for as keywords 
and the corresponding value being a vector of all metabolites missing that database name.
"""

function unannotated_metabolites(model::MetabolicModel,
    annotation_keywords = ["pubchem.compound", "kegg.compound", "seed.compound", "inchi_key", "inchi", 
    "chebi", "hmdb", "reactome", "metanetx.chemical", "bigg.metabolite", "biocyc"])
   missing_annos = Dict{String, Vector{String}}()
   for keyword in annotation_keywords
       missing_annos[keyword] = []
       for metabolite_id in metabolites(model)
           metabolite_annos = metabolite_annotations(model, metabolite_id)
           !haskey(metabolite_annos, keyword) && push!(missing_annos[keyword], metabolite_id)
       end
   end
   return missing_annos
end

"""
$(TYPEDSIGNATURES)

Metabolite Annotation Conformity Per Database
Uses the following regex dictionary to check if the annotations conform to patterns 
defined according to the MIRIAM guidelines,
i.e. matching those that are defined at https://identifiers.org/.
Returns a dictionary with a key for every database that has annotations 
with the corresponding value being an array of all metabolites whose annotations
do not match the given regex pattern.
"""

metabolites_regex = Dict("pubchem.compound" => r"^\d+$",
    "kegg.compound" => r"^C\d+$",
    "seed.compound" => r"^cpd\d+$",
    "inchi_key" => r"^[A-Z]{14}\-[A-Z]{10}(\-[A-Z])?",
    "inchi" => r"^InChI\=1S?\/[A-Za-z0-9\.]+(\+[0-9]+)?(\/[cnpqbtmsih][A-Za-z0-9\-\+\(\)\,\/\?\;\.]+)*$",
    "chebi" => r"^CHEBI:\d+$",
    "hmdb" => r"^HMDB\d{5}$",
    "reactome" => r"(^R-[A-Z]{3}-[0-9]+(-[0-9]+)?$)|(^REACT_\d+(\.\d+)?$)",
    "metanetx.chemical" => r"^MNXM\d+$",
    "bigg.metabolite" => r"^[a-z_A-Z0-9]+$",
    "biocyc" => r"^[A-Z-0-9]+(?<!CHEBI)(\:)?[A-Za-z0-9+_.%-]+$")
    

function metabolite_annotation_conformity(model::MetabolicModel, annotation_standards = metabolites_regex)
    nonconform_annos = Dict{String, Vector{String}}()
    for metabolite_id in metabolites(model)
        metabolite_annos = metabolite_annotations(model, metabolite_id)
        for key in keys(annotation_standards)
            if haskey(metabolite_annos, key) == true
                !in(key, keys(nonconform_annos)) && (nonconform_annos[key] = [])
                in.(nothing, Ref(match.(metabolites_regex[key],  metabolite_annos[key]))) && push!(nonconform_annos[key], metabolite_id)
            end
        end
    end
    return nonconform_annos
end
