"""
    module Metabolite

A module testing various metabolite properties.
"""
module Metabolite

using COBREXA
using DocStringExtensions

import ..Config

"""
$(TYPEDSIGNATURES)

Return a list of all boundary reactions that allow flux into the model and
create a metabolite. Assume that boundary reactions only create a single
metabolite. Use the testing config, `config.metabolite.only_imported =
false`, to also return metabolites that can be produced by the model under
default conditions.
"""
function metabolites_medium_components(model::MetabolicModel; config = Config.memote_config)
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
metabolites_no_formula(model::MetabolicModel; config = Config.memote_config) = [
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
metabolites_no_charge(model::MetabolicModel; config = Config.memote_config) = [
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
function metabolites_are_duplicated(
    model::MetabolicModel,
    m1,
    m2;
    config = Config.memote_config,
)
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
function metabolites_unique(model::MetabolicModel; config = Config.memote_config)
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
function metabolites_duplicated_in_compartment(
    model::MetabolicModel;
    config = Config.memote_config,
)
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


end # module
