"""
    module Metabolites

A module testing various metabolite properties.
"""
module Metabolites

using COBREXA
using DocStringExtensions
using SparseArrays

import ..Config
import ..Utils: parse_annotations

"""
$(TYPEDSIGNATURES)

Test if metabolites `m1` and `m2` are different by comparing their
`config.metabolite.test_annotation` field in the annotations of each
metabolite. Note, if no annotations are present for one or both of the
metabolites, then return `false`.
"""
function metabolites_are_duplicated(model::MetabolicModel, m1, m2, test_annotation)
    k1s =
        get(parse_annotations(metabolite_annotations(model, m1)), test_annotation, nothing)
    isnothing(k1s) && return false
    k2s =
        get(parse_annotations(metabolite_annotations(model, m2)), test_annotation, nothing)
    isnothing(k2s) && return false
    any(in.(k1s, Ref(k2s)))
end

"""
$(TYPEDSIGNATURES)

Return a dictionary of metabolites that are duplicated in their compartment. If
any of the test annotations, stored in `config.metabolite.test_annotations`, are
repeated, then the metabolite is counted as duplicated. Missing annotations are ignored.
"""
function metabolites_duplicated_in_compartment(
    model::MetabolicModel;
    config = Config.memote_config,
)
    unique_metabolites = Set{String}()
    for m1 in metabolites(model)
        c1 = metabolite_compartment(model, m1)
        for m2 in metabolites(model)
            c2 = metabolite_compartment(model, m2)
            if c1 == c2 &&
               m1 != m2 &&
               any(
                   metabolites_are_duplicated(model, m1, m2, test_annotation) for
                   test_annotation in config.metabolite.test_annotations
               )
                push!(unique_metabolites, m1)
            end
        end
    end
    return unique_metabolites
end

"""
$(TYPEDSIGNATURES)

Helper function to find orphan or deadend metabolites. Specify `consumed=true`
to consider orphan metabolites or `false` to consider deadend metabolites. Set
`complete_medium=true` to open all boundary reactions to simulate a complete
medium.
"""
function _find_orphan_or_deadend_metabolites(model::MetabolicModel; consumed = true)
    mids = metabolites(model)
    mets = String[]
    S = stoichiometry(model)
    lbs, ubs = bounds(model)
    for idx in axes(S, 1)
        rids, vals = findnz(S[idx, :])
        if length(vals) == 1 && (consumed ? first(vals) < 0 : first(vals) > 0)
            ridx = first(rids)
            rid = reactions(model)[ridx]
            met = mids[idx]
            v = reaction_stoichiometry(model, rid)[met]
            # check if reaction can actually make or consume the metabolite
            lbs[ridx] < 0 < ubs[ridx] && continue # ignore reversible reactions
            consumed && lbs[ridx] * v <= 0 && push!(mets, met)
            !consumed && ubs[ridx] * v >= 0 && push!(mets, met)
        end
    end
    return mets
end

"""
$(TYPEDSIGNATURES)

Find all metabolites that can only (excludes reversible reactions) be consumed
in the `model` by inspecting the stoichiometric matrix.
"""
find_orphan_metabolites(model::MetabolicModel) =
    _find_orphan_or_deadend_metabolites(model, consumed = true)

"""
$(TYPEDSIGNATURES)

Find all metabolites that can only (excludes reversible reactions) be produced
in the `model` by inspecting the stoichiometric matrix.
"""
find_deadend_metabolites(model::MetabolicModel) =
    _find_orphan_or_deadend_metabolites(model, consumed = false)


"""
$(TYPEDSIGNATURES)

Returns a list of all metabolites that aren't part of any reactions.
"""
find_disconnected_metabolites(model::MetabolicModel) =
    metabolites(model)[all(stoichiometry(model) .== 0, dims = 2)[:]]


end # module
