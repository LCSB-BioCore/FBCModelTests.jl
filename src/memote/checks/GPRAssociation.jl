"""
    module GPRAssociation

A module testing various facets of the gene reaction associations.
"""
module GPRAssociation

using COBREXA
using DocStringExtensions

import ..Config

"""
$(TYPEDSIGNATURES)

Check if a reaction has a gene reaction rule, and that each gene in the rule is
contained in the model.
"""
function reaction_has_sensible_gpr(model::MetabolicModel, rid)
    grrs = reaction_gene_association(model, rid)
    isnothing(grrs) && return false
    isempty(grrs) && return false
    any(isempty.(grrs)) && return false
    any("" in grr for grr in grrs) && return false

    gids = Set(reduce(vcat, grrs))

    return all(in.(gids, Ref(genes(model))))
end

"""
$(TYPEDSIGNATURES)

Return a list of reaction ids that have protein complexes assigned to them.
"""
reactions_with_complexes(model::MetabolicModel) = [
    rid for rid in reactions(model) if reaction_has_sensible_gpr(model, rid) &&
    any(length(grr) > 1 for grr in reaction_gene_association(model, rid))
]

end # module
