"""
    module GPRAssociation

A module testing various facets of the gene reaction associations.
"""
module GPRAssociation

using DocStringExtensions
using COBREXA
import ..Config
import ..Utils: _has_sensible_gpr, _probably_transport_reaction

"""
$(TYPEDSIGNATURES)

Return a list of reaction ids that do not have gene reaction rules (aka gene
protein reaction associations).
"""
reactions_without_gpr(model) = [
    rid for rid in reactions(model) if
    !_has_sensible_gpr(model, rid) && !is_boundary(model, rid)
]

"""
$(TYPEDSIGNATURES)

Return a list of reaction ids that have protein complexes assigned to them.
"""
reactions_with_complexes(model) = [
    rid for rid in reactions(model) if _has_sensible_gpr(model, rid) &&
    any(length(grr) > 1 for grr in reaction_gene_association(model, rid))
]

"""
$(TYPEDSIGNATURES)

Return a list of transport reactions that do not have gene reaction rules
assigned to them.
"""
reactions_transport_no_gpr(model; config = Config.memote_config) = [
    rid for rid in reactions(model) if !_has_sensible_gpr(model, rid) &&
    _probably_transport_reaction(model, rid, config.metabolite.test_annotation)
]


end # module
