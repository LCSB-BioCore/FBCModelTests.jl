module GPRAssociation

using ..DocStringExtensions
using ..ModuleTools
using ..COBREXA
import ..Config.memote_config
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

reactions_with_complexes(model) = [
    rid for rid in reactions(model) if _has_sensible_gpr(model, rid) &&
    any(length(grr) > 1 for grr in reaction_gene_association(model, rid))
]

reactions_transport_no_gpr(model; config = memote_config) = [
    rid for rid in reactions(model) if !_has_sensible_gpr(model, rid) &&
    _probably_transport_reaction(model, rid, config.metabolite.test_annotation)
]

@export_locals

end # module