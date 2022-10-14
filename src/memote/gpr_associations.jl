#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Return a list of reaction ids that do not have gene reaction rules (aka gene
protein reaction associations).
"""
n_reactions_without_gpr(model) = [rid for rid in reactions(model) if !_has_sensible_gpr(model, rid)]

n_reactions_with_complexes(model) = [rid for rid in reactions(model) if _has_sensible_gpr(model, rid) && any(length(grr) > 1 for grr in reaction_gene_association(model, rid))]

n_reactions_transport_no_gpr(model) = []