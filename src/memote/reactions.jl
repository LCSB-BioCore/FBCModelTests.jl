#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Test if a reaction is constrained.
"""
_is_constrained(lb, ub, default) = lb ∉ [-default, 0, default] || ub ∉ [-default, 0, default] 

"""
$(TYPEDSIGNATURES)

Find all purely metabolic reactions and indicate which are constrained beyond
directionality (if value of dictionary is true, then constrained). For a
constraint to be not purely directional, the lower or upper bound needs to be
different from `[-config.reactions.bound_default, 0,
config.reactions.bound_default]`. Metabolic reactions exclude transport,
boundary and biomass reactions.
"""
function find_all_purely_metabolic_reactions(model; config = memote_config)
    biomass_rxns = model_biomass_reactions(model; config)
    metabolic_reactions = Dict{String, Bool}()
    for (lb, ub, rid) in zip(bounds(model)..., reactions(model))
        is_boundary(model, rid) && continue
        rid in biomass_rxns && continue
        _probably_transport_reaction(model, rid, config.reaction.test_annotation) && continue
        metabolic_reactions[rid] = _is_constrained(lb, ub, config.reaction.bound_default)
    end
    return metabolic_reactions
end

"""
$(TYPEDSIGNATURES)

Find all transport reactions and indicate which are constrained beyond
directionality. For a constraint to be not purely directional, the lower or
upper bound needs to be different from `[-config.reactions.bound_default, 0,
config.reactions.bound_default]`. Transport reactions are heuristically
identified, see [`!_probably_transport_reaction`](@ref), which partially uses
reaction annotations. Set the annotation field to use via
`config.reactions.rest_annotation`. 
"""
function find_all_transport_reactions(model; config = memote_config)
    transport_reactions = Dict{String, Bool}()
    for (lb, ub, rid) in zip(bounds(model)..., reactions(model))
        !_probably_transport_reaction(model, rid, config.reaction.test_annotation) && continue
        transport_reactions[rid] = _is_constrained(lb, ub, config.reaction.bound_default)
    end
    return transport_reactions
end

"""
$(TYPEDSIGNATURES)

Find all reactions with overlapping annotation information. Internally calls
[`annotation_index`](@ref). Some annotations, like sbo terms will necessarily be
non-unique, ignore annotations like this by editing `config.reaction.ignore_annotations`.
"""
function reactions_with_partially_identical_annotations(model; config = memote_config)
    stdmodel = convert(StandardModel, model)
    idx = annotation_index(stdmodel.reactions)
    for anno in config.reaction.ignore_annotations
        delete!(idx, anno)
    end
    ambiguously_identified_items(idx)
end

"""
$(TYPEDSIGNATURES)

Return a list of all reactions that are duplicated.
"""
function duplicate_reactions(model)
    stdmodel = convert(StandardModel, model)
    duplicated_rxns = Set{String}()
    for rid in reactions(stdmodel) 
        dup = check_duplicate_reaction(
            stdmodel.reactions[rid],
            stdmodel.reactions;
            only_metabolites = false,
        )
        isnothing(dup) && continue
        push!(duplicated_rxns, dup)
    end
    return duplicated_rxns
end

"""
$(TYPEDSIGNATURES)

Identify reactions that have the same gene reaction rules. Does not take
directionality or compartments into account.
"""
function reactions_with_identical_genes(model)
    grr_rids = Dict{Vector{String}, Set{String}}()
    rids = reactions(model)
    for rid1 in rids
        !_has_sensible_gpr(model, rid1) && continue
        grrs1 = reaction_gene_association(model, rid1)
        for rid2 in rids
            rid1 == rid2 && continue
            !_has_sensible_gpr(model, rid2) && continue
            grrs2 = reaction_gene_association(model, rid2)
            
            for grr1 in grrs1
                for grr2 in grrs2
                    issetequal(grr1, grr2) && push!(get!(grr_rids, sort(grr1), Set{Vector{String}}()), rid1)
                end
            end
        end 
    end
    return grr_rids
end