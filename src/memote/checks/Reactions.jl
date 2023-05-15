"""
    module Reactions

A module testing reaction properties.
"""
module Reactions

using COBREXA
using DocStringExtensions

import ..Config
import ..Biomass: findall_biomass_reactions
import ..GPRAssociation: reaction_has_sensible_gpr

"""
$(TYPEDSIGNATURES)

Return a list of all reactions that are duplicated.
"""
function findall_duplicated_reactions(model::MetabolicModel)
    stdmodel = convert(StandardModel, model) # not modified
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

Check if the charges across a reaction balance.
"""
function reaction_is_charge_balanced(model::MetabolicModel, rid::String)
    rs = reaction_stoichiometry(model, rid)
    isempty(rs) && return false

    _rbal = 0
    for (mid, stoich) in rs
        mc = metabolite_charge(model, mid)
        _rbal += isnothing(mc) ? Inf : mc * stoich
    end
    isapprox(_rbal, 0)
end

"""
$(TYPEDSIGNATURES)

Check if the mass across a reaction balances.
"""
function reaction_is_mass_balanced(model::MetabolicModel, rid::String)
    try
        _rbal = reaction_atom_balance(model, rid) # this throws an error if no formula for metabolite
        return all(values(_rbal) .== 0)
    catch
        return false
    end
end

"""
$(TYPEDSIGNATURES)

Check if model has an ATP maintenance reaction built in (also called a
non-growth associated maintenance cost). Looks for reaction annotations
corresponding to the sbo maintenance term.
"""
model_has_atpm_reaction(model::MetabolicModel) =
    any(is_atp_maintenance_reaction(model, rid) for rid in reactions(model))

end # module
