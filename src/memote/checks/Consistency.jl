"""
    module Consistency

This module checks if the metabolic model is overall consistent.
"""
module Consistency

using DocStringExtensions
using COBREXA
using JuMP
import ..Config

"""
$(TYPEDSIGNATURES)

Iterates through all the reactions in `model` and checks if the charges across
each reaction balance. Returns a list of reaction IDs that are charge
unbalanced.
Optionally, use `config.consistency.mass_ignored_reactions` to pass a vector
of reaction ids to ignore in this process. Internally biomass and exchang
reactions are ignored.
"""
function reactions_charge_unbalanced(model; config = Config.memote_config)
    unbalanced_rxns = String[]
    ignored_reactions = [
        find_biomass_reaction_ids(model)
        find_exchange_reaction_ids(model)
        config.consistency.mass_ignored_reactions
    ]

    for rid in reactions(model)
        rid in ignored_reactions && continue
        _rbal = 0
        for (mid, stoich) in reaction_stoichiometry(model, rid)
            try
                mc = metabolite_charge(model, mid)
                _rbal += isnothing(mc) ? Inf : mc * stoich
            catch
                throw(error("Something is wrong with reaction $rid."))
            end
        end
        !isapprox(_rbal, 0) && push!(unbalanced_rxns, rid)
    end

    return unbalanced_rxns
end

"""
$(TYPEDSIGNATURES)

Iterates through all the reactions in `model` and checks if the mass across each
reaction balances. Returns a list of reaction IDs that are mass unbalanced.
Optionally, use `config.consistency.charge_ignored_reactions` to pass a vector
of reaction ids to ignore in this process. Internally biomass and exchang
reactions are ignored.
"""
function reactions_mass_unbalanced(model; config = Config.memote_config)
    unbalanced_rxns = String[]
    ignored_reactions = [
        find_biomass_reaction_ids(model)
        find_exchange_reaction_ids(model)
        config.consistency.charge_ignored_reactions
    ]

    for rid in reactions(model)
        rid in ignored_reactions && continue
        try
            _rbal = reaction_atom_balance(model, rid)
            all(values(_rbal) .== 0) || push!(unbalanced_rxns, rid)
        catch
            throw(error("Something is wrong with reaction $rid."))
        end
    end

    return unbalanced_rxns
end

"""
$(TYPEDSIGNATURES)

Determines if the model is stoichiometrically consistent. Note, stoichiometric
consistency does not guarantee that mass balances must hold in the model. A more
robust check is [`reactions_mass_unbalanced`](@ref), but this works if not all
metabolites have mass assigned to them.
Based on Gevorgyan, Albert, Mark G. Poolman, and David A. Fell. "Detection of
stoichiometric inconsistencies in biomolecular models." Bioinformatics (2008).
Optionally ignore some reactions in this analysis by adding reaction IDs to
`config.consistency.consistency_ignored_reactions`.
"""
function model_is_consistent(model, optimizer; config = Config.memote_config)
    #=
    Note, there is a MILP method that can be used to find the unconserved metabolite,
    but the problem is a MILP (and probably why the original MEMOTE takes so long to run).
    Note, it may be better to add additional constraints on the model to ensure that mass
    cannot be create (through lower and upper bounds on m). This is to prevent things like:
    A -> x*B -> C where x can be anything. This test will not catch these kinds of errors.
    =#

    # need to remove reactions, only implemented for Core and StdModel
    if model isa StandardModel
        _model = deepcopy(model) # copy because will add stuff to it
    else
        _model = convert(StandardModel, deepcopy(model))
    end

    remove_reactions!(
        _model,
        [
            [rid for rid in reactions(model) if is_boundary(model, rid)]
            find_biomass_reaction_ids(model)
            config.consistency.consistency_ignored_reactions
        ],
    )

    N = stoichiometry(_model)
    n_mets, _ = size(N)

    opt_model = Model(optimizer)
    m = @variable(opt_model, 1 <= m[1:n_mets])
    @constraint(opt_model, N' * m .== 0)
    @objective(opt_model, Min, sum(m))
    optimize!(opt_model)
    value.(m)
    termination_status(opt_model) == OPTIMAL
end

end # module
