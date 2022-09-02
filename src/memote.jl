#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Iterates through all the reactions in `model` and checks if the charges across
each reaction balance. Returns a list of reaction IDs that are charge unbalanced, 
which is empty if the test passes. 
"""
function is_model_charge_balanced(model::COBREXA.MetabolicModel)
    unbalanced_rxns = String[]

    for rid in reactions(model)
        _rbal = 0
        for (mid, stoich) in reacton_stoichiometry(model, rid)
            try
                _rbal += metabolite_charge(model, mid) * stoich
            catch 
                throw(error("Something is wrong with reaction $rid."))
            end
        end
        !approx(_rbal, 0) && push!(unbalanced_rxns, rid) 
    end

    return unbalanced_rxns
end

"""
$(TYPEDSIGNATURES)

Iterates through all the reactions in `model` and checks if the mass across each
reaction balances. Returns a list of reaction IDs that are mass unbalanced, which
is empty if the test passes.
"""
function is_model_mass_balanced(model::COBREXA.MetabolicMode)
    unbalanced_rxns = String[]
    
    for rid in reactions(model)
        try 
            _rbal = reaction_atom_balance(model, rid)
            all(values(_rbal) .== 0) || push!(unbalanced_rxns, rid)
        catch
            throw(error("Something is wrong with reaction $rid."))
        end
    end

    return unbalanced_rxns
end