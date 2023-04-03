"""
    module Consistency

This module checks if the metabolic model is overall consistent.
"""
module Consistency

using COBREXA
using DocStringExtensions
using JuMP

import ..Config
import ..Utils

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
function model_is_consistent(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
)
    #=
    Note, there is a MILP method that can be used to find the unconserved metabolite,
    but the problem is a MILP (and probably why the original MEMOTE takes so long to run).
    Note, it may be better to add additional constraints on the model to ensure that mass
    cannot be created (through lower and upper bounds on m). This is to prevent things like:
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
            [rid for rid in reactions(_model) if is_boundary(_model, rid)]
            filter(x -> is_biomass_reaction(_model, x), reactions(_model))
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
    termination_status(opt_model) == OPTIMAL
end

end # module
