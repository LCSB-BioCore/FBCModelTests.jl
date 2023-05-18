"""
    module Network

A module testing the network and topology properties of the model.
"""
module Network

using COBREXA
using DocStringExtensions
using JuMP
using SparseArrays
using Distributed

import ..Config

"""
$(TYPEDSIGNATURES)

Return the ratio of the absolute maximum and minimum value of the nonzero
coefficients in the stoichiometric matrix of `model`.
"""
stoichiometric_max_min_ratio(model::MetabolicModel) =
    /(reverse(extrema(abs, x for x in stoichiometry(model) if x != 0.0))...)

"""
$(TYPEDSIGNATURES)

Find all reactions that participate in stoichiometrically balanced cycles by
closing all boundary reactions and running fva on the resultant model.
"""
function find_cycle_reactions(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
    workers = [myid()],
)
    if model isa StandardModel
        stdmodel = deepcopy(model) # copy because will add stuff to it
    else
        stdmodel = convert(StandardModel, deepcopy(model))
    end

    for rid in reactions(stdmodel)
        if is_boundary(stdmodel, rid)
            change_bound!(stdmodel, rid, lower = 0, upper = 0)
        else # remove fixed constraints
            if stdmodel.reactions[rid].lb < 0
                stdmodel.reactions[rid].lb = -1000.0
            else
                stdmodel.reactions[rid].lb = 0
            end
            if stdmodel.reactions[rid].ub > 0
                stdmodel.reactions[rid].ub = 1000.0
            else
                stdmodel.reactions[rid].ub = 0
            end
        end
    end

    idx_id = Dict(
        idx => id for (idx, id) in enumerate(reactions(model)) if !is_boundary(model, id)
    )
    fvas = flux_variability_analysis(
        stdmodel,
        collect(keys(idx_id)),
        optimizer;
        modifications = config.network.optimizer_modifications,
        workers,
        optimal_objective_value = -Inf,
    )

    cycle_reactions = Int64[]
    for (i, vs) in zip(collect(keys(idx_id)), eachrow(fvas))
        for v in vs
            isnothing(v) && continue
            abs(v) > config.network.cycle_tol && (push!(cycle_reactions, i); break)
        end
    end

    [idx_id[idx] for idx in cycle_reactions]
end

"""
$(TYPEDSIGNATURES)

Make all boundary reactions reversible and run FVA on the model to find all
reactions that are universally blocked. Optimizer modifications can be passed
through `config.network.optimizer_modifications`
"""
function find_all_universally_blocked_reactions(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
    workers = [myid()],
)

    if model isa StandardModel
        stdmodel = deepcopy(model) # copy because will add stuff to it
    else
        stdmodel = convert(StandardModel, deepcopy(model))
    end

    for rid in reactions(stdmodel)
        if is_boundary(stdmodel, rid)
            change_bound!(stdmodel, rid, lower = -1000, upper = 1000)
        end
    end

    fvas = flux_variability_analysis(
        stdmodel,
        collect(1:n_reactions(stdmodel)),
        optimizer;
        modifications = config.network.optimizer_modifications,
        workers,
        optimal_objective_value = -Inf,
    )

    blocked_rxns = String[]
    for (vs, rid) in zip(eachrow(fvas), reactions(model))
        isnothing(vs[1]) && continue
        isnothing(vs[2]) && continue
        abs(vs[1]) <= config.network.blocked_tol &&
            abs(vs[2]) <= config.network.blocked_tol && push!(blocked_rxns, rid)
    end

    blocked_rxns
end

end # module
