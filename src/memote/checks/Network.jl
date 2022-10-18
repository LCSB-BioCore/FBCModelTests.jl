"""
module Network

A module testing the network and topology properties of the model.
"""
module Network

using DocStringExtensions
using COBREXA
using SparseArrays
import ..Config

"""
$(TYPEDSIGNATURES)

Return the ratio of the absolute maximum and minimum value of the nonzero
coefficients in the stoichiometric matrix of `model`.
"""
stoichiometric_max_min_ratio(model) =
    /(reverse(extrema(abs, x for x in stoichiometry(model) if x != 0.0))...)

"""
$(TYPEDSIGNATURES)

Test if the stoichiometric matrix is well conditioned by determining if
[`stoichiometric_max_min_ratio`](@ref) is less than 10⁹ (which can be set in `config.network.condition_number`).
"""
stoichiometric_matrix_is_well_conditioned(model; config = Config.memote_config) =
    stoichiometric_max_min_ratio(model) < config.network.condition_number

"""
$(TYPEDSIGNATURES)

Make all boundary reactions reversible and run FVA on the model to find all
reactions that are universally blocked. Optimizer modifications can be passed
through `config.network.optimizer_modifications`
"""
function find_all_universally_blocked_reactions(model, optimizer; config = Config.memote_config)
    stdmodel = convert(StandardModel, model)
    for rid in reactions(stdmodel)
        if is_boundary(stdmodel, rid)
            change_bound!(stdmodel, rid, lower = -1000, upper = 1000)
        end
    end
    mins, maxs = flux_variability_analysis_dict(
        stdmodel,
        optimizer;
        bounds = objective_bounds(config.network.fva_bound),
        modifications = config.network.optimizer_modifications,
    )
    return [
        rid for rid = reactions(stdmodel) if
        isapprox(abs(mins[rid][rid]), 0.0) && isapprox(abs(maxs[rid][rid]), 0.0)
    ]
end

"""
$(TYPEDSIGNATURES)

Helper function to find orphan or deadend metabolites. Specify `consumed=true`
to consider orphan metabolites or `false` to consider deadend metabolites. Set
`complete_medium=true` to open all boundary reactions to simulate a complete
medium.
"""
function _find_orphan_or_deadend_metabolites(model; consumed = true)
    mids = metabolites(model)
    mets = String[]
    S = stoichiometry(model)
    lbs, ubs = bounds(model)
    for idx in axes(S, 1)
        rids, vals = findnz(S[idx, :])
        if length(vals) == 1 && (consumed ? first(vals) < 0 : first(vals) > 0)
            ridx = first(rids)
            rid = reactions(model)[ridx]
            met = mids[idx]
            v = reaction_stoichiometry(model, rid)[met]
            # check if reaction can actually make or consume the metabolite
            lbs[ridx] < 0 < ubs[ridx] && continue # ignore reversible reactions
            consumed && lbs[ridx] * v <= 0 && push!(mets, met)
            !consumed && ubs[ridx] * v >= 0 && push!(mets, met)
        end
    end
    return mets
end

"""
$(TYPEDSIGNATURES)

Find all metabolites that can only (excludes reversible reactions) be consumed
in the `model` by inspecting the stoichiometric matrix.
"""
find_orphan_metabolites(model) = _find_orphan_or_deadend_metabolites(model, consumed = true)

"""
$(TYPEDSIGNATURES)

Find all metabolites that can only (excludes reversible reactions) be produced
in the `model` by inspecting the stoichiometric matrix.
"""
find_deadend_metabolites(model) =
    _find_orphan_or_deadend_metabolites(model, consumed = false)

"""
$(TYPEDSIGNATURES)

Find all reactions that participate in stoichiometrically balanced cycles by
closing all boundary reactions and running fva on the resultant model.
"""
function find_cycle_reactions(model, optimizer; config = Config.memote_config)
    stdmodel = convert(StandardModel, model)
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
    cycle_reactions = Set(String[])
    for rid in filter(x -> !is_boundary(model, x), reactions(model))
        for sense in [COBREXA.MIN_SENSE, COBREXA.MAX_SENSE]
            mu = solved_objective_value(
                flux_balance_analysis(
                    stdmodel,
                    optimizer;
                    modifications = [
                        config.network.optimizer_modifications
                        change_objective(rid)
                        change_sense(sense)
                    ],
                ),
            )
            isnothing(mu) && continue
            abs(mu) > config.network.cycle_tol && (push!(cycle_reactions, rid); break)
        end
    end
    return cycle_reactions
end

"""
$(TYPEDSIGNATURES)

Return a dictionary mapping orphan and deadend metabolites that occur in the
model in complete medium. Complete medium is modeled by opening all the bounday
reactions. For each metabolite FBA is run with a temporary reaction either
consuming or producing the metabolite in question. At minimum a flux of
`config.network.minimum_metabolite_flux` must be attained for the metabolite to
pass the test.
"""
function find_complete_medium_orphans_and_deadends(model, optimizer; config = Config.memote_config)
    stdmodel = convert(StandardModel, model)
    for rid in reactions(stdmodel)
        change_bound!(stdmodel, rid, lower = -1000.0, upper = 1000.0)
    end

    found_mets = Dict{Symbol,Vector{String}}()
    for mid in metabolites(model)
        for (k, v) in [(:produce, -1), (:consume, 1)]
            temp_rid = "MEMOTE_TEMP_RXN"
            add_reaction!(stdmodel, Reaction(temp_rid, Dict(mid => v), :forward))
            mu = solved_objective_value(
                flux_balance_analysis(
                    stdmodel,
                    optimizer;
                    modifications = [
                        config.network.optimizer_modifications
                        change_objective(temp_rid)
                    ],
                ),
            )
            remove_reaction!(stdmodel, temp_rid)
            if isnothing(mu) || mu < config.network.minimum_metabolite_flux
                push!(get!(found_mets, k, String[]), mid)
            end
        end
    end
    return found_mets
end

end # module