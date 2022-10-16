#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Return the ratio of the absolute maximum and minimum value of the nonzero
coefficients in the stoichiometric matrix of `model`.
"""
stoichiometric_max_min_ratio(model) = /(reverse(extrema(abs, x for x in stoichiometry(model) if x != 0.0))...)

"""
$(TYPEDSIGNATURES)

Test if the stoichiometric matrix is well conditioned by determining if
[`stoichiometric_max_min_value`](@ref) is less than 10⁹ (which can be set in `config.network.condition_number`).
"""
stoichiometric_matrix_is_well_conditioned(model; config = memote_config) = stoichiometric_max_min_ratio(model) < config.network.condition_number

"""
$(TYPEDSIGNATURES)

Make all boundary reactions reversible and run FVA on the model to find all
reactions that are universally blocked. Optimizer modifications can be passed
through `config.network.optimizer_modifications`
"""
function find_all_universally_blocked_reactions(model, optimizer; config = memote_config)
    stdmodel = convert(StandardModel, model)
    for rid in reactions(stdmodel)
        if is_boundary(stdmodel, rid)
            change_bound!(stdmodel, rid, lower=-1000, upper=1000)
        end
    end
    mins, maxs = flux_variability_analysis_dict(
        stdmodel, 
        optimizer;
        bounds = objective_bounds(config.network.fva_bound),
        modifications = config.network.optimizer_modifications
    )
    blocked_rids = String[]
    for rid in reactions(stdmodel)
        isapprox(abs(mins[rid][rid]), 0.0) && isapprox(abs(maxs[rid][rid]), 0.0) && push!(blocked_rids, rid)
    end
    return blocked_rids
end

"""
$(TYPEDSIGNATURES)

Helper function to find orphan or deadend metabolites. Specify `consumed=true`
to consider orphan metabolites or `false` to consider deadend metabolites. Set
`complete_medium=true` to open all boundary reactions to simulate a complete
medium.
"""
function _find_orphan_or_deadend_metabolites(model; consumed=true)
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
find_orphan_metabolites(model) = _find_orphan_or_deadend_metabolites(model, consumed=true)

"""
$(TYPEDSIGNATURES)

Find all metabolites that can only (excludes reversible reactions) be produced
in the `model` by inspecting the stoichiometric matrix.
"""
find_deadend_metabolites(model) = _find_orphan_or_deadend_metabolites(model, consumed=false)


"""
$(TYPEDSIGNATURES)

Find all reactions that participate in stoichiometrically balanced cycles by
closing all boundary reactions and running fva on the resultant model.
"""
function find_cycle_reactions(model, optimizer; config = memote_config)
    stdmodel = convert(StandardModel, model)
    for rid in reactions(stdmodel)
        if is_boundary(stdmodel, rid)
            change_bound!(stdmodel, rid, lower=0, upper=0)
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
        for sense in [MIN_SENSE, MAX_SENSE]
            mu = solved_objective_value(
                flux_balance_analysis(
                    stdmodel,
                    optimizer;
                    modifications = [
                        config.network.optimizer_modifications
                        change_objective(rid)
                        change_sense(sense)
                    ]
                )
            )
            isnothing(mu) && continue
            abs(mu) > config.network.cycle_tol && (push!(cycle_reactions, rid); break)
        end
    end
    return cycle_reactions
end