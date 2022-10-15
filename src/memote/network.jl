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
stoichiometric_max_min_value(model) = /(reverse(extrema(abs, x for x in stoichiometry(model) if x != 0.0))...)

"""
$(TYPEDSIGNATURES)

Test if the stoichiometric matrix is well conditioned by determining if
[`stoichiometric_max_min_value`](@ref) is less than 10⁹.
"""
stoichiometric_well_conditioned(model) = stoichiometric_max_min_value(model) < 10^9

"""
$(TYPEDSIGNATURES)

Make all boundary reactions reversible and run FVA on the model to find all
reactions that are universally blocked.
"""
function find_all_universally_blocked_reactions(model, optimizer; config = memote_config)
    stdmodel = convert(StandardModel, model)
    for rid in reactions(stdmodel)
        if is_boundary(stdmodel, rid)
            change_bound!(stdmodel, rid, lower=-1000, upper=1000)
        end
    end
    fva = flux_variability_analysis_dict(stdmodel, optimizer; modifications = config.network.)
end