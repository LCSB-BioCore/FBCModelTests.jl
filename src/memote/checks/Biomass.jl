"""
    module Biomass

This module contains tests that check the consistency of the biomass reaction.
"""
module Biomass

using DocStringExtensions
using COBREXA

import ..Config
import ..Utils

"""
$(TYPEDSIGNATURES)

Identify all the biomass reactions in the `model` using both sbo annotations, as
well biomass strings typically contained in their reaction IDs. Use
`config.biomass.biomass_strings` to update the list of strings to look for.
"""
model_biomass_reactions(model::MetabolicModel; config = Config.memote_config) = Set(
    [
        [rid for rid in reactions(model) if is_biomass_reaction(model, rid)]
        find_biomass_reaction_ids(model; biomass_strings = config.biomass.biomass_strings)
    ],
)

"""
$(TYPEDSIGNATURES)

Check if the biomass reaction consumes ATP and H₂O, and produces ADP, HO₄P, and
H⁺. Each of these metabolites have a lookup table mapping them to the name space
of the model, defined in `config.biomass.growth_metabolites`. These need to be
set if you use anything other than the BiGG namespace.
"""
function atp_present_in_biomass(model::MetabolicModel; config = Config.memote_config)
    biomass_rxns = model_biomass_reactions(model; config)
    x = Dict(biomass_rxns .=> true)
    for rid in biomass_rxns
        d = reaction_stoichiometry(model, rid)
        for (m1, m2) in config.biomass.growth_metabolites
            if m1 in ["atp", "h2o"] # must consume
                x[rid] &= get(d, m2, Inf) <= 0 # thanks mirek!
            else # must produce
                x[rid] &= get(d, m2, -Inf) >= 0
            end
        end
    end
    return x
end

"""
$(TYPEDSIGNATURES)

For each biomass reaction, identified by [`model_biomass_reactions`](@ref),
calculate the molar weight of the reaction by summing the products of the
associated metabolite coefficients with their molar masses.
"""
function model_biomass_molar_mass(model::MetabolicModel; config = Config.memote_config)
    biomass_rxns = model_biomass_reactions(model; config)
    get_molar_mass(mid) = begin
        rs = metabolite_formula(model, mid)
        sum(v * Utils.to_element(k).atomic_mass for (k, v) in rs).val
    end

    x = Dict(biomass_rxns .=> 0.0)
    for rid in biomass_rxns
        d = reaction_stoichiometry(model, rid)
        x[rid] += sum(v * get_molar_mass(k) for (k, v) in d)
    end
    for (k, v) in x
        x[k] /= -1000.0
    end
    return x
end

"""
$(TYPEDSIGNATURES)

Check that the molar mass of each biomass reactions falls within `[1 - 1e-3, 1 +
1e-6]` by calling [`model_biomass_molar_mass`](@ref) internally.
"""
model_biomass_is_consistent(model::MetabolicModel; config = Config.memote_config) =
    all(-1e-3 .< (collect(values(model_biomass_molar_mass(model; config))) .- 1) .< 1e-6)

"""
$(TYPEDSIGNATURES)

Check if the model can synthesize all of the biomass precursors in all of the
biomass functions, except those listed in `config.biomass.ignored_precursors` in
the default medium. Set any optimizer modifications with
`config.biomass.optimizer_modifications`.
"""
function find_blocked_biomass_precursors(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
)
    stdmodel = convert(StandardModel, model) # convert to stdmodel so that reactions can be added/removed
    biomass_rxns = model_biomass_reactions(stdmodel; config)
    blocked_precursors = Dict{String,Vector{String}}()

    for rid in biomass_rxns
        mids = String[]
        for mid in [k for (k, v) in reaction_stoichiometry(stdmodel, rid) if v < 0]
            push!(mids, mid)
        end

        for mid in mids
            mid in config.biomass.ignored_precursors && continue
            temp_rid = "MEMOTE_TEMP_RXN"
            add_reaction!(stdmodel, Reaction(temp_rid, Dict(mid => -1), :forward))
            mu = solved_objective_value(
                flux_balance_analysis(
                    stdmodel,
                    optimizer;
                    modifications = [
                        config.biomass.optimizer_modifications
                        change_objective(temp_rid)
                    ],
                ),
            )
            remove_reaction!(stdmodel, temp_rid)
            config.biomass.minimum_growth_rate < mu && continue # no max bound required
            push!(get!(blocked_precursors, rid, String[]), mid)
        end
    end

    return blocked_precursors
end

"""
$(TYPEDSIGNATURES)

Tests if each biomass reaction contains a set of essential precursors, listed in
`config.biomass.essential_precursors`. Note, this function only works on a
lumped biomass function.
"""
function biomass_missing_essential_precursors(
    model::MetabolicModel;
    config = Config.memote_config,
)
    biomass_rxns = model_biomass_reactions(model; config)
    num_missing_essential_precursors = Dict{String,Vector{String}}() # for some reason can't do biomass_rxns .=> String[]
    for rid in biomass_rxns
        b = reaction_stoichiometry(model, rid)
        for mid in values(config.biomass.essential_precursors)
            haskey(b, mid) && b[mid] < 0 && continue # is a precursor
            push!(get!(num_missing_essential_precursors, rid, String[]), mid)
        end
    end
    return num_missing_essential_precursors
end

end # module
