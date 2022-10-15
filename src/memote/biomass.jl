#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Identify all the biomass reactions in the `model` using both sbo annotations, as
well biomass strings typically contained in their reaction IDs. Use
`config.biomass.biomass_strings` to update the list of strings to look for.
"""
model_biomass_reactions(model; config = memote_config) = Set(
    [
        [rid for rid in reactions(model) if is_biomass_reaction(model, rid)]
        find_biomass_reaction_ids(model; biomass_strings = config.biomass.biomass_strings)
    ],
)

"""
$(TYPEDSIGNATURES)

Check if model has an ATP maintenance reaction built in (also called a
non-growth associated maintenance cost). Looks for reaction annotations
corresponding to the sbo maintenance term, or looks for reaction ids that
contain the strings listed in `config.biomass.atpm_strings`.
"""
model_has_atpm_reaction(model; config = memote_config) =
    any(is_atp_maintenance_reaction(model, rid) for rid in reactions(model)) ||
    any(any(occursin.(config.biomass.atpm_strings, Ref(rid))) for rid in reactions(model))

"""
$(TYPEDSIGNATURES)

Check if the biomass reaction consumes ATP and H₂O, and produces ADP, HO₄P, and
H⁺. Each of these metabolites have a lookup table mapping them to the name space
of the model, defined in `config.biomass.growth_metabolites`. These need to be
set if you use anything other than the BiGG namespace.
"""
function model_has_growth_atp_in_biomass(model; config = memote_config)
    biomass_rxns = model_biomass_reactions(model; config)
    x = fill(true, length(biomass_rxns))
    for (idx, rid) in enumerate(biomass_rxns)
        d = reaction_stoichiometry(model, rid)
        for (m1, m2) in config.biomass.growth_metabolites
            if m1 in ["atp", "h2o"] # must consume
                !haskey(d, m2) && !(d[m2] < 0) && (x[idx] = false)
            else # must produce
                !haskey(d, m2) && !(d[m2] > 0) && (x[idx] = false)
            end
        end
    end
    return x
end

"""
$(TYPEDSIGNATURES)

For each biomass reaction, identified by [`model_biomass_reactions`](@ref),
calculate the molar weight of the reaction by summing the products of the
associated metabolite coefficients with their molar masses. This sum should fall
within [1 - 1e-3, 1 + 1e-6].
"""
function model_biomass_is_consistent(model; config = memote_config)
    biomass_rxns = model_biomass_reactions(model; config)
    to_symbol(x::String) =
        length(x) > 1 ? Symbol(uppercase(first(x)) * x[2:end]) : Symbol(uppercase(first(x)))
    get_molar_mass(mid) = begin
        rs = metabolite_formula(model, mid)
        sum(v * elements[to_symbol(k)].atomic_mass for (k, v) in rs).val
    end

    x = fill(0.0, length(biomass_rxns))
    for (idx, rid) in enumerate(biomass_rxns)
        d = reaction_stoichiometry(model, rid)
        x[idx] = sum(v * get_molar_mass(k) for (k, v) in d)
    end
    return -x ./ 1000.0
end

"""
$(TYPEDSIGNATURES)

Check if the model can be solved under default conditions and yield a reasonable
growth rate. Here reasonable is set via `config.biomass.minimum_growth_rate` and
`config.biomass.maximum_growth_rate`. Optionally, pass optimization
modifications to the solver through `config.biomass.optimizer_modifications`.
"""
function model_solves_default(model, optimizer; config = memote_config)
    mu = solved_objective_value(
        flux_balance_analysis(
            model,
            optimizer;
            modifications = config.biomass.optimizer_modifications,
        ),
    )
    config.biomass.minimum_growth_rate < mu < config.biomass.maximum_growth_rate
end


"""
$(TYPEDSIGNATURES)

Check if the model can synthesize all of the biomass precursors in all of the
biomass functions, except those listed in `config.biomass.ignored_precursors` in
the default medium. Set any optimizer modifications with
`config.biomass.optimizer_modifications`.
"""
function model_can_produce_biomass_precursors(model, optimizer; config = memote_config)
    stdmodel = convert(StandardModel, model) # convert to stdmodel so that reactions can be added/removed
    biomass_rxns = model_biomass_reactions(stdmodel; config)
    mids = Set(String[])
    for rid in biomass_rxns
        for mid in [k for (k, v) in reaction_stoichiometry(stdmodel, rid) if v < 0]
            push!(mids, mid)
        end
    end
    blocked_precursors = String[]
    for mid in filter(x -> x ∉ config.biomass.ignored_precursors, mids)
        rid = "MEMOTE_TEMP_RXN"
        add_reaction!(stdmodel, Reaction(rid, Dict(mid => -1), :forward))
        mu = solved_objective_value(
            flux_balance_analysis(
                stdmodel,
                optimizer;
                modifications = [
                    config.biomass.optimizer_modifications
                    change_objective(rid)
                ],
            ),
        )
        remove_reaction!(stdmodel, rid)
        config.biomass.minimum_growth_rate < mu && continue # no max bound required
        push!(blocked_precursors, mid)
    end
    return blocked_precursors
end

"""
$(TYPEDSIGNATURES)

Tests if each biomass reaction contains a set of essential precursors, listed in
`config.biomass.essential_precursors`. Note, this function only works on a
lumped biomass function.
"""
function biomass_contains_essential_precursors(model; config = memote_config)
    biomass_rxns = model_biomass_reactions(model; config)
    num_missing_essential_precursors = fill(String[], length(biomass_rxns))
    for (idx, rid) in enumerate(biomass_rxns)
        mids = [k for (k, v) in reaction_stoichiometry(model, rid) if v < 0]
        for mid in mids
            mid ∉ keys(config.biomass.essential_precursors) &&
                push!(num_missing_essential_precursors[idx], mid)
        end
    end
    return num_missing_essential_precursors
end
