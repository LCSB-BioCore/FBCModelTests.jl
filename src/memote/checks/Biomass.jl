"""
    module Biomass

This module contains tests that check the consistency of the biomass reaction.
"""
module Biomass

using DocStringExtensions
using COBREXA

import ..Config
import ..Utils: get_molar_mass, parse_annotations

"""
$(TYPEDSIGNATURES)

Identify all the biomass reactions in the `model` using only sbo annotations.
"""
findall_biomass_reactions(model::MetabolicModel) =
    [rid for rid in reactions(model) if is_biomass_reaction(model, rid)]

"""
$(TYPEDSIGNATURES)

Check if the biomass reaction consumes ATP and H₂O, and produces ADP, HO₄P, and
H⁺. Annotations are parsed, and the BiGG namespace is used to identify these
metabolites in the underlying model.
"""
function atp_present_in_biomass_reaction(model::MetabolicModel, rid::String;)

    atp_found = false
    adp_found = false
    phos_found = false
    water_found = true
    proton_found = false

    for (mid, v) in reaction_stoichiometry(model, rid)
        annos = get(
            parse_annotations(metabolite_annotations(model, mid)),
            "bigg.metabolite",
            String[],
        ) # return safe fallback

        # consume
        if any(in.(Config.atp, annos)) && v < 0
            atp_found = true
        end
        if any(in.(Config.h2o, annos)) && v < 0
            water_found = true
        end
        # produce
        if any(in.(Config.adp, annos)) && v > 0
            adp_found = true
        end
        if any(in.(Config.phosphate, annos)) && v > 0
            phos_found = true
        end
        if any(in.(Config.H, annos)) && v > 0
            proton_found = true
        end
    end

    return atp_found && adp_found && phos_found && water_found && proton_found
end

"""
$(TYPEDSIGNATURES)

For a biomass reaction `rid`, calculate the molar weight of the reaction by
summing the products of the associated metabolite coefficients with their molar
masses (g/mol).
"""
function biomass_reaction_molar_mass(model::MetabolicModel, rid::String;)
    d = reaction_stoichiometry(model, rid)
    isempty(d) && return NaN
    x = sum(v * get_molar_mass(model, k) for (k, v) in d) / -1000.0
    return x
end

"""
$(TYPEDSIGNATURES)

Check that the molar mass of a biomass reactions falls within `[1 - 1e-3, 1 +
1e-6]` by calling [`biomass_reaction_molar_mass`](@ref) internally.
"""
biomass_reaction_is_consistent(model::MetabolicModel, rid::String;) =
    -1e-3 < (biomass_reaction_molar_mass(model, rid) - 1) < 1e-6

"""
$(TYPEDSIGNATURES)

Tests if the biomass reaction contains a set of essential precursors, listed in
`config.biomass.essential_precursors`. Uses the BiGG namespace as a reference
and the metabolite annotations.
"""
function biomass_missing_essential_precursors(
    model::MetabolicModel,
    rid::String;
    config = Config.memote_config,
)

    biomass = Dict{String,Float64}()
    for (mid, coeff) in reaction_stoichiometry(model, rid)
        ids = get(
            parse_annotations(metabolite_annotations(model, mid)),
            "bigg.metabolite",
            ["missing"],
        )
        for id in ids
            biomass[id] = coeff
        end
    end

    for mid in config.biomass.essential_precursors
        (haskey(biomass, mid) && biomass[mid] < 0) || return true
    end
    return false
end

end # module
