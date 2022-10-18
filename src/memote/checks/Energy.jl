"""
    module Energy

This module checks if the model is energetically sensible.
"""
module Energy

using COBREXA
using DocStringExtensions
using JuMP

import ..Config

"""
$(TYPEDSIGNATURES)

Check if model has an ATP maintenance reaction built in (also called a
non-growth associated maintenance cost). Looks for reaction annotations
corresponding to the sbo maintenance term, or looks for reaction ids that
contain the strings listed in `config.energy.atpm_strings`.
"""
model_has_atpm_reaction(model::MetabolicModel; config = Config.memote_config) =
    any(is_atp_maintenance_reaction(model, rid) for rid in reactions(model)) ||
    any(any(occursin.(config.energy.atpm_strings, Ref(rid))) for rid in reactions(model))

"""
$(TYPEDSIGNATURES)

Attempts to detect if the `model` contains any erroneous energy generating
cycles by closing all exchange reactions and using flux balance analysis to
maximize the sum of fluxes through a set of energy dissipating reactions. The
flux sum should be zero if the model is free of energy generating reactions.
This function is based on Fritzemeier, Claus Jonathan, et al. "Erroneous
energy-generating cycles in published genome scale metabolic networks:
Identification and removal." PLoS computational biology (2017).
The energy dissipating reactions are based on the source paper, and include:
```
ATP + H2O --> ADP + H + Phosphate
CTP + H2O --> CDP + H + Phosphate
GTP + H2O --> GDP + H + Phosphate
UTP + H2O --> UDP + H + Phosphate
ITP + H2O --> IDP + H + Phosphate
NADH --> H + NAD
NADPH --> H + NADP
FADH2 --> 2 H + FAD
FMNH2 --> 2 H + FMN
Ubiquinol-8 --> 2 H + Ubiquinone-8
Menaquinol-8 --> 2 H + Menaquinone-8
2-Demethylmenaquinol-8 --> 2 H + 2-Demethylmenaquinone-8
H2O + ACCOA --> H + Acetate + COA
L-Glutamate + H2O --> 2-Oxoglutarate + Ammonium + 2 H
H[external] --> H
```
Additional energy dissipating reactions can be directly specified through
`config.energy.additional_energy_generating_reactions`, which should be
vector of COBREXA `Reaction`s using the same metabolite name space as the
`model`. Internally, the `model` is converted to a COBREXA `StandardModel`, so
ensure that the appropriate accessors are defined for it.
Since models use different name spaces,
`config.energy.energy_dissipating_metabolites` is used to create the
energy dissipating reactions. By default it uses the BiGG name space, but this
will be changed to ChEBI in due course. If your model uses a different name
space, then you have to change the values (NOT the keys) of
`config.energy.energy_dissipating_metabolites`. Each energy dissipating
reaction is added to the test only if all its associated metabolites are
present. Any `config.energy.optimizer_modifications` to the solver are
passed directly through to COBREXA's `flux_balance_analysis` function. All
`config.energy.boundary_reactions` and
`config.energy.ignored_energy_reactions` are deleted from an internal
copy of `model`; this internal copy is used for analysis.
Returns `true` if the model has no energy generating cycles.
"""
function model_has_no_erroneous_energy_generating_cycles(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
)
    # need to add reactions, only implemented for Core and StdModel
    if model isa StandardModel
        _model = deepcopy(model) # copy because will add stuff to it
    else
        _model = convert(StandardModel, deepcopy(model))
    end

    objective_ids = String[]
    """
    Helper function, adds the energy dissipating reaction only if all the
    metabolites are in the model.
    """
    function maybe_add_energy_reaction(mets, id)
        if all(
            in.(
                [config.energy.energy_dissipating_metabolites[x] for x in keys(mets)],
                Ref(metabolites(_model)),
            ),
        )
            _mets = Dict(
                config.energy.energy_dissipating_metabolites[k] => v for (k, v) in mets
            )
            rid = "MEMOTE_TEMP_RXN_$id"
            add_reaction!(_model, Reaction(rid, _mets, :forward))
            push!(objective_ids, rid)
        end
    end

    # add basic energy dissipating reactions
    maybe_add_energy_reaction(
        Dict("ATP" => -1, "H" => 1, "H2O" => -1, "ADP" => 1, "Phosphate" => 1),
        "ATP",
    )
    maybe_add_energy_reaction(
        Dict("CTP" => -1, "H" => 1, "H2O" => -1, "CDP" => 1, "Phosphate" => 1),
        "CTP",
    )
    maybe_add_energy_reaction(
        Dict("GTP" => -1, "H" => 1, "H2O" => -1, "GDP" => 1, "Phosphate" => 1),
        "GTP",
    )
    maybe_add_energy_reaction(
        Dict("UTP" => -1, "H" => 1, "H2O" => -1, "UDP" => 1, "Phosphate" => 1),
        "UTP",
    )
    maybe_add_energy_reaction(
        Dict("ITP" => -1, "H" => 1, "H2O" => -1, "IDP" => 1, "Phosphate" => 1),
        "ITP",
    )

    maybe_add_energy_reaction(Dict("NADH" => -1, "H" => 1, "NAD" => 1), "NADH")
    maybe_add_energy_reaction(Dict("NADPH" => -1, "H" => 1, "NADP" => 1), "NADPH")

    maybe_add_energy_reaction(Dict("FADH2" => -1, "H" => 2, "FAD" => 1), "FADH2")
    maybe_add_energy_reaction(Dict("FMNH2" => -1, "H" => 2, "FMN" => 1), "FMNH2")
    maybe_add_energy_reaction(
        Dict("Ubiquinol-8" => -1, "H" => 2, "Ubiquinone-8" => 1),
        "Q8H2",
    )
    maybe_add_energy_reaction(
        Dict("Menaquinol-8" => -1, "H" => 2, "Menaquinone-8" => 1),
        "MQL8",
    )
    maybe_add_energy_reaction(
        Dict("2-Demethylmenaquinol-8" => -1, "H" => 2, "2-Demethylmenaquinone-8" => 1),
        "DMMQL8",
    )

    maybe_add_energy_reaction(
        Dict("ACCOA" => -1, "H2O" => -1, "H" => 1, "Acetate" => 1, "COA" => 1),
        "ACCOA",
    )

    maybe_add_energy_reaction(
        Dict(
            "L-Glutamate" => -1,
            "H2O" => -1,
            "2-Oxoglutarate" => 1,
            "Ammonium" => 1,
            "H" => 2,
        ),
        "GLU",
    )

    maybe_add_energy_reaction(Dict("H[external]" => -1, "H" => 1), "PROTON")

    # add user specified reactions
    for rxn in config.energy.additional_energy_generating_reactions
        push!(objective_ids, rxn.id)
        add_reaction!(_model, rxn)
    end

    # ignore reactions by just removing them
    for rid in [
        config.energy.ignored_energy_reactions
        [rid for rid in reactions(_model) if is_boundary(_model, rid)]
    ]
        rid in reactions(_model) && remove_reaction!(_model, rid)
    end

    objval = solved_objective_value(
        flux_balance_analysis(
            _model,
            optimizer;
            modifications = [
                config.energy.optimizer_modifications
                change_objective(objective_ids)
            ],
        ),
    )
    # if problem does not solve then also fail
    isnothing(objval) && return false
    # if objective is approximately 0 then no energy generating cycles present
    isapprox(objval, 0; atol = 1e-6)
end

end # module
