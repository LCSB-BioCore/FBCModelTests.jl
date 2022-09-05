#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Iterates through all the reactions in `model` and checks if the charges across
each reaction balance. Returns a list of reaction IDs that are charge
unbalanced, which is empty if the test passes.

Optionally, pass a list of reactions to ignore in this process through
`ignored_reactions`. It makes sense to include the biomass and
exchange reactions in this list (default).
"""
function is_model_charge_balanced(
    model::COBREXA.MetabolicModel;
    ignored_reactions = [
        find_biomass_reaction_ids(model)
        find_exchange_reaction_ids(model)
    ],
)
    unbalanced_rxns = String[]

    for rid in reactions(model)
        rid in ignored_reactions && continue
        _rbal = 0
        for (mid, stoich) in reaction_stoichiometry(model, rid)
            try
                _rbal += metabolite_charge(model, mid) * stoich
            catch
                throw(error("Something is wrong with reaction $rid."))
            end
        end
        !isapprox(_rbal, 0) && push!(unbalanced_rxns, rid)
    end

    return unbalanced_rxns
end

"""
$(TYPEDSIGNATURES)

Iterates through all the reactions in `model` and checks if the mass across each
reaction balances. Returns a list of reaction IDs that are mass unbalanced, which
is empty if the test passes.

Optionally, pass a list of reactions to ignore in this process through
`ignored_reactions`. It makes sense to include the biomass and
exchange reactions in this list (default).
"""
function is_model_mass_balanced(
    model::COBREXA.MetabolicModel;
    ignored_reactions = [
        find_biomass_reaction_ids(model)
        find_exchange_reaction_ids(model)
    ],
)
    unbalanced_rxns = String[]

    for rid in reactions(model)
        rid in ignored_reactions && continue
        try
            _rbal = reaction_atom_balance(model, rid)
            all(values(_rbal) .== 0) || push!(unbalanced_rxns, rid)
        catch
            throw(error("Something is wrong with reaction $rid."))
        end
    end

    return unbalanced_rxns
end

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
`additional_energy_generating_reactions`, which should be vector of COBREXA
`Reaction`s using the same metabolite name space as the `model`. Internally, the
`model` is converted to a COBREXA `StandardModel`, so ensure that the
appropriate accessors are defined for it.

Since models use different name spaces, `energy_dissipating_metabolites` is used
to create the energy dissipating reactions. By default it uses the BiGG name
space, but this will be changed to ChEBI in due course. Any `modifications` to
the solver are passed directly through to COBREXA's `flux_balance_analysis`
function. All `exchange_reactions` and `ignore_reactions` are deleted from an
internal copy of `model`; this internal copy is used for analysis.

Returns `true` if the model has energy generating cycles.
"""
function has_erroneous_energy_generating_cycles(
    model,
    optimizer;
    energy_dissipating_metabolites = Dict(
        "ATP" => "atp_c",
        "CTP" => "ctp_c",
        "GTP" => "gtp_c",
        "UTP" => "utp_c",
        "ITP" => "itp_c",
        "ADP" => "adp_c",
        "CDP" => "cdp_c",
        "GDP" => "gdp_c",
        "UDP" => "udp_c",
        "IDP" => "idp_c",
        "NADH" => "nadh_c",
        "NAD" => "nad_c",
        "NADPH" => "nadph_c",
        "NADP" => "nadp_c",
        "FADH2" => "fadh2_c",
        "FAD" => "fad_c",
        "FMNH2" => "fmn_c",
        "FMN" => "fmnh2_c",
        "Ubiquinol-8" => "q8h2_c",
        "Ubiquinone-8" => "q8_c",
        "Menaquinol-8" => "mql8_c",
        "Menaquinone-8" => "mqn8_c",
        "2-Demethylmenaquinol-8" => "2dmmql8_c",
        "2-Demethylmenaquinone-8" => "2dmmq8_c",
        "ACCOA" => "accoa_c",
        "COA" => "coa_c",
        "L-Glutamate" => "glu__L_c",
        "2-Oxoglutarate" => "akg_c",
        "Ammonium" => "nh4_c",
        "H" => "h_c",
        "H[external]" => "h_e",
        "H2O" => "h2o_c",
        "Phosphate" => "pi_c",
        "Acetate" => "ac_c",
    ),
    additional_energy_generating_reactions = [],
    ignore_reactions = ["ATPM"],
    modifications = [],
    exchange_reactions = find_exchange_reaction_ids(model),
)
    # need to add reactions, only implemented for Core and StdModel
    if model isa StandardModel
        _model = deepcopy(model) # copy because will add stuff to it
    else
        _model = convert(StandardModel, model)
    end

    objective_ids = String[]
    """
    Helper function, adds the energy dissipating reaction only if all the
    metabolites are in the model.
    """
    function maybe_add_energy_reaction(mets, id)
        if all(
            in.(
                [energy_dissipating_metabolites[x] for x in keys(mets)],
                Ref(metabolites(_model)),
            ),
        )
            _mets = Dict(energy_dissipating_metabolites[k] => v for (k, v) in mets)
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
    for rxn in additional_energy_generating_reactions
        push!(objective_ids, rxn.id)
        add_reaction!(_model, rxn)
    end

    # ignore reactions by just removing them
    for rxn in [ignore_reactions; exchange_reactions]
        remove_reaction!(_model, rxn)
    end

    # if objective is approximately 0 then no energy generating cycles present
    !isapprox(
        solved_objective_value(
            flux_balance_analysis(
                _model,
                optimizer;
                modifications = [modifications; change_objective(objective_ids)],
            ),
        ),
        0;
        atol = 1e-6,
    )
end
