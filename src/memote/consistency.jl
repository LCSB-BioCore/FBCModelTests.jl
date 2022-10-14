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

Optionally, use `config.consistency.mass_ignored_reactions` to pass a vector
of reaction ids to ignore in this process. Internally biomass and exchang
reactions are ignored.
"""
function reactions_charge_unbalanced(model; config = memote_config)
    unbalanced_rxns = String[]
    ignored_reactions = [
        find_biomass_reaction_ids(model)
        find_exchange_reaction_ids(model)
    ]

    for rid in reactions(model)
        rid in [ignored_reactions; config.consistency.mass_ignored_reactions] && continue
        _rbal = 0
        for (mid, stoich) in reaction_stoichiometry(model, rid)
            try
                mc = metabolite_charge(model, mid)
                _rbal += isnothing(mc) ? Inf : mc * stoich
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
reaction balances. Returns a list of reaction IDs that are mass unbalanced,
which is empty if the test passes.

Optionally, use `config.consistency.charge_ignored_reactions` to pass a vector
of reaction ids to ignore in this process. Internally biomass and exchang
reactions are ignored.
"""
function reactions_mass_unbalanced(model; config = memote_config)
    unbalanced_rxns = String[]
    ignored_reactions = [
        find_biomass_reaction_ids(model)
        find_exchange_reaction_ids(model)
    ]

    for rid in reactions(model)
        rid in [ignored_reactions; config.consistency.charge_ignored_reactions] &&
            continue
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
`config.consistency.additional_energy_generating_reactions`, which should be
vector of COBREXA `Reaction`s using the same metabolite name space as the
`model`. Internally, the `model` is converted to a COBREXA `StandardModel`, so
ensure that the appropriate accessors are defined for it.

Since models use different name spaces,
`config.consistency.energy_dissipating_metabolites` is used to create the
energy dissipating reactions. By default it uses the BiGG name space, but this
will be changed to ChEBI in due course. If your model uses a different name
space, then you have to change the values (NOT the keys) of
`config.consistency.energy_dissipating_metabolites`. Each energy dissipating
reaction is added to the test only if all its associated metabolites are
present. Any `config.consistency.optimizer_modifications` to the solver are
passed directly through to COBREXA's `flux_balance_analysis` function. All
`config.consistency.boundary_reactions` and
`config.consistency.ignored_energy_reactions` are deleted from an internal
copy of `model`; this internal copy is used for analysis.

Returns `true` if the model has no energy generating cycles.
"""
function model_has_no_erroneous_energy_generating_cycles(
    model,
    optimizer;
    config = memote_config,
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
                [
                    config.consistency.energy_dissipating_metabolites[x] for
                    x in keys(mets)
                ],
                Ref(metabolites(_model)),
            ),
        )
            _mets = Dict(
                config.consistency.energy_dissipating_metabolites[k] => v for
                (k, v) in mets
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
    for rxn in config.consistency.additional_energy_generating_reactions
        push!(objective_ids, rxn.id)
        add_reaction!(_model, rxn)
    end

    # ignore reactions by just removing them
    for rid in [
        config.consistency.ignored_energy_reactions
        [rid for rid in reactions(_model) if is_boundary(_model, rid)]
    ]
        rid in reactions(_model) && remove_reaction!(_model, rid)
    end

    objval = solved_objective_value(
        flux_balance_analysis(
            _model,
            optimizer;
            modifications = [
                # config.consistency.optimizer_modifications
                change_objective(objective_ids)
            ],
        ),
    )
    # if problem does not solve then also fail
    isnothing(objval) && return false
    # if objective is approximately 0 then no energy generating cycles present
    isapprox(
        objval,
        0;
        atol = 1e-6,
    )
end

"""
$(TYPEDSIGNATURES)

Determines if the model is stoichiometrically consistent. Note, stoichiometric
consistency does not guarantee that mass balances must hold in the model. A more
robust check is [`is_model_mass_balanced`](@ref), but this works if not all
metabolites have mass assigned to them.

Based on Gevorgyan, Albert, Mark G. Poolman, and David A. Fell. "Detection of
stoichiometric inconsistencies in biomolecular models." Bioinformatics (2008).

Optionally ignore some reactions in this analysis by adding reaction IDs to
`config.consistency.consistency_ignored_reactions`.
"""
function model_is_consistent(model, optimizer; config = memote_config)
    #=
    Note, there is a MILP method that can be used to find the unconserved metabolite,
    but the problem is a MILP (and probably why the original MEMOTE takes so long to run).

    Note, it may be better to add additional constraints on the model to ensure that mass
    cannot be create (through lower and upper bounds on m). This is to prevent things like:
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
            [rid for rid in reactions(model) if is_boundary(model, rid)]
            find_biomass_reaction_ids(model)
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
    value.(m)
    termination_status(opt_model) == OPTIMAL
end

"""
$(TYPEDSIGNATURES)

Test if model is consistent by checking that:
1. the model is stoichiometrically consistent, tested with
   [`model_is_consistent`](@ref)
2. there are no energy generating cycles, tested with
   [`model_has_no_erroneous_energy_generating_cycles`](@ref)
3. the model is both mass and charge balanced, tested with
   [`reaction_charge_unbalanced`](@ref) and [`reaction_mass_unbalanced`]

Each function called in this test function can be called individually. The
kwargs are forwarded as indicated by the prefix.
"""
function test_consistency(model, optimizer; config = memote_config)
    @testset "Consistency" begin
        @test model_is_consistent(model, optimizer; config)
        # @test model_has_no_erroneous_energy_generating_cycles(model, optimizer; config)
        @test isempty(reactions_mass_unbalanced(model; config))
        @test isempty(reactions_charge_unbalanced(model; config))
    end
end
