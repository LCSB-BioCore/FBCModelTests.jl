"""
$(TYPEDSIGNATURES)

Run a MEMOTE-like test suite on `model` using `optimizer` to solve optimization
problems. Forwards arguments to [`run_tests`](@ref), but wraps all of the test
sets in a `FBCModelTests.QuietTestSet` to generate more concise user
facing results.
"""
function run_tests_toplevel(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
    filename = nothing,
    workers = [myid()],
)
    @testset QuietTestSet "Memote.jl" begin
        run_tests(model, optimizer; config, workers, filename)
    end
end

"""
$(TYPEDSIGNATURES)

Overload of [`run_tests_toplevel`](@ref) that works directly with a file.
"""
run_tests_toplevel(filename::String, optimizer; kwargs...) = run_tests_toplevel(
    load_model(StandardModel, filename),
    optimizer;
    filename = basename(filename),
    kwargs...,
)

"""
$(TYPEDSIGNATURES)

Overload of [`run_tests`](@ref) that works directly with a file.
"""
run_tests(filename::String, optimizer; kwargs...) = run_tests(
    load_model(StandardModel, filename),
    optimizer;
    filename = basename(filename),
    kwargs...,
)

"""
$(TYPEDSIGNATURES)

Run a MEMOTE-like test suite on `model` using `optimizer` to solve optimization
problems; some basic parameters and thresholds are taken from `config`.

Some of the tests internally convert the input model to `StandardModel` of
COBREXA; supplying a `StandardModel` may thus increase efficiency of the whole
process.
"""
function run_tests(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
    filename = nothing,
    workers = [myid()],
)
    @testset "Metabolic Model Tests$(isnothing(filename) ? "" : ": $filename")" begin
        @info "Testing basic information about the model..."
        @testset "Basic information" begin
            @atest Basic.model_has_reactions(model) "The model has reactions"
            @atest Basic.model_has_metabolites(model) "The model has metabolites"
            @atest Basic.model_has_genes(model) "The model has genes"
            @atest Basic.model_has_compartments(model) "The model has at least 1 comparment"
            metcov = Basic.model_metabolic_coverage(model)
            @atest config.basic.minimum_metabolic_coverage <= metcov "The metabolic coverage of the model is $(round(metcov, digits=2)), but it should be bigger than $(round(config.basic.minimum_metabolic_coverage, digits=2)) to pass"
            @atest Basic.model_solves_in_default_medium(model, optimizer; config) "The model solves in the default medium with an objective value within [$(round(config.basic.minimum_growth_rate, digits=2)), $(round(config.basic.maximum_growth_rate,digits=2))]"
        end

        @info "Testing the model annotations..."
        @testset "Annotations" begin
            @testset "SBO annotations" begin
                @atest any(is_metabolic_reaction(model, rid) for rid in reactions(model)) "At least 1 reaction is SBO annotated as a metabolic reaction"
                @atest any(is_transport_reaction(model, rid) for rid in reactions(model)) "At least 1 reaction is SBO annotated as a transport reaction"
                @atest any(is_biomass_reaction(model, rid) for rid in reactions(model)) "At least 1 reaction is SBO annotated as a biomass reaction"
                @atest any(is_exchange_reaction(model, rid) for rid in reactions(model)) "At least 1 reaction is SBO annotated as an exchange reaction"
                @atest any(
                    is_atp_maintenance_reaction(model, rid) for rid in reactions(model)
                ) "At least 1 reaction is SBO annotated as an ATP maintenance reaction"
                for rid in reactions(model)
                    @atest any(
                        in.(
                            keys(
                                Memote.Utils.parse_annotations(
                                    reaction_annotations(model, rid),
                                ),
                            ),
                            Ref(["sbo", "SBO"]),
                        ),
                    ) "$rid has an SBO annotation"
                end
            end

            @testset "Reactions" begin
                @testset "Any annotations present" begin
                    for rid in reactions(model)
                        @atest Annotation.is_annotated_reaction(model, rid) "$rid has at least 1 annotation"
                    end
                end

                @testset "At least $(config.annotation.minimum_crossreferences) cross-references per entry" begin
                    for rid in reactions(model)
                        @atest config.annotation.minimum_crossreferences <= length(
                            Annotation.findall_annotated_reaction_databases(
                                model,
                                rid;
                                config,
                            ),
                        ) "$rid needs at least $(config.annotation.minimum_crossreferences) cross-references"
                    end
                end
                @testset "At least $(config.annotation.minimum_conformal_crossreferences) recognizable cross-references per entry" begin
                    for rid in reactions(model)
                        @atest config.annotation.minimum_conformal_crossreferences <=
                               length(
                            Annotation.findall_conformal_reaction_annotations(
                                model,
                                rid;
                                config,
                            ),
                        ) "$rid needs at least $(config.annotation.minimum_conformal_crossreferences) recognizable cross-references"
                    end
                end
            end

            @testset "Metabolites" begin
                @testset "Any annotations present" begin
                    for mid in metabolites(model)
                        @atest Annotation.is_annotated_metabolite(model, mid) "$mid has at least 1 annotation"
                    end
                end

                @testset "At least $(config.annotation.minimum_crossreferences) cross-references per entry" begin
                    for mid in metabolites(model)
                        @atest config.annotation.minimum_crossreferences <= length(
                            Annotation.findall_annotated_metabolite_databases(
                                model,
                                mid;
                                config,
                            ),
                        ) "$mid dneeds at least $(config.annotation.minimum_crossreferences) cross-references"
                    end
                end
                @testset "At least $(config.annotation.minimum_conformal_crossreferences) recognizable cross-references per entry" begin
                    for mid in metabolites(model)
                        @atest config.annotation.minimum_conformal_crossreferences <=
                               length(
                            Annotation.findall_conformal_metabolite_annotations(
                                model,
                                mid;
                                config,
                            ),
                        ) "$mid needs at least $(config.annotation.minimum_conformal_crossreferences) recognizable cross-references"
                    end
                end
            end

            @testset "Genes" begin
                @testset "Any annotations present" begin
                    for gid in genes(model)
                        @atest Annotation.is_annotated_gene(model, gid) "$gid has at least 1 annotation"
                    end
                end

                @testset "At least $(config.annotation.minimum_crossreferences) cross-references per entry" begin
                    for gid in genes(model)
                        @atest config.annotation.minimum_crossreferences <= length(
                            Annotation.findall_annotated_gene_databases(model, gid; config),
                        ) "$gid needs at least $(config.annotation.minimum_crossreferences) cross-references"
                    end
                end
                @testset "At least $(config.annotation.minimum_conformal_crossreferences) recognizable cross-references per entry" begin
                    for gid in genes(model)
                        @atest config.annotation.minimum_conformal_crossreferences <=
                               length(
                            Annotation.findall_conformal_gene_annotations(
                                model,
                                gid;
                                config,
                            ),
                        ) "$gid needs at least $(config.annotation.minimum_conformal_crossreferences) recognizable cross-references"
                    end
                end
            end
        end

        @info "Testing the biomass function(s)..."
        @testset "Biomass" begin
            @atest length(Biomass.findall_biomass_reactions(model)) > 0 "At least 1 biomass reaction is in the model"

            @testset "ATP + H2O -> ADP + Pi + H included" begin
                for rid in Biomass.findall_biomass_reactions(model)
                    @atest Biomass.atp_present_in_biomass_reaction(model, rid) "$rid (biomass reaction) hydrolyses ATP into ADP"
                end
            end
            @testset "Molar mass consistent" begin
                for rid in Biomass.findall_biomass_reactions(model)
                    @atest Biomass.biomass_reaction_is_consistent(model, rid) "$rid (biomass reaction) is consistent"
                end
            end
            @testset "Has common precursors" begin
                for rid in Biomass.findall_biomass_reactions(model)
                    @atest Biomass.biomass_missing_essential_precursors(model, rid; config) "$rid (biomass reaction) contains all the common precursors"
                end
            end
        end

        @info "Testing the stoichiometric consistency of the model..."
        @testset "Stoichiometric consistency" begin
            @atest Consistency.model_is_consistent(model, optimizer; config) "The model is stoichiometrically consistent"
        end

        @info "Testing the energy metabolism..."
        @testset "Energy metabolism" begin # TODO make this play nicely with compartments
            @testset "Erroneous energy cycles" begin
                @test_skip Energy.model_has_no_erroneous_energy_generating_cycles(
                    model,
                    optimizer;
                    config,
                )
            end
        end

        @info "Testing the gene-protein-reaction associations..."
        @testset "GPR associations" begin
            @testset "Metabolic reactions have GPRs" begin
                for rid in filter(x -> is_metabolic_reaction(model, x), reactions(model))
                    @atest GPRAssociation.reaction_has_sensible_gpr(model, rid) "$rid (metabolic reaction) has at least 1 GPRs"
                end
            end

            @testset "Transport reactions have GPRs" begin
                for rid in filter(x -> is_transport_reaction(model, x), reactions(model))
                    @atest GPRAssociation.reaction_has_sensible_gpr(model, rid) "$rid (transport reaction) has at least 1 GPRs"
                end
            end

            @atest length(GPRAssociation.reactions_with_complexes(model)) != 0 "At least 1 enzyme complex exists in the model"
        end

        @info "Testing the reaction information..."
        @testset "Reaction information" begin
            tms_rxns = [
                rid for rid in reactions(model) if is_transport_reaction(model, rid) ||
                is_metabolic_reaction(model, rid) ||
                is_spontaneous_reaction(model, rid)
            ]
            @testset "Charge balanced (transport, metabolic, spontaneous)" begin
                for rid in tms_rxns
                    if rid in config.reaction.charge_ignored_reactions
                        @test_skip Reactions.reaction_is_charge_balanced(model, rid)
                    else
                        @atest Reactions.reaction_is_charge_balanced(model, rid) "$rid is charge balanced"
                    end
                end
            end
            @testset "Mass balanced (transport, metabolic, spontaneous)" begin
                for rid in tms_rxns
                    if rid in config.reaction.mass_ignored_reactions
                        @test_skip Reactions.reaction_is_mass_balanced(model, rid)
                    else
                        @atest Reactions.reaction_is_mass_balanced(model, rid) "$rid is mass balanced"
                    end
                end
            end
            @testset "No duplicated reactions" begin
                dup_rxns = Reactions.findall_duplicated_reactions(model)
                for rid in reactions(model)
                    @atest rid ∉ dup_rxns "$rid is not duplicated"
                end
            end

            @atest Reactions.model_has_atpm_reaction(model) "There is an ATP maintenance reaction in the model"
        end

        @info "Testing the metabolite information..."
        @testset "Metabolite information" begin
            @testset "No missing formulas" begin
                for mid in metabolites(model)
                    @atest !isnothing(metabolite_formula(model, mid)) "$mid has a formula"
                end
            end
            @testset "No missing charges" begin
                for mid in metabolites(model)
                    @atest !isnothing(metabolite_charge(model, mid)) "$mid has a charge"
                end
            end
            @testset "No duplicated metabolites" begin
                dup_mids = Metabolites.metabolites_duplicated_in_compartment(model; config)
                for mid in metabolites(model)
                    @atest mid ∉ dup_mids "$mid is not duplicated"
                end
            end
            @testset "No orphan metabolites" begin
                orphans = Metabolites.find_orphan_metabolites(model)
                for mid in metabolites(model)
                    @atest mid ∉ orphans "$mid is not an orphan metabolite"
                end
            end
            @testset "No deadend metabolites" begin
                deadends = Metabolites.find_deadend_metabolites(model)
                for mid in metabolites(model)
                    @atest mid ∉ deadends "$mid is not a deadend metabolite"
                end
            end
            @testset "No disconnected metabolites" begin
                disconnected = Metabolites.find_disconnected_metabolites(model)
                for mid in metabolites(model)
                    @atest mid ∉ disconnected "$mid is not disconnected"
                end
            end
        end

        @info "Testing the network information..."
        @testset "Network information" begin
            @atest Network.stoichiometric_max_min_ratio(model) <
                   config.network.condition_number "The stoichiometric matrix is well conditioned"

            @testset "No stoichiometrically balanced cycles (reactions)" begin
                rxns_in_cycles =
                    Network.find_cycle_reactions(model, optimizer; config, workers)
                for rid in reactions(model)
                    @atest rid ∉ rxns_in_cycles "$rid does not take part in a stoichiometrically balanced cycle"
                end
            end
            @testset "No universally blocked reactions" begin
                blocked_rxns = Network.find_all_universally_blocked_reactions(
                    model,
                    optimizer;
                    config,
                    workers,
                )
                for rid in reactions(model)
                    @atest rid ∉ blocked_rxns "$rid is not a universally blocked reaction"
                end
            end
        end
    end
end
