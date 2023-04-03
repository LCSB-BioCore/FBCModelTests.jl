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
    workers = workers(),
)
    @testset "Metabolic Model Tests$(isnothing(filename) ? "" : ": $filename")" begin

        @testset "Basic information" begin
            @testset "Has reactions" begin
                @test Basic.model_has_reactions(model)
            end
            @testset "Has metabolites" begin
                @test Basic.model_has_metabolites(model)
            end
            @testset "Has genes" begin
                @test Basic.model_has_genes(model)
            end
            @testset "Has compartments" begin
                @test Basic.model_has_compartments(model)
            end
            @testset "Metabolic coverage exceeds threshold" begin
                @test config.basic.minimum_metabolic_coverage <=
                      Basic.model_metabolic_coverage(model)
            end
            @testset "Solves in default medium" begin
                @test Basic.model_solves_in_default_medium(model, optimizer; config)
            end
        end

        @testset "Annotations" begin
            @testset "SBO annotations" begin
                @testset "Has SBO annotated metabolic reactions" begin
                    @test any(is_metabolic_reaction(model, rid) for rid in reactions(model))
                end
                @testset "Has SBO annotationed transport reactions" begin
                    @test any(is_transport_reaction(model, rid) for rid in reactions(model))
                end
            end
            @testset "Reactions" begin
                @testset "Any annotations" begin
                    for rid in reactions(model)
                        @test Annotation.is_annotated_reaction(model, rid)
                    end
                end

                @testset "At most $(config.annotation.maximum_missing_databases) missing cross-references" begin
                    for rid in reactions(model)
                        @test length(
                            Annotation.findall_unannotated_reaction_databases(
                                model,
                                rid;
                                config,
                            ),
                        ) <= config.annotation.maximum_missing_databases
                    end
                end
                @testset "At most $(config.annotation.maximum_nonconformal_references) recognizable cross-references" begin
                    for rid in reactions(model)
                        @test length(
                            Annotation.findall_nonconformal_reaction_annotations(
                                model,
                                rid;
                                config,
                            ),
                        ) <= config.annotation.maximum_nonconformal_references
                    end
                end
            end

            @testset "Metabolites" begin
                @testset "Any annotations" begin
                    for mid in metabolites(model)
                        @test Annotation.is_annotated_metabolite(model, mid)
                    end
                end

                @testset "At most $(config.annotation.maximum_missing_databases) missing cross-references" begin
                    for mid in metabolites(model)
                        @test length(
                            Annotation.findall_unannotated_metabolite_databases(
                                model,
                                mid;
                                config,
                            ),
                        ) <= config.annotation.maximum_missing_databases
                    end
                end
                @testset "At most $(config.annotation.maximum_nonconformal_references) recognizable cross-references" begin
                    for mid in metabolites(model)
                        @test length(
                            Annotation.findall_nonconformal_metabolite_annotations(
                                model,
                                mid;
                                config,
                            ),
                        ) <= config.annotation.maximum_nonconformal_references
                    end
                end
            end


            @testset "Genes" begin
                @testset "Any annotations" begin
                    for gid in genes(model)
                        @test Annotation.is_annotated_gene(model, gid)
                    end
                end

                @testset "At most $(config.annotation.maximum_missing_databases) missing cross-references" begin
                    for gid in genes(model)
                        @test length(
                            Annotation.findall_unannotated_gene_databases(
                                model,
                                gid;
                                config,
                            ),
                        ) <= config.annotation.maximum_missing_databases
                    end
                end
                @testset "At most $(config.annotation.maximum_nonconformal_references) recognizable cross-references" begin
                    for gid in genes(model)
                        @test length(
                            Annotation.findall_nonconformal_gene_annotations(
                                model,
                                gid;
                                config,
                            ),
                        ) <= config.annotation.maximum_nonconformal_references
                    end
                end
            end
        end

        @testset "Biomass" begin
            @testset "Is present" begin
                @test length(Biomass.findall_biomass_reactions(model)) > 0
            end
            @testset "ATP burn included" begin
                for rid in Biomass.findall_biomass_reactions(model)
                    @test Biomass.atp_present_in_biomass_reaction(model, rid)
                end
            end
            @testset "Molar mass consistent" begin
                for rid in Biomass.findall_biomass_reactions(model)
                    @test Biomass.biomass_reaction_is_consistent(model, rid)
                end
            end
            @testset "Missing common precursors" begin
                for rid in Biomass.findall_biomass_reactions(model)
                    @test Biomass.biomass_missing_essential_precursors(model, rid; config)
                end
            end
        end

        @testset "Stoichiometrically consistent" begin
            @test Consistency.model_is_consistent(model, optimizer; config)
        end

        @testset "Energy metabolism" begin # TODO make this play nicely with compartments
            @testset "Erroneous energy cycles" begin
                @test_skip Energy.model_has_no_erroneous_energy_generating_cycles(
                    model,
                    optimizer;
                    config,
                )
            end
        end

        @testset "GPR associations" begin
            @testset "Metabolic reactions have GPRs" begin
                for rid in filter(x -> is_metabolic_reaction(model, x), reactions(model))
                    @test GPRAssociation.reaction_has_sensible_gpr(model, rid)
                end
            end

            @testset "Transport reactions have GPRs" begin
                for rid in filter(x -> is_transport_reaction(model, x), reactions(model))
                    @test GPRAssociation.reaction_has_sensible_gpr(model, rid)
                end
            end

            @testset "At least one enzyme complex exists" begin
                @test length(GPRAssociation.reactions_with_complexes(model)) != 0
            end
        end

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
                        @test Reactions.reaction_is_charge_balanced(model, rid)
                    end
                end
            end
            @testset "Mass balanced (transport, metabolic, spontaneous)" begin
                for rid in tms_rxns
                    if rid in config.reaction.mass_ignored_reactions
                        @test_skip Reactions.reaction_is_mass_balanced(model, rid)
                    else
                        @test Reactions.reaction_is_mass_balanced(model, rid)
                    end
                end
            end
            @testset "No duplicated reactions" begin
                dup_rxns = Reactions.duplicate_reactions(model)
                for rid in reactions(model)
                    @test rid ∉ dup_rxns
                end
            end
            @testset "ATP maintenance reaction is present" begin
                @test Reactions.model_has_atpm_reaction(model)
            end
        end

        @testset "Metabolite information" begin
            @testset "No missing formulas" begin
                for mid in metabolites(model)
                    @test !isnothing(metabolite_formula(model, mid))
                end
            end
            @testset "No missing charges" begin
                for mid in metabolites(model)
                    @test !isnothing(metabolite_charge(model, mid))
                end
            end
            @testset "Duplicated metabolites" begin
                dup_mids = Metabolites.metabolites_duplicated_in_compartment(model; config)
                for mid in metabolites(model)
                    @test mid ∉ dup_mids
                end
            end
            @testset "No orphan metabolites" begin
                orphans = Metabolites.find_orphan_metabolites(model)
                for mid in metabolites(model)
                    @test mid ∉ orphans
                end
            end
            @testset "No deadend metabolites" begin
                deadends = Metabolites.find_deadend_metabolites(model)
                for mid in metabolites(model)
                    @test mid ∉ deadends
                end
            end
            @testset "No disconnected metabolites" begin
                disconnected = Metabolites.find_disconnected_metabolites(model)
                for mid in metabolites(model)
                    @test mid ∉ disconnected
                end
            end
        end

        @testset "Network information" begin
            @testset "Stoichiometric matrix is well conditioned" begin
                @test Network.stoichiometric_max_min_ratio(model) <
                      config.network.condition_number
            end
            @testset "No stoichiometrically balanced cycles (reactions)" begin
                rxns_in_cycles =
                    Network.find_cycle_reactions(model, optimizer; config, workers)
                for rid in reactions(model)
                    @test rid ∉ rxns_in_cycles
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
                    @test rid ∉ blocked_rxns
                end
            end
        end
    end
end
