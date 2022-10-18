"""
$(TYPEDSIGNATURES)

Run the standard memote test suite on `model` using `optimizer` to solve
optimization problems and configure the test parameters with `config`.

Note, for best results, make sure the model can be converted into a
`COBREXA.StandardModel`.
"""
function run_tests(model, optimizer; config = Config.memote_config)
    @testset "Metabolic model tests" begin

        @testset "Basic information" begin
            @test Basic.model_has_name(model)
            @test Basic.model_has_reactions(model)
            @test Basic.model_has_metabolites(model)
            @test Basic.model_has_genes(model)
            @test Basic.model_metabolic_coverage_exceeds_minimum(model; config)
            @test Basic.model_has_compartments(model)
            @test Basic.model_solves_in_default_medium(model, optimizer; config)
        end

        @testset "Annotations" begin
            @testset "Reaction annotations" begin
                @test isempty(Annotation.find_all_unannotated_reactions(model))

                num_annotated = reduce(
                    max,
                    length.(
                        values(
                            Annotation.find_database_unannotated_reactions(model; config),
                        )
                    ),
                )
                @test config.annotation.minimum_fraction_database_annotations *
                      n_reactions(model) < num_annotated

                num_conform = reduce(
                    max,
                    length.(
                        values(
                            Annotation.find_nonconformal_reaction_annotations(
                                model;
                                config,
                            ),
                        )
                    ),
                )
                @test config.annotation.minimum_fraction_database_conformity *
                      n_reactions(model) < num_conform
            end

            @testset "Metabolite annotations" begin
                @test isempty(Annotation.find_all_unannotated_metabolites(model))

                num_annotated = reduce(
                    max,
                    length.(
                        values(
                            Annotation.find_database_unannotated_metabolites(model; config),
                        )
                    ),
                )
                @test config.annotation.minimum_fraction_database_annotations *
                      n_metabolites(model) < num_annotated

                num_conform = reduce(
                    max,
                    length.(
                        values(
                            Annotation.find_nonconformal_metabolite_annotations(
                                model;
                                config,
                            ),
                        )
                    ),
                )
                @test config.annotation.minimum_fraction_database_conformity *
                      n_metabolites(model) < num_conform
            end

            @testset "Gene annotations" begin
                @test isempty(Annotation.find_all_unannotated_genes(model))

                num_annotated = reduce(
                    max,
                    length.(
                        values(Annotation.find_database_unannotated_genes(model; config))
                    ),
                )
                @test config.annotation.minimum_fraction_database_annotations *
                      n_genes(model) < num_annotated

                num_conform = reduce(
                    max,
                    length.(
                        values(
                            Annotation.find_nonconformal_gene_annotations(model; config),
                        )
                    ),
                )
                @test config.annotation.minimum_fraction_database_conformity *
                      n_genes(model) < num_conform
            end
        end

        @testset "Biomass" begin
            @test !isempty(Biomass.model_biomass_reactions(model; config))
            @test all(values(Biomass.atp_present_in_biomass(model; config)))
            @test Biomass.model_biomass_is_consistent(model; config)
        end

        @testset "Consistency" begin
            @testset "Mass and charge balances" begin
                @test isempty(Consistency.reactions_charge_unbalanced(model; config))
                @test isempty(Consistency.reactions_mass_unbalanced(model; config))
            end

            @testset "Stoichiometric consistency" begin
                @test Consistency.model_is_consistent(model, optimizer; config)
            end
        end

        @testset "Energy metabolism" begin
            @testset "Erroneous energy cycles" begin
                @test Energy.model_has_no_erroneous_energy_generating_cycles(
                    model,
                    optimizer;
                    config,
                )
            end
            @test Energy.model_has_atpm_reaction(model; config)
        end

        @testset "GPR associations" begin
            @test length(GPRAssociation.reactions_without_gpr(model)) != n_reactions(model)
            @test length(GPRAssociation.reactions_with_complexes(model)) != 0
        end

        @testset "Metabolite information" begin
            @test isempty(Metabolite.metabolites_no_formula(model; config))
            @test isempty(Metabolite.metabolites_no_charge(model; config))
            @test isempty(Metabolite.metabolites_duplicated_in_compartment(model; config))
        end

        @testset "Network topology" begin
            @test isempty(
                Network.find_all_universally_blocked_reactions(model, optimizer; config),
            )
            @test isempty(Network.find_orphan_metabolites(model))
            @test isempty(Network.find_deadend_metabolites(model))
            @test isempty(
                Network.find_complete_medium_orphans_and_deadends(model, optimizer; config),
            )
        end

        @testset "Matrix conditioning" begin
            @test Network.stoichiometric_matrix_is_well_conditioned(model; config)
        end

        @testset "Reaction information" begin
            @test isempty(Reaction.duplicate_reactions(model))
        end

    end
end

"""
$(TYPEDSIGNATURES)

Generate a report of model characteristics that are typically important measures
of the scope of the model.
"""
function generate_memote_report(model, optimizer; config = Config.memote_config)
    result = Dict()

    # Basic information
    result["basic"] = Dict(
        "number_reactions" => n_reactions(model),
        "number_metabolites" => n_metabolites(model),
        "number_genes" => n_genes(model),
        "metabolic_coverage" => Basic.model_metabolic_coverage(model),
        "compartments" => Basic.model_compartments(model),
    )

    # Reaction annotations
    result["reaction_annotations"] = Dict(
        "all_unannotated" => Annotation.find_all_unannotated_reactions(model),
        "missing_databases" =>
            Annotation.find_database_unannotated_reactions(model; config),
        "conformity" =>
            Annotation.find_nonconformal_reaction_annotations(model; config),
    )

    # Metabolite annotations
    result["metabolite_annotations"] = Dict(
        "all_unannotated" => Annotation.find_all_unannotated_metabolites(model),
        "missing_databases" =>
            Annotation.find_database_unannotated_metabolites(model; config),
        "conformity" =>
            Annotation.find_nonconformal_metabolite_annotations(model; config),
    )

    # Gene annotations
    result["gene_annotations"] = Dict(
        "all_unannotated" => Annotation.find_all_unannotated_genes(model),
        "missing_databases" =>
            Annotation.find_database_unannotated_genes(model; config),
        "conformity" => Annotation.find_nonconformal_gene_annotations(model; config),
    )

    # Biomass
    result["biomass"] = Dict(
        "biomas_reactions" => Biomass.model_biomass_reactions(model; config),
        "biomass_molar_masses" => Biomass.model_biomass_molar_mass(model; config),
        "blocked_biomass_precursors" =>
            Biomass.find_blocked_biomass_precursors(model, optimizer; config),
        "missing_essential_precursors_in_biomass_reaction" =>
            Biomass.biomass_missing_essential_precursors(model; config),
    )

    # Consistency
    result["consistency"] = Dict(
        "mass_unbalanced_reactions" =>
            Consistency.reactions_mass_unbalanced(model; config),
        "charge_unbalanced_reactions" =>
            Consistency.reactions_charge_unbalanced(model; config),
    )

    # Gene protein reaction associations
    result["gpr_associations"] = Dict(
        "reactions_no_gpr" => GPRAssociation.reactions_without_gpr(model),
        "reactions_with_complexes" => GPRAssociation.reactions_with_complexes(model),
        "transporters_without_gpr" =>
            GPRAssociation.reactions_transport_no_gpr(model; config),
    )

    # Metabolite information
    result["metabolite_information"] = Dict(
        "number_unique_metablites" => Metabolite.metabolites_unique(model; config),
        "metabolite_only_imported" =>
            Metabolite.metabolites_medium_components(model; config),
    )

    # Network topology
    result["network_topology"] = Dict(
        "universally_blocked_reactions" =>
            Network.find_all_universally_blocked_reactions(model, optimizer; config),
        "orphan_metabolites" => Network.find_orphan_metabolites(model),
        "deadend_metabolites" => Network.find_deadend_metabolites(model),
        "reaction_in_stoichiometrically_balanced_cycles" =>
            Network.find_cycle_reactions(model, optimizer; config),
        "metabolites_consumed_produced_complete_medium" =>
            Network.find_complete_medium_orphans_and_deadends(model, optimizer; config),
    )

    # Matrix conditioning
    result["conditioning"] = Dict(
        "stoichiometric_matrix_conditioning" =>
            Network.stoichiometric_max_min_ratio(model),
    )

    # Reaction information
    uncon_met, con_met = Reaction.find_all_purely_metabolic_reactions(model; config)
    uncon_trans, con_trans = Reaction.find_all_transport_reactions(model; config)
    result["reaction_information"] = Dict(
        "unconstrained_metabolic_reactions" => uncon_met,
        "constrained_metabolic_reactions" => con_met,
        "unconstrained_transporters" => uncon_trans,
        "constrained_transporters" => con_trans,
        "reactions_identical_genes" => Reaction.reactions_with_identical_genes(model),
        "duplicated_reactions" => Reaction.duplicate_reactions(model),
        "reactions_partially_identical_annotations" =>
            Reaction.reactions_with_partially_identical_annotations(model; config),
    )

    return result
end
