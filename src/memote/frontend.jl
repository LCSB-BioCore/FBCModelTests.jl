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
)
    @testset "Metabolic Model Tests$(isnothing(filename) ? "" : ": $filename")" begin

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

Overload of [`generate_report`](@ref) that works directly on a given
filename.
"""
generate_report(filename::String, optimizer; kwargs...) =
    generate_memote_report(load_model(filename), optimizer; kwargs...)

"""
$(TYPEDSIGNATURES)

Generate a machine-readable report of model properties that can be used as a
baseline for checking the quality of the model. Suitable for being exported as
JSON.
"""
generate_report(model::MetabolicModel, optimizer; config = Config.memote_config) = Dict(
    "basic" => Dict(
        "number_reactions" => n_reactions(model),
        "number_metabolites" => n_metabolites(model),
        "number_genes" => n_genes(model),
        "metabolic_coverage" => Basic.model_metabolic_coverage(model),
        "compartments" => Basic.model_compartments(model),
    ),
    "reaction_annotations" => Dict(
        "all_unannotated" => Annotation.find_all_unannotated_reactions(model),
        "missing_databases" =>
            Annotation.find_database_unannotated_reactions(model; config),
        "conformity" =>
            Annotation.find_nonconformal_reaction_annotations(model; config),
    ),
    "metabolite_annotations" => Dict(
        "all_unannotated" => Annotation.find_all_unannotated_metabolites(model),
        "missing_databases" =>
            Annotation.find_database_unannotated_metabolites(model; config),
        "conformity" =>
            Annotation.find_nonconformal_metabolite_annotations(model; config),
    ),
    "gene_annotations" => Dict(
        "all_unannotated" => Annotation.find_all_unannotated_genes(model),
        "missing_databases" =>
            Annotation.find_database_unannotated_genes(model; config),
        "conformity" => Annotation.find_nonconformal_gene_annotations(model; config),
    ),
    "biomass" => Dict(
        "biomas_reactions" => Biomass.model_biomass_reactions(model; config),
        "biomass_molar_masses" => Biomass.model_biomass_molar_mass(model; config),
        "blocked_biomass_precursors" =>
            Biomass.find_blocked_biomass_precursors(model, optimizer; config),
        "missing_essential_precursors_in_biomass_reaction" =>
            Biomass.biomass_missing_essential_precursors(model; config),
    ),
    "consistency" => Dict(
        "mass_unbalanced_reactions" =>
            Consistency.reactions_mass_unbalanced(model; config),
        "charge_unbalanced_reactions" =>
            Consistency.reactions_charge_unbalanced(model; config),
    ),
    "gpr_associations" => Dict(
        "reactions_no_gpr" => GPRAssociation.reactions_without_gpr(model),
        "reactions_with_complexes" => GPRAssociation.reactions_with_complexes(model),
        "transporters_without_gpr" =>
            GPRAssociation.reactions_transport_no_gpr(model; config),
    ),
    "metabolite_information" => Dict(
        "number_unique_metablites" => Metabolite.metabolites_unique(model; config),
        "metabolite_only_imported" =>
            Metabolite.metabolites_medium_components(model; config),
    ),
    "network_topology" => Dict(
        "universally_blocked_reactions" =>
            Network.find_all_universally_blocked_reactions(model, optimizer; config),
        "orphan_metabolites" => Network.find_orphan_metabolites(model),
        "deadend_metabolites" => Network.find_deadend_metabolites(model),
        "reaction_in_stoichiometrically_balanced_cycles" =>
            Network.find_cycle_reactions(model, optimizer; config),
        "metabolites_consumed_produced_complete_medium" =>
            Network.find_complete_medium_orphans_and_deadends(model, optimizer; config),
    ),
    "conditioning" => Dict(
        "stoichiometric_matrix_conditioning" =>
            Network.stoichiometric_max_min_ratio(model),
    ),
    "reaction_information" => let
        (uncon_met, con_met) = Reaction.find_all_purely_metabolic_reactions(model; config)
        (uncon_trans, con_trans) = Reaction.find_all_transport_reactions(model; config)
        Dict(
            "unconstrained_metabolic_reactions" => uncon_met,
            "constrained_metabolic_reactions" => con_met,
            "unconstrained_transporters" => uncon_trans,
            "constrained_transporters" => con_trans,
            "reactions_identical_genes" =>
                Reaction.reactions_with_identical_genes(model),
            "duplicated_reactions" => Reaction.duplicate_reactions(model),
            "reactions_partially_identical_annotations" =>
                Reaction.reactions_with_partially_identical_annotations(model; config),
        )
    end,
)
