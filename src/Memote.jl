"""
# Memote

Main entry point into metabolic model quality tests.
"""
module Memote

using ..DocStringExtensions
using ..COBREXA
using ..Test
using ..JuMP
using ..Utils

include(joinpath("memote_tests", "Config.jl"))
include(joinpath("memote_tests", "Annotation.jl"))
include(joinpath("memote_tests", "Basic.jl"))
include(joinpath("memote_tests", "Biomass.jl"))
include(joinpath("memote_tests", "Consistency.jl"))
include(joinpath("memote_tests", "Energy.jl"))
include(joinpath("memote_tests", "GPRAssociation.jl"))
include(joinpath("memote_tests", "Metabolite.jl"))
include(joinpath("memote_tests", "Network.jl"))
# include(joinpath("memote_tests", "Reaction.jl"))

using .Config
using .Annotation
using .Basic
using .Biomass
using .Consistency
using .Energy
using .GPRAssociation
using .Metabolite
using .Network
# using .Reaction

"""
$(TYPEDSIGNATURES)

Run the standard memote test suite on `model` using `optimizer` to solve
optimization problems and configure the test parameters with `config`.

Note, for best results, make sure the model can be converted into a
`COBREXA.StandardModel`.
"""
function run_tests(model, optimizer; config = memote_config)
    @testset "Metabolic model tests" begin

        @testset "Basic information" begin
            @test model_has_name(model)
            @test model_has_reactions(model)
            @test model_has_metabolites(model)
            @test model_has_genes(model)
            @test model_metabolic_coverage_exceeds_minimum(model; config)
            @test model_has_compartments(model)
        end

#         @testset "Annotations" begin
#             @testset "Reaction annotations" begin
#                 @test isempty(find_all_unannotated_reactions(model))

#                 frac_annotated =
#                     1 -
#                     reduce(
#                         min,
#                         length.(values(find_database_unannotated_reactions(model; config))),
#                     ) / n_reactions(model)
#                 @test config.annotation.minimum_fraction_database_annotations <
#                       frac_annotated

#                 frac_conform =
#                     1 -
#                     reduce(
#                         min,
#                         length.(
#                             values(find_nonconformal_reaction_annotations(model; config))
#                         ),
#                     ) / n_reactions(model)
#                 @test config.annotation.minimum_fraction_database_conformity < frac_conform
#             end

#             @testset "Metabolite annotations" begin
#                 @test isempty(find_all_unannotated_metabolites(model))

#                 frac_annotated =
#                     1 -
#                     reduce(
#                         min,
#                         length.(
#                             values(find_database_unannotated_metabolites(model; config))
#                         ),
#                     ) / n_metabolites(model)
#                 @test config.annotation.minimum_fraction_database_annotations <
#                       frac_annotated

#                 frac_conform =
#                     1 -
#                     reduce(
#                         min,
#                         length.(
#                             values(find_nonconformal_metabolite_annotations(model; config))
#                         ),
#                     ) / n_metabolites(model)
#                 @test config.annotation.minimum_fraction_database_conformity < frac_conform
#             end

#             @testset "Gene annotations" begin
#                 @test isempty(find_all_unannotated_genes(model))

#                 frac_annotated =
#                     1 -
#                     reduce(
#                         min,
#                         length.(values(find_database_unannotated_genes(model; config))),
#                     ) / n_genes(model)
#                 @test config.annotation.minimum_fraction_database_annotations <
#                       frac_annotated

#                 frac_conform =
#                     1 -
#                     reduce(
#                         min,
#                         length.(values(find_nonconformal_gene_annotations(model; config))),
#                     ) / n_genes(model)
#                 @test config.annotation.minimum_fraction_database_conformity < frac_conform
#             end
#         end

#         @testset "Reaction information" begin
#             @test isempty(duplicate_reactions(model))
#         end

#         @testset "Metabolite information" begin
#             @test isempty(metabolites_no_formula(model; config))
#             @test isempty(metabolites_no_charge(model; config))
#             @test isempty(metabolites_duplicated_in_compartment(model; config))
#         end

#         @testset "GPR associations" begin
#             @test length(reactions_without_gpr(model)) != n_reactions(model)
#             @test length(reactions_with_complexes(model)) != 0
#         end

#         @testset "Consistency" begin
#             @testset "Mass and charge balances" begin
#                 @test isempty(reactions_charge_unbalanced(model; config))
#                 @test isempty(reactions_mass_unbalanced(model; config))
#             end

#             @testset "Stoichiometric consistency" begin
#                 @test model_is_consistent(model, optimizer; config)
#             end
#         end

#         @testset "Biomass" begin
#             @test !isempty(model_biomass_reactions(model; config))
#             @test model_has_atpm_reaction(model; config)
#             @test all(values(atp_present_in_biomass(model; config)))
#             @test model_biomass_is_consistent(model; config)
#             @test model_solves_in_default_medium(model, optimizer; config)
#         end

#         @testset "Energy metabolism" begin
#             @testset "Erroneous energy cycles" begin
#                 @test model_has_no_erroneous_energy_generating_cycles(
#                     model,
#                     optimizer;
#                     config,
#                 )
#             end
#         end

#         @testset "Network topology" begin
#             @test isempty(find_all_universally_blocked_reactions(model, optimizer; config))
#             @test isempty(find_orphan_metabolites(model))
#             @test isempty(find_deadend_metabolites(model))
#             @test isempty(
#                 find_complete_medium_orphans_and_deadends(model, optimizer; config),
#             )
#         end

#         @testset "Matrix conditioning" begin
#             @test stoichiometric_matrix_is_well_conditioned(model; config)
#         end

    end
end

# """
# $(TYPEDSIGNATURES)

# Generate a report of model characteristics that are typically important measures
# of the scope of the model.
# """
# function generate_memote_report(model, optimizer; config = memote_config)
#     result = Dict()

#     # Basic information
#     result["basic"] = Dict(
#         "number_reactions" => n_reactions(model),
#         "number_metabolites" => n_metabolites(model),
#         "number_genes" => n_genes(model),
#         "metabolic_coverage" => model_metabolic_coverage(model),
#         "compartments" => model_compartments(model),
#     )

#     # Reaction annotations
#     result["reaction_annotations"] = Dict(
#         "all_unannotated" => find_all_unannotated_reactions(model),
#         "missing_databases" => find_database_unannotated_reactions(model; config),
#         "conformity" => find_nonconformal_reaction_annotations(model; config),
#     )

#     # Metabolite annotations
#     result["metabolite_annotations"] = Dict(
#         "all_unannotated" => find_all_unannotated_metabolites(model),
#         "missing_databases" => find_database_unannotated_metabolites(model; config),
#         "conformity" => find_nonconformal_metabolite_annotations(model; config),
#     )

#     # Gene annotations
#     result["gene_annotations"] = Dict(
#         "all_unannotated" => find_all_unannotated_genes(model),
#         "missing_databases" => find_database_unannotated_genes(model; config),
#         "conformity" => find_nonconformal_gene_annotations(model; config),
#     )

#     # Reaction information
#     uncon_met, con_met = find_all_purely_metabolic_reactions(model; config)
#     uncon_trans, con_trans = find_all_transport_reactions(model; config)
#     result["reaction_information"] => Dict(
#         "unconstrained_metabolic_reactions" => uncon_met,
#         "constrained_metabolic_reactions" => con_met,
#         "unconstrained_transporters" => uncon_trans,
#         "constrained_transporters" => con_trans,
#         "reactions_identical_genes" => reactions_with_identical_genes(model),
#         "duplicated_reactions" => duplicate_reactions(model),
#         "reactions_partially_identical_annotations" =>
#             reactions_with_partially_identical_annotations(model; config),
#     )

#     # Metabolite information
#     result["metabolite_information"] = Dict(
#         "number_unique_metablites" => metabolites_unique(model; config),
#         "metabolite_only_imported" => metabolites_medium_components(model; config),
#     )

#     # Gene protein reaction associations
#     result["gpr_associations"] = Dict(
#         "reactions_no_gpr" => reactions_without_gpr(model),
#         "reactions_with_complexes" => reactions_with_complexes(model),
#         "transporters_without_gpr" => reactions_transport_no_gpr(model; config),
#     )

#     # Consistency
#     result["consistency"] = Dict(
#         "mass_unbalanced_reactions" => reactions_mass_unbalanced(model; config),
#         "charge_unbalanced_reactions" => reactions_charge_unbalanced(model; config),
#     )

#     # Biomass
#     result["biomass"] = Dict(
#         "biomas_reactions" => model_biomass_reactions(model; config),
#         "biomass_molar_masses" => model_biomass_molar_mass(model; config),
#         "blocked_biomass_precursors" =>
#             find_blocked_biomass_precursors(model, optimizer; config),
#         "missing_essential_precursors_in_biomass_reaction" =>
#             biomass_missing_essential_precursors(model; config),
#     )

#     # Network topology
#     result["network_topology"] = Dict(
#         "universally_blocked_reactions" =>
#             find_all_universally_blocked_reactions(model, optimizer; config),
#         "orphan_metabolites" => find_orphan_metabolites(model),
#         "deadend_metabolites" => find_deadend_metabolites(model),
#         "reaction_in_stoichiometrically_balanced_cycles" =>
#             find_cycle_reactions(model, optimizer; config),
#         "metabolites_consumed_produced_complete_medium" =>
#             find_complete_medium_orphans_and_deadends(model, optimizer; config),
#     )

#     # Matrix conditioning
#     result["conditioning"] =
#         Dict("stoichiometric_matrix_conditioning" => stoichiometric_max_min_ratio(model))

#     return result
# end

@export_locals

end # module