@testset "Basic" begin
    @test !Basic.model_has_name(model) # TODO without accessors to JSONModel, this should fail
    @test Basic.model_has_metabolites(model)
    @test Basic.model_has_reactions(model)
    @test Basic.model_has_genes(model)
    @test Basic.model_metabolic_coverage_exceeds_minimum(model)
    @test length(Basic.model_compartments(model)) == 2
    @test Basic.model_has_compartments(model)
    @test Basic.model_solves_in_default_medium(model, Tulip.Optimizer)
end

@testset "Annotations" begin
    # identify unannotated components
    all_mets_no_anno = Annotation.find_all_unannotated_metabolites(model)
    @test isempty(all_mets_no_anno)

    all_rxns_no_anno = Annotation.find_all_unannotated_reactions(model)
    @test isempty(all_rxns_no_anno)

    all_gene_no_anno = Annotation.find_all_unannotated_genes(model)
    @test isempty(all_gene_no_anno)

    # find databases missing in annotations
    gene_databases = Annotation.find_database_unannotated_genes(model)
    @test length(gene_databases["uniprot"]) == 1
    @test length(gene_databases["ncbiprotein"]) == 137

    met_databases = Annotation.find_database_unannotated_metabolites(model)
    @test length(met_databases["kegg.compound"]) == 1
    @test length(met_databases["biocyc"]) == 2

    rxn_databases = Annotation.find_database_unannotated_reactions(model)
    @test length(rxn_databases["rhea"]) == 33
    @test length(rxn_databases["ec-code"]) == 44

    # find nonconformal annotations
    gene_nonconformal = Annotation.find_nonconformal_gene_annotations(model)
    # TODO fix regexes

    met_nonconformal = Annotation.find_nonconformal_metabolite_annotations(model)
    # TODO fix regexes

    rxn_nonconformal = Annotation.find_nonconformal_reaction_annotations(model)
    # TODO fix regexes
end

@testset "Biomass" begin
    @test all(values(Biomass.atp_present_in_biomass(model)))
    @test "BIOMASS_Ecoli_core_w_GAM" in Biomass.model_biomass_reactions(model)

    @test Biomass.model_biomass_molar_mass(model)["BIOMASS_Ecoli_core_w_GAM"] ==
          1.5407660614638816
    @test !Biomass.model_biomass_is_consistent(model)

    @test length(
        Biomass.find_blocked_biomass_precursors(model, Tulip.Optimizer)["BIOMASS_Ecoli_core_w_GAM"],
    ) == 3

    @test length(
        Biomass.biomass_missing_essential_precursors(model)["BIOMASS_Ecoli_core_w_GAM"],
    ) == 32
end

@testset "Consistency" begin
    # consistency
    @test Consistency.model_is_consistent(model, Tulip.Optimizer)
    temp_model = convert(StandardModel, model)
    temp_model.reactions["PFK"].metabolites["fdp_c"] = 2
    @test !Consistency.model_is_consistent(temp_model, Tulip.Optimizer)

    # use default conditions to exclude biomass and exchanges
    @test isempty(Consistency.reactions_charge_unbalanced(model))
    @test isempty(Consistency.reactions_mass_unbalanced(model))

    # test if biomass and exchanges are identified
    wrong_model = convert(StandardModel, model)
    wrong_model.metabolites["pyr_c"].charge = nothing
    wrong_model.metabolites["pyr_c"].formula = "C2H3X"
    @test !isempty(Consistency.reactions_charge_unbalanced(wrong_model))
    @test !isempty(Consistency.reactions_mass_unbalanced(wrong_model))

    # test metabolite connectivity
    dm = Consistency.find_disconnected_metabolites(model)
    @test isempty(dm)

    # test unbounded flux
    fva_result = flux_variability_analysis_dict(
        model,
        Tulip.Optimizer;
        bounds = objective_bounds(0.99),
    )
    mb = median_bounds(model)
    low_unlimited_flux, high_unlimited_flux =
        Consistency.unbounded_flux_in_default_medium(model, fva_result)
    @test mb == (-1000.0, 1000.0)
    @test isempty(low_unlimited_flux)
    @test isapprox(high_unlimited_flux["FRD7"][2], 1000.0, atol = 1e-07)
end

@testset "Energy metabolism" begin
    @test Energy.model_has_atpm_reaction(model)
    wrong_model = convert(StandardModel, model)
    remove_reaction!(wrong_model, "ATPM")
    @test !Energy.model_has_atpm_reaction(wrong_model)

    # energy cycles
    @test Energy.model_has_no_erroneous_energy_generating_cycles(model, Tulip.Optimizer)
    memote_config.energy.ignored_energy_reactions = ["BIOMASS_KT_TEMP", "ATPM"]
    @test !Energy.model_has_no_erroneous_energy_generating_cycles(iJN746, Tulip.Optimizer)
end

@testset "GPR" begin
    @test length(GPRAssociation.reactions_without_gpr(model)) == 6
    @test length(GPRAssociation.reactions_with_complexes(model)) == 15
    @test length(GPRAssociation.reactions_transport_no_gpr(model)) == 4
end

@testset "Metabolite" begin
    memote_config.metabolite.medium_only_imported = false
    @test "glc__D_e" in FBCModelTests.Memote.Metabolite.metabolites_medium_components(
        model;
        config = memote_config,
    )

    wrong_model = convert(StandardModel, model)
    wrong_model.reactions["EX_h2o_e"].ub = 0
    @test "h2o_e" in
          FBCModelTests.Memote.Metabolite.metabolites_medium_components(wrong_model)

    @test isempty(FBCModelTests.Memote.Metabolite.metabolites_no_formula(model))
    wrong_model.metabolites["pyr_c"].formula = ""
    @test !isempty(FBCModelTests.Memote.Metabolite.metabolites_no_formula(wrong_model))
    wrong_model.metabolites["pyr_c"].formula = "C2X"
    @test !isempty(FBCModelTests.Memote.Metabolite.metabolites_no_formula(wrong_model))

    @test isempty(FBCModelTests.Memote.Metabolite.metabolites_no_charge(model))
    wrong_model.metabolites["pyr_c"].charge = nothing
    @test !isempty(FBCModelTests.Memote.Metabolite.metabolites_no_charge(wrong_model))

    @test length(FBCModelTests.Memote.Metabolite.metabolites_unique(model)) == 54
    wrong_model.metabolites["pyr_c"].annotations["inchi_key"] =
        wrong_model.metabolites["etoh_c"].annotations["inchi_key"]
    wrong_model.metabolites["pyr_e"].annotations["inchi_key"] =
        wrong_model.metabolites["etoh_c"].annotations["inchi_key"]
    @test length(FBCModelTests.Memote.Metabolite.metabolites_unique(wrong_model)) == 53

    @test isempty(
        FBCModelTests.Memote.Metabolite.metabolites_duplicated_in_compartment(model),
    )
    @test !isempty(
        FBCModelTests.Memote.Metabolite.metabolites_duplicated_in_compartment(wrong_model),
    )
end

@testset "Network" begin
    @test Network.stoichiometric_matrix_is_well_conditioned(model)

    @test isempty(Network.find_all_universally_blocked_reactions(model, Tulip.Optimizer))
    wrong_model = convert(StandardModel, model)
    change_bound!(wrong_model, "MDH"; lower = 0.0, upper = 0.0)
    @test first(
        Network.find_all_universally_blocked_reactions(wrong_model, Tulip.Optimizer),
    ) == "MDH"

    @test isempty(Network.find_orphan_metabolites(model))
    @test length(Network.find_orphan_metabolites(iJN746)) == 40

    @test isempty(Network.find_deadend_metabolites(model))
    @test length(Network.find_deadend_metabolites(iJN746)) == 51

    crs = Network.find_cycle_reactions(model, Tulip.Optimizer)
    @test "FRD7" in crs && "SUCDi" in crs

    d = Network.find_complete_medium_orphans_and_deadends(model, Tulip.Optimizer)
    @test length(d[:consume]) == 12
    @test length(d[:produce]) == 12
end

@testset "Reactions" begin
    ident_grrs = FBCModelTests.Memote.Reaction.reactions_with_identical_genes(model)
    @test length(ident_grrs) == 11
    @test issetequal(ident_grrs[["b1602", "b1603"]], ["NADTRHD", "THD2"])

    metabolic_reactions_unconstrained, metabolic_reactions_constrained =
        FBCModelTests.Memote.Reaction.find_all_purely_metabolic_reactions(model)
    @test length(metabolic_reactions_unconstrained) == 50
    @test length(metabolic_reactions_constrained) == 1

    transport_unconstrained, transport_constrained =
        FBCModelTests.Memote.Reaction.find_all_transport_reactions(model)
    @test length(transport_unconstrained) == 23
    @test length(transport_constrained) == 0

    @test length(
        FBCModelTests.Memote.Reaction.reactions_with_partially_identical_annotations(model),
    ) == 14

    @test issetequal(
        FBCModelTests.Memote.Reaction.duplicate_reactions(model),
        ["FRD7", "SUCDi"],
    )
end
