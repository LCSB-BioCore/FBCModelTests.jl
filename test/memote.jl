@testset "Consistency" begin
    # consistency
    @test model_is_consistent(model, Tulip.Optimizer)
    temp_model = convert(StandardModel, model)
    temp_model.reactions["PFK"].metabolites["fdp_c"] = 2
    @test !model_is_consistent(temp_model, Tulip.Optimizer)

    # energy cycles
    @test model_has_no_erroneous_energy_generating_cycles(model, Tulip.Optimizer)
    memote_config.consistency.ignored_energy_reactions = ["BIOMASS_KT_TEMP", "ATPM"]
    @test !model_has_no_erroneous_energy_generating_cycles(iJN746, Tulip.Optimizer)

    # use default conditions to exclude biomass and exchanges
    @test isempty(reactions_charge_unbalanced(model))
    @test isempty(reactions_mass_unbalanced(model))

    # test if biomass and exchanges are identified
    wrong_model = convert(StandardModel, model)
    wrong_model.metabolites["pyr_c"].charge = nothing
    wrong_model.metabolites["pyr_c"].formula = "C2H3X"
    @test !isempty(reactions_charge_unbalanced(wrong_model))
    @test !isempty(reactions_mass_unbalanced(wrong_model))

    # test all
    test_consistency(model, Tulip.Optimizer)
end

@testset "Metabolite" begin
    memote_config.metabolite.medium_only_imported = false
    @test "glc__D_e" in metabolites_medium_components(model)

    memote_config.metabolite.medium_only_imported = true
    wrong_model = convert(StandardModel, model)
    wrong_model.reactions["EX_h2o_e"].ub = 0
    @test first(metabolites_medium_components(wrong_model)) == "h2o_e"

    @test isempty(metabolites_no_formula(model))
    wrong_model.metabolites["pyr_c"].formula = ""
    @test !isempty(metabolites_no_formula(wrong_model))
    wrong_model.metabolites["pyr_c"].formula = "C2X"
    @test !isempty(metabolites_no_formula(wrong_model))

    @test isempty(metabolites_no_charge(model))
    wrong_model.metabolites["pyr_c"].charge = nothing
    @test !isempty(metabolites_no_charge(wrong_model))

    @test length(metabolites_unique(model)) == 54
    wrong_model.metabolites["pyr_c"].annotations["inchi_key"] =
        wrong_model.metabolites["etoh_c"].annotations["inchi_key"]
    wrong_model.metabolites["pyr_e"].annotations["inchi_key"] =
        wrong_model.metabolites["etoh_c"].annotations["inchi_key"]
    @test length(metabolites_unique(wrong_model)) == 53

    @test isempty(metabolites_duplicated_in_compartment(model))
    @test !isempty(metabolites_duplicated_in_compartment(wrong_model))

    memote_config.metabolite.medium_only_imported = false
    test_metabolites(model)
end

@testset "Basic" begin
    # these tests are too basic to split out into multiple subtests
    test_basic(model)
end