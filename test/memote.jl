@testset "Consistency" begin
    # consistency
    @test is_consistent(model, Tulip.Optimizer)
    temp_model = convert(StandardModel, model)
    temp_model.reactions["PFK"].metabolites["fdp_c"] = 2
    @test !is_consistent(temp_model, Tulip.Optimizer)

    # energy cycles
    @test !has_erroneous_energy_generating_cycles(model, Tulip.Optimizer)
    @test has_erroneous_energy_generating_cycles(
        iJN746,
        Tulip.Optimizer;
        ignored_reactions = ["BIOMASS_KT_TEMP"],
    )

    # use default conditions to exclude biomass and exchanges
    cbal = reaction_charge_unbalanced(model)
    mbal = reaction_mass_unbalanced(model)
    @test isempty(cbal)
    @test isempty(mbal)

    # test if biomass and exchanges are identified
    cbal = reaction_charge_unbalanced(model; ignored_reactions = [])
    mbal = reaction_mass_unbalanced(model; ignored_reactions = [])
    @test "BIOMASS_Ecoli_core_w_GAM" in cbal
    @test "BIOMASS_Ecoli_core_w_GAM" in mbal
    @test all(startswith(rid, "EX") for rid in cbal if rid != "BIOMASS_Ecoli_core_w_GAM")
    @test all(startswith(rid, "EX") for rid in mbal if rid != "BIOMASS_Ecoli_core_w_GAM")

    # test all 
    test_consistency(model, Tulip.Optimizer)
end

@testset "Metabolite" begin
    test_metabolites(model; medium_only_into=false)
end
