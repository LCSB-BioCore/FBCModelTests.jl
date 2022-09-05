@testset "Basic tests" begin

    # use default conditions to exclude biomass and exchanges
    cbal = is_model_charge_balanced(model)
    mbal = is_model_mass_balanced(model)
    @test isempty(cbal)
    @test isempty(mbal)

    # test if biomass and exchanges are identified
    cbal = is_model_charge_balanced(model; ignored_reactions = [])
    mbal = is_model_mass_balanced(model; ignored_reactions = [])
    @test "BIOMASS_Ecoli_core_w_GAM" in cbal
    @test "BIOMASS_Ecoli_core_w_GAM" in mbal
    @test all(startswith(rid, "EX") for rid in cbal if rid != "BIOMASS_Ecoli_core_w_GAM")
    @test all(startswith(rid, "EX") for rid in mbal if rid != "BIOMASS_Ecoli_core_w_GAM")
end

@testset "Erroneous energy generation" begin
    @test !has_erroneous_energy_generating_cycles(model, Tulip.Optimizer)
    @test has_erroneous_energy_generating_cycles(iJN746, Tulip.Optimizer; ignore_reactions=["BIOMASS_KT_TEMP",])
end
