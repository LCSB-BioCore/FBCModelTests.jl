
using JSON
using GLPK

using FBCModelTests.Memote
using FBCModelTests.Memote.Annotation
using FBCModelTests.Memote.Basic
using FBCModelTests.Memote.Biomass
using FBCModelTests.Memote.Consistency
using FBCModelTests.Memote.Energy
using FBCModelTests.Memote.GPRAssociation
using FBCModelTests.Memote.Network
using FBCModelTests.Memote.Reactions
using FBCModelTests.Memote.Metabolites

@testset "Front-end" begin
    result = @testset CountTests "Testing a model that _should_ be OK" begin
        @test Memote.run_tests(
            load_model(StandardModel, model_file["e_coli_core.json"]),
            GLPK.Optimizer,
        )
    end

    @test result.passes == 1783
    @test result.fails == 12
    @test result.errs == 1 # broken counts here I think...
end

@testset "Basic" begin
    @test Basic.model_has_metabolites(model)
    @test Basic.model_has_reactions(model)
    @test Basic.model_has_genes(model)
    @test isapprox(
        Basic.model_metabolic_coverage(model),
        0.6934306569343066;
        atol = TEST_TOL,
    )
    @test length(Basic.model_compartments(model)) == 2
    @test Basic.model_has_compartments(model)
    @test Basic.model_solves_in_default_medium(model, GLPK.Optimizer)

    # Negative tests
    empty_model = StandardModel("test") # should fail all of these tests as it has no reactions, genes or metabolites
    @test !Basic.model_has_metabolites(empty_model)
    @test !Basic.model_has_reactions(empty_model)
    @test !Basic.model_has_genes(empty_model)
    @test isnan(Basic.model_metabolic_coverage(empty_model))
    @test isempty(Basic.model_compartments(empty_model))
    @test !Basic.model_has_compartments(empty_model)
    @test !Basic.model_solves_in_default_medium(empty_model, GLPK.Optimizer)
end

@testset "Annotations" begin
    # identify unannotated components
    mid = "akg_c"
    rid = "ENO"
    gid = "b1478"

    @test Annotation.is_annotated_reaction(model, rid)
    @test Annotation.is_annotated_metabolite(model, mid)
    @test Annotation.is_annotated_gene(model, gid)

    @test "ncbiprotein" in Annotation.findall_unannotated_gene_databases(model, gid)
    @test "inchi" in Annotation.findall_unannotated_metabolite_databases(model, mid)
    @test "brenda" in Annotation.findall_unannotated_reaction_databases(model, rid)

    @test "ncbigi" in Annotation.findall_nonconformal_gene_annotations(model, gid)
    @test "reactome.compound" in
          Annotation.findall_nonconformal_metabolite_annotations(model, mid)
    @test isempty(Annotation.findall_nonconformal_reaction_annotations(model, rid))
end

@testset "Biomass" begin
    @test Biomass.atp_present_in_biomass_reaction(model, "BIOMASS_Ecoli_core_w_GAM")
    @test "BIOMASS_Ecoli_core_w_GAM" in Biomass.findall_biomass_reactions(model)

    @test isapprox(
        Biomass.biomass_reaction_molar_mass(model, "BIOMASS_Ecoli_core_w_GAM"),
        1.5407660614638816;
        atol = TEST_TOL,
    )
    @test !Biomass.biomass_reaction_is_consistent(model, "BIOMASS_Ecoli_core_w_GAM")

    @test Biomass.biomass_missing_essential_precursors(model, "BIOMASS_Ecoli_core_w_GAM")

    # Negative tests
    broken_model = convert(StandardModel, deepcopy(model))

    delete!(broken_model.reactions["BIOMASS_Ecoli_core_w_GAM"].metabolites, "atp_c")
    @test !Biomass.atp_present_in_biomass_reaction(broken_model, "BIOMASS_Ecoli_core_w_GAM")

    delete!(broken_model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    @test isempty(Biomass.findall_biomass_reactions(broken_model))

    # load a model that is better specified: iML1515
    @test isapprox(
        Biomass.biomass_reaction_molar_mass(iML1515, "BIOMASS_Ec_iML1515_WT_75p37M"),
        0.9992745719378807;
        atol = TEST_TOL,
    )
    @test !Biomass.biomass_missing_essential_precursors(
        iML1515,
        "BIOMASS_Ec_iML1515_core_75p37M",
    )
end

@testset "Consistency" begin
    @test Consistency.model_is_consistent(model, GLPK.Optimizer)
    temp_model = convert(StandardModel, deepcopy(model))
    temp_model.reactions["PFK"].metabolites["fdp_c"] = 2
    @test !Consistency.model_is_consistent(temp_model, GLPK.Optimizer)
end

@testset "Energy metabolism" begin
    @test_skip Energy.model_has_no_erroneous_energy_generating_cycles(model, GLPK.Optimizer)
    @test_skip !Energy.model_has_no_erroneous_energy_generating_cycles(
        iJN746,
        GLPK.Optimizer,
    )
end

@testset "GPR" begin
    @test length(GPRAssociation.reactions_with_complexes(model)) == 15
    # fix the model to improve stats
    new_model = convert(StandardModel, deepcopy(model))
    new_model.reactions["GLUSy"].grr = [["b2097"]]
    @test length(GPRAssociation.reactions_with_complexes(new_model)) == 14
end

@testset "Metabolite" begin
    @test isempty(Metabolites.metabolites_duplicated_in_compartment(model))
    @test isempty(Metabolites.find_orphan_metabolites(model))
    @test isempty(Metabolites.find_deadend_metabolites(model))
    @test isempty(Metabolites.find_disconnected_metabolites(model))

    negative_model = convert(StandardModel, deepcopy(model))
    dup = deepcopy(negative_model.metabolites["atp_c"])
    dup.id = "atp_dup"
    add_metabolite!(negative_model, dup)
    @test "atp_dup" in Metabolites.metabolites_duplicated_in_compartment(negative_model)
    @test !isempty(Metabolites.find_orphan_metabolites(iJN746))
    @test !isempty(Metabolites.find_deadend_metabolites(iJN746))
    @test !isempty(Metabolites.find_disconnected_metabolites(negative_model))
end

@testset "Network" begin
    @test isapprox(
        Network.stoichiometric_max_min_ratio(model),
        843.5825105782792;
        atol = TEST_TOL,
    )
    @test isempty(Network.find_all_universally_blocked_reactions(model, GLPK.Optimizer))
    crs = Network.find_cycle_reactions(model, GLPK.Optimizer)
    @test "FRD7" in crs && "SUCDi" in crs

    # negative tests
    negative_model = convert(StandardModel, deepcopy(model))

    change_bound!(negative_model, "MDH"; lower = 0.0, upper = 0.0)
    @test first(
        Network.find_all_universally_blocked_reactions(negative_model, GLPK.Optimizer),
    ) == "MDH"

    negative_model.reactions["FRD7"].ub = 0
    negative_model.reactions["FRD7"].lb = -1000
    crs = Network.find_cycle_reactions(negative_model, GLPK.Optimizer)
    @test isempty(crs)
end

@testset "Reactions" begin
    @test issetequal(Reactions.findall_duplicated_reactions(model), ["FRD7", "SUCDi"]) # TODO should be empty, will be broken once #785 in COBREXA gets fixed

    @test !Reactions.model_has_atpm_reaction(model)
    @test Reactions.reaction_is_charge_balanced(model, "ENO")
    @test Reactions.reaction_is_mass_balanced(model, "ENO")

    negative_model = convert(StandardModel, deepcopy(model))
    push!(negative_model.reactions["ATPM"].annotations["sbo"], "SBO:0000630")
    @test Reactions.model_has_atpm_reaction(negative_model)
    negative_model.metabolites["2pg_c"].charge = nothing
    negative_model.metabolites["2pg_c"].formula = "C2H3X"
    @test !Reactions.reaction_is_charge_balanced(negative_model, "ENO")
    @test !Reactions.reaction_is_mass_balanced(negative_model, "ENO")
end
