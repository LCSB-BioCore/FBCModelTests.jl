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

@testset "Metabolite Annotations" begin
    #test all_unannotated_metabolites()
    all_m = all_unannotated_metabolites(model)
    @test isempty(all_m)

    #test unannotated_metabolites()
    u_m = unannotated_metabolites(model)
    @test u_m["kegg.compound"] == ["q8h2_c"]
    @test u_m["biocyc"] == ["icit_c", "fdp_c"]
    @test u_m["hmdb"] == ["q8_c", "r5p_c", "fdp_c"]
    for db in
        ["seed.compound", "inchi_key", "chebi", "metanetx.chemical", "bigg.metabolite"]
        @test isempty(u_m[db])
    end
    for db2 in ["pubchem.compound", "inchi", "reactome"]
        @test length(u_m[db2]) == 72
    end

    #test metabolite_annotation_conformity()
    c_m = metabolite_annotation_conformity(model)
    @test length(c_m) == 8
    for db in [
        "chebi",
        "metanetx.chemical",
        "inchi_key",
        "hmdb",
        "bigg.metabolite",
        "biocyc",
        "seed.compound",
    ]
        @test isempty(c_m[db])
    end
end

@testset "Basic" begin
    # these tests are too basic to split out into multiple subtests
    test_basic(model)
end

@testset "GPR" begin
    @test length(reactions_without_gpr(model)) == 6

    @test length(reactions_with_complexes(model)) == 15

    @test length(reactions_transport_no_gpr(model; config = memote_config)) == 4
end

@testset "Biomass" begin
    @test model_has_atpm_reaction(model)
    wrong_model = convert(StandardModel, model)
    remove_reaction!(wrong_model, "ATPM")
    @test !model_has_atpm_reaction(wrong_model)

    @test all(values(atp_present_in_biomass(model)))

    @test "BIOMASS_Ecoli_core_w_GAM" in model_biomass_reactions(model)

    @test model_biomass_molar_mass(model)["BIOMASS_Ecoli_core_w_GAM"] == 1.5407660614638816
    @test !model_biomass_is_consistent(model)

    @test model_solves_in_default_medium(model, Tulip.Optimizer)

    @test length(
        find_blocked_biomass_precursors(model, Tulip.Optimizer)["BIOMASS_Ecoli_core_w_GAM"],
    ) == 3

    @test length(biomass_missing_essential_precursors(model)["BIOMASS_Ecoli_core_w_GAM"]) ==
          32
end

@testset "Gene Annotations" begin

    # use default conditions to exclude biomass and exchanges
    missing_gene_anno = all_unannotated_genes(model)
    missing_gene_key = unannotated_genes(model)
    gene_anno_confi = gene_annotation_conformity(model)

    #Test1 all_unannotated_genes()
    @test isempty(missing_gene_anno)

    #Test2 unannotated_genes()
    for db in ["ncbiprotein", "ccds", "kegg.genes", "hprd", "refseq"]
        @test length(missing_gene_key[db]) == 137
    end
    for db in ["uniprot", "ecogene", "ncbigi", "ncbigene", "asap"]
        @test missing_gene_key[db] == ["s0001"]
    end

    #Test3 gene_annotation_conformity()
    for db in ["uniprot", "ecogene", "ncbigene", "asap"]
        @test isempty(gene_anno_confi[db])
    end

    @test length(gene_anno_confi["ncbigi"]) == 136

end

@testset "Network" begin
    @test stoichiometric_matrix_is_well_conditioned(model)

    @test isempty(find_all_universally_blocked_reactions(model, Tulip.Optimizer))
    wrong_model = convert(StandardModel, model)
    change_bound!(wrong_model, "MDH"; lower = 0.0, upper = 0.0)
    @test first(find_all_universally_blocked_reactions(wrong_model, Tulip.Optimizer)) ==
          "MDH"

    @test isempty(find_orphan_metabolites(model))
    @test length(find_orphan_metabolites(iJN746)) == 40

    @test isempty(find_deadend_metabolites(model))
    @test length(find_deadend_metabolites(iJN746)) == 51

    crs = find_cycle_reactions(model, Tulip.Optimizer)
    @test "FRD7" in crs && "SUCDi" in crs

    d = find_complete_medium_orphans_and_deadends(model, Tulip.Optimizer)
    @test length(d[:consume]) == 12
    @test length(d[:produce]) == 12
end
