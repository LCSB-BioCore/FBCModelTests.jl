var documenterSearchIndex = {"docs":
[{"location":"functions/#FROG","page":"Reference","title":"FROG","text":"","category":"section"},{"location":"functions/","page":"Reference","title":"Reference","text":"Modules = [FBCModelTests]\nPages = [\"frog.jl\"]","category":"page"},{"location":"functions/#FBCModelTests.ResetObjective","page":"Reference","title":"FBCModelTests.ResetObjective","text":"struct ResetObjective <: COBREXA.ModelWrapper\n\nFields\n\nmodel::COBREXA.MetabolicModel\nobjective::SparseArrays.SparseVector{Float64, Int64}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.frog_compare_reports-Tuple{String, String}","page":"Reference","title":"FBCModelTests.frog_compare_reports","text":"frog_compare_reports(\n    report_dir_a::String,\n    report_dir_b::String\n) -> Union{Test.FallbackTestSet, Test.DefaultTestSet}\n\n\nA simple wrapper for comparing 2 previously generated reports in their respective directories.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_generate_report-Tuple{String}","page":"Reference","title":"FBCModelTests.frog_generate_report","text":"frog_generate_report(\n    filename::String;\n    report_dir,\n    optimizer,\n    workers,\n    basefilename\n)\n\n\nA complete function for one-shot generation of FROG reports. Use frog_model_report, frog_metadata and frog_write_to_directory for finer control of the process.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_metadata-Tuple{String}","page":"Reference","title":"FBCModelTests.frog_metadata","text":"frog_metadata(filename::String; optimizer, basefilename)\n\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_model_report-Tuple{COBREXA.SBMLModel}","page":"Reference","title":"FBCModelTests.frog_model_report","text":"frog_model_report(\n    model::COBREXA.SBMLModel;\n    optimizer,\n    workers\n)\n\n\nGenerate FROGReportData for a model.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_objective_report-Tuple{COBREXA.SBMLModel, String}","page":"Reference","title":"FBCModelTests.frog_objective_report","text":"frog_objective_report(\n    sbml_model::COBREXA.SBMLModel,\n    objective::String;\n    optimizer,\n    workers\n)\n\n\nGenerate a FROGObjectiveReport containing the reproducibility data for a single objective in the SBML model.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_read_from_directory-Tuple{String}","page":"Reference","title":"FBCModelTests.frog_read_from_directory","text":"frog_read_from_directory(\n    report_dir::String\n) -> NamedTuple{(:metadata, :report), _A} where _A<:Tuple{Dict{String, String}, Dict}\n\n\nReverse of frog_write_to_directory.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_test_metadata_compatibility-Tuple{Dict{String, String}, Dict{String, String}}","page":"Reference","title":"FBCModelTests.frog_test_metadata_compatibility","text":"frog_test_metadata_compatibility(\n    a::Dict{String, String},\n    b::Dict{String, String}\n)\n\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_test_report_equality-Tuple{Dict{String, FROGObjectiveReport}, Dict{String, FROGObjectiveReport}}","page":"Reference","title":"FBCModelTests.frog_test_report_equality","text":"frog_test_report_equality(\n    a::Dict{String, FROGObjectiveReport},\n    b::Dict{String, FROGObjectiveReport};\n    absolute_tolerance,\n    relative_tolerance\n) -> Union{Test.FallbackTestSet, Test.DefaultTestSet}\n\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_write_to_directory-Tuple{Dict{String, FROGObjectiveReport}, Dict{String, String}}","page":"Reference","title":"FBCModelTests.frog_write_to_directory","text":"frog_write_to_directory(\n    r::Dict{String, FROGObjectiveReport},\n    metadata::Dict{String, String};\n    report_dir,\n    basefilename\n)\n\n\nWrite the contents of FROGReportData to the 4 TSV files as specified by FROG standard, and additionally write the metadata into the JSON file.\n\n\n\n\n\n","category":"method"},{"location":"functions/#MEMOTE","page":"Reference","title":"MEMOTE","text":"","category":"section"},{"location":"functions/","page":"Reference","title":"Reference","text":"Modules = [FBCModelTests]\nPages = [\"src/memote/basic.jl\", \"src/memote/consistency.jl\", \"src/memote/metabolites.jl\", \"src/memote/reactions.jl\", \"src/memote/gpr_associations.jl\"]\n","category":"page"},{"location":"functions/#FBCModelTests.model_has_no_erroneous_energy_generating_cycles-Tuple{Any, Any}","page":"Reference","title":"FBCModelTests.model_has_no_erroneous_energy_generating_cycles","text":"model_has_no_erroneous_energy_generating_cycles(\n    model,\n    optimizer;\n    config\n) -> Bool\n\n\nAttempts to detect if the model contains any erroneous energy generating cycles by closing all exchange reactions and using flux balance analysis to maximize the sum of fluxes through a set of energy dissipating reactions. The flux sum should be zero if the model is free of energy generating reactions.\n\nThis function is based on Fritzemeier, Claus Jonathan, et al. \"Erroneous energy-generating cycles in published genome scale metabolic networks: Identification and removal.\" PLoS computational biology (2017).\n\nThe energy dissipating reactions are based on the source paper, and include:\n\nATP + H2O --> ADP + H + Phosphate\nCTP + H2O --> CDP + H + Phosphate\nGTP + H2O --> GDP + H + Phosphate\nUTP + H2O --> UDP + H + Phosphate\nITP + H2O --> IDP + H + Phosphate\nNADH --> H + NAD\nNADPH --> H + NADP\nFADH2 --> 2 H + FAD\nFMNH2 --> 2 H + FMN\nUbiquinol-8 --> 2 H + Ubiquinone-8\nMenaquinol-8 --> 2 H + Menaquinone-8\n2-Demethylmenaquinol-8 --> 2 H + 2-Demethylmenaquinone-8\nH2O + ACCOA --> H + Acetate + COA\nL-Glutamate + H2O --> 2-Oxoglutarate + Ammonium + 2 H\nH[external] --> H\n\nAdditional energy dissipating reactions can be directly specified through config.consistency.additional_energy_generating_reactions, which should be vector of COBREXA Reactions using the same metabolite name space as the model. Internally, the model is converted to a COBREXA StandardModel, so ensure that the appropriate accessors are defined for it.\n\nSince models use different name spaces, config.consistency.energy_dissipating_metabolites is used to create the energy dissipating reactions. By default it uses the BiGG name space, but this will be changed to ChEBI in due course. If your model uses a different name space, then you have to change the values (NOT the keys) of config.consistency.energy_dissipating_metabolites. Each energy dissipating reaction is added to the test only if all its associated metabolites are present. Any config.consistency.optimizer_modifications to the solver are passed directly through to COBREXA's flux_balance_analysis function. All config.consistency.boundary_reactions and config.consistency.ignored_energy_reactions are deleted from an internal copy of model; this internal copy is used for analysis.\n\nReturns true if the model has no energy generating cycles.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.model_is_consistent-Tuple{Any, Any}","page":"Reference","title":"FBCModelTests.model_is_consistent","text":"model_is_consistent(model, optimizer; config) -> Bool\n\n\nDetermines if the model is stoichiometrically consistent. Note, stoichiometric consistency does not guarantee that mass balances must hold in the model. A more robust check is reactions_mass_unbalanced, but this works if not all metabolites have mass assigned to them.\n\nBased on Gevorgyan, Albert, Mark G. Poolman, and David A. Fell. \"Detection of stoichiometric inconsistencies in biomolecular models.\" Bioinformatics (2008).\n\nOptionally ignore some reactions in this analysis by adding reaction IDs to config.consistency.consistency_ignored_reactions.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.reactions_charge_unbalanced-Tuple{Any}","page":"Reference","title":"FBCModelTests.reactions_charge_unbalanced","text":"reactions_charge_unbalanced(model; config) -> Vector{String}\n\n\nIterates through all the reactions in model and checks if the charges across each reaction balance. Returns a list of reaction IDs that are charge unbalanced, which is empty if the test passes.\n\nOptionally, use config.consistency.mass_ignored_reactions to pass a vector of reaction ids to ignore in this process. Internally biomass and exchang reactions are ignored.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.reactions_mass_unbalanced-Tuple{Any}","page":"Reference","title":"FBCModelTests.reactions_mass_unbalanced","text":"reactions_mass_unbalanced(model; config) -> Vector{String}\n\n\nIterates through all the reactions in model and checks if the mass across each reaction balances. Returns a list of reaction IDs that are mass unbalanced, which is empty if the test passes.\n\nOptionally, use config.consistency.charge_ignored_reactions to pass a vector of reaction ids to ignore in this process. Internally biomass and exchang reactions are ignored.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.test_consistency-Tuple{Any, Any}","page":"Reference","title":"FBCModelTests.test_consistency","text":"test_consistency(\n    model,\n    optimizer;\n    config\n) -> Union{Test.FallbackTestSet, Test.DefaultTestSet}\n\n\nTest if model is consistent by checking that:\n\nthe model is stoichiometrically consistent, tested with model_is_consistent\nthere are no energy generating cycles, tested with model_has_no_erroneous_energy_generating_cycles\nthe model is both mass and charge balanced, tested with reactions_charge_unbalanced and [reactions_mass_unbalanced]\n\nEach function called in this test function can be called individually. The kwargs are forwarded as indicated by the prefix.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.metabolites_are_duplicated-Tuple{Any, Any, Any}","page":"Reference","title":"FBCModelTests.metabolites_are_duplicated","text":"metabolites_are_duplicated(model, m1, m2; config) -> Any\n\n\nTest if metabolites m1 and m2 are different by comparing their config.metabolite.test_annotation field in the annotations of each metabolite. Note, if no annotations are present for one or both of the metabolites, then return true.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.metabolites_duplicated_in_compartment-Tuple{Any}","page":"Reference","title":"FBCModelTests.metabolites_duplicated_in_compartment","text":"metabolites_duplicated_in_compartment(\n    model;\n    config\n) -> Dict{String, Set{String}}\n\n\nReturn a dictionary of metabolites that are duplicated in their compartment.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.metabolites_medium_components-Tuple{Any}","page":"Reference","title":"FBCModelTests.metabolites_medium_components","text":"metabolites_medium_components(\n    model;\n    config\n) -> Vector{String}\n\n\nReturn a list of all boundary reactions that allow flux into the model and create a metabolite. Assume that boundary reactions only create a single metabolite. Use the testing config, default.metabolite.only_imported = false, to also return metabolites that can be produced by the model under default conditions.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.metabolites_no_charge-Tuple{Any}","page":"Reference","title":"FBCModelTests.metabolites_no_charge","text":"metabolites_no_charge(model; config) -> Any\n\n\nList all metabolites without a charge. Use config.metabolite.charge_corner_cases to specify an extra case to check for charge's that are not properly assigned.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.metabolites_no_formula-Tuple{Any}","page":"Reference","title":"FBCModelTests.metabolites_no_formula","text":"metabolites_no_formula(model; config) -> Any\n\n\nList all metabolites without a formula. Use config.metabolite.formula_corner_cases to specify an extra case to check for formula's that are not properly assigned.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.metabolites_unique-Tuple{Any}","page":"Reference","title":"FBCModelTests.metabolites_unique","text":"metabolites_unique(model; config) -> Set{String}\n\n\nReturn a list of unique metabolites in model. Uses metabolites_are_duplicated internally and forwards test_annotation to it. The latter argument is used to determine if two metabolites are the same by checking for any correspondence.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.test_metabolites-Tuple{Any}","page":"Reference","title":"FBCModelTests.test_metabolites","text":"test_metabolites(\n    model;\n    config\n) -> Union{Test.FallbackTestSet, Test.DefaultTestSet}\n\n\nTest if the metabolites contained in the model:\n\nare not duplicated, tested with metabolites_duplicated_in_compartment\nall have a formula and charge associated with them, tested with metabolites_no_formula and metabolites_no_charge\nthe default medium of the cell is not empty, tested with metabolites_medium_components.\n\nEach of the basic functions can be run independently. THe kwargs are forwarded as indicated by the prefix.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.reactions_without_gpr-Tuple{Any}","page":"Reference","title":"FBCModelTests.reactions_without_gpr","text":"reactions_without_gpr(model) -> Any\n\n\nReturn a list of reaction ids that do not have gene reaction rules (aka gene protein reaction associations).\n\n\n\n\n\n","category":"method"},{"location":"functions/#Types-and-utilities","page":"Reference","title":"Types and utilities","text":"","category":"section"},{"location":"functions/","page":"Reference","title":"Reference","text":"Modules = [FBCModelTests]\nPages = [\"structs.jl\", \"common.jl\"]","category":"page"},{"location":"functions/#FBCModelTests.ObjectiveValue","page":"Reference","title":"FBCModelTests.ObjectiveValue","text":"\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.ConsistencyConfig","page":"Reference","title":"FBCModelTests.ConsistencyConfig","text":"mutable struct ConsistencyConfig\n\nParameters used by the consistency tests.\n\nFields\n\nmass_ignored_reactions::Vector{String}\ncharge_ignored_reactions::Vector{String}\nconsistency_ignored_reactions::Vector{String}\nenergy_dissipating_metabolites::Dict{String, String}\nadditional_energy_generating_reactions::Vector{COBREXA.Reaction}\nignored_energy_reactions::Vector{String}\noptimizer_modifications::Vector{Function}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.FROGMetadata","page":"Reference","title":"FBCModelTests.FROGMetadata","text":"mutable struct Dict{String, String} <: AbstractDict{String, String}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.FROGObjectiveReport","page":"Reference","title":"FBCModelTests.FROGObjectiveReport","text":"struct FROGObjectiveReport\n\nFields\n\noptimum::Union{Nothing, Float64}\nreactions::Dict{String, FROGReactionReport}\ngene_deletions::Dict{String, Union{Nothing, Float64}}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.FROGReactionReport","page":"Reference","title":"FBCModelTests.FROGReactionReport","text":"struct FROGReactionReport\n\nFields\n\nflux::Union{Nothing, Float64}\nvariability_min::Union{Nothing, Float64}\nvariability_max::Union{Nothing, Float64}\ndeletion::Union{Nothing, Float64}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.FROGReportData","page":"Reference","title":"FBCModelTests.FROGReportData","text":"mutable struct Dict{String, FROGObjectiveReport} <: AbstractDict{String, FROGObjectiveReport}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.MemoteConfig","page":"Reference","title":"FBCModelTests.MemoteConfig","text":"mutable struct MemoteConfig\n\nA grouping of parameters used by the metabolic testing infrastructure.\n\nFields\n\nmetabolite::MetaboliteConfig\nconsistency::ConsistencyConfig\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.MetaboliteConfig","page":"Reference","title":"FBCModelTests.MetaboliteConfig","text":"mutable struct MetaboliteConfig\n\nParameters used by the metabolite tests.\n\nFields\n\nformula_corner_cases::Vector{String}\ncharge_corner_cases::Vector{Int64}\nmedium_only_imported::Bool\ntest_annotation::String\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests._has_sensible_gpr-Tuple{Any, Any}","page":"Reference","title":"FBCModelTests._has_sensible_gpr","text":"_has_sensible_gpr(model, rid) -> Any\n\n\nInternal helper function that determines if a reaction has a gene reaction rule and that each gene in the rule is contained in the model.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests._probably_transport_reaction-Tuple{Any, Any, Any}","page":"Reference","title":"FBCModelTests._probably_transport_reaction","text":"_probably_transport_reaction(\n    model,\n    rid,\n    test_annotation\n) -> Bool\n\n\nDetermine if a reaction is probably a transport reaction by checking if:\n\nit has sbo annotations corresponding to a transport reaction\nthe reaction contains metabolites from at least 2 different compartments\nif at least 1 metabolite does not undergo a chemical transformation (via formula or annotation checks)\n\nNote, PTS type transport reactions will be missed if they do not have sbo annotations. This test may yield false negatives.\n\n\n\n\n\n","category":"method"},{"location":"#FBCModelTests.jl-—-testing-and-reproducibility-of-constraint-based-metabolic-modeling","page":"Home","title":"FBCModelTests.jl — testing and reproducibility of constraint-based metabolic modeling","text":"","category":"section"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"functions.md\"]\nDepth = 2","category":"page"}]
}
