var documenterSearchIndex = {"docs":
[{"location":"functions/#FROG","page":"Reference","title":"FROG","text":"","category":"section"},{"location":"functions/","page":"Reference","title":"Reference","text":"Modules = [FBCModelTests]\nPages = [\"frog.jl\"]","category":"page"},{"location":"functions/#FBCModelTests.ResetObjective","page":"Reference","title":"FBCModelTests.ResetObjective","text":"struct ResetObjective <: COBREXA.ModelWrapper\n\nFields\n\nmodel::COBREXA.MetabolicModel\nobjective::SparseArrays.SparseVector{Float64, Int64}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.frog_compare_reports-Tuple{String, String}","page":"Reference","title":"FBCModelTests.frog_compare_reports","text":"frog_compare_reports(\n    report_dir_a::String,\n    report_dir_b::String\n) -> Union{Test.FallbackTestSet, Test.DefaultTestSet}\n\n\nA simple wrapper for comparing 2 previously generated reports in their respective directories.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_generate_report-Tuple{String}","page":"Reference","title":"FBCModelTests.frog_generate_report","text":"frog_generate_report(\n    filename::String;\n    report_dir,\n    optimizer,\n    workers,\n    basefilename\n)\n\n\nA complete function for one-shot generation of FROG reports. Use frog_model_report, frog_metadata and frog_write_to_directory for finer control of the process.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_metadata-Tuple{String}","page":"Reference","title":"FBCModelTests.frog_metadata","text":"frog_metadata(filename::String; optimizer, basefilename)\n\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_model_report-Tuple{COBREXA.SBMLModel}","page":"Reference","title":"FBCModelTests.frog_model_report","text":"frog_model_report(\n    model::COBREXA.SBMLModel;\n    optimizer,\n    workers\n)\n\n\nGenerate FROGReportData for a model.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_objective_report-Tuple{COBREXA.SBMLModel, String}","page":"Reference","title":"FBCModelTests.frog_objective_report","text":"frog_objective_report(\n    sbml_model::COBREXA.SBMLModel,\n    objective::String;\n    optimizer,\n    workers\n)\n\n\nGenerate a FROGObjectiveReport containing the reproducibility data for a single objective in the SBML model.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_read_from_directory-Tuple{String}","page":"Reference","title":"FBCModelTests.frog_read_from_directory","text":"frog_read_from_directory(\n    report_dir::String\n) -> NamedTuple{(:metadata, :report), _A} where _A<:Tuple{Dict{String, String}, Dict}\n\n\nReverse of frog_write_to_directory.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_test_metadata_compatibility-Tuple{Dict{String, String}, Dict{String, String}}","page":"Reference","title":"FBCModelTests.frog_test_metadata_compatibility","text":"frog_test_metadata_compatibility(\n    a::Dict{String, String},\n    b::Dict{String, String}\n)\n\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_test_report_equality-Tuple{Dict{String, FBCModelTests.FROGObjectiveReport}, Dict{String, FBCModelTests.FROGObjectiveReport}}","page":"Reference","title":"FBCModelTests.frog_test_report_equality","text":"frog_test_report_equality(\n    a::Dict{String, FBCModelTests.FROGObjectiveReport},\n    b::Dict{String, FBCModelTests.FROGObjectiveReport};\n    absolute_tolerance,\n    relative_tolerance\n) -> Union{Test.FallbackTestSet, Test.DefaultTestSet}\n\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.frog_write_to_directory-Tuple{Dict{String, FBCModelTests.FROGObjectiveReport}, Dict{String, String}}","page":"Reference","title":"FBCModelTests.frog_write_to_directory","text":"frog_write_to_directory(\n    r::Dict{String, FBCModelTests.FROGObjectiveReport},\n    metadata::Dict{String, String};\n    report_dir,\n    basefilename\n)\n\n\nWrite the contents of FROGReportData to the 4 TSV files as specified by FROG standard, and additionally write the metadata into the JSON file.\n\n\n\n\n\n","category":"method"},{"location":"functions/#MEMOTE","page":"Reference","title":"MEMOTE","text":"","category":"section"},{"location":"functions/","page":"Reference","title":"Reference","text":"Modules = [FBCModelTests]\nPages = [\"memote.jl\"]","category":"page"},{"location":"functions/#FBCModelTests.has_erroneous_energy_generating_cycles-Tuple{Any, Any}","page":"Reference","title":"FBCModelTests.has_erroneous_energy_generating_cycles","text":"has_erroneous_energy_generating_cycles(\n    model,\n    optimizer;\n    energy_dissipating_metabolites,\n    additional_energy_generating_reactions,\n    ignore_reactions,\n    modifications,\n    boundary_reactions\n) -> Bool\n\n\nAttempts to detect if the model contains any erroneous energy generating cycles by closing all exchange reactions and using flux balance analysis to maximize the sum of fluxes through a set of energy dissipating reactions. The flux sum should be zero if the model is free of energy generating reactions.\n\nThis function is based on Fritzemeier, Claus Jonathan, et al. \"Erroneous energy-generating cycles in published genome scale metabolic networks: Identification and removal.\" PLoS computational biology (2017).\n\nThe energy dissipating reactions are based on the source paper, and include:\n\nATP + H2O --> ADP + H + Phosphate\nCTP + H2O --> CDP + H + Phosphate\nGTP + H2O --> GDP + H + Phosphate\nUTP + H2O --> UDP + H + Phosphate\nITP + H2O --> IDP + H + Phosphate\nNADH --> H + NAD\nNADPH --> H + NADP\nFADH2 --> 2 H + FAD\nFMNH2 --> 2 H + FMN\nUbiquinol-8 --> 2 H + Ubiquinone-8\nMenaquinol-8 --> 2 H + Menaquinone-8\n2-Demethylmenaquinol-8 --> 2 H + 2-Demethylmenaquinone-8\nH2O + ACCOA --> H + Acetate + COA\nL-Glutamate + H2O --> 2-Oxoglutarate + Ammonium + 2 H\nH[external] --> H\n\nAdditional energy dissipating reactions can be directly specified through additional_energy_generating_reactions, which should be vector of COBREXA Reactions using the same metabolite name space as the model. Internally, the model is converted to a COBREXA StandardModel, so ensure that the appropriate accessors are defined for it.\n\nSince models use different name spaces, energy_dissipating_metabolites is used to create the energy dissipating reactions. By default it uses the BiGG name space, but this will be changed to ChEBI in due course. If your model uses a different name space, then you have to change the values (NOT the keys) of energy_dissipating_metabolites. Each energy dissipating reaction is added to the test only if all its associated metabolites are present. Any modifications to the solver are passed directly through to COBREXA's flux_balance_analysis function. All boundary_reactions and ignore_reactions are deleted from an internal copy of model; this internal copy is used for analysis.\n\nReturns true if the model has energy generating cycles.\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.is_consistent-Tuple{Any, Any}","page":"Reference","title":"FBCModelTests.is_consistent","text":"is_consistent(\n    model,\n    optimizer;\n    boundary_reactions,\n    ignored_reactions\n) -> Bool\n\n\nDetermines if the model is stiochiometrically consistent. Note, stoichiometric consistency does not guarantee that mass balances must hold in the model. A more robust check is is_model_mass_balanced, but this works if not all metabolites have mass assigned to them.\n\nBased on Gevorgyan, Albert, Mark G. Poolman, and David A. Fell. \"Detection of stoichiometric inconsistencies in biomolecular models.\" Bioinformatics (2008).\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.is_model_charge_balanced-Tuple{COBREXA.MetabolicModel}","page":"Reference","title":"FBCModelTests.is_model_charge_balanced","text":"is_model_charge_balanced(\n    model::COBREXA.MetabolicModel;\n    ignored_reactions\n) -> Vector{String}\n\n\nIterates through all the reactions in model and checks if the charges across each reaction balance. Returns a list of reaction IDs that are charge unbalanced, which is empty if the test passes.\n\nOptionally, pass a list of reactions to ignore in this process through ignored_reactions. It makes sense to include the biomass and exchange reactions in this list (default).\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.is_model_mass_balanced-Tuple{COBREXA.MetabolicModel}","page":"Reference","title":"FBCModelTests.is_model_mass_balanced","text":"is_model_mass_balanced(\n    model::COBREXA.MetabolicModel;\n    ignored_reactions\n) -> Vector{String}\n\n\nIterates through all the reactions in model and checks if the mass across each reaction balances. Returns a list of reaction IDs that are mass unbalanced, which is empty if the test passes.\n\nOptionally, pass a list of reactions to ignore in this process through ignored_reactions. It makes sense to include the biomass and exchange reactions in this list (default).\n\n\n\n\n\n","category":"method"},{"location":"functions/#Types-and-utilities","page":"Reference","title":"Types and utilities","text":"","category":"section"},{"location":"functions/","page":"Reference","title":"Reference","text":"Modules = [FBCModelTests]\nPages = [\"structs.jl\", \"common.jl\"]","category":"page"},{"location":"functions/#FBCModelTests.ObjectiveValue","page":"Reference","title":"FBCModelTests.ObjectiveValue","text":"\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.FROGMetadata","page":"Reference","title":"FBCModelTests.FROGMetadata","text":"mutable struct Dict{String, String} <: AbstractDict{String, String}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.FROGObjectiveReport","page":"Reference","title":"FBCModelTests.FROGObjectiveReport","text":"struct FROGObjectiveReport\n\nFields\n\noptimum::Union{Nothing, Float64}\nreactions::Dict{String, FBCModelTests.FROGReactionReport}\ngene_deletions::Dict{String, Union{Nothing, Float64}}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.FROGReactionReport","page":"Reference","title":"FBCModelTests.FROGReactionReport","text":"struct FROGReactionReport\n\nFields\n\nflux::Union{Nothing, Float64}\nvariability_min::Union{Nothing, Float64}\nvariability_max::Union{Nothing, Float64}\ndeletion::Union{Nothing, Float64}\n\n\n\n\n\n","category":"type"},{"location":"functions/#FBCModelTests.FROGReportData","page":"Reference","title":"FBCModelTests.FROGReportData","text":"mutable struct Dict{String, FBCModelTests.FROGObjectiveReport} <: AbstractDict{String, FBCModelTests.FROGObjectiveReport}\n\n\n\n\n\n","category":"type"},{"location":"#FBCModelTests.jl-—-testing-and-reproducibility-of-constraint-based-metabolic-modeling","page":"Home","title":"FBCModelTests.jl — testing and reproducibility of constraint-based metabolic modeling","text":"","category":"section"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"functions.md\"]\nDepth = 2","category":"page"}]
}
