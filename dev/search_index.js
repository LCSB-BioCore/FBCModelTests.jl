var documenterSearchIndex = {"docs":
[{"location":"functions/#FROG","page":"Reference","title":"FROG","text":"","category":"section"},{"location":"functions/","page":"Reference","title":"Reference","text":"Modules = [FBCModelTests]\nPages = [\"frog.jl\"]","category":"page"},{"location":"functions/#FBCModelTests.generate_frog_report-Tuple{COBREXA.SBMLModel, String}","page":"Reference","title":"FBCModelTests.generate_frog_report","text":"generate_frog_report(\n    model::COBREXA.SBMLModel,\n    filename::String;\n    report_dir,\n    optimizer,\n    workers,\n    basefilename\n)\n\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.generate_frog_report-Tuple{String}","page":"Reference","title":"FBCModelTests.generate_frog_report","text":"generate_frog_report(filename::String; kwargs...)\n\n\n\n\n\n\n","category":"method"},{"location":"functions/#MEMOTE","page":"Reference","title":"MEMOTE","text":"","category":"section"},{"location":"functions/","page":"Reference","title":"Reference","text":"Modules = [FBCModelTests]\nPages = [\"memote.jl\"]","category":"page"},{"location":"functions/#FBCModelTests.is_model_charge_balanced-Tuple{COBREXA.MetabolicModel}","page":"Reference","title":"FBCModelTests.is_model_charge_balanced","text":"is_model_charge_balanced(\n    model::COBREXA.MetabolicModel;\n    ignored_reactions\n) -> Vector{String}\n\n\nIterates through all the reactions in model and checks if the charges across each reaction balance. Returns a list of reaction IDs that are charge unbalanced, which is empty if the test passes.\n\nOptionally, pass a list of reactions to ignore in this process through ignored_reactions. It makes sense to include the biomass and exchange reactions in this list (default).\n\n\n\n\n\n","category":"method"},{"location":"functions/#FBCModelTests.is_model_mass_balanced-Tuple{COBREXA.MetabolicModel}","page":"Reference","title":"FBCModelTests.is_model_mass_balanced","text":"is_model_mass_balanced(\n    model::COBREXA.MetabolicModel;\n    ignored_reactions\n) -> Vector{String}\n\n\nIterates through all the reactions in model and checks if the mass across each reaction balances. Returns a list of reaction IDs that are mass unbalanced, which is empty if the test passes.\n\nOptionally, pass a list of reactions to ignore in this process through ignored_reactions. It makes sense to include the biomass and exchange reactions in this list (default).\n\n\n\n\n\n","category":"method"},{"location":"#FBCModelTests.jl-—-testing-and-reproducibility-of-constraint-based-metabolic-modeling","page":"Home","title":"FBCModelTests.jl — testing and reproducibility of constraint-based metabolic modeling","text":"","category":"section"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"functions.md\"]\nDepth = 2","category":"page"}]
}
