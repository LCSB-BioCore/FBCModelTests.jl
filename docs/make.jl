using Documenter, FBCModelTests

makedocs(
    modules = [FBCModelTests],
    clean = false,
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://lcsb-biocore.github.io/FBCModelTests.jl/stable/",
    ),
    sitename = "FBCModelTests.jl",
    authors = "The developers of FBCModelTests.jl",
    linkcheck = !("skiplinks" in ARGS),
    pages = ["Home" => "index.md", "Reference" => "functions.md"],
    strict = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/LCSB-BioCore/FBCModelTests.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
