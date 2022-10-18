
@testset "FROG" begin

    import FBCModelTests: FROG
    using GLPK

    report_path = joinpath(datadir, "report-generated")
    rm(report_path, recursive = true, force = true)

    FROG.frog_generate_report(
        model_file["e_coli_core.xml"],
        report_dir = report_path,
        optimizer = GLPK.Optimizer,
    )

    @test isdir(report_path)
    @test isfile(joinpath(report_path, "00_metadata.json"))
    @test isfile(joinpath(report_path, "01_objective.tsv"))
    @test isfile(joinpath(report_path, "02_fva.tsv"))
    @test isfile(joinpath(report_path, "03_gene_deletion.tsv"))
    @test isfile(joinpath(report_path, "04_reaction_deletion.tsv"))
end
