
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

    @testset "Report generates expected files" begin
        @test isdir(report_path)
        @test isfile(joinpath(report_path, "00_metadata.json"))
        @test isfile(joinpath(report_path, "01_objective.tsv"))
        @test isfile(joinpath(report_path, "02_fva.tsv"))
        @test isfile(joinpath(report_path, "03_gene_deletion.tsv"))
        @test isfile(joinpath(report_path, "04_reaction_deletion.tsv"))
    end

    @testset "Report equality is reflexive" begin
        result = @testset CountTests "Same reports" begin
            FROG.frog_compare_reports(report_path, report_path)
        end

        @test result.passes == 758
        @test result.fails == 0
        @test result.errs == 0
    end

    @testset "Generated report looks as expected" begin
        result = @testset CountTests "Expected report" begin
            FROG.frog_compare_reports(report_path, joinpath(datadir, "report-ecoli-good"))
        end

        @test result.passes == 758
        @test result.fails == 0
        @test result.errs == 0
    end

    @testset "Broken report does not pass comparison" begin
        result = @testset CountTests "Different reports" begin
            FROG.frog_compare_reports(report_path, joinpath(datadir, "report-ecoli-damaged"))
        end

        @test result.passes == 745
        @test result.fails == 13
        @test result.errs == 0
    end
end
