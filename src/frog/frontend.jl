
"""
$(TYPEDSIGNATURES)

A complete function for one-shot generation of FROG reports. Use
[`frog_model_report`](@ref), [`frog_metadata`](@ref) and
[`frog_write_to_directory`](@ref) for finer control of the process.
"""
function frog_generate_report(
    filename::String;
    report_dir::String,
    optimizer,
    workers = [Distributed.myid()],
    basefilename::String = basename(filename),
)
    m = load_sbml_model(filename)
    @info "Generating FROG report for $filename ..."
    r = ReportGenerators.frog_model_report(m; optimizer, workers)

    @info "Generating metadata for $filename ..."
    metadata = ReportGenerators.frog_metadata(filename; optimizer, basefilename)

    @info "Writing output to $report_dir ..."
    ReportIO.frog_write_to_directory(r, metadata; report_dir, basefilename)

    @info "FROG report done for $filename."
end

"""
$(TYPEDSIGNATURES)

A simple wrapper for comparing 2 previously generated reports in their
respective directories.
"""
function frog_compare_reports(report_dir_a::String, report_dir_b::String)
    a = ReportIO.frog_read_from_directory(report_dir_a)
    b = ReportIO.frog_read_from_directory(report_dir_b)

    @testset "Comparing FROG reports in $report_dir_a and $report_dir_b" begin
        @testset "Metadata" begin
            ReportTests.frog_test_metadata_compatibility(a.metadata, b.metadata)
        end
        @testset "Objectives and solution values" begin
            ReportTests.frog_test_report_equality(a.report, b.report)
        end
    end
end
