
"""
$(TYPEDSIGNATURES)

A complete function for one-shot generation of FROG reports. Use
[`ReportGenerators.generate_report_data`](@ref),
[`ReportGenerators.generate_metadata`](@ref) and
[`ReportIO.save_report`](@ref) for finer control of the process.
"""
function generate_report(
    filename::String;
    report_dir::String,
    optimizer,
    modifications = [],
    workers = [Distributed.myid()],
    basefilename::String = basename(filename),
)
    m = SBMLModel(readSBML(filename, throw_severities = ["Fatal"]))
    @info "Generating FROG report for $filename ..."
    r = ReportGenerators.generate_report_data(m; optimizer, modifications, workers)

    @info "Generating metadata for $filename ..."
    metadata = ReportGenerators.generate_metadata(filename; optimizer, basefilename)

    @info "Writing output to $report_dir ..."
    ReportIO.save_report(r, metadata; report_dir, basefilename)

    @info "FROG report done for $filename."
end

"""
$(TYPEDSIGNATURES)

A simple wrapper for comparing 2 previously generated reports in their
respective directories. Additional arguments are fowarded to
[`ReportTests.test_report_compatibility`](@ref).
"""
function compare_reports(report_dir_a::String, report_dir_b::String; kwargs...)
    a = ReportIO.load_report(report_dir_a)
    b = ReportIO.load_report(report_dir_b)

    @testset "Metadata" begin
        ReportTests.test_metadata_compatibility(a.metadata, b.metadata)
    end
    @testset "Objectives and solution values" begin
        ReportTests.test_report_compatibility(a.report, b.report; kwargs...)
    end
    nothing
end

"""
$(TYPEDSIGNATURES)

Like [`compare_reports`](@ref), but with a specific testset settings that
provide better user-facing reports.
"""
function compare_reports_toplevel(report_dir_a::String, report_dir_b::String; kwargs...)
    @testset QuietTestSet "Comparison of FROG reports $report_dir_a and $report_dir_b" begin
        compare_reports(report_dir_a, report_dir_b; kwargs...)
    end
end
