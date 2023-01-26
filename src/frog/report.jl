
"""
    module ReportGenerators

Functions for generating FROG report data and metadata file contents.
"""
module ReportGenerators

using ..FROG: FROGReactionReport, FROGObjectiveReport, FROGMetadata, FROGReportData

using ...FBCModelTests: FBCMT_VERSION

using COBREXA
using JuMP
using Distributed
using DocStringExtensions
using MD5
using SBML
using SHA
using SparseArrays

import InteractiveUtils

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
struct ResetObjective <: ModelWrapper
    model::MetabolicModel
    objective::SparseVec
end

COBREXA.unwrap_model(x::ResetObjective) = x.model
COBREXA.objective(x::ResetObjective) = x.objective

get_objective(m::SBMLModel, objective::String) =
    let ds = Dict(keys(m.sbml.reactions) .=> SBML.fbc_flux_objective(m.sbml, objective))
    sparse([ds[rid] for rid=reactions(m)])
end

"""
$(TYPEDSIGNATURES)

Generate a [`FROGObjectiveReport`](@ref) containing the reproducibility data
for a single objective in the SBML model.
"""
function frog_objective_report(
    sbml_model::SBMLModel,
    objective::String;
    optimizer,
    modifications = [],
    workers = [Distributed.myid()],
    fraction_optimum = 1.0,
)::FROGObjectiveReport
    @info "Creating report for objective $objective ..."
    # this prevents the default SBMLModel fireworks in case there's multiple objectives
    model = ResetObjective(sbml_model, get_objective(sbml_model, objective))

    # run the first FBA
    @info "Finding model objective value ..."
    solved_model = flux_balance_analysis(model, optimizer; modifications)
    optobj = solved_objective_value(solved_model)

    fvas = if isnothing(optobj)
        @warn "Model does not have a feasible solution, skipping FVA."
        zeros(0, 2)
    else
        @info "Optimal solution found." optobj
        @info "Calculating model variability ..."
        flux_variability_analysis(
            model,
            optimizer;
            bounds = objective_bounds(fraction_optimum),
            optimal_objective_value = optobj,
            workers = workers,
            modifications,
            ret = m -> (JuMP.objective_value(m), JuMP.value.(m[:x])' * model.objective),
        )
    end

    @info "Calculating gene knockouts ..."
    gids = genes(model)
    gs = if isempty(gids)
        @info "Model has no genes"
        zeros(0)
    else
        screen(
            model,
            args = tuple.(gids),
            analysis = (m, gene) -> solved_objective_value(
                flux_balance_analysis(
                    m,
                    optimizer,
                    modifications = vcat(modifications, knockout(gene)),
                ),
            ),
            workers = workers,
        )
    end

    @info "Calculating reaction knockouts ..."
    rids = reactions(model)
    rs = screen(
        model,
        args = tuple.(rids),
        analysis = (m, rid) -> solved_objective_value(
            flux_balance_analysis(
                m,
                optimizer,
                modifications = vcat(
                    modifications,
                    change_constraint(rid, lb = 0.0, ub = 0.0),
                ),
            ),
        ),
        workers = workers,
    )

    @info "Objective $objective done."
    return FROGObjectiveReport(
        optimum = optobj,
        reactions = Dict(
            rid => FROGReactionReport(
                objective_flux = isnothing(fvarow[1]) ? nothing : last(fvarow[1]),
                fraction_optimum = fraction_optimum,
                variability_min = isnothing(fvarow[1]) ? nothing : first(fvarow[1]),
                variability_max = isnothing(fvarow[2]) ? nothing : first(fvarow[2]),
                deletion = ko,
            ) for (rid, fvarow, ko) in
            zip(rids, eachrow(fvas), rs)
        ),
        gene_deletions = Dict(gids .=> gs),
    )
end

"""
$(TYPEDSIGNATURES)

Generate [`FROGReportData`](@ref) for a model.
"""
generate_report_data(model::SBMLModel; kwargs...) = Dict([
    objective => frog_objective_report(model, objective; kwargs...) for
    (objective, _) in model.sbml.objectives
])

"""
$(TYPEDSIGNATURES)
"""
generate_metadata(filename::String; optimizer, basefilename::String = basename(filename)) =
    FROGMetadata(
        "software.name" => "FBCModelTests.jl",
        "software.version" => string(FBCMT_VERSION),
        "software.url" => "https://github.com/LCSB-BioCore/FBCModelTests.jl/",
        "environment" => begin
            x = IOBuffer()
            InteractiveUtils.versioninfo(x)
            replace(String(take!(x)), r"\n *" => " ")
        end,
        "model.filename" => basefilename,
        "model.md5" => bytes2hex(open(f -> md5(f), filename, "r")),
        "model.sha256" => bytes2hex(open(f -> sha256(f), filename, "r")),
        "solver.name" => "COBREXA.jl $COBREXA_VERSION ($(COBREXA.JuMP.MOI.get(optimizer(), COBREXA.JuMP.MOI.SolverName())))",
    )

end

"""
    module ReportTests

Function for testing the compatibility of FROG report data.
"""
module ReportTests

using ..FROG: FROGReportData, FROGMetadata
using ...FBCModelTests: test_dicts

using Test, DocStringExtensions

"""
$(TYPEDSIGNATURES)
"""
function test_report_compatibility(
    a::FROGReportData,
    b::FROGReportData;
    absolute_tolerance = 1e-6,
    relative_tolerance = 1e-4,
)
    intol(a, b) =
        (isnothing(a) && isnothing(b)) || (
            !isnothing(a) &&
            !isnothing(b) &&
            abs(a - b) < absolute_tolerance + relative_tolerance * max(abs(a), abs(b))
        )

    @testset "Comparing objectives" begin
        test_dicts(
            (_, a, b) -> begin
                @test intol(a.optimum, b.optimum)
                @testset "Reactions" begin
                    test_dicts(
                        (_, a, b) -> begin
                            @test intol(a.objective_flux, b.objective_flux)
                            @test intol(a.variability_min, b.variability_min)
                            @test intol(a.variability_max, b.variability_max)
                            @test intol(a.fraction_optimum, b.fraction_optimum)
                            @test intol(a.deletion, b.deletion)
                        end,
                        a.reactions,
                        b.reactions,
                    )
                end
                @testset "Gene deletions" begin
                    test_dicts(
                        (_, a, b) -> begin
                            @test intol(a, b)
                        end,
                        a.gene_deletions,
                        b.gene_deletions,
                    )
                end
            end,
            a,
            b,
        )
    end
end

"""
$(TYPEDSIGNATURES)
"""
function test_metadata_compatibility(a::FROGMetadata, b::FROGMetadata)
    for k in ["model.filename", "model.md5"]
        @testset "$k is present" begin
            @test haskey(a, k)
            @test haskey(b, k)
        end
    end

    for k in ["model.filename", "model.md5", "model.sha256"]
        if haskey(a, k) && haskey(b, k)
            @testset "$k matches" begin
                @test a[k] == b[k]
            end
        end
    end
end

end
