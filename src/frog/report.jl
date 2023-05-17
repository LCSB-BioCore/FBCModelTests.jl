
"""
    module ReportGenerators

Functions for generating FROG report data and metadata file contents.
"""
module ReportGenerators

using ..FROG: FROGReactionReport, FROGObjectiveReport, FROGMetadata, FROGReportData

using ...FBCModelTests: FBCMT_VERSION

using COBREXA
using JuMP
using Dates
using Distributed
using DocStringExtensions
using MD5
using SBML
using SHA
using SparseArrays
using Test

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
        sparse([ds[rid] for rid in reactions(m)])
    end

"""
$(TYPEDSIGNATURES)

Generate a [`FROGObjectiveReport`](@ref) containing the reproducibility data
for a single objective in the SBML model.
"""
function frog_objective_report(
    sbml_model::SBMLModel,
    objective::Maybe{String} = nothing;
    optimizer,
    modifications = [],
    workers = [Distributed.myid()],
    fraction_optimum = 1.0,
)::FROGObjectiveReport
    @info "Creating report for objective $objective ..."
    # this prevents the default SBMLModel fireworks in case there's multiple objectives
    model =
        isnothing(objective) ? sbml_model :
        ResetObjective(sbml_model, get_objective(sbml_model, objective))

    # run the first FBA
    @info "Finding model objective value ..."
    solved_model = flux_balance_analysis(model, optimizer; modifications)
    fbas = flux_vector(model, solved_model)
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
                objective_flux = fba,
                fraction_optimum = fraction_optimum,
                variability_min = fvarow[1],
                variability_max = fvarow[2],
                deletion = ko,
            ) for (rid, fba, fvarow, ko) in zip(rids, fbas, eachrow(fvas), rs)
        ),
        gene_deletions = Dict(gids .=> gs),
    )
end

"""
$(TYPEDSIGNATURES)

Generate [`FROGReportData`](@ref) for a model.
"""
generate_report_data(model::SBMLModel; kwargs...) =
    isempty(model.sbml.objectives) ?
    Dict("obj" => frog_objective_report(model; kwargs...)) :
    Dict([
        objective => frog_objective_report(model, objective; kwargs...) for
        (objective, _) in model.sbml.objectives
    ])

"""
$(TYPEDSIGNATURES)
"""
generate_metadata(filename::String; optimizer, basefilename::String = basename(filename)) =
    FROGMetadata(
        "frog_version" => "0.1.4",
        "frog_date" => string(Dates.today()),
        "model_filename" => basefilename,
        "model_md5" => bytes2hex(open(md5, filename, "r")),
        "model_sha256" => bytes2hex(open(sha256, filename, "r")),
        "environment" => begin
            x = IOBuffer()
            InteractiveUtils.versioninfo(x)
            replace(String(take!(x)), r"\n *" => " ")
        end,
        "software" => Dict(
            "frog" => Dict(
                "name" => "FBCModelTests.jl",
                "version" => string(FBCMT_VERSION),
                "url" => "https://github.com/LCSB-BioCore/FBCModelTests.jl",
            ),
            "toolbox" => Dict(
                "name" => "COBREXA.jl",
                "version" => string(COBREXA_VERSION),
                "url" => "https://github.com/LCSB-BioCore/COBREXA.jl",
            ),
            "solver" => Dict(
                "name" => COBREXA.JuMP.MOI.get(optimizer(), COBREXA.JuMP.MOI.SolverName()),
                "version" => COBREXA.JuMP.MOI.get(
                    optimizer(),
                    COBREXA.JuMP.MOI.SolverVersion(),
                ),
            ),
        ),
    )

end

"""
    module ReportTests

Function for testing the compatibility of FROG report data.
"""
module ReportTests

using ..FROG: FROGReportData, FROGMetadata
using ...FBCModelTests: test_dicts, @atest

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

    invar(a, min, max) =
        all(isnothing.((a, min, max))) || (
            !any(isnothing.((a, min, max))) &&
            (intol(a, min) || intol(a, max) || (min <= a && a <= max))
        )

    test_dicts(
        (objid, a, b) -> begin
            @testset "Objective $objid" begin
                @atest intol(a.optimum, b.optimum) "objective values are the same"
                @testset "Reactions" begin
                    test_dicts(
                        (rid, a, b) -> begin
                            @atest intol(a.variability_min, b.variability_min) "lower variability bound of $rid matches"
                            @atest intol(a.variability_max, b.variability_max) "upper variability bound of $rid matches"
                            @atest invar(
                                a.objective_flux,
                                b.variability_min,
                                b.variability_max,
                            ) "flux through $rid ($(a.objective_flux)) is in variability range $((b.variability_min, b.variability_max))"
                            @atest invar(
                                b.objective_flux,
                                a.variability_min,
                                a.variability_max,
                            ) "flux through $rid ($(b.objective_flux)) is in variability range $((a.variability_min,a.variability_max))"
                            @atest intol(a.fraction_optimum, b.fraction_optimum) "optimum fraction at $rid is the same"
                            @atest intol(a.deletion, b.deletion) "deletion of reaction $rid gives the same objective value"
                        end,
                        a.reactions,
                        b.reactions,
                    )
                end
                @testset "Gene deletions" begin
                    test_dicts(
                        (gid, a, b) -> begin
                            @atest intol(a, b) "deletion of gene $gid gives the same objective value"
                        end,
                        a.gene_deletions,
                        b.gene_deletions,
                    )
                end
            end
        end,
        a,
        b,
    )
end

"""
$(TYPEDSIGNATURES)
"""
function test_metadata_compatibility(a::FROGMetadata, b::FROGMetadata)
    for k in ["model_filename", "model_md5"]
        @atest haskey(a, k) "$k is present in metadata"
        @atest haskey(b, k) "$k is present in metadata"
    end

    for k in ["model_filename", "model_md5", "model_sha256"]
        if haskey(a, k) && haskey(b, k)
            @atest a[k] == b[k] "$k matches $((a[k], b[k]))"
        end
    end
end

end
