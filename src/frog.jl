
#TODO any other kinds of suboptimal solutions?
objstatus(::Nothing) = "infeasible"
objstatus(::Any) = "optimal"

objvalue(::Nothing) = ""
objvalue(x::Number) = string(x)

parse_objvalue(x::String) = x == "" ? nothing : parse(Float64, x)

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

"""
$(TYPEDSIGNATURES)

Generate a [`FROGObjectiveReport`](@ref) containing the reproducibility data
for a single objective in the SBML model.
"""
function frog_objective_report(
    sbml_model::SBMLModel,
    objective::String;
    optimizer,
    workers,
)::FROGObjectiveReport
    @info "Creating report for objective $objective ..."
    # this prevents the default SBMLModel fireworks in case there's multiple objectives
    model = ResetObjective(sbml_model, SBML.fbc_flux_objective(sbml_model.sbml, objective))

    # run the first FBA
    @info "Finding model objective value ..."
    solved_model = flux_balance_analysis(model, optimizer)
    obj = solved_objective_value(solved_model)

    fvas = if isnothing(obj)
        @warn "Model does not have a feasible solution, skipping FVA."
        zeros(0, 2)
    else
        @info "Optimal solution found." obj
        @info "Calculating model variability ..."
        flux_variability_analysis(
            model,
            optimizer;
            bounds = objective_bounds(1.0),
            optimal_objective_value = obj,
            workers = workers,
        )
    end

    @info "Calculating gene knockouts ..."
    gids = genes(model)
    gs = screen(
        model,
        args = tuple.(gids),
        analysis = (m, gene) -> solved_objective_value(
            flux_balance_analysis(m, optimizer, modifications = [knockout(gene)]),
        ),
        workers = workers,
    )

    @info "Calculating reaction knockouts ..."
    rids = reactions(model)
    rs = screen(
        model,
        args = tuple.(rids),
        analysis = (m, rid) -> solved_objective_value(
            flux_balance_analysis(
                m,
                optimizer,
                modifications = [change_constraint(rid, lb = 0.0, ub = 0.0)],
            ),
        ),
        workers = workers,
    )

    @info "Objective $objective done."
    return FROGObjectiveReport(
        optimum = obj,
        reactions = Dict(
            rid => FROGReactionReport(
                flux = flx,
                variability_min = vmin,
                variability_max = vmax,
                deletion = ko,
            ) for (rid, flx, vmin, vmax, ko) in
            zip(rids, flux_vector(model, solved_model), fvas[:, 1], fvas[:, 2], rs)
        ),
        gene_deletions = Dict(gids .=> gs),
    )
end

"""
$(TYPEDSIGNATURES)

Generate [`FROGReportData`](@ref) for a model.
"""
frog_model_report(model::SBMLModel; optimizer, workers = [Distributed.myid()]) = Dict([
    objective => frog_objective_report(model, objective; optimizer, workers) for
    (objective, _) in model.sbml.objectives
])

"""
$(TYPEDSIGNATURES)
"""
function frog_test_report_equality(
    a::FROGReportData,
    b::FROGReportData;
    absolute_tolerance = 1e-6,
    relative_tolerance = 1e-4,
)
    test_dicts(match::Function, a::Dict, b::Dict) =
        for k in union(keys(a), keys(b))
            @testset "$k" begin
                @test haskey(a, k) && haskey(b, k)
                if haskey(a, k) && haskey(b, k)
                    match(k, a[k], b[k])
                end
            end
        end

    intol(a, b) =
        (isnothing(a) && isnothing(b)) || (
            !isnothing(a) &&
            !isnothing(b) &&
            abs(a - b) <= absolute_tolerance &&
            (
                (a * b > 0) && abs(a * (1 + relative_tolerance)) >= abs(b) ||
                abs(b * (1 + relative_tolerance)) >= abs(a)
            )
        )

    @testset "Comparing objectives" begin
        test_dicts(
            (k, a, b) -> begin
                @test intol(a.optimum, b.optimum)
                @testset "Reactions" begin
                    test_dicts(
                        (k, a, b) -> begin
                            @test intol(a.flux, b.flux)
                            @test intol(a.variability_min, b.variability_min)
                            @test intol(a.variability_max, b.variability_max)
                            @test intol(a.deletion, b.deletion)
                        end,
                        a.reactions,
                        b.reactions,
                    )
                end
                @testset "Gene deletions" begin
                    test_dicts(
                        (k, a, b) -> begin
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
frog_metadata(filename::String; optimizer, basefilename::String = basename(filename)) =
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

"""
$(TYPEDSIGNATURES)
"""
function frog_test_metadata_compatibility(a::FROGMetadata, b::FROGMetadata)
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

const frog_headers = Dict(
    :objective => ["model" "objective" "status" "value"],
    :fva => ["model" "objective" "reaction" "flux" "status" "minimum" "maximum"],
    :gene_deletion => ["model" "objective" "gene" "status" "value"],
    :reaction_deletion => ["model" "objective" "reaction" "status" "value"],
)

"""
$(TYPEDSIGNATURES)

Write the contents of [`FROGReportData`](@ref) to the 4 TSV files as specified
by FROG standard, and additionally write the metadata into the JSON file.
"""
function frog_write_to_directory(
    r::FROGReportData,
    metadata::FROGMetadata;
    report_dir::String,
    basefilename::String,
)
    outname(x) = joinpath(report_dir, x)
    writeto(x::Function, fn) = open(x, outname(fn), "w")
    squash(x) = permutedims(hcat(x...), (2, 1))

    mkpath(report_dir)

    writeto("01_objective.tsv") do f
        writedlm(
            f,
            [
                frog_headers[:objective]
                squash(
                    [basefilename, obj, objstatus(o.optimum), objvalue(o.optimum)] for
                    (obj, o) in r
                )
            ],
            '\t',
        )
    end

    writeto("02_fva.tsv") do f
        writedlm(
            f,
            [
                frog_headers[:fva]
                squash(
                    [
                        basefilename,
                        obj,
                        rxn,
                        objvalue(r.flux),
                        objstatus(r.variability_min),
                        objvalue(r.variability_min),
                        objvalue(r.variability_max),
                    ] for (obj, o) = r if !isnothing(o.optimum) for (rxn, r) in o.reactions
                )
            ],
            '\t',
        )
    end

    writeto("03_gene_deletion.tsv") do f
        writedlm(
            f,
            [
                frog_headers[:gene_deletion]
                squash(
                    [basefilename, obj, gene, objstatus(geneval), objvalue(geneval)] for
                    (obj, o) in r for (gene, geneval) in o.gene_deletions
                )
            ],
            '\t',
        )
    end

    writeto("04_reaction_deletion.tsv") do f
        writedlm(
            f,
            [
                frog_headers[:reaction_deletion]
                squash(
                    [basefilename, obj, rxn, objstatus(r.deletion), objvalue(r.deletion)]
                    for (obj, o) in r for (rxn, r) in o.reactions
                )
            ],
            '\t',
        )
    end

    writeto("00_metadata.json") do f
        JSON.print(f, metadata, 2)
    end

    nothing
end

"""
$(TYPEDSIGNATURES)

Reverse of [`frog_write_to_directory`](@ref).
"""
function frog_read_from_directory(report_dir::String)
    outname(x) = joinpath(report_dir, x)
    readfrom(x::Function, fn) = open(x, outname(fn), "r")

    metadata = readfrom("00_metadata.json") do f
        FROGMetadata(JSON.parse(f))
    end

    basefilename = metadata["model.filename"]

    # objectives
    objdata = readfrom("01_objective.tsv") do f
        readdlm(f, '\t', String)
    end
    @assert objdata[1, :] == frog_headers[:objective][1, :]
    obj_vals = Dict(
        o => parse_objvalue(v) for
        (m, o, _, v) = eachrow(objdata[2:end, :]) if m == basefilename
    )

    # FVA
    fvadata = readfrom("02_fva.tsv") do f
        readdlm(f, '\t', String)
    end
    @assert fvadata[1, :] == frog_headers[:fva][1, :]
    fva_vals = Dict(
        obj => Dict(
            r => parse_objvalue.((flx, min, max)) for
            (m, o, r, flx, _, min, max) = eachrow(fvadata[2:end, :]) if
            m == basefilename && o == obj
        ) for obj in keys(obj_vals)
    )

    # gene deletions
    gdata = readfrom("03_gene_deletion.tsv") do f
        readdlm(f, '\t', String)
    end
    @assert gdata[1, :] == frog_headers[:gene_deletion][1, :]
    gene_vals = Dict(
        obj => Dict(
            g => parse_objvalue(v) for (m, o, g, _, v) = eachrow(gdata[2:end, :]) if
            m == basefilename && o == obj
        ) for obj in keys(obj_vals)
    )

    # reaction deletions
    rdata = readfrom("04_reaction_deletion.tsv") do f
        readdlm(f, '\t', String)
    end
    @assert rdata[1, :] == frog_headers[:reaction_deletion][1, :]
    rxn_vals = Dict(
        obj => Dict(
            r => parse_objvalue(v) for (m, o, r, _, v) = eachrow(rdata[2:end, :]) if
            m == basefilename && o == obj
        ) for obj in keys(obj_vals)
    )

    all_rxns = unique(
        vcat(
            [rxn for (_, o) in fva_vals for (rxn, _) in o],
            [rxn for (_, o) in rxn_vals for (rxn, _) in o],
        ),
    )

    return (
        metadata = metadata,
        report = Dict(
            obj => FROGObjectiveReport(
                optimum = opt,
                reactions = Dict(
                    rxn => FROGReactionReport(
                        flux = gets(fva_vals, nothing, obj, rxn, 1),
                        variability_min = gets(fva_vals, nothing, obj, rxn, 2),
                        variability_max = gets(fva_vals, nothing, obj, rxn, 3),
                        deletion = gets(rxn_vals, nothing, obj, rxn),
                    ) for rxn in all_rxns
                ),
                gene_deletions = gets(gene_vals, Dict(), obj),
            ) for (obj, opt) in obj_vals
        ),
    )
end

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
    r = frog_model_report(m; optimizer, workers)

    @info "Generating metadata for $filename ..."
    metadata = frog_metadata(filename; optimizer, basefilename)

    @info "Writing output to $report_dir ..."
    frog_write_to_directory(r, metadata; report_dir, basefilename)

    @info "FROG report done for $filename."
end

"""
$(TYPEDSIGNATURES)

A simple wrapper for comparing 2 previously generated reports in their
respective directories.
"""
function frog_compare_reports(report_dir_a::String, report_dir_b::String)
    a = frog_read_from_directory(report_dir_a)
    b = frog_read_from_directory(report_dir_b)

    @testset "Comparing FROG reports in $report_dir_a and $report_dir_b" begin
        @testset "Metadata" begin
            frog_test_metadata_compatibility(a.metadata, b.metadata)
        end
        @testset "Objectives and solution values" begin
            frog_test_report_equality(a.report, b.report)
        end
    end
end
