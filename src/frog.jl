
#TODO any other kinds of suboptimal solutions?
objstatus(::Nothing) = "infeasible"
objstatus(::Any) = "optimal"

objvalue(::Nothing) = ""
objvalue(x::Number) = string(x)

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
    gs = screen(
        model,
        args = tuple.(genes(model)),
        analysis = (m, gene) -> solved_objective_value(
            flux_balance_analysis(m, optimizer, modifications = [knockout(gene)]),
        ),
        workers = workers,
    )

    @info "Calculating reaction knockouts ..."
    rs = screen(
        model,
        args = tuple.(reactions(model)),
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
        flux = flux_vector(model, solved_model),
        variabilities_min = fvas[:, 1],
        variabilities_max = fvas[:, 2],
        reaction_knockouts = rs,
        gene_knockouts = gs,
    )
end

"""
$(TYPEDSIGNATURES)

Generate [`FROGReportData`](@ref) for a model.
"""
frog_model_report(model::SBMLModel; optimizer, workers = [Distributed.myid()]) =
    FROGReportData(
        reactions = reactions(model),
        genes = genes(model),
        objectives = Dict(
            objective => frog_objective_report(model, objective; optimizer, workers) for
            (objective, _) in model.sbml.objectives
        ),
    )

"""
$(TYPEDSIGNATURES)

Write the contents of [`FROGReportData`](@ref) to the 4 TSV files as specified
by FROG standard, and additionally write the metadata into the JSON file.
"""
function frog_write_to_directory(
    r::FROGReportData,
    metadata::Dict;
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
                "model" "objective" "status" "value"
                squash(
                    [basefilename, obj, objstatus(o.optimum), objvalue(o.optimum)] for
                    (obj, o) in r.objectives
                )
            ],
        )
    end

    writeto("02_fva.tsv") do f
        writedlm(
            f,
            [
                "model" "objective" "reaction" "flux" "status" "minimum" "maximum"
                squash(
                    [
                        basefilename,
                        obj,
                        rxn,
                        objvalue(flx),
                        objstatus(varmin),
                        objvalue(varmin),
                        objvalue(varmax),
                    ] for (obj, o) = r.objectives if !isnothing(o.optimum) for
                    (rxn, flx, varmin, varmax) in
                    zip(r.reactions, o.flux, o.variabilities_min, o.variabilities_max)
                )
            ],
        )
    end

    writeto("03_gene_deletion.tsv") do f
        writedlm(
            f,
            [
                "model" "objective" "gene" "status" "value"
                squash(
                    [basefilename, obj, gene, objstatus(geneval), objvalue(geneval)] for
                    (obj, o) in r.objectives for
                    (gene, geneval) in zip(r.genes, o.gene_knockouts)
                )
            ],
        )
    end

    writeto("04_reaction_deletion.tsv") do f
        writedlm(
            f,
            [
                "model" "objective" "reaction" "status" "value"
                squash(
                    [basefilename, obj, rxn, objstatus(rxnval), objvalue(rxnval)] for
                    (obj, o) in r.objectives for
                    (rxn, rxnval) in zip(r.reactions, o.reaction_knockouts)
                )
            ],
        )
    end

    writeto("00_metadata.json") do f
        JSON.print(f, metadata, 2)
    end

    nothing
end

"""
$(TYPEDSIGNATURES)
"""
frog_metadata(filename::String; optimizer, basefilename::String = basename(filename)) =
    Dict(
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
