"""
$(TYPEDSIGNATURES)

"""
function generate_frog_report(filename::String; kwargs...)
    generate_frog_report(load_sbml_model(filename), filename; kwargs...)
end

"""
$(TYPEDSIGNATURES)

"""
function generate_frog_report(
    model::SBMLModel,
    filename::String;
    report_dir::String,
    optimizer,
    workers = [Distributed.myid()],
    basename::String = filename,
)
    @info "Generating FROG report for $filename..."
    mkpath(report_dir)

    outname(x) = joinpath(report_dir, x)

    @info "Finding model objective value..."
    n = n_reactions(model)
    solved_model = flux_balance_analysis(model, optimizer)
    obj = solved_objective_value(solved_model)
    open(outname("01_objective.tsv"), "w") do f
        writedlm(
            f,
            [
                "model" "objective" "status" "value"
                filename "obj" objstatus(obj) objvalue(obj)
            ],
        )
    end

    isnothing(obj) && throw("Model does not have an objective value, can't continue.")

    @info "Calculating model variability..."
    fvas = flux_variability_analysis(
        model,
        optimizer;
        bounds = objective_bounds(1.0),
        optimal_objective_value = obj,
        workers = workers,
    )
    open(outname("02_fva.tsv"), "w") do f
        writedlm(
            f,
            [
                "model" "objective" "reaction" "flux" "status" "minimum" "maximum"
                hcat(
                    fill(filename, n),
                    fill("obj", n),
                    reactions(model),
                    objvalue.(flux_vector(model, solved_model)),
                    objstatus.(fvas[:, 1]),
                    objvalue.(fvas[:, 1]),
                    objvalue.(fvas[:, 2]),
                )
            ],
        )
    end

    @info "Calculating gene knockouts..."
    gs = screen(
        model,
        args = tuple.(genes(model)),
        analysis = (m, gene) -> solved_objective_value(
            flux_balance_analysis(m, optimizer, modifications = [knockout(gene)]),
        ),
        workers = workers,
    )
    ng = n_genes(model)

    open(outname("03_gene_deletion.tsv"), "w") do f
        writedlm(
            f,
            [
                "model" "objective" "gene" "status" "value"
                hcat(
                    fill(filename, ng),
                    fill("obj", ng),
                    genes(model),
                    objstatus.(gs),
                    objvalue.(gs),
                )
            ],
        )
    end

    @info "Calculating reaction knockouts..."
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

    open(outname("04_reaction_deletion.tsv"), "w") do f
        writedlm(
            f,
            [
                "model" "objective" "reaction" "status" "value"
                hcat(
                    fill(filename, n),
                    fill("obj", n),
                    reactions(model),
                    objstatus.(rs),
                    objvalue.(rs),
                )
            ],
        )
    end

    @info "Writing the JSON metadata..."
    open(outname("00_metadata.json"), "w") do f
        JSON.print(
            f,
            Dict(
                "software.name" => "FBCModelTests.jl",
                "software.version" => string(FBCMT_VERSION),
                "software.url" => "https://github.com/LCSB-BioCore/FBCModelTests.jl/",
                "environment" => begin
                    x = IOBuffer()
                    InteractiveUtils.versioninfo(x)
                    replace(String(take!(x)), r"\n *" => " ")
                end,
                "model.filename" => filename,
                "model.md5" => bytes2hex(open(f -> md5(f), filename, "r")),
                "solver.name" => "COBREXA.jl $COBREXA_VERSION ($(COBREXA.JuMP.MOI.get(optimizer(), COBREXA.JuMP.MOI.SolverName())))",
            ),
            2,
        )
    end
    @info "Generated FROG report for $filename into $report_dir"
end

#TODO properly mark various other kinds of suboptimal solutions
objstatus(::Nothing) = "infeasible"
objstatus(::Any) = "optimal"

objvalue(::Nothing) = ""
objvalue(x::Number) = string(x)
