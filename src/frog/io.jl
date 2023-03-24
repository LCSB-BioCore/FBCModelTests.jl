
"""
    module ReportIO

Functions for reading and writing FROG reports.
"""
module ReportIO

using ..FROG: FROGReactionReport, FROGObjectiveReport, FROGMetadata, FROGReportData

using ...FBCModelTests: gets

using DelimitedFiles
using DocStringExtensions
using JSON

#TODO any other kinds of suboptimal solutions?
objstatus(::Nothing) = "infeasible"
objstatus(::Any) = "optimal"

objvalue(::Nothing) = ""
objvalue(x::Number) = string(x)

parse_objvalue(x::String) = x == "" ? nothing : parse(Float64, x)

const frog_report_tsv_headers = Dict(
    :objective => ["model" "objective" "status" "value"],
    :fva =>
        ["model" "objective" "reaction" "flux" "status" "minimum" "maximum" "fraction_optimum"],
    :gene_deletion => ["model" "objective" "gene" "status" "value"],
    :reaction_deletion => ["model" "objective" "reaction" "status" "value"],
)

"""
$(TYPEDSIGNATURES)

Write the contents of [`FROGReportData`](@ref) to the 4 TSV files as specified
by FROG standard, and additionally write the metadata into the JSON file.
"""
function save_report(
    r::FROGReportData,
    metadata::FROGMetadata;
    report_dir::String,
    basefilename::String,
)
    outname(x) = joinpath(report_dir, x)
    writeto(x::Function, fn) = open(x, outname(fn), "w")
    squash(cols, x) =
        isempty(x) ? [cols; fill("", (0, length(cols)))] :
        [cols; permutedims(hcat(x...), (2, 1))]

    mkpath(report_dir)

    writeto("01_objective.tsv") do f
        writedlm(
            f,
            squash(
                frog_report_tsv_headers[:objective],
                [basefilename, obj, objstatus(o.optimum), objvalue(o.optimum)] for
                (obj, o) in r
            ),
            '\t',
        )
    end

    writeto("02_fva.tsv") do f
        writedlm(
            f,
            squash(
                frog_report_tsv_headers[:fva],
                [
                    basefilename,
                    obj,
                    rxn,
                    objvalue(r.objective_flux),
                    objstatus(r.variability_min),
                    objvalue(r.variability_min),
                    objvalue(r.variability_max),
                    string(r.fraction_optimum),
                ] for (obj, o) in r if !isnothing(o.optimum) for (rxn, r) in o.reactions
            ),
            '\t',
        )
    end

    writeto("03_gene_deletion.tsv") do f
        writedlm(
            f,
            squash(
                frog_report_tsv_headers[:gene_deletion],
                [basefilename, obj, gene, objstatus(geneval), objvalue(geneval)] for
                (obj, o) in r for (gene, geneval) in o.gene_deletions
            ),
            '\t',
        )
    end

    writeto("04_reaction_deletion.tsv") do f
        writedlm(
            f,
            squash(
                frog_report_tsv_headers[:reaction_deletion],
                [basefilename, obj, rxn, objstatus(r.deletion), objvalue(r.deletion)] for
                (obj, o) in r for (rxn, r) in o.reactions
            ),
            '\t',
        )
    end

    writeto("metadata.json") do f
        JSON.print(f, metadata, 2)
    end

    nothing
end

"""
$(TYPEDSIGNATURES)

Reverse of [`save_report`](@ref).
"""
function load_report(
    report_dir::String,
)::NamedTuple{(:metadata, :report),Tuple{FROGMetadata,FROGReportData}}
    outname(x) = joinpath(report_dir, x)
    readfrom(x::Function, fn) = open(x, outname(fn), "r")

    metadata = readfrom("metadata.json") do f
        FROGMetadata(JSON.parse(f))
    end

    basefilename = get(metadata, "model_filename", get(metadata, "model.filename", nothing))
    isnothing(basefilename) && error("report does not list its filename")

    # objectives
    objdata = readfrom("01_objective.tsv") do f
        readdlm(f, '\t', String)
    end
    @assert objdata[1, :] == frog_report_tsv_headers[:objective][1, :]
    obj_vals = Dict(
        o => parse_objvalue(v) for
        (m, o, _, v) in eachrow(objdata[2:end, :]) if m == basefilename
    )

    # FVA
    fvadata = readfrom("02_fva.tsv") do f
        readdlm(f, '\t', String)
    end
    @assert fvadata[1, :] == frog_report_tsv_headers[:fva][1, :]
    fva_vals = Dict(
        obj => Dict(
            r => parse_objvalue.((flx, min, max, frac)) for
            (m, o, r, flx, _, min, max, frac) in eachrow(fvadata[2:end, :]) if
            m == basefilename && o == obj
        ) for obj in keys(obj_vals)
    )

    # gene deletions
    gdata = readfrom("03_gene_deletion.tsv") do f
        readdlm(f, '\t', String)
    end
    @assert gdata[1, :] == frog_report_tsv_headers[:gene_deletion][1, :]
    gene_vals = Dict(
        obj => Dict(
            g => parse_objvalue(v) for
            (m, o, g, _, v) in eachrow(gdata[2:end, :]) if m == basefilename && o == obj
        ) for obj in keys(obj_vals)
    )

    # reaction deletions
    rdata = readfrom("04_reaction_deletion.tsv") do f
        readdlm(f, '\t', String)
    end
    @assert rdata[1, :] == frog_report_tsv_headers[:reaction_deletion][1, :]
    rxn_vals = Dict(
        obj => Dict(
            r => parse_objvalue(v) for
            (m, o, r, _, v) in eachrow(rdata[2:end, :]) if m == basefilename && o == obj
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
                        objective_flux = gets(fva_vals, nothing, obj, rxn, 1),
                        variability_min = gets(fva_vals, nothing, obj, rxn, 2),
                        variability_max = gets(fva_vals, nothing, obj, rxn, 3),
                        fraction_optimum = gets(fva_vals, NaN, obj, rxn, 4),
                        deletion = gets(rxn_vals, nothing, obj, rxn),
                    ) for rxn in all_rxns
                ),
                gene_deletions = gets(gene_vals, Dict(), obj),
            ) for (obj, opt) in obj_vals
        ),
    )
end

end
