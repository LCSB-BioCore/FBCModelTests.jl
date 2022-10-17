#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

#ToDo: change all instances of ::MetabolicModel to ::AbstractMetabolicModel once COBREXA 2.0 releases


"""
$(TYPEDSIGNATURES)

Checks if every gene has an annotation and returns an array of genes which do not have annotatons.
"""
function all_unannotated_genes(model::MetabolicModel)
    [gid for gid in genes(model) if isempty(gene_annotations(model, gid))]
end

"""
$(TYPEDSIGNATURES)

Checks if the annotation_keywords are present in the annotations and returns an array of genes which do not include the keyword.
"""

function unannotated_genes(
    model::MetabolicModel;
    annotation_keywords = [
        "kegg.genes",
        "refseq",
        "uniprot",
        "ecogene",
        "ncbigi",
        "ncbigene",
        "ncbiprotein",
        "ccds",
        "hprd",
        "asap",
    ],
)
    missing_annos = Dict{String,Vector{String}}()
    for anno_kw in annotation_keywords
        missing_annos[anno_kw] = []
        for gid in genes(model) # accessor: reactions
            annos = gene_annotations(model, gid) # accessor: reaction_annotations
            !haskey(annos, anno_kw) && push!(missing_annos[anno_kw], gid)

        end
    end
    return missing_annos
end

"""
$(TYPEDSIGNATURES)

Checks if the entry of the genes which have the annotation keyword are formally correct and returns a Dict with the keywords as keys and
an array of the formally incorrect genes as values.
"""

function gene_annotation_conformity(
    model::MetabolicModel;
    GENE_PRODUCT_ANNOTATIONS = Dict(
        "refseq" =>
            r"^((AC|AP|NC|NG|NM|NP|NR|NT|NW|XM|XP|XR|YP|ZP)_\d+|(NZ\_[A-Z]{4}\d+))(\.\d+)?$",
        "uniprot" =>
            r"^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$",
        "ecogene" => r"^EG\d+$",
        "kegg.genes" => r"^\w+:[\w\d\.-]*$",
        "ncbigi" => r"^(GI|gi)\:\d+$",
        "ncbigene" => r"^\d+$",
        "ncbiprotein" => r"^(\w+\d+(\.\d+)?)|(NP_\d+)$",
        "ccds" => r"^CCDS\d+\.\d+$",
        "hprd" => r"^\d+$",
        "asap" => r"^[A-Za-z0-9-]+$",
    ),
)
    no_annos = Dict{String,Vector{String}}()
    for gid in genes(model)
        annos = gene_annotations(model, gid)
        for anno_kw in keys(GENE_PRODUCT_ANNOTATIONS)
            if haskey(annos, anno_kw) == true
                !in(anno_kw, keys(no_annos)) && (no_annos[anno_kw] = [])
                in.(
                    nothing,
                    Ref(match.(GENE_PRODUCT_ANNOTATIONS[anno_kw], annos[anno_kw])),
                ) && push!(no_annos[anno_kw], gid)
            end
        end
    end
    return no_annos
end
