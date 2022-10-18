"""
$(TYPEDEF)

This module contains tests used to check the coverage and conformance of
reaction, gene, and metabolite annotations.
"""
module Annotation

using ..COBREXA

using ..DocStringExtensions
import ..Config.memote_config

"""
$(TYPEDSIGNATURES)

Checks if every gene has an annotation and returns an array of genes that do not
have annotatons.
"""
find_all_unannotated_genes(model) =
    [gid for gid in genes(model) if isempty(gene_annotations(model, gid))]

"""
$(TYPEDSIGNATURES)

Checks if every metabolite has an annotation and returns an array of metabolites
that do not have annotatons.
"""
find_all_unannotated_metabolites(model) =
    [mid for mid in metabolites(model) if isempty(metabolite_annotations(model, mid))]

"""
$(TYPEDSIGNATURES)

Checks if every reaction has an annotation and returns an array of reactions
that do not have annotatons.
"""
find_all_unannotated_reactions(model) =
    [rid for rid in reactions(model) if isempty(reaction_annotations(model, rid))]

"""
$(TYPEDSIGNATURES)

Helper function to find all the unannotated components in the model.
"""
function _find_unannotated_components(
    model,
    id_accessor,
    annotation_accessor,
    annotation_kws,
)
    missing_annos = Dict{String,Vector{String}}()
    for anno_kw in annotation_kws
        missing_annos[anno_kw] = String[]
        for id in id_accessor(model)
            annos = annotation_accessor(model, id)
            !haskey(annos, anno_kw) && push!(missing_annos[anno_kw], id)
        end
    end
    return missing_annos
end

"""
$(TYPEDSIGNATURES)

Checks if the databases listed in `config.annotation.gene_annotation_keywords`
are present in the gene annotations. Returns a dictionary of annotation keywords
mapped to a list of genes that do not include the keyword.
"""
find_database_unannotated_genes(model; config = memote_config) =
    _find_unannotated_components(
        model,
        genes,
        gene_annotations,
        config.annotation.gene_annotation_keywords,
    )

"""
$(TYPEDSIGNATURES)

Checks if the databases listed in
`config.annotation.metabolite_annotation_keywords` are present in the metabolite
annotations. Returns a dictionary of annotation keywords mapped to a list of
metabolites that do not include the keyword.
"""
find_database_unannotated_metabolites(model; config = memote_config) =
    _find_unannotated_components(
        model,
        metabolites,
        metabolite_annotations,
        config.annotation.metabolite_annotation_keywords,
    )

"""
$(TYPEDSIGNATURES)

Checks if the databases listed in
`config.annotation.reaction_annotation_keywords` are present in the reaction
annotations. Returns a dictionary of annotation keywords mapped to a list of
reactions that do not include the keyword.
"""
find_database_unannotated_reactions(model; config = memote_config) =
    _find_unannotated_components(
        model,
        reactions,
        reaction_annotations,
        config.annotation.reaction_annotation_keywords,
    )

"""
$(TYPEDSIGNATURES)

Helper function to find all the annotations that do not conform in the model.
"""
function _find_nonconformal_components(
    model,
    id_accessor,
    annotation_accessor,
    annotation_regex,
)
    no_annos = Dict{String,Vector{String}}()
    for (anno_kw, anno_regex) in annotation_regex
        no_annos[anno_kw] = String[]
        for id in id_accessor(model)
            annos = annotation_accessor(model, id)
            if haskey(annos, anno_kw)
                in.(nothing, Ref(match.(anno_regex, annos[anno_kw]))) &&
                    push!(no_annos[anno_kw], id)
            end
        end
    end
    return no_annos
end

"""
$(TYPEDSIGNATURES)

Check if the gene annotation entry conforms to commonly recognized formats of
annotation database using regex patterns. Uses the database formats listed in
`config.annotation.gene_annotation_regexes` to test the conformity. Returns a
dictionary mapping the database id to a list of genes that do not conform.
"""
find_nonconformal_gene_annotations(model; config = memote_config) =
    _find_nonconformal_components(
        model,
        genes,
        gene_annotations,
        config.annotation.gene_annotation_regexes,
    )

"""
$(TYPEDSIGNATURES)

Check if the metabolite annotation entry conforms to commonly recognized formats
of annotation database using regex patterns. Uses the database formats listed in
`config.annotation.metabolite_annotation_regexes` to test the conformity.
Returns a dictionary mapping the database id to a list of metabolites that do
not conform.
"""
find_nonconformal_metabolite_annotations(model; config = memote_config) =
    _find_nonconformal_components(
        model,
        genes,
        gene_annotations,
        config.annotation.metabolite_annotation_regexes,
    )

"""
$(TYPEDSIGNATURES)

Check if the reaction annotation entry conforms to commonly recognized formats
of annotation database using regex patterns. Uses the database formats listed in
`config.annotation.reaction_annotation_regexes` to test the conformity. Returns
a dictionary mapping the database id to a list of reactions that do not conform.
"""
find_nonconformal_reaction_annotations(model; config = memote_config) =
    _find_nonconformal_components(
        model,
        reactions,
        reaction_annotations,
        config.annotation.reaction_annotation_regexes,
    )

end # module