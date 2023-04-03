"""
    module Annotation

This module contains tests used to check the coverage and conformance of
reaction, gene, and metabolite annotations.
"""
module Annotation

using DocStringExtensions
using COBREXA

import ..Config
import ..Utils: parse_annotations

is_annotated_reaction(model, rid) = !isempty(reaction_annotations(model, rid))
is_annotated_metabolite(model, mid) = !isempty(metabolite_annotations(model, mid))
is_annotated_gene(model, gid) = !isempty(gene_annotations(model, gid))

"""
$(TYPEDSIGNATURES)

Helper function to find all the annotations that are missing for a component in a model.
"""
function _find_unannotated_components(
    model::MetabolicModel,
    id_accessor::String,
    annotation_accessor,
    annotation_kws,
)
    missing_annos = String[]
    for anno_kw in annotation_kws
        annos = parse_annotations(annotation_accessor(model, id_accessor))
        !haskey(annos, anno_kw) && push!(missing_annos, anno_kw)
    end
    return missing_annos
end

"""
$(TYPEDSIGNATURES)

Checks if the databases listed in `config.annotation.gene_annotation_keywords`
are present in the gene annotations. Returns a vector of annotation keywords
that were not found.
"""
findall_unannotated_gene_databases(model, gid::String; config = Config.memote_config) =
    _find_unannotated_components(
        model,
        gid,
        gene_annotations,
        config.annotation.gene_annotation_keywords,
    )

"""
$(TYPEDSIGNATURES)

Checks if the databases listed in `config.annotation.metabolite_annotation_keywords`
are present in the metabolite annotations. Returns a vector of annotation keywords
that were not found.
"""
findall_unannotated_metabolite_databases(
    model,
    mid::String;
    config = Config.memote_config,
) = _find_unannotated_components(
    model,
    mid,
    metabolite_annotations,
    config.annotation.metabolite_annotation_keywords,
)

"""
$(TYPEDSIGNATURES)

Checks if the databases listed in `config.annotation.reaction_annotation_keywords`
are present in the reaction annotations. Returns a vector of annotation keywords
that were not found.
"""
findall_unannotated_reaction_databases(model, rid::String; config = Config.memote_config) =
    _find_unannotated_components(
        model,
        rid,
        reaction_annotations,
        config.annotation.reaction_annotation_keywords,
    )

"""
$(TYPEDSIGNATURES)

Helper function to find all the annotations for a component that do not conform in the model.
"""
function _find_nonconformal_annotations(
    model::MetabolicModel,
    id_accessor::String,
    annotation_accessor,
    annotation_regex,
)
    no_annos = String[]
    for (anno_kw, anno_regex) in annotation_regex
        annos = parse_annotations(annotation_accessor(model, id_accessor))
        if haskey(annos, anno_kw)
            in.(nothing, Ref(match.(anno_regex, annos[anno_kw]))) &&
                push!(no_annos, anno_kw)
        end
    end
    return no_annos
end

"""
$(TYPEDSIGNATURES)

Check if the gene annotation entry conforms to commonly recognized formats of
annotation database using regex patterns. Uses the database formats listed in
`config.annotation.gene_annotation_regexes` to test the conformity. Returns a
string vector of database ids do not conform.
"""
findall_nonconformal_gene_annotations(model, gid::String; config = Config.memote_config) =
    _find_nonconformal_annotations(
        model,
        gid,
        gene_annotations,
        config.annotation.gene_annotation_regexes,
    )

"""
$(TYPEDSIGNATURES)

Check if the metabolite annotation entry conforms to commonly recognized formats
of annotation database using regex patterns. Uses the database formats listed in
`config.annotation.metabolite_annotation_regexes` to test the conformity.
Returns a string vector of database ids do not conform.
"""
findall_nonconformal_metabolite_annotations(
    model,
    mid::String;
    config = Config.memote_config,
) = _find_nonconformal_annotations(
    model,
    mid,
    metabolite_annotations,
    config.annotation.metabolite_annotation_regexes,
)

"""
$(TYPEDSIGNATURES)

Check if the reaction annotation entry conforms to commonly recognized formats
of annotation database using regex patterns. Uses the database formats listed in
`config.annotation.reaction_annotation_regexes` to test the conformity. Returns
a string vector of database ids do not conform.
"""
findall_nonconformal_reaction_annotations(
    model,
    rid::String;
    config = Config.memote_config,
) = _find_nonconformal_annotations(
    model,
    rid,
    reaction_annotations,
    config.annotation.reaction_annotation_regexes,
)

end # module
