#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

#TODO: change all instances of ::MetabolicModel to ::AbstractMetabolicModel once COBREXA 2.0 releases

"""
$(TYPEDSIGNATURES)

Checks for the presence of reaction annotations.
Returns a vector of all reactions without any annotations.
"""
function all_unannotated_reactions(model::MetabolicModel)
    return [reaction_id for reaction_id in reactions(model) if isempty(reaction_annotations(model, reaction_id))]
end

"""
$(TYPEDSIGNATURES)

Iterates through a model's reactions to check if any common biochemical databases appear in the annotations field.
Returns a dictionary with the biochemical database names being checked for as keywords 
and the corresponding value being a vector of all reactions missing that database name.
"""
function unannotated_reactions(model::MetabolicModel, annotation_keywords = ["rhea", "kegg.reaction", 
    "seed.reaction", "metanetx.reaction", "bigg.reaction", "reactome", "ec-code", "brenda", "biocyc"])
    missing_annos = Dict{String, Vector{String}}()
    for keyword in annotation_keywords
        missing_annos[keyword] = []
        for reaction_id in reactions(model)
            reaction_annos = reaction_annotations(model, reaction_id)
            !haskey(reaction_annos, keyword) && push!(missing_annos[keyword], reaction_id)
        end
    end
    return missing_annos
end

"""
$(TYPEDSIGNATURES)

Checks if the entry of the reaction which have the annotation keyword are formally correct 
and returns a Dict with the keywords as keys and 
an array of the formally incorrect reactions as values.
"""

function reactions_annotation_conformity(
    model::MetabolicModel;
    REACTION_ANNOTATIONS = Dict("rhea" => r"^\d{5}$", 
    "kegg.reaction" => r"^R\d+$",
    "seed.reaction" => r"^rxn\d+$",
    "metanetx.reaction" => r"^MNXR\d+$",
    "bigg.reaction" => r"^[a-z_A-Z0-9]+$",
    "reactome" => r"(^R-[A-Z]{3}-[0-9]+(-[0-9]+)?$)|(^REACT_\d+(\.\d+)?$)",
    "ec-code" => r"^\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$",
    "brenda" => r"^\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$",
    "biocyc" => r"^[A-Z-0-9]+(?<!CHEBI)(\:)?[A-Za-z0-9+_.%-]+$"
    )
    
  )
  no_annos = Dict{String, Vector{String}}()
  for rid in reactions(model) 
    annos = reaction_annotations(model, rid)
    for anno_kw in keys(REACTION_ANNOTATIONS)
      if haskey(annos, anno_kw) == true 
        !in(anno_kw, keys(no_annos)) && (no_annos[anno_kw] = [])
        in.(nothing, Ref(match.(REACTION_ANNOTATIONS[anno_kw], annos[anno_kw]))) && push!(no_annos[anno_kw], rid)
      end  
    end
  end
  return no_annos
end
