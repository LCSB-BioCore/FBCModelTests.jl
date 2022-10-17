function model_has_name(model)
    model isa StandardModel && model.id != "" && return true
    return false # TODO need accessors for the model name
end

model_has_metabolites(model) = n_metabolites(model) > 0

model_has_reactions(model) = n_reactions(model) > 0

model_has_genes(model) = n_genes(model) > 0

model_metabolic_coverage(model) = n_reactions(model) / n_genes(model)

model_compartments(model) =
    Set(metabolite_compartment(model, mid) for mid in metabolites(model))

model_has_compartments(model) = length(model_compartments(model)) > 0

model_metabolic_coverage_exceeds_minimum(model; config = memote_config) =
    model_metabolic_coverage(model) > config.basic.minimum_metabolic_coverage
