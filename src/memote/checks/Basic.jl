"""
    module Basic

Basic metabolic model tests.
"""
module Basic

using DocStringExtensions
using COBREXA

import ..Config

"""
$(TYPEDSIGNATURES)

Test if the model has any metabolites.
"""
model_has_metabolites(model::MetabolicModel) = n_metabolites(model) > 0

"""
$(TYPEDSIGNATURES)

Test if the model has any reactions.
"""
model_has_reactions(model::MetabolicModel) = n_reactions(model) > 0

"""
$(TYPEDSIGNATURES)

Test if the model has any genes.
"""
model_has_genes(model::MetabolicModel) = n_genes(model) > 0

"""
$(TYPEDSIGNATURES)

Calculate the metabolic coverage by dividing the number of reactions by the
number of genes.
"""
model_metabolic_coverage(model::MetabolicModel) = n_reactions(model) / n_genes(model)

"""
$(TYPEDSIGNATURES)

Return the number of unique compartments in the model.
"""
model_compartments(model::MetabolicModel) = Set(
    metabolite_compartment(model, mid) for
    mid in metabolites(model) if !isnothing(metabolite_compartment(model, mid))
)

"""
$(TYPEDSIGNATURES)

Test if the model has one or more compartments.
"""
model_has_compartments(model::MetabolicModel) = length(model_compartments(model)) > 0

"""
$(TYPEDSIGNATURES)

Check if the model can be solved under default conditions and yield a reasonable
growth rate. Here reasonable is set via `config.basic.minimum_growth_rate` and
`config.basic.maximum_growth_rate`. Optionally, pass optimization
modifications to the solver through `config.basic.optimizer_modifications`.
"""
function model_solves_in_default_medium(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
)
    mu = solved_objective_value(
        flux_balance_analysis(
            model,
            optimizer;
            modifications = config.basic.optimizer_modifications,
        ),
    )
    isnothing(mu) && return false
    config.basic.minimum_growth_rate < mu < config.basic.maximum_growth_rate
end

end
