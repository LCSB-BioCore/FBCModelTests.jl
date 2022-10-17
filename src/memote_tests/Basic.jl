"""
$(TYPEDEF)

Basic metabolic model tests.
"""
module Basic

using ..DocStringExtensions
using ..ModuleTools
using ..COBREXA
import ..Config.memote_config

"""
$(TYPEDSIGNATURES)

Test if the model has an id. Note, this test will fail unless you use a
`StandardModel` and its id field is assigned.
"""
model_has_name(model) = begin
    model isa StandardModel && model.id != "" && return true
    return false # TODO need accessors for the model name
end

"""
$(TYPEDSIGNATURES)

Test if the model has any metabolites.
"""
model_has_metabolites(model) = n_metabolites(model) > 0

"""
$(TYPEDSIGNATURES)

Test if the model has any reactions.
"""
model_has_reactions(model) = n_reactions(model) > 0

"""
$(TYPEDSIGNATURES)

Test if the model has any genes.
"""
model_has_genes(model) = n_genes(model) > 0

"""
$(TYPEDSIGNATURES)

Calculate the metabolic coverage by dividing the number of reactions by the
number of genes.
"""
model_metabolic_coverage(model) = n_reactions(model) / n_genes(model)

"""
$(TYPEDSIGNATURES)

Return the number of unique compartments in the model.
"""
model_compartments(model) =
    Set(metabolite_compartment(model, mid) for mid in metabolites(model))

"""
$(TYPEDSIGNATURES)

Test if the model has one or more compartments.
"""
model_has_compartments(model) = length(model_compartments(model)) > 0

"""
$(TYPEDSIGNATURES)

Test if the metabolic coverage, calculated with
[`model_metabolic_coverage`](@ref), exceeds
`config.basic.minimum_metabolic_coverage`.
"""
model_metabolic_coverage_exceeds_minimum(model; config = memote_config) =
    model_metabolic_coverage(model) > config.basic.minimum_metabolic_coverage

"""
$(TYPEDSIGNATURES)

Check if the model can be solved under default conditions and yield a reasonable
growth rate. Here reasonable is set via `config.basic.minimum_growth_rate` and
`config.basic.maximum_growth_rate`. Optionally, pass optimization
modifications to the solver through `config.basic.optimizer_modifications`.
"""
function model_solves_in_default_medium(model, optimizer; config = memote_config)
    mu = solved_objective_value(
        flux_balance_analysis(
            model,
            optimizer;
            modifications = config.biomass.optimizer_modifications,
        ),
    )
    config.basic.minimum_growth_rate < mu < config.basic.maximum_growth_rate
end

@export_locals

end