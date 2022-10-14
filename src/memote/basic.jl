#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

function model_has_name(model)
    model isa StandardModel && model.id != "" && return true
    return true # TODO need accessors for the model name
end

model_has_metabolites(model) = n_metabolites(model) > 0

model_has_reactions(model) = n_reactions(model) > 0

model_has_genes(model) = n_genes(model) > 0

model_metabolite_coverage(model) = n_reactions(model) / n_genes(model)

model_compartments(model) = Set(metabolite_compartment(model, mid) for mid in metabolites(model))

model_has_compartments(model) = length(model_compartments(model)) > 0

function test_basic(model)
    @testset "Basic information" begin
        @test model_has_name(model)
        @test model_has_metabolites(model)
        @test model_has_reactions(model)
        @test model_has_genes(model)
        @test model_has_compartments(model)
    end
end