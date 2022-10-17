"""
$(TYPEDSIGNATURES)

Run the standard memote test suite on `model` using `optimizer` to solve
optimization problems and configure the test parameters with `config`.

Note, for best results, make sure the model can be converted into a
`COBREXA.StandardModel`.
"""
function run_tests(model, optimizer; config = memote_config)
    @testset "Metabolic model tests" begin
        
        @testset "Basic information" begin
            @test model_has_name(model)
            @test model_has_reactions(model)
            @test model_has_metabolites(model)
            @test model_has_genes(model)
            @test model_metabolic_coverage_exceeds_minimum(model; config)
            @test model_has_compartments(model)
        end

        @testset "Reaction annotations" begin
            
        end

        @testset "Metabolite annotations" begin
            
        end

        @testset "Gene annotations" begin
            
        end

        @testset "SBO term annotations" begin
            
        end

        @testset "Reaction information" begin
            
        end

        @testset "Metabolite information" begin
            
        end

        @testset "GPR associations" begin
            
        end

        @testset "Consistency" begin
            
        end

        @testset "Biomass" begin
            
        end

        @testset "Energy metabolism" begin
            
        end

        @testset "Network topology" begin
            
        end

        @testset "Matrix conditioning" begin
            
        end

    end
end

"""
$(TYPEDSIGNATURES)

Generate a report of model characteristics that are typically important measures
of the scope of the model.
"""
function generate_memote_report(model, optimizer; config = memote_config)
    # Basic information

    # Reaction annotations

    # Metabolite annotations

    # Gene annotations

    # SBO term annotations

    # Reaction information

    # Metabolite information

    # Gene protein reaction associations

    # Consistency

    # Biomass

    # Energy metabolism

    # Network topology

    # Matrix conditioning
end