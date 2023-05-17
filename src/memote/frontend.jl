"""
$(TYPEDSIGNATURES)

Run a MEMOTE-like test suite on `model` using `optimizer` to solve optimization
problems; some basic parameters and thresholds are taken from `config`.

Some of the tests internally convert the input model to `StandardModel` of
COBREXA; supplying a `StandardModel` may thus increase efficiency of the whole
process.
"""
function run_tests(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
    filename = nothing,
    workers = [myid()],
)
    @info "Starting MEMOTE tests..."
    @testset QuietTestSet "Metabolic Model Tests$(isnothing(filename) ? "" : ": $filename")" begin
        run_testsets(model, optimizer; config, workers)
    end
end

"""
$(TYPEDSIGNATURES)

Overload of [`run_tests`](@ref) that works directly with a file.
"""
run_tests(filename::String, optimizer; kwargs...) = run_tests(
    load_model(StandardModel, filename),
    optimizer;
    filename = basename(filename),
    kwargs...,
)

"""
$(TYPEDSIGNATURES)

Execute all tests.
"""
function run_testsets(
    model::MetabolicModel,
    optimizer;
    config = Config.memote_config,
    workers = [myid()],
)

    @info "Testing basic information about the model..."
    @testset "Basic information" begin
        @testset "Has reactions" begin
            @atest Basic.model_has_reactions(model) "The model has no reactions."
        end
        @testset "Has metabolites" begin
            @atest Basic.model_has_metabolites(model) "The model has no metabolites."
        end
        @testset "Has genes" begin
            @atest Basic.model_has_genes(model) "The model has no genes."
        end
        @testset "Has compartments" begin
            @atest Basic.model_has_compartments(model) "The model has no comparments."
        end
        @testset "Metabolic coverage exceeds threshold" begin
            @atest config.basic.minimum_metabolic_coverage <=
                   Basic.model_metabolic_coverage(model) "The metabolic coverage of the model is too small."
        end
        @testset "Solves in default medium" begin
            @atest Basic.model_solves_in_default_medium(model, optimizer; config) "The model does not solve in the default medium."
        end
    end

    @info "Testing the model annotations..."
    @testset "Annotations" begin
        @testset "SBO annotations" begin
            @testset "Has SBO annotated metabolic reactions" begin
                @atest any(is_metabolic_reaction(model, rid) for rid in reactions(model)) "No reactions are SBO annotated as metabolic reactions."
            end
            @testset "Has SBO annotated transport reactions" begin
                @atest any(is_transport_reaction(model, rid) for rid in reactions(model)) "No reactions are SBO annotated as transport reactions."
            end
            @testset "Has SBO annotated biomass reactions" begin
                @atest any(is_biomass_reaction(model, rid) for rid in reactions(model)) "No reactions are SBO  annotated as biomass reactions."
            end
            @testset "Has SBO annotated exchange reactions" begin
                @atest any(is_exchange_reaction(model, rid) for rid in reactions(model)) "No reactions are SBO annotated as exchange reactions."
            end
            @testset "Has SBO annotated atp maintenance reactions" begin
                @atest any(
                    is_atp_maintenance_reaction(model, rid) for rid in reactions(model)
                ) "No reactions are SBO annotated as ATP maintenance reactions."
            end
        end

        @testset "Reactions" begin
            @testset "Any annotations present" begin
                for rid in reactions(model)
                    @atest Annotation.is_annotated_reaction(model, rid) "$rid does not have any annotations."
                end
            end

            @testset "At most $(config.annotation.maximum_missing_databases) missing cross-references per entry." begin
                for rid in reactions(model)
                    @atest length(
                        Annotation.findall_unannotated_reaction_databases(
                            model,
                            rid;
                            config,
                        ),
                    ) <= config.annotation.maximum_missing_databases "$rid does not have enough cross-references."
                end
            end
            @testset "At most $(config.annotation.maximum_nonconformal_references) unrecognizable cross-references per entry" begin
                for rid in reactions(model)
                    @atest length(
                        Annotation.findall_nonconformal_reaction_annotations(
                            model,
                            rid;
                            config,
                        ),
                    ) <= config.annotation.maximum_nonconformal_references "$rid does not have enough recognizable cross-references."
                end
            end
        end

        @testset "Metabolites" begin
            @testset "Any annotations present" begin
                for mid in metabolites(model)
                    @atest Annotation.is_annotated_metabolite(model, mid) "$mid does not have any annotations."
                end
            end

            @testset "At most $(config.annotation.maximum_missing_databases) missing cross-references per entry" begin
                for mid in metabolites(model)
                    @atest length(
                        Annotation.findall_unannotated_metabolite_databases(
                            model,
                            mid;
                            config,
                        ),
                    ) <= config.annotation.maximum_missing_databases "$mid does not have enough cross-references."
                end
            end
            @testset "At most $(config.annotation.maximum_nonconformal_references) unrecognizable cross-references per entry" begin
                for mid in metabolites(model)
                    @atest length(
                        Annotation.findall_nonconformal_metabolite_annotations(
                            model,
                            mid;
                            config,
                        ),
                    ) <= config.annotation.maximum_nonconformal_references "$mid does not have enough recognizable cross-references."
                end
            end
        end


        @testset "Genes" begin
            @testset "Any annotations present" begin
                for gid in genes(model)
                    @atest Annotation.is_annotated_gene(model, gid) "$gid does not have any annotations."
                end
            end

            @testset "At most $(config.annotation.maximum_missing_databases) missing cross-references per entry" begin
                for gid in genes(model)
                    @atest length(
                        Annotation.findall_unannotated_gene_databases(model, gid; config),
                    ) <= config.annotation.maximum_missing_databases "$gid does not have enough cross-references."
                end
            end
            @testset "At most $(config.annotation.maximum_nonconformal_references) unrecognizable cross-references per entry" begin
                for gid in genes(model)
                    @atest length(
                        Annotation.findall_nonconformal_gene_annotations(
                            model,
                            gid;
                            config,
                        ),
                    ) <= config.annotation.maximum_nonconformal_references "$gid does not have enough recognizable cross-references."
                end
            end
        end
    end

    @info "Testing the biomass function(s)..."
    @testset "Biomass" begin
        @testset "Is present" begin
            @atest length(Biomass.findall_biomass_reactions(model)) > 0 "No biomass reactions present."
        end
        @testset "ATP burn included" begin
            for rid in Biomass.findall_biomass_reactions(model)
                @atest Biomass.atp_present_in_biomass_reaction(model, rid) "$rid (biomass reaction) does not consume ATP."
            end
        end
        @testset "Molar mass consistent" begin
            for rid in Biomass.findall_biomass_reactions(model)
                @atest Biomass.biomass_reaction_is_consistent(model, rid) "$rid (biomass reaction) is not consistent."
            end
        end
        @testset "Missing common precursors" begin
            for rid in Biomass.findall_biomass_reactions(model)
                @atest Biomass.biomass_missing_essential_precursors(model, rid; config) "$rid (biomass reaction) does not contain all the common precursors."
            end
        end
    end

    @info "Testing the stoichiometric consistency of the model..."
    @testset "Stoichiometrically consistent" begin
        @atest Consistency.model_is_consistent(model, optimizer; config) "The model is not stoichiometrically consistent."
    end

    @info "Testing the energy metabolism..."
    @testset "Energy metabolism" begin # TODO make this play nicely with compartments
        @testset "Erroneous energy cycles" begin
            @test_skip Energy.model_has_no_erroneous_energy_generating_cycles(
                model,
                optimizer;
                config,
            )
        end
    end

    @info "Testing the gene-protein-reaction associations..."
    @testset "GPR associations" begin
        @testset "Metabolic reactions have GPRs" begin
            for rid in filter(x -> is_metabolic_reaction(model, x), reactions(model))
                @atest GPRAssociation.reaction_has_sensible_gpr(model, rid) "$rid (metabolic reaction) does not have any GPRs."
            end
        end

        @testset "Transport reactions have GPRs" begin
            for rid in filter(x -> is_transport_reaction(model, x), reactions(model))
                @atest GPRAssociation.reaction_has_sensible_gpr(model, rid) "$rid (transport reaction) does not have any GPRs."
            end
        end

        @testset "At least one enzyme complex exists" begin
            @atest length(GPRAssociation.reactions_with_complexes(model)) != 0 "No enzyme complexes exist in the model."
        end
    end

    @info "Testing the reaction information..."
    @testset "Reaction information" begin
        tms_rxns = [
            rid for rid in reactions(model) if is_transport_reaction(model, rid) ||
            is_metabolic_reaction(model, rid) ||
            is_spontaneous_reaction(model, rid)
        ]
        @testset "Charge balanced (transport, metabolic, spontaneous)" begin
            for rid in tms_rxns
                if rid in config.reaction.charge_ignored_reactions
                    @test_skip Reactions.reaction_is_charge_balanced(model, rid)
                else
                    @atest Reactions.reaction_is_charge_balanced(model, rid) "$rid is not charge balanced."
                end
            end
        end
        @testset "Mass balanced (transport, metabolic, spontaneous)" begin
            for rid in tms_rxns
                if rid in config.reaction.mass_ignored_reactions
                    @test_skip Reactions.reaction_is_mass_balanced(model, rid)
                else
                    @atest Reactions.reaction_is_mass_balanced(model, rid) "$rid is not mass balanced."
                end
            end
        end
        @testset "No duplicated reactions" begin
            dup_rxns = Reactions.findall_duplicated_reactions(model)
            for rid in reactions(model)
                @atest rid ∉ dup_rxns "$rid is duplicated"
            end
        end
        @testset "ATP maintenance reaction is present" begin
            @atest Reactions.model_has_atpm_reaction(model) "There is no atp maintenance reaction in the model."
        end
    end

    @info "Testing the metabolite information..."
    @testset "Metabolite information" begin
        @testset "No missing formulas" begin
            for mid in metabolites(model)
                @atest !isnothing(metabolite_formula(model, mid)) "$mid has no formula."
            end
        end
        @testset "No missing charges" begin
            for mid in metabolites(model)
                @atest !isnothing(metabolite_charge(model, mid)) "$mid has no charge."
            end
        end
        @testset "Duplicated metabolites" begin
            dup_mids = Metabolites.metabolites_duplicated_in_compartment(model; config)
            for mid in metabolites(model)
                @atest mid ∉ dup_mids "$mid is duplicated."
            end
        end
        @testset "No orphan metabolites" begin
            orphans = Metabolites.find_orphan_metabolites(model)
            for mid in metabolites(model)
                @atest mid ∉ orphans "$mid is an orphan metabolite."
            end
        end
        @testset "No deadend metabolites" begin
            deadends = Metabolites.find_deadend_metabolites(model)
            for mid in metabolites(model)
                @atest mid ∉ deadends "$mid is a deadend metabolite."
            end
        end
        @testset "No disconnected metabolites" begin
            disconnected = Metabolites.find_disconnected_metabolites(model)
            for mid in metabolites(model)
                @atest mid ∉ disconnected "$mid is disconnected."
            end
        end
    end

    @info "Testing the network information..."
    @testset "Network information" begin
        @testset "Stoichiometric matrix is well conditioned" begin
            @atest Network.stoichiometric_max_min_ratio(model) <
                   config.network.condition_number "The stoichiometric matrix is poorly conditioned."
        end
        @testset "No stoichiometrically balanced cycles (reactions)" begin
            rxns_in_cycles = Network.find_cycle_reactions(model, optimizer; config, workers)
            for rid in reactions(model)
                @atest rid ∉ rxns_in_cycles "$rid takes part in a balanced cycle."
            end
        end
        @testset "No universally blocked reactions" begin
            blocked_rxns = Network.find_all_universally_blocked_reactions(
                model,
                optimizer;
                config,
                workers,
            )
            for rid in reactions(model)
                @atest rid ∉ blocked_rxns "$rid is a blocked reaction."
            end
        end
    end
end
