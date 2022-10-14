#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Return a list of all boundary reactions that allow flux into the model and
create a metabolite. Assume that boundary reactions only create a single
metabolite. Use the testing defaults, `default.metabolite.only_imported =
false`, to also return metabolites that can be produced by the model under
default conditions.
"""
function metabolites_medium_components(model; defaults = memote_defaults)
    mets = String[]
    for (rid, lb, ub) in zip(reactions(model), bounds(model)...)
        if is_boundary(model, rid)
            mid = first(keys(reaction_stoichiometry(model, rid)))
            defaults.metabolite.medium_only_imported &&
                lb < 0 &&
                ub <= 0 &&
                push!(mets, mid)
            !defaults.metabolite.medium_only_imported && lb < 0 && push!(mets, mid)
        end
    end
    mets
end

"""
$(TYPEDSIGNATURES)

List all metabolites without a formula. Use
`defaults.metabolite.formula_corner_cases` to specify an extra case to check for
formula's that are not properly assigned.
"""
metabolites_no_formula(model; defaults = memote_defaults) = [
    mid for mid in metabolites(model) if isnothing(metabolite_formula(model, mid)) || isempty(metabolite_formula(model, mid)) ||
    any(in.(keys(metabolite_formula(model, mid)), Ref(defaults.metabolite.formula_corner_cases)))
]

"""
$(TYPEDSIGNATURES)

List all metabolites without a charge. Use
`defaults.metabolite.charge_corner_cases` to specify an extra case to check for
charge's that are not properly assigned.
"""
metabolites_no_charge(model; defaults = memote_defaults) = [
    mid for mid in metabolites(model) if isnothing(metabolite_charge(model, mid)) ||
    metabolite_charge(model, mid) in defaults.metabolite.charge_corner_cases
]

"""
$(TYPEDSIGNATURES)

Test if metabolites `m1` and `m2` are different by comparing their
`defaults.metabolite.test_annotation` field in the annotations of each
metabolite. Note, if no annotations are present for one or both of the
metabolites, then return `true`. 
"""
function metabolites_are_duplicated(model, m1, m2; defaults = memote_defaults)
    k1s =
        get(metabolite_annotations(model, m1), defaults.metabolite.test_annotation, nothing)
    isnothing(k1s) && return true
    k2s =
        get(metabolite_annotations(model, m2), defaults.metabolite.test_annotation, nothing)
    isnothing(k2s) && return true
    any(in.(k1s, Ref(k2s)))
end

"""
$(TYPEDSIGNATURES)

Return a list of unique metabolites in model. Uses
[`metabolites_are_duplicated`](@ref) internally and forwards `test_annotation`
to it. The latter argument is used to determine if two metabolites are the same
by checking for any correspondence.
"""
function metabolites_unique(model; defaults = memote_defaults)
    unique_metabolites = Set[]
    for m1 in metabolites(model)
        duplicate = false
        for m2 in unique_metabolites
            duplicate = metabolites_are_duplicated(model, m1, m2; defaults)
            duplicate && break
        end
        !duplicate && push!(unique_metabolites, m1)
    end
    return unique_metabolites
end

"""
$(TYPEDSIGNATURES)

Return a dictionary of metabolites that are duplicated in their compartment. 
"""
function metabolites_duplicated_in_compartment(model; defaults = memote_defaults)
    unique_metabolites = Dict{String,Set{String}}()
    for m1 in metabolites(model)
        c1 = metabolite_compartment(model, m1)
        for m2 in metabolites(model)
            c2 = metabolite_compartment(model, m2)
            if c1 == c2 && m1 != m2 && metabolites_are_duplicated(model, m1, m2; defaults)
                if haskey(unique_metabolites, c1)
                    push!(unique_metabolites[c1], m1)
                else
                    unique_metabolites[c1] = [m1]
                end
            end
        end
    end
    return unique_metabolites
end

"""
$(TYPEDSIGNATURES)

Test if the metabolites contained in the `model`:
1. are not duplicated, tested with [`metabolites_duplicated_in_compartment`](@ref)
2. all have a formula and charge associated with them, tested with [`metabolites_no_formula`](ref) and [`metabolites_no_charge`](@ref)
3. the default medium of the cell is not empty, tested with [`metabolites_medium_components`](ref).

Each of the basic functions can be run independently. THe kwargs are forwarded as indicated by the prefix.
"""
function test_metabolites(model; defaults = memote_defaults)
    @testset "Metabolite Information" begin
        @test isempty(metabolites_duplicated_in_compartment(model; defaults))
        @test isempty(metabolites_no_formula(model; defaults))
        @test isempty(metabolites_no_charge(model; defaults))
        @test !isempty(metabolites_medium_components(model; defaults))
    end
end
