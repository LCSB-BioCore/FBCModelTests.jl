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
metabolite. Set the kwarg `only_into=false` to also return metabolites that can
be produced by the model under default conditions.
"""
function metabolites_medium_components(model; only_into=true)
    mets = String[]
    for (rid, lb, ub) in zip(reactions(model), bounds(model)...)
        if is_boundary(model, rid)
            mid = first(keys(reaction_stoichiometry(model, rid))) 
            only_into && lb < 0 && ub <= 0 && push!(mets, mid)
            !only_into && lb < 0 && push!(mets, mid)
        end
    end
    mets
end

"""
$(TYPEDSIGNATURES)

List all metabolites without a formula. Use `corner_cases` to specify an extra
case to check for formula's that are not properly assigned.
"""
metabolites_no_formula(model; corner_cases = ["X", "x", ""]) = [mid for mid in metabolites(model) if isnothing(metabolite_formula(model, mid)) || metabolite_formula(model, mid) in corner_cases]

"""
$(TYPEDSIGNATURES)

List all metabolites without a charge. Use `corner_cases` to specify an extra
case to check for charge's that are not properly assigned.
"""
metabolites_no_charge(model; corner_cases = []) = [mid for mid in metabolites(model) if isnothing(metabolite_charge(model, mid)) || metabolite_charge(model, mid) in corner_cases]

"""
$(TYPEDSIGNATURES)

Test if metabolites `m1` and `m2` are different by comparing their
`test_annotation` fields in the annotations of each metabolite. Note, if no
annotations are present for one or both of the metabolites, then return `true`. 
"""
function metabolites_are_duplicated(model, m1, m2; test_annotation = "inchi_key")
    k1s = get(metabolite_annotations(model, m1), test_annotation, nothing)
    isnothing(k1s) && return true
    k2s = get(metabolite_annotations(model, m2), test_annotation, nothing)
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
function metabolites_unique(model; test_annotation = "inchi_key")
    unique_metabolites = String[]
    for m1 in metabolites(model)
        duplicate = false
        for m2 in unique_metabolites
            duplicate = metabolites_are_duplicated(model, m1, m2; test_annotation)
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
function metabolites_duplicated_in_compartment(model; test_annotation = "inchi_key")
    unique_metabolites = Dict{String, Vector{String}}()
    for m1 in metabolites(model)
        c1 = metabolite_compartment(model, m1)
        for m2 in metabolites(model)
            c2 = metabolite_compartment(model, m2)
            if c1 == c2 && m1 !=  m2 && metabolites_are_duplicated(model, m1, m2; test_annotation) 
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
function test_metabolites(model; duplicates_test_annotation = "inchi_key", formula_corner_cases = ["X", "x", ""], charge_corner_cases=[], medium_only_into=true)
    @testset "Metabolite Information" begin
        @test isempty(metabolites_duplicated_in_compartment(model; test_annotation = duplicates_test_annotation))
        @test isempty(metabolites_no_formula(model; corner_cases = formula_corner_cases))
        @test isempty(metabolites_no_charge(model; corner_cases = charge_corner_cases))
        @test !isempty(metabolites_medium_components(model; only_into=medium_only_into))
    end
end