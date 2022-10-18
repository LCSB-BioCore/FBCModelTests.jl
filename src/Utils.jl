"""
# Utils

Miscellaneous functions, not typically user facing.
"""
module Utils

using ..DocStringExtensions
using ..COBREXA
using ..Test
using PeriodicTable

const Maybe{T} = Union{Nothing,T}

gets(c::Dict, def, k, ks...) = haskey(c, k) ? gets(c[k], def, ks...) : def
gets(c::Tuple, def, k, ks...) = k in 1:length(c) ? gets(c[k], def, ks...) : def
gets(c, _) = c

test_dicts(match::Function, a::Dict, b::Dict) =
    for k in union(keys(a), keys(b))
        @testset "$k" begin
            @test haskey(a, k) && haskey(b, k)
            if haskey(a, k) && haskey(b, k)
                match(k, a[k], b[k])
            end
        end
    end

"""
$(TYPEDSIGNATURES)

Internal helper function that determines if a reaction has a gene reaction rule
and that each gene in the rule is contained in the model.
"""
function _has_sensible_gpr(model, rid)
    grrs = reaction_gene_association(model, rid)
    isnothing(grrs) && return false
    isempty(grrs) && return false
    any(isempty.(grrs)) && return false
    any("" in grr for grr in grrs) && return false

    gids = Set(reduce(vcat, grrs))
    all(in.(gids, Ref(genes(model))))
end

"""
$(TYPEDSIGNATURES)

Determine if a reaction is probably a transport reaction by checking if:
1. it has sbo annotations corresponding to a transport reaction
2. the reaction contains metabolites from at least 2 different compartments
3. if at least 1 metabolite does not undergo a chemical transformation (via
   formula or annotation checks)

Note, PTS type transport reactions will be missed if they do not have sbo
annotations. This test may yield false negatives.
"""
function _probably_transport_reaction(model, rid, test_annotation)
    is_transport_reaction(model, rid) && return true
    allequal(x) = all(isequal(first(x), x)) #  TODO remove when Julia LTS is > v1.8
    allequal(
        metabolite_compartment(model, mid) for
        mid in keys(reaction_stoichiometry(model, rid))
    ) && return false

    comp_mid = Dict{String,Set{String}}()
    for mid in keys(reaction_stoichiometry(model, rid))
        push!(get!(comp_mid, metabolite_compartment(model, mid), Set{String}()), mid)
    end

    # must have annotations for all metabolites
    any(
        !haskey(metabolite_annotations(model, mid), test_annotation) for
        mid in keys(reaction_stoichiometry(model, rid))
    ) && return false
    # must have formula for all metabolites
    any(
        isnothing(metabolite_formula(model, mid)) for
        mid in keys(reaction_stoichiometry(model, rid))
    ) && return false

    get_annotation(mid) = metabolite_annotations(model, mid)[test_annotation]
    get_formula(mid) = begin
        d = metabolite_formula(model, mid)
        ks = sort(collect(keys(d)))
        [join(k * string(d[k]) for k in ks)]
    end

    #=
    Compare the formulas and annotations of metabolites in different
    compartments. If the any metabolite has the same formula or same annotation
    but occurs in different comparments, then assume that it is probably a
    transported.
    =#
    for (k1, v1) in comp_mid
        for (k2, v2) in comp_mid
            k1 == k2 && continue
            any(
                in.(
                    reduce(vcat, get_formula(x) for x in v1),
                    Ref(reduce(vcat, get_formula(x) for x in v2)),
                ),
            ) && return true
            any(
                in.(
                    reduce(vcat, get_annotation(x) for x in v1),
                    Ref(reduce(vcat, get_annotation(x) for x in v2)),
                ),
            ) && return true
        end
    end

    return false
end

"""
$(TYPEDSIGNATURES)

Return the chemical element of `x`.
"""
to_element(x::String) = begin
    sym =
        length(x) > 1 ? Symbol(uppercase(first(x)) * x[2:end]) : Symbol(uppercase(first(x)))
    elements[sym]
end

end # module