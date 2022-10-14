
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
    any("" in grr for grr in grrs) && return false

    gids = Set(reduce(vcat, grrs))
    all(in.(gids, Ref(genes(model))))
end

"""
$(TYPEDSIGNATURES)

Determine if a reaction is probably a transport reaction by checking if:
1. the reaction contains metabolites from at least 2 different compartments
2. if at least 1 metabolite does not undergo a chemical transformation
"""
function _resembles_transport_reaction(model, rid)
    
end