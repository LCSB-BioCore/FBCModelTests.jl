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
