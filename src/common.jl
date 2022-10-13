
const Maybe{T} = Union{Nothing,T}

gets(c::Dict, def, k, ks...) = haskey(c, k) ? gets(c[k], def, ks...) : def
gets(c::Tuple, def, k, ks...) = k in 1:length(c) ? gets(c[k], def, ks...) : def
gets(c, _) = c
