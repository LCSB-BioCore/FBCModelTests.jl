"""
    module  Utils

Miscellaneous functions used by memote style tests, not typically user facing.
"""
module Utils

using COBREXA
using DocStringExtensions
using PeriodicTable
using Statistics

"""
$(TYPEDSIGNATURES)

Return the chemical element of `x`.
"""
to_element(x::String) = begin
    sym =
        length(x) > 1 ? Symbol(uppercase(first(x)) * x[2:end]) : Symbol(uppercase(first(x)))
    elements[sym]
end

"""
$(TYPEDSIGNATURES)

Return molar mass of `mid`. Return `NaN` if formula does not exist.
"""
get_molar_mass(model, mid) = begin
    rs = metabolite_formula(model, mid)
    isnothing(rs) && return NaN # if metabolite has no molar mass
    sum(v * to_element(k).atomic_mass for (k, v) in rs).val
end

"""
$(TYPEDSIGNATURES)

A helper function that returns the median upper and lower bounds in a tuple. If none can be calculated,
constants from COBREXA are used as the default values.
"""
function median_bounds(model::MetabolicModel)
    default = COBREXA._constants.default_reaction_bound
    lower_bound, upper_bound = bounds(model)
    lb_list = [element for element in lower_bound if !isapprox(element, 0.0)]
    ub_list = [element for element in upper_bound if !isapprox(element, 0.0)]
    isempty(lb_list) ? m_lower_bound = -default : m_lower_bound = median(lb_list)
    isempty(ub_list) ? m_upper_bound = default : m_upper_bound = median(ub_list)
    return m_lower_bound, m_upper_bound
end

"""
$(TYPEDSIGNATURES)

Internal helper function to compare fluxes with specific bounds.
"""
function _compare_flux_bounds(fluxes, bound, tol, comparison_operator)
    unlimited_flux = Dict{String,Tuple{String,Float64}}()
    for (rid, d) in fluxes, (frid, flux) in d
        if comparison_operator(flux, bound) || isapprox(flux, bound; atol = tol)
            unlimited_flux[rid] = (frid, flux)
        end
    end
    return unlimited_flux
end

"""
$(TYPEDSIGNATURES)

Ensure annotations are in a standardized format. Some models represent
annotations like:
```
Dict{String, Vector{String}} with 2 entries:
  "sbo"          => ["SBO:0000176"]
  "RESOURCE_URI" => ["https://identifiers.org/ec-code/4.1.99.12", "https://identifiers.org/bigg.reaction/DB4PS", ...]
```
but the key used to map to URI annotations is not standarized. This helper
function constructs a new annotation dictionary to ensure that all annotations
looks like:
```
Dict{String, Vector{String}} with 7 entries:
  "bigg.reaction"     => ["DB4PS"]
  "pubmed"            => ["12595523"]
  "sbo"               => ["SBO:0000176"]
  "kegg.pathway"      => ["sce00740", "sce01110"]
  "metanetx.reaction" => ["MNXR97178"]
  "kegg.reaction"     => ["R07281"]
  "ec-code"           => ["4.1.99.12"]
```
"""
function parse_annotations(annos)
    parsed_annos = Dict{String,Vector{String}}()
    for (k, vs) in annos
        for v in vs
            if contains(v, "https://identifiers.org")
                s = split(v, "/")
                new_k = s[end-1]
                new_v = last(s)
            else
                new_k = k
                new_v = v
            end
            push!(get!(parsed_annos, new_k, String[]), new_v)
        end
    end
    return parsed_annos
end

end # module
