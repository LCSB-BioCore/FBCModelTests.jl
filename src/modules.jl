"""
ModuleTools

Helper functions used to structure the package.
"""
module ModuleTools # thanks Mirek!

# export everything from the local namespace that seems exportable
# (inspired by JuMP.jl, thanks!)
macro export_locals()
    quote
        for sym in names(@__MODULE__; all = true, imported = true)
            sym in [Symbol(@__MODULE__), :eval, :include] && continue
            startswith(string(sym), ['_', '#']) && continue
            sym == :Internal && continue
            @eval export $(Expr(:$, :sym))
        end
    end
end

@export_locals

end
