#=
This file contains a collection of tests based on Memote. See Lieven, C., Beber,
M.E., Olivier, B.G. et al. MEMOTE for standardized genome-scale metabolic model
testing. Nat Biotechnol 38, 272–276 (2020).
https://doi.org/10.1038/s41587-020-0446-y for details.
=#

"""
$(TYPEDSIGNATURES)

Check if model has an ATP maintenance reaction built in (also called a
non-growth associated maintenance cost). Looks for reaction annotations
corresponding to the sbo maintenance term, or looks for reaction ids that
contain the strings listed in `config.biomass.atpm_strings`.
"""
model_has_atpm_reaction(model; config = memote_config) = any(is_atp_maintenance_reaction(model, rid) for rid in reactions(model)) || any(any(occursin.(config.biomass.atpm_strings, Ref(rid))) for rid in reactions(model))

"""
$(TYPEDSIGNATURES)

Check if the biomass reaction consumes ATP and H₂O, and produces ADP, HO₄P, and
H⁺. Each of these metabolites have a lookup table mapping them to the name space
of the model, defined in `config.biomass.growth_metabolites`. These need to be
set if you use anything other than the BiGG namespace.
"""
function model_has_growth_atp_in_biomass(model; config = memote_config)
    
end