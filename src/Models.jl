"""
    AbstractModel

Base type for all PPP models.
"""
abstract type AbstractModel end

"""
    t(model::AbstractModel, system, i, j) -> Float64

Calculate the hopping matrix element between sites i and j.
Dispatches to t_ii for diagonal elements and t_ij for off-diagonal elements.
"""
function t(model::AbstractModel, system, i, j)
    if i == j
        return t_ii(model, system, i)
    else
        return t_ij(model, system, i, j)
    end
end

"""
    γ(model::AbstractModel, system, i, j) -> Float64

Calculate the Ohno parameter between sites i and j.
Dispatches to γ_ii for diagonal elements and γ_ij for off-diagonal elements.
"""
function γ(model::AbstractModel, system, i, j)
    if i == j
        return γ_ii(model, system, i)
    else
        return γ_ij(model, system, i, j)
    end
end

# Internal dispatch functions - these will be implemented by each model
function t_ii(model::AbstractModel, system, i)
    warn("t_ii not implemented for $(typeof(model)). Returning 0.0.")
    return 0.0
end

function t_ij(model::AbstractModel, system, i, j)
    warn("t_ij not implemented for $(typeof(model)). Returning 0.0.")
    return 0.0
end

function γ_ii(model::AbstractModel, system, i)
    warn("γ_ii not implemented for $(typeof(model)). Returning 0.0.")
    return 0.0
end

function γ_ij(model::AbstractModel, system, i, j)
    warn("γ_ij not implemented for $(typeof(model)). Returning 0.0.")
    return 0.0
end

# Include Models
include("Models/Bedogni2024.jl")
include("Models/Jorner2024.jl")
include("Models/Zhang2011.jl")
