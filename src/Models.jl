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

# Functions for getting pi electrons and updating atom parameters
function get_pi_electrons(model::AbstractModel, atom)
    return 1 # very simple Huckel style models assume one pi electron per atom
end

function update_atom_params(model::AbstractModel, atom)
    return atom 
end

# Utility function to determine specific atom type (e.g., :N -> :N1 or :N2)
function get_atomtype(atom)
    if atom.symbol == :N
        # N1: Pyridine-like (2 bonds, 1 pi electron)
        # N2: Pyrrole-like (3 bonds, 2 pi electrons)
        return atom.n_bonds == 3 ? :N2 : :N1
    elseif atom.symbol == :O
        # O1: Carbonyl-like (1 bond?, 1 pi electron)
        # O2: Ether/Furan-like (2 bonds, 2 pi electrons)
        return atom.n_bonds == 2 ? :O2 : :O1
    elseif atom.symbol == :S
        # S1: Thio-carbonyl?
        # S2: Thiophene-like
        return atom.n_bonds == 2 ? :S2 : :S1
    elseif atom.symbol == :P
        # P1: Phosphinine?
        # P2: Phosphole?
        return atom.n_bonds == 3 ? :P2 : :P1 
    else
        return atom.symbol
    end
end

# Include Models
include("Models/Bedogni2024.jl")
include("Models/Jorner2024.jl")
include("Models/Zhang2011.jl")
