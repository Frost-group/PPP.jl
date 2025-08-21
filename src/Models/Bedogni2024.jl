using Base: @kwdef
using LinearAlgebra

# Fundamental physical constants
const q = 1.602176634e-19    # Elementary charge (C)
const ε₀ = 8.8541878128e-12  # Vacuum permittivity (F/m)

"""
    Bedogni2024Model <: AbstractModel

    Simple Carbon / Nitrogen PPP model, bare Ohno exchange.
    Bedogni, M., Giavazzi, D., Di Maiolo, F., Painelli, A., 2024. 
    Shining Light on Inverted Singlet-Triplet Emitters. 
    J. Chem. Theory Comput. 20, 902-913. https://doi.org/10.1021/acs.jctc.3c01112

All energies are in eV, distances in Angstroms.

Fields:
- `ohnoconstant`: e²/4πε₀ in eV⋅Å (classical interaction energy of two point charges at 1 Å)
- `T`: Hopping integral (eV)
- `UC`: Carbon atom Hubbard U (eV)
- `UN`: Pyrrole nitrogen Hubbard U (eV)
- `UN_AZA`: Aza nitrogen Hubbard U (eV)
- `ES_C`: Carbon site energy (eV)
- `ES_NPY`: Pyrrole nitrogen site energy (eV)
- `ES_NAZA`: Aza nitrogen site energy (eV)
"""
@kwdef struct Bedogni2024Model <: AbstractModel
    # Fundamental physical constant (directly from SI constants)
    ohnoconstant = q / (4π * ε₀ * 1e-10)  # e²/4πε₀ in eV⋅Å
    
    # Model parameters (eV) - using same names as original Constants
    T = -2.4        # hopping integral

    UC = 11.26      # carbon atom Hubbard U
    UN = 15.0       # pyrrole nitrogen Hubbard U
    UN_AZA = 15.5   # aza nitrogen Hubbard U

    ES_C = 0.0      # carbon site energy
    ES_NPY = -13.0  # pyrrole nitrogen site energy
    ES_NAZA = -5.0  # aza nitrogen site energy
    
    # Slater orbital parameters (needed for atom creation)
    Z_EFF_C = 3.25 # carbon Z effective
    Z_EFF_NPY = 3.90 # pyrrole nitrogen Z effective
    Z_EFF_NAZA = 4.25 # aza nitrogen Z effective

    N_C = 2.0 # carbon principal quantum number
    N_NPY = 2.0 # pyrrole nitrogen principal quantum number
    N_NAZA = 2.0 # aza nitrogen principal quantum number
end

# Implement internal dispatch for Bedogni model
function t_ii(model::Bedogni2024Model, system, i)
    atom = system.atoms[i]
    if atom.symbol == :C
        return model.ES_C
    elseif atom.symbol == :N
        return atom.n_bonds == 3 ? model.ES_NPY : model.ES_NAZA
    else
        error("Unsupported atom type: $(atom.symbol)")
    end
end

function t_ij(model::Bedogni2024Model, system, i, j)
    # Fixed hopping integral for Bedogni model
    return model.T
end

function γ_ii(model::Bedogni2024Model, system, i)
    atom = system.atoms[i]
    if atom.symbol == :C
        return model.UC
    elseif atom.symbol == :N
        return atom.n_bonds == 3 ? model.UN : model.UN_AZA
    else
        error("Unsupported atom type: $(atom.symbol)")
    end
end

function γ_ij(model::Bedogni2024Model, system, i, j)
    # Ohno formula for off-diagonal elements
    atom_i = system.atoms[i]
    atom_j = system.atoms[j]
    r_ij = norm(atom_i.position - atom_j.position)
    U_i = γ_ii(model, system, i)  # Get diagonal element
    U_j = γ_ii(model, system, j)  # Get diagonal element
    return model.ohnoconstant / sqrt(r_ij^2 + (2*model.ohnoconstant/(U_i + U_j))^2)
end
