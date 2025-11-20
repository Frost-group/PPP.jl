using Base: @kwdef
using LinearAlgebra

# Fundamental physical constants
const q = 1.602176634e-19    # Elementary charge (C)
const ε₀ = 8.8541878128e-12  # Vacuum permittivity (F/m)

# ============================================================================ #
# Data Tables for Bedogni2024 Model
# ============================================================================ #

const BEDOGNI_SITE_ENERGY = Dict{Symbol, Float64}(
    :C => 0.0,
    :N1 => -5.0,   # Aza nitrogen (2 bonds)
    :N2 => -13.0   # Pyrrole nitrogen (3 bonds)
)

const BEDOGNI_HUBBARD_U = Dict{Symbol, Float64}(
    :C => 11.26,
    :N1 => 15.5,   # Aza nitrogen
    :N2 => 15.0    # Pyrrole nitrogen
)

const BEDOGNI_Z_EFF = Dict{Symbol, Float64}(
    :C => 3.25,
    :N1 => 4.25,   # Aza nitrogen
    :N2 => 3.90    # Pyrrole nitrogen
)

const BEDOGNI_N_PRINCIPAL = Dict{Symbol, Int}(
    :C => 2,
    :N1 => 2,
    :N2 => 2
)

const BEDOGNI_PI_ELECTRONS = Dict{Symbol, Int}(
    :C => 1,
    :N1 => 1,  # Aza nitrogen
    :N2 => 2   # Pyrrole nitrogen
)

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
- `SITE_ENERGY`: Dictionary of site energies (eV)
- `HUBBARD_U`: Dictionary of Hubbard U parameters (eV)
- `Z_EFF`: Dictionary of Effective Nuclear Charges
- `N_PRINCIPAL`: Dictionary of Principal Quantum Numbers
- `PI_ELECTRONS`: Dictionary of pi electron counts
"""
@kwdef struct Bedogni2024Model <: AbstractModel
    # Fundamental physical constant (directly from SI constants)
    ohnoconstant::Float64 = q / (4π * ε₀ * 1e-10)  # e²/4πε₀ in eV⋅Å
    cutoff::Float64 = 1.4 # Cut-off for connectivity matrix, in Angstroms
    
    # Model parameters (eV) - using dictionaries
    T::Float64 = -2.4  # hopping integral

    SITE_ENERGY::Dict{Symbol, Float64} = BEDOGNI_SITE_ENERGY
    HUBBARD_U::Dict{Symbol, Float64} = BEDOGNI_HUBBARD_U
    Z_EFF::Dict{Symbol, Float64} = BEDOGNI_Z_EFF
    N_PRINCIPAL::Dict{Symbol, Int} = BEDOGNI_N_PRINCIPAL
    PI_ELECTRONS::Dict{Symbol, Int} = BEDOGNI_PI_ELECTRONS
end

# Implement internal dispatch for Bedogni model
function t_ii(model::Bedogni2024Model, system, i)
    atom = system.atoms[i]
    atomtype = get_atomtype(atom)
    
    if !haskey(model.SITE_ENERGY, atomtype)
        error("Unsupported atom type in Bedogni2024Model: $atomtype (symbol: $(atom.symbol))")
    end
    
    return model.SITE_ENERGY[atomtype]
end

function t_ij(model::Bedogni2024Model, system, i, j)
    # Fixed hopping integral for Bedogni model
    return model.T
end

function γ_ii(model::Bedogni2024Model, system, i)
    atom = system.atoms[i]
    atomtype = get_atomtype(atom)
    
    if !haskey(model.HUBBARD_U, atomtype)
        error("Unsupported atom type in Bedogni2024Model: $atomtype (symbol: $(atom.symbol))")
    end
    
    return model.HUBBARD_U[atomtype]
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

# ============================================================================ #
# Interface Implementations
# ============================================================================ #

function get_pi_electrons(model::Bedogni2024Model, atom)
    atomtype = get_atomtype(atom)
    if !haskey(model.PI_ELECTRONS, atomtype)
        error("Unknown atom type for PI electrons: $atomtype")
    end
    return model.PI_ELECTRONS[atomtype]
end

function update_atom_params(model::Bedogni2024Model, atom)
    atomtype = get_atomtype(atom)
    
    # Update nz (pi electrons)
    nz = get_pi_electrons(model, atom)
    
    Z_eff = get(model.Z_EFF, atomtype, 0.0)
    n_int = get(model.N_PRINCIPAL, atomtype, 2)
    
    return Atom(atom.symbol, atom.position, atom.n_bonds, nz, atom.site_energy, atom.Hubbard_U, Z_eff, Float64(n_int))
end
