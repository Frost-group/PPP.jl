using Base: @kwdef
using LinearAlgebra

# ============================================================================ #
# Data Tables (from coulson/parameters.py 'CRC' / 'MODERN')
# ============================================================================ #

const IP_CRC = Dict{Symbol, Float64}(
    :B => 1.011,
    :C => 11.164,
    :N1 => 14.093, # Pyridine-like (1 e)
    :N2 => 28.717, # Pyrrole-like (2 e)
    :O1 => 17.701,
    :O2 => 34.122,
    :F => 40.697,
    :Si => 9.177,
    :P1 => 11.146,
    :P2 => 20.732,
    :S1 => 12.706,
    :S2 => 23.740,
    :Cl => 27.296
)

const EA_CRC = Dict{Symbol, Float64}(
    :B => 0.115,
    :C => 0.168,
    :N1 => 1.659,
    :N2 => 11.956,
    :O1 => 2.456,
    :O2 => 15.305,
    :F => 18.519,
    :Si => 1.925,
    :P1 => 1.776,
    :P2 => 10.266,
    :S1 => 2.764,
    :S2 => 11.648,
    :Cl => 14.501
)

const Z_EFF_SLATER = Dict{Symbol, Float64}(
    :B => 2.25,
    :C => 3.25,
    :N1 => 3.90,
    :N2 => 4.25,
    :O1 => 4.55,
    :O2 => 4.9,
    :F => 5.55,
    :Si => 4.15,
    :P1 => 4.80,
    :P2 => 5.15,
    :S1 => 5.45,
    :S2 => 5.80,
    :Cl => 6.45,
    :Br => 7.95,
    :I => 7.95
)

const N_PRINCIPAL = Dict{Symbol, Int}(
    :B => 2,
    :C => 2,
    :N1 => 2,
    :N2 => 2,
    :O1 => 2,
    :O2 => 2,
    :F => 2,
    :Si => 2,
    :P1 => 3,
    :P2 => 3,
    :S1 => 3,
    :S2 => 3,
    :Cl => 3,
    :Br => 4,
    :I => 5
)

# N_STAR map from n (int) to n_star (float)
const N_STAR_MAP = Dict{Int, Float64}(
    1 => 1.0,
    2 => 2.0,
    3 => 3.0,
    4 => 3.7,
    5 => 4.0,
    6 => 4.2
)

const PI_ELECTRONS = Dict{Symbol, Int}(
    :B => 0,
    :C => 1,
    :N1 => 1,
    :N2 => 2,
    :O1 => 1,
    :O2 => 2,
    :F => 2,
    :Si => 1,
    :P1 => 1,
    :P2 => 2,
    :S1 => 1,
    :S2 => 2,
    :Cl => 2,
    :Br => 2, 
    :I => 2   
)

"""
    Jorner, K., Pollice, R., Lavigne, C., Aspuru-Guzik, A., 2024. 
    Ultrafast Computational Screening of Molecules with Inverted Singlet-Triplet Energy Gaps 
    Using the Pariser-Parr-Pople Semiempirical Quantum Chemistry Method.
    J. Phys. Chem. A 128, 2445-2456. https://doi.org/10.1021/acs.jpca.3c06357

All energies are in eV, distances in Angstroms.

Fields:
- `ohnoconstant`: e²/4πε₀ in eV⋅Å (classical interaction energy of two point charges at 1 Å)
- `IP`: Dictionary of Ionization Potentials (eV)
- `EA`: Dictionary of Electron Affinities (eV)
- `Z_EFF`: Dictionary of Effective Nuclear Charges
- `N_PRINCIPAL`: Dictionary of Principal Quantum Numbers
- `N_STAR_MAP`: Dictionary mapping n to n*
"""
@kwdef struct Jorner2024Model <: AbstractModel
    # Fundamental physical constants (computed from SI constants)
    ohnoconstant::Float64 = q / (4π * ε₀ * 1e-10)  # e²/4πε₀ in eV⋅Å
    cutoff::Float64 = 1.4 # Cut-off for connectivity matrix, in Angstroms
    
    # Model parameters (eV) - using dictionaries
    IP::Dict{Symbol, Float64} = IP_CRC
    EA::Dict{Symbol, Float64} = EA_CRC
    Z_EFF::Dict{Symbol, Float64} = Z_EFF_SLATER
    N_PRINCIPAL::Dict{Symbol, Int} = N_PRINCIPAL
    N_STAR_MAP::Dict{Int, Float64} = N_STAR_MAP
    PI_ELECTRONS::Dict{Symbol, Int} = PI_ELECTRONS
end

# Implement internal dispatch for Jorner model
function t_ii(model::Jorner2024Model, system, i)
    atom = system.atoms[i]
    type = get_atomtype(atom)
    
    if !haskey(model.IP, type)
        error("Unsupported atom type in Jorner2024Model: $type (symbol: $(atom.symbol))")
    end
    
    # Site energy alpha ~= -IP
    return -model.IP[type]
end

function t_ij(model::Jorner2024Model, system, i, j)
    # Distance-dependent hopping using Slater orbital overlap
    atom_i = system.atoms[i]
    atom_j = system.atoms[j]
    distance = norm(atom_i.position - atom_j.position)
    
    # Get Slater parameters
    type_i = get_atomtype(atom_i)
    type_j = get_atomtype(atom_j)
    
    if !haskey(model.Z_EFF, type_i) || !haskey(model.N_PRINCIPAL, type_i)
        error("Missing Slater parameters for atom type: $type_i")
    end
    if !haskey(model.Z_EFF, type_j) || !haskey(model.N_PRINCIPAL, type_j)
        error("Missing Slater parameters for atom type: $type_j")
    end

    Z_eff_i = model.Z_EFF[type_i]
    n_i_int = model.N_PRINCIPAL[type_i]
    n_star_i = model.N_STAR_MAP[n_i_int]
    
    Z_eff_j = model.Z_EFF[type_j]
    n_j_int = model.N_PRINCIPAL[type_j]
    n_star_j = model.N_STAR_MAP[n_j_int]
    
    exp_i = Z_eff_i / n_star_i
    exp_j = Z_eff_j / n_star_j
    
    # Calculate overlap gradient and convert to hopping integral
    r_bohr = Angstrom2Bohr * distance 
    # Convert integer n to float for the overlap function which expects floats
    overlap_grad = slater_grad(r_bohr, Float64(n_i_int), Float64(n_j_int), exp_i, exp_j, 0.01)
    beta = overlap_grad / r_bohr  
    
    return beta * Ha2eV
end

function γ_ii(model::Jorner2024Model, system, i)
    atom = system.atoms[i]
    type = get_atomtype(atom)
    
    if !haskey(model.IP, type) || !haskey(model.EA, type)
        error("Missing IP/EA parameters for atom type: $type")
    end
    
    # Gamma = IP - EA
    return model.IP[type] - model.EA[type]
end

function γ_ij(model::Jorner2024Model, system, i, j)
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

function get_pi_electrons(model::Jorner2024Model, atom)
    type = get_atomtype(atom)
    if !haskey(model.PI_ELECTRONS, type)
        error("Unknown atom type for PI electrons: $type")
    end
    return model.PI_ELECTRONS[type]
end

function update_atom_params(model::Jorner2024Model, atom)
    type = get_atomtype(atom)
    
    # Update nz (pi electrons)
    nz = get_pi_electrons(model, atom)
    
    Z_eff = get(model.Z_EFF, type, 0.0)
    n_int = get(model.N_PRINCIPAL, type, 2)
    
    return Atom(atom.symbol, atom.position, atom.n_bonds, nz, atom.site_energy, atom.Hubbard_U, Z_eff, Float64(n_int))
end

# ============================================================================ #
# Slater orbital functions (from Coulson code)
# ============================================================================ #
"""
Helper functions A & B for Slater orbital overlap
"""
function A(k::Int, p::Float64)
    summed = 0.0
    for i in 1:(k+1)
        summed += factorial(k) / (p^i*factorial(k-i+1))
    end
    return exp(-p)*summed
end

function B(k::Int, pt::Float64)
    summed_1  = 0.0
    summed_2 = 0.0
    for i in 1:(k+1)
        term = factorial(k) / ((pt)^i *factorial(k-i+1))
        summed_1 += term
        summed_2 += (-1)^(k-i)*term
    end
    return -exp(-pt)*summed_1 - exp(pt)*summed_2
end

"""
Sort orbitals to use with Mulliken's STO orbital overlap formulas:  https://doi.org/10.1063/1.1747150
"""
function sort_orbitals(n_1::Float64, n_2::Float64, exp_1::Float64, exp_2::Float64)
    if n_1 > n_2
        return n_2, n_1, exp_2, exp_1
    elseif (n_1==n_2) && (exp_1<exp_2)
        return n_1, n_2, exp_2, exp_1
    else
        return n_1, n_2, exp_1, exp_2
    end
end

"""
Calculates STO overlap based on Mulliken formulas:  https://doi.org/10.1063/1.1747150
"""
function calculate_overlap(n_1::Float64, n_2::Float64, p::Float64, t::Float64)
    overlap = 0.0
    if p == 0.0
        if n_1 == n_2 == 2
            overlap = (1-t^2)^(5/2)
        elseif n_1 == n_2 == 3
            overlap = (1-t^2)^(7/2)
        elseif (n_1 ==2) && (n_2 ==3)
            overlap = ((5/6)*(1+t)^5*(1-t)^7)^(1/2)
        end
    end

    if t == 0.0 && (n_1 == n_2)
        # 2p-2p
        if n_1 == n_2 == 2
            overlap = exp(-p)*(1+p+(2/5)*p^2+(1/15)*p^3)
        end
        # 3p-3p
        if n_1 == n_2 ==3
            overlap = exp(-p)*(1+p+(34/75)*p^2+(3/25)*p^3+(31/1575)*p^4+(1/525)*p^5)
        end
    else
        pt = p*t
        # 2p-2p
        if n_1 == n_2 == 2
            overlap = ((1 / 32)* p^5*(1 - t^2)^(5/2)*(A(4, p) * (B(0, pt) - B(2, pt))+ A(2, p) * (B(4, pt) - B(0, pt))+ A(0, p) * (B(2, pt) - B(4, pt))))
        # 2p-3p
        elseif (n_1 == 2) && (n_2 == 3)
            if t == 0.0
                overlap = (1/(120*sqrt(30))*p^6*(5 * A(5, p) - 6 * A(3, p) + A(1, p)))
            else
                overlap = (1/(32*sqrt(30)*p^6*(1+t)^(5/2)*(1-t)^(7/2)*(A(5, p) * (B(0, p * t) - B(2, pt))+ A(4, p) * (B(3, pt) - B(1, pt))+ A(3, p) * (B(4, pt) - B(0, pt))
                + A(2, p) * (B(1, pt) - B(5, pt))+ A(1, p) * (B(2, pt) - B(4, pt)) + A(0, p) * (B(5, pt) - B(3, pt)))))
            end
        # 3p-3p
        elseif n_1 == n_2 == 3
            overlap = (1/960*p^7*(1-t^2)^(7/2)*(A(6, p) * (B(0, pt) - B(2, pt))+ A(4, p) * (2 * B(4, pt) - B(0, pt) - B(2, pt))+ A(2, p) * (2 * B(2, pt) - B(4, pt) - B(6, pt))
                    + A(0, p) * (B(6, pt) - B(4, pt))))
        end
    end
    return overlap
end

"""
Returns Slater overlap between two Slater p orbitals
"""
function slater_overlap(r::Float64, n_1::Float64, n_2::Float64, exp_1::Float64, exp_2::Float64)
    # Sort orbitals to match Mulliken formulas
    n_1, n_2, exp_1, exp_2 = sort_orbitals(n_1,n_2,exp_1,exp_2)
    # Calculate p and t parameters for Mulliken formulas: https://doi.org/10.1063/1.1747150
    p = r*(exp_1+exp_2)/2
    t = (exp_1-exp_2)/(exp_1+exp_2)
    # Calculate overlap
    overlap = calculate_overlap(n_1, n_2, p, t)
    return overlap
end

"""
Calculates gradient of orbital overlap integral by numerical differentiation
"""
function slater_grad(r::Float64, n_1::Float64, n_2::Float64, exp_1::Float64, exp_2::Float64, dx::Float64)
    f = x -> slater_overlap(x, n_1, n_2, exp_1, exp_2)
    grad = (f(r + dx) - f(r - dx)) / (2 * dx)
    return grad
end
