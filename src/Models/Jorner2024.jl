using Base: @kwdef
using LinearAlgebra


"""
    Jorner, K., Pollice, R., Lavigne, C., Aspuru-Guzik, A., 2024. 
    Ultrafast Computational Screening of Molecules with Inverted Singlet-Triplet Energy Gaps 
    Using the Pariser-Parr-Pople Semiempirical Quantum Chemistry Method.
    J. Phys. Chem. A 128, 2445-2456. https://doi.org/10.1021/acs.jpca.3c06357

All energies are in eV, distances in Angstroms.

Fields:
- `ohnoconstant`: e²/4πε₀ in eV⋅Å (classical interaction energy of two point charges at 1 Å)
- `UC`: Carbon atom ionisation potential (eV)
- `UN`: Pyrrole nitrogen ionisation potential (eV)
- `UN_AZA`: Aza nitrogen ionisation potential (eV)
- `ES_C`: Carbon electron affinity (eV)
- `ES_NPY`: Pyrrole nitrogen electron affinity (eV)
- `ES_NAZA`: Aza nitrogen electron affinity (eV)
- `Z_EFF`: Effective nuclear charge
- `N_C`: Carbon atom principal quantum number
- `N_NPY`: Pyrrole nitrogen principal quantum number
- `N_NAZA`: Aza nitrogen principal quantum number
"""
@kwdef struct Jorner2024Model <: AbstractModel
    # Fundamental physical constants (computed from SI constants)
    ohnoconstant = q / (4π * ε₀ * 1e-10)  # e²/4πε₀ in eV⋅Å
    cutoff=1.4 # Cut-off for connectivity matrix, in Angstroms
    
    # Model parameters (eV) - using same names as original Constants

    # FIXME: Hubbard U is ionisation potential - electron affinity from Coulson; site energies ripped from Bedogni
    UC = 10.992      # carbon atom ionisation potential - electron affinity (11.16-0.168)
    UN = 12.461    # pyrrole nitrogen ionisation potential - electron affinity (14.12-1.659)
    UN_AZA = 16.754   # aza nitrogen ionisation potential - electron affinity (28.71-11.956)

    ES_C = 0.0      # carbon site energy
    ES_NPY = -13.0  # pyrrole nitrogen site energy
    ES_NAZA = -5.0  # aza nitrogen site energy

    Z_EFF_C = 3.25 # carbon Z effective
    Z_EFF_NPY = 3.90 # pyrrole nitrogen Z effective
    Z_EFF_NAZA = 4.25 # aza nitrogen Z effective

    N_C = 2.0 # carbon principal quantum number
    N_NPY = 2.0 # pyrrole nitrogen principal quantum number
    N_NAZA = 2.0 # aza nitrogen principal quantum number
end

# Implement internal dispatch for Jorner model
function t_ii(model::Jorner2024Model, system, i)
    atom = system.atoms[i]
    if atom.symbol == :C
        return model.ES_C
    elseif atom.symbol == :N
        return atom.n_bonds == 3 ? model.ES_NPY : model.ES_NAZA
    else
        error("Unsupported atom type: $(atom.symbol)")
    end
end

function t_ij(model::Jorner2024Model, system, i, j)
    # Distance-dependent hopping using Slater orbital overlap
    atom_i = system.atoms[i]
    atom_j = system.atoms[j]
    distance = norm(atom_i.position - atom_j.position)
    
    # Get Slater parameters
    Z_eff_i, n_i = get_slater_parameters(model, atom_i)
    Z_eff_j, n_j = get_slater_parameters(model, atom_j)
    
    exp_i = Z_eff_i / n_i
    exp_j = Z_eff_j / n_j
    
    # Calculate overlap gradient and convert to hopping integral
    r_bohr = Angstrom2Bohr * distance 
    overlap_grad = slater_grad(r_bohr, n_i, n_j, exp_i, exp_j, 0.01)
    beta = overlap_grad / r_bohr  
    
    return beta * Ha2eV
end

function γ_ii(model::Jorner2024Model, system, i)
    atom = system.atoms[i]
    if atom.symbol == :C
        return model.UC
    elseif atom.symbol == :N
        return atom.n_bonds == 3 ? model.UN : model.UN_AZA
    else
        error("Unsupported atom type: $(atom.symbol)")
    end
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

# Helper function for Slater parameters
function get_slater_parameters(model::Jorner2024Model, atom)
    if atom.symbol == :C
        return (model.Z_EFF_C, model.N_C)
    elseif atom.symbol == :N
        if atom.n_bonds == 3
            return (model.Z_EFF_NPY, model.N_NPY)
        else
            return (model.Z_EFF_NAZA, model.N_NAZA)
        end
    else
        error("Unsupported atom type: $(atom.symbol)")
    end
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
