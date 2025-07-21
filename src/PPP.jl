module PPP

using LinearAlgebra
using StaticArrays
using Printf


# ============================================================================ #
# Constants
# ============================================================================ #
"""
Physical and model parameters for PPP calculations.
All energies are in eV, distances in Angstroms.
"""
module Constants
    # Fundamental physical constants
    const q = 1.602176634e-19    # Elementary charge (C)
    const ε₀ = 8.8541878128e-12  # Vacuum permittivity (F/m)
    # Derived constant for PPP model (e²/4πε₀ in eV⋅Å)
    const ohnoconstant = q²_to_ev = q / (4π * ε₀ * 1e-10)  # 1 Å = 1e-10 m
    #  i.e. classical interaction energy of two point charges at 1 Å, in eV
  #  @assert e²_to_ev ≈ 14.397 "Calculated e²/Å differs from expected value"

    # Model parameters (eV)
    const T = -2.4        # hopping integral

    const UC = 11.26      # carbon atom Hubbard U
    const UN = 15.0       # pyrrole nitrogen Hubbard U
    const UN_AZA = 15.5   # aza nitrogen Hubbard U
    
    const ES_C = 0.0      # carbon site energy
    const ES_NPY = -13.0  # pyrrole nitrogen site energy
    const ES_NAZA = -5.0  # aza nitrogen site energy
    const cutoff = 1.4    # bond length cutoff (Å), used to definbe connectivity
end

# ============================================================================ #
# Types
# ============================================================================ #
"""
    Atom

Represents an atom in the π-conjugated system.

Fields:
- `symbol::Symbol`: Atomic symbol (e.g., :C, :N)
- `position::SVector{2,Float64}`: 2D coordinates (x,y) in Angstroms
- `n_bonds::Int`: Number of bonds in π structure
- `nz::Int`: Number of valence electrons; terrible hack for Nitrogen
- `site_energy::Float64`: Site energy in eV
- `Hubbard_U::Float64`: Hubbard U parameter in eV
"""
struct Atom
    symbol::Symbol
    position::SVector{2,Float64}
    n_bonds::Int
    nz::Int
    site_energy::Float64
    Hubbard_U::Float64
end

"""
    MolecularSystem

Represents a molecular system for PPP calculations.

Fields:
- `atoms::Vector{Atom}`: List of atoms in the system
- `connectivity::Matrix{Int}`: Connectivity matrix (0/1)
- `n_electrons::Int`: Number of π electrons
"""
struct MolecularSystem
    atoms::Vector{Atom}
    connectivity::Matrix{Int}
    n_electrons::Int
end

"""
    HuckelResult
"""
struct HuckelResult
    hamiltonian::Matrix{Float64}
    energies::Vector{Float64}
    eigenvectors::Matrix{Float64}
    n_occupied::Int
    density_matrix::Matrix{Float64}
    total_energy::Float64
end

"""
    SCFResult
"""
struct SCFResult
    F::Matrix{Float64}
    K::Matrix{Float64}
    J::Matrix{Float64}
    V_raw::Matrix{Float64}
    energies::Vector{Float64}
    eigenvectors::Matrix{Float64}
    n_occupied::Int
    density_matrix::Matrix{Float64}
    total_energy::Float64
    iterations::Int
    converged::Bool
end

# ============================================================================ #
# Geometry and System Setup
# ============================================================================ #

function create_atom(symbol::Symbol, x::Float64, y::Float64)
    position = SVector{2,Float64}(x, y)
   # treatment of nz factor is an absolute mess, because I didn't understand what it was as I was coding along
   # FIXME: The horror that remains works, but only for C and N, and only in the one specific geometry tested 
    if symbol == :C
        return Atom(symbol, position, 2, 1, Constants.ES_C, Constants.UC)
    elseif symbol == :N
        # n_bonds will be updated later based on connectivity
        return Atom(symbol, position, 0, 0,Constants.ES_NPY, Constants.UN)
    else
        throw(ArgumentError("Unsupported atomic symbol: $symbol"))
    end
end

"""
    parse_xyz_line(line::String) -> Tuple{Symbol,Float64,Float64}

Parse a line from XYZ file format, returning atomic symbol and coordinates.
"""
function parse_xyz_line(line::String)
    parts = split(strip(line))
    symbol = Symbol(parts[1])
    x = parse(Float64, parts[2])
    y = parse(Float64, parts[3])
    return symbol, x, y
end

"""
    calculate_connectivity(positions::Vector{SVector{2,Float64}})::Matrix{Int}

Calculate the connectivity matrix based on atomic positions.
Returns a symmetric matrix where 1 indicates bonded atoms and 0 indicates non-bonded atoms.
"""
function calculate_connectivity(positions::Vector{SVector{2,Float64}})::Matrix{Int}
    n_atoms = length(positions)
    connectivity = zeros(Int, n_atoms, n_atoms)
    
    for i in 1:n_atoms
        for j in (i+1):n_atoms
            distance = norm(positions[i] - positions[j])
            if distance < Constants.cutoff
                connectivity[i,j] = connectivity[j,i] = 1
            end
        end
    end
    
    return connectivity
end

"""
    count_bonds(connectivity::Matrix{Int}, idx::Int)::Int

Count the number of bonds for a given atom index.
"""
function count_bonds(connectivity::Matrix{Int}, idx::Int)::Int
    sum(connectivity[idx, :])
end

"""
    calculate_n_electrons(atoms::Vector{Atom})::Int

Calculate the total number of π electrons in the system.
"""
function calculate_n_electrons(atoms::Vector{Atom})::Int
    # Absolute terrible hack; only works for C and N
    # FIXME: Some kind of periodic lookup table, with valency calculation c.f. group ?
    sum(atom.symbol == :C ? 1 : (atom.n_bonds == 3 ? 2 : 1) for atom in atoms)
end

# ============================================================================ #
# Hückel Calculations
# ============================================================================ #

"""
    calculate_density_matrix(eigenvectors::Matrix{Float64}, n_occupied::Int)::Matrix{Float64}

Calculate the density matrix from molecular orbital coefficients.
"""
function calculate_density_matrix(eigenvectors::Matrix{Float64}, n_occupied::Int)::Matrix{Float64}
    n_sites = size(eigenvectors, 1)
    P = zeros(Float64, n_sites, n_sites)
    
    for i in 1:n_sites
        for j in 1:n_sites
            P[i,j] = 2.0 * sum(eigenvectors[i,k] * eigenvectors[j,k] for k in 1:n_occupied)
        end
    end
    
    return P
end

"""
    calculate_total_energy(energies::Vector{Float64}, n_occupied::Int)::Float64

Calculate the total electronic energy from orbital energies.
"""
function calculate_total_energy(energies::Vector{Float64}, n_occupied::Int)::Float64
    2.0 * sum(energies[1:n_occupied])
end

"""
    calculate_Huckel_Hamiltonian(system::MolecularSystem)::HuckelResult

Calculate the Hückel Hamiltonian and its eigenvalues/eigenvectors.
"""
function calculate_Huckel_Hamiltonian(system::MolecularSystem)::HuckelResult
    n_sites = length(system.atoms)
    H = zeros(Float64, n_sites, n_sites)
    
    # Diagonal elements (site energies)
    for i in 1:n_sites
        H[i,i] = system.atoms[i].site_energy
    end
    
    # Off-diagonal elements (hopping)
    for i in 1:n_sites
        for j in (i+1):n_sites
            if system.connectivity[i,j] == 1
                H[i,j] = H[j,i] = Constants.T
            end
        end
    end
    
    # Solve eigenvalue problem
    F = eigen(Symmetric(H))
    n_occupied = system.n_electrons ÷ 2
    
    density_matrix = calculate_density_matrix(F.vectors, n_occupied)
    total_energy = calculate_total_energy(F.values, n_occupied)
    
    return HuckelResult(H, F.values, F.vectors, n_occupied, density_matrix, total_energy)
end

# ============================================================================ #
# SCF with PPP 'Ohno' Exchange
# ============================================================================ #

"""
    calculate_Ohno_parameters(system::MolecularSystem)::Matrix{Float64}

Calculate the Ohno-Klopman matrix elements for electron-electron interactions.
"""
function calculate_Ohno_parameters(system::MolecularSystem)::Matrix{Float64}
    n_sites = length(system.atoms)
    V = zeros(Float64, n_sites, n_sites)
    
    for i in 1:n_sites
        V[i,i] = system.atoms[i].Hubbard_U
        
        for j in (i+1):n_sites
            r_ij = norm(system.atoms[i].position - system.atoms[j].position)
            U_i = system.atoms[i].Hubbard_U
            U_j = system.atoms[j].Hubbard_U
            V[i,j] = V[j,i] = Constants.ohnoconstant / sqrt(r_ij^2 + (2*Constants.ohnoconstant/(U_i + U_j))^2)
        end
    end
    return V
end

"""
    calculate_coulomb_matrix(P::Matrix{Float64}, V::Matrix{Float64}, 
                           system::MolecularSystem)::Matrix{Float64}

Calculate the Coulomb (J) matrix.
"""
function calculate_coulomb_matrix(P::Matrix{Float64}, V::Matrix{Float64}, 
                                system::MolecularSystem)::Matrix{Float64}
    n_sites = length(system.atoms)
    J = zeros(Float64, n_sites, n_sites)
    
    for i in 1:n_sites
        J[i,i] = sum((P[j,j] - Float64(system.atoms[j].nz)) * V[j,i] 
                     for j in 1:n_sites)
    end

    return J
end


"""
    calculate_exchange_matrix(P::Matrix{Float64}, V::Matrix{Float64})::Matrix{Float64}

Calculate the exchange (K) matrix.
"""
function calculate_exchange_matrix(P::Matrix{Float64}, V::Matrix{Float64}, system::MolecularSystem)::Matrix{Float64}
    n_sites = size(P, 1)
    K = zeros(Float64, n_sites, n_sites)
    
    for i in 1:n_sites
        for j in 1:n_sites
            K[i,j] = 0.5 * P[i,j] * V[i,j]
        end
        K[i,i] = (0.5 * P[i,i] - Float64(system.atoms[i].nz)) * V[i,i]
        # special case for diagonal element : picks up Hubbard U * extra electron density vs. valence
    end
    
    return K
end

"""
    run_SCF(system::MolecularSystem, huckel_result::HuckelResult;
            max_iterations::Int=10000, threshold::Float64=1e-12)::SCFResult

Perform SCF iterations to obtain the final electronic structure.
"""
function run_SCF(system::MolecularSystem, huckel_result::HuckelResult;
                max_iterations::Int=1000, threshold::Float64=1e-12)::SCFResult
    n_sites = length(system.atoms)
    n_occupied = system.n_electrons ÷ 2
    
    P_old = copy(huckel_result.density_matrix)
    H_core = copy(huckel_result.hamiltonian)
    V = calculate_Ohno_parameters(system)
    V_raw = V
    
    for iter in 1:max_iterations
        J = calculate_coulomb_matrix(P_old, V, system)
        K = calculate_exchange_matrix(P_old, V, system)
        F = H_core + J - K
        
        F_eig = eigen(Symmetric(F))
        P_new = calculate_density_matrix(F_eig.vectors, n_occupied)

        println("Iteration $iter:")

        # Debug output of matrices
#        println("\nCoulomb (J) Matrix:")
#        display(J)
#        println("\nExchange (K) Matrix:") 
#        display(K)
#        println("\nFock (F) Matrix:")
#        display(F)
#        println("\nNew Density Matrix:")
#        display(P_new)
        
        # HF ground state energy: sum of occupied orbital eigenvalues
        E_hf = calculate_total_energy(F_eig.values, n_occupied)
        println("HF ground state energy: $E_hf eV")


        if maximum(abs.(P_new - P_old)) < threshold
            return SCFResult(
                F, K, J, V_raw, F_eig.values, F_eig.vectors, n_occupied, P_new,
                calculate_total_energy(F_eig.values, n_occupied),
                iter, true
            )
        end
        
        P_old = copy(P_new)
    end
    
    @warn "SCF failed to converge after $max_iterations iterations"
    return SCFResult(
        F, K, J, V_raw, F_eig.values, F_eig.vectors, n_occupied, P_old,
        calculate_total_energy(F_eig.values, n_occupied),
        max_iterations, false
    )
end

"""
    calculate_HOMO_LUMO_exchange(scf_result::SCFResult, system::MolecularSystem)::Float64

Calculate the HOMO-LUMO exchange integral from the converged SCF results.
This represents the exchange interaction between HOMO and LUMO states.

Returns the exchange integral in eV.
"""
function calculate_HOMO_LUMO_exchange(scf_result::SCFResult, 
                                    system::MolecularSystem)::Float64
    K = scf_result.K
    C = scf_result.eigenvectors  # Using eigenvectors instead of Fock matrix
    
    n_sites = length(system.atoms)
    Krot = zeros(Float64, n_sites, n_sites)

    # Transform K matrix to MO basis
    for i in 1:n_sites
        for j in 1:n_sites
            for α in 1:n_sites
                for β in 1:n_sites
                    Krot[i,j] += C[α,i] * C[α,j] * K[α,β] * C[β,j] * C[β,i]
                end
            end
        end
    end

    homo_idx = scf_result.n_occupied
    lumo_idx = scf_result.n_occupied + 1
    
    # Check for degeneracy
    if abs(scf_result.energies[homo_idx] - scf_result.energies[homo_idx-1]) < 1e-5
        @warn "HOMO state shows degeneracy! Results may be unreliable."
    end

    return Krot[homo_idx,lumo_idx]
end

# ============================================================================ #
# Main Interface
# ============================================================================ #

"""
    read_geometry(filename::String)::MolecularSystem

Read molecular geometry from XYZ file and create a MolecularSystem instance.
"""
function read_geometry(filename::String)::MolecularSystem
    lines = readlines(filename)
    n_atoms_total = parse(Int, lines[1])
    
    temp_atoms = Atom[]
    
    for i in 2:(n_atoms_total+1)
        symbol, x, y = parse_xyz_line(lines[i])
        symbol == :H && continue  # skip hydrogen
        push!(temp_atoms, create_atom(symbol, x, y)) # Nb: only works for C and N!
    end
    
    positions = [atom.position for atom in temp_atoms]
    connectivity = calculate_connectivity(positions)
    

# now for the horrific hack to get the nz factor right
    atoms = map(enumerate(temp_atoms)) do (i, atom)
        n_bonds = count_bonds(connectivity, i)
        nz=1

        if atom.symbol == :N
            nz =  n_bonds == 3 ? 2 : 1
            site_energy = n_bonds == 3 ? Constants.ES_NPY : Constants.ES_NAZA
            Hubbard_U = n_bonds == 3 ? Constants.UN : Constants.UN_AZA
            return Atom(atom.symbol, atom.position, n_bonds, nz, site_energy, Hubbard_U)
        else
            return Atom(atom.symbol, atom.position, n_bonds, nz, atom.site_energy, atom.Hubbard_U)
        end
    end
    
    n_electrons = calculate_n_electrons(atoms)
    
    return MolecularSystem(atoms, connectivity, n_electrons)
end

"""
    run_ppp_calculation(xyz_file::String)::Tuple{MolecularSystem,HuckelResult,SCFResult}

Run a complete PPP calculation for a molecule specified in an XYZ file.
"""
function run_ppp_calculation(xyz_file::String)::Tuple{MolecularSystem,HuckelResult,SCFResult}
    system = read_geometry(xyz_file)
    display(system)
    
    Huckel_result = calculate_Huckel_Hamiltonian(system)
    display(Huckel_result)
  
    # TODO: Check for degeneracy
# @warn "HOMO/LUMO states show degeneracy! Results may be unreliable."
    
    SCF_result = run_SCF(system, Huckel_result)
    display(SCF_result)
    
    return system, Huckel_result, SCF_result
end

# ============================================================================ #
# Display Methods
# ============================================================================ #

"""
    Base.show(io::IO, mime::MIME"text/plain", system::MolecularSystem)

Custom display method for MolecularSystem type.
Provides detailed information about the molecular system when displayed in REPL or notebooks.
"""
function Base.show(io::IO, mime::MIME"text/plain", system::MolecularSystem)
    println(io, "MolecularSystem with $(length(system.atoms)) atoms and $(system.n_electrons) π electrons")
    
    # Print atom information in table format
    println(io, "\nAtom Properties:")
    println(io, " #  Atom    X (Å)    Y (Å)   Bonds  Site E (eV)  Hubbard U (eV)")
    println(io, "─"^65)  # Unicode box drawing character for line
    for (i, atom) in enumerate(system.atoms)
        @printf(io, "%2d  %2s   %7.3f  %7.3f    %d     %8.3f    %8.3f\n",
                i, atom.symbol, atom.position[1], atom.position[2],
                atom.n_bonds, atom.site_energy, atom.Hubbard_U)
    end
    
    # Print connectivity matrix
    println(io, "\nConnectivity Matrix:")
    println(io, "─"^(2 * size(system.connectivity, 1) + 1))
    for row in eachrow(system.connectivity)
        print(io, "│")
        for val in row
            print(io, " $(val)")
        end
        println(io)
    end
    println(io, "─"^(2 * size(system.connectivity, 1) + 1))
end

"""
    Base.show(io::IO, system::MolecularSystem)

Compact display method for MolecularSystem type.
Used when the object is shown in arrays or other containers.
"""
function Base.show(io::IO, system::MolecularSystem)
    print(io, "MolecularSystem($(length(system.atoms)) atoms, $(system.n_electrons) π electrons)")
end

"""
    Base.show(io::IO, mime::MIME"text/plain", result::HuckelResult)

Custom display method for HuckelResult type.
Shows detailed Hückel calculation results including energies, orbital occupations,
and key electronic properties.
"""
function Base.show(io::IO, mime::MIME"text/plain", result::HuckelResult)
    n_sites = size(result.hamiltonian, 1)
    n_occupied = result.n_occupied
 
    println(io, "Hückel Calculation Results")
    println(io, "========================")
    
    println(io, "\nHamiltonian Matrix (eV):")
    for i in 1:n_sites
        for j in 1:n_sites
            @printf(io, " %5.1f", result.hamiltonian[i,j])
        end
        println(io)
    end
    
    println(io, "\nOrbital Energies (eV):")
    for (i, energy) in enumerate(result.energies)
        occupation = i <= n_occupied ? 2 : 0
        @printf(io, "MO %3d: %8.3f eV (occupation: %d)\n", i, energy, occupation)
    end
    
    # Print HOMO-LUMO information
    @printf(io, "\nHOMO energy:     %8.3f eV", result.energies[n_occupied])
    @printf(io, "\nLUMO energy:     %8.3f eV", result.energies[n_occupied + 1])
    @printf(io, "\nHOMO-LUMO gap:   %8.3f eV", 
            result.energies[n_occupied + 1] - result.energies[n_occupied])
    
    println(io, "\n\nDensity Matrix:")
    for i in 1:n_sites
        for j in 1:n_sites
            @printf(io, " %4.1f", result.density_matrix[i,j])
        end
        println(io)
    end
    
    @printf(io, "\nTotal Electronic Energy: %8.3f eV", result.total_energy)
    println(io, "\n========================")
end

"""
    Base.show(io::IO, result::HuckelResult)

Compact display method for HuckelResult type.
"""
function Base.show(io::IO, result::HuckelResult)
    n_sites = size(result.hamiltonian, 1)
    @printf(io, "HuckelResult(%d atoms, E = %.3f eV)", n_sites, result.total_energy)
end

"""
    Base.show(io::IO, mime::MIME"text/plain", result::SCFResult)

Custom display method for SCFResult type.
Shows detailed SCF calculation results including convergence information,
energies, and matrices with 3-decimal formatting.
"""
function Base.show(io::IO, mime::MIME"text/plain", result::SCFResult)
    n_sites = size(result.F, 1)
    n_occupied = result.n_occupied
    
    println(io, "SCF Calculation Results")
    println(io, "=====================")
    
    # Convergence information
    println(io, "\nConvergence Status:")
    println(io, "  Converged:   ", result.converged ? "Yes" : "No")
    println(io, "  Iterations:  ", result.iterations)
    @printf(io, "  Final Energy: %.6f eV\n", result.total_energy)
    
    # Orbital energies
    println(io, "\nOrbital Energies (eV):")
    println(io, "  #    Energy    Occupation")
    println(io, "  ----------------------")
    for (i, energy) in enumerate(result.energies)
        occupation = i <= n_occupied ? 2 : 0
        @printf(io, "  %2d  %8.3f      %d\n", i, energy, occupation)
    end
    
    # HOMO-LUMO information
    println(io, "\nHOMO-LUMO Properties:")
    @printf(io, "\n  HOMO energy:     %8.3f eV", result.energies[n_occupied])
    @printf(io, "\n  LUMO energy:     %8.3f eV", result.energies[n_occupied + 1])
    @printf(io, "\n  HOMO-LUMO gap:   %8.3f eV", 
            result.energies[n_occupied + 1] - result.energies[n_occupied])
    @printf(io,"\n")
    
    # Matrix outputs
    println(io, "\nFinal Coulomb (J) Matrix Diagonal:")
    for i in 1:n_sites
        @printf(io, " %8.3f\n", result.J[i,i])
    end
    
    println(io, "\nFinal Exchange (K) Matrix:")
    for i in 1:n_sites
        for j in 1:n_sites
            @printf(io, " %5.2f", result.K[i,j])
        end
        println(io)
    end
    
    println(io, "\n\nFinal Fock Matrix:")
    for i in 1:n_sites
        for j in 1:n_sites
            @printf(io, " %5.2f", result.F[i,j])
        end
        println(io)
    end
    
    println(io, "\nFinal Density Matrix:")
    for i in 1:n_sites
        for j in 1:n_sites
            @printf(io, " %5.2f", result.density_matrix[i,j])
        end
        println(io)
    end
end

"""
    Base.show(io::IO, result::SCFResult)

Compact display method for SCFResult type.
"""
function Base.show(io::IO, result::SCFResult)
    status = result.converged ? "converged" : "not converged"
    @printf(io, "SCFResult(%d iterations, %s, E = %.6f eV)", 
            result.iterations, status, result.total_energy)
end

include("CI.jl")

export read_geometry, run_ppp_calculation


end # module 
