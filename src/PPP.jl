module PPP

using LinearAlgebra
using StaticArrays
using Printf

# Include model definitions
include("Models.jl")

const Ha2eV = 27.21138505
const eV2Ha = 1 / Ha2eV # ~= 0.0367

# Re-export model types and functions
export AbstractModel, Bedogni2024Model, Jorner2024Model, t, γ, t_ii, t_ij, γ_ii, γ_i# ============================================================================ #
# Types
# ============================================================================ #
"""
    Atom

Represents an atom in the π-conjugated system.

Fields:
- `symbol::Symbol`: Atomic symbol (e.g., :C, :N)
- `position::SVector{3,Float64}`: 3D coordinates (x,y,z) in Angstroms
- `n_bonds::Int`: Number of bonds in π structure
- `nz::Int`: Number of valence electrons; terrible hack for Nitrogen
- `site_energy::Float64`: Site energy in eV
- `Hubbard_U::Float64`: Hubbard U parameter in eV
"""
struct Atom
    symbol::Symbol
    position::SVector{3,Float64}
    n_bonds::Int
    nz::Int
    site_energy::Float64
    Hubbard_U::Float64
    Z_eff::Float64
    n_number::Float64
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
    Hamiltonian::Matrix{Float64}
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

function create_atom(symbol::Symbol, x::Float64, y::Float64, z::Float64, model::AbstractModel)
    position = SVector{3,Float64}(x, y, z)
    
    if symbol == :C
        # Carbon: 2 bonds, nz=1, will be updated later based on connectivity
        return Atom(symbol, position, 2, 1, 0.0, 0.0, model.Z_EFF_C, model.N_C)
    elseif symbol == :N
        # Nitrogen: 0 bonds initially, nz=0, will be updated later based on connectivity
        return Atom(symbol, position, 0, 0, 0.0, 0.0, model.Z_EFF_NPY, model.N_NPY)
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
    symbol = Symbol(parts[1]) # Julia symbol, like :C, :N, etc.
    x = parse(Float64, parts[2])
    y = parse(Float64, parts[3])
    z = parse(Float64, parts[4])
    return symbol, x, y, z
end

"""
    calculate_connectivity(positions::Vector{SVector{3,Float64}})::Matrix{Int}

Calculate the connectivity matrix based on atomic positions.
Returns a symmetric matrix where 1 indicates bonded atoms and 0 indicates non-bonded atoms.
"""
function calculate_connectivity(positions::Vector{SVector{3,Float64}}, model::AbstractModel)::Matrix{Int}
    n_atoms = length(positions)
    connectivity = zeros(Int, n_atoms, n_atoms)
    
    cutoff = 1.4  # Gone back to hard coding this for now! 
    
    for i in 1:n_atoms
        for j in (i+1):n_atoms
            distance = norm(positions[i] - positions[j])
            if distance < cutoff
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
# Matrix Operations
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
function calculate_Huckel_Hamiltonian(system::MolecularSystem, model::AbstractModel)::HuckelResult
    n_sites = length(system.atoms)
    H = zeros(Float64, n_sites, n_sites)
    
    # Fill matrix using unified interface
    for i in 1:n_sites
        for j in 1:n_sites
            if system.connectivity[i,j] == 1 || i == j
                H[i,j] = t(model, system, i, j)
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
function calculate_Ohno_parameters(system::MolecularSystem, model::AbstractModel)::Matrix{Float64}
    n_sites = length(system.atoms)
    V = zeros(Float64, n_sites, n_sites)
    
    # Fill matrix using unified interface
    for i in 1:n_sites
        for j in 1:n_sites
            V[i,j] = γ(model, system, i, j)
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
            max_iterations::Int=1000, threshold::Float64=1e-12)::SCFResult

Perform SCF iterations to obtain the final electronic structure.
"""
function run_SCF(system::MolecularSystem, huckel_result::HuckelResult, model::AbstractModel;
                max_iterations::Int=1000, threshold::Float64=1e-12)::SCFResult
    n_sites = length(system.atoms)
    n_occupied = system.n_electrons ÷ 2

    P_old = copy(huckel_result.density_matrix)
    H_core = copy(huckel_result.Hamiltonian)
            V = calculate_Ohno_parameters(system, model)
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

# Check for convergence
# FIXME: Make a bit more general?
        if maximum(abs.(P_new - P_old)) < threshold
            return SCFResult(
                F, K, J, V_raw, F_eig.values, F_eig.vectors, n_occupied, P_new,
                calculate_total_energy(F_eig.values, n_occupied),
                iter, true
            )
        end
        
        # Mixing parameter α, to prevent charge sloshing. Arnoldi mixing or something?
#        α = 0.6
#        P_old = α * P_new + (1-α) * P_old
        P_old = copy(P_new)
    end
    
    @warn "SCF failed to converge after $max_iterations iterations. Returning the last iteration. BEWARE! DRAGONS!"
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
function read_geometry(filename::String, model::AbstractModel)::MolecularSystem
    lines = readlines(filename)
    n_atoms_total = parse(Int, lines[1])
    
    temp_atoms = Atom[]
    
    for i in 2:(n_atoms_total+1)
        symbol, x, y, z = parse_xyz_line(lines[i])
        symbol == :H && continue  # skip hydrogen
        push!(temp_atoms, create_atom(symbol, x, y, z, model)) # Nb: only works for C and N!
    end
    
    positions = [atom.position for atom in temp_atoms]
    connectivity = calculate_connectivity(positions, model)
    

    # Update atoms with correct bond counts and parameters
    atoms = map(enumerate(temp_atoms)) do (i, atom)
        n_bonds = count_bonds(connectivity, i)
        nz = atom.symbol == :N && n_bonds == 3 ? 2 : 1
        
        # Update Slater parameters for nitrogen based on bond count
        if atom.symbol == :N
            Z_eff = n_bonds == 3 ? model.Z_EFF_NPY : model.Z_EFF_NAZA
            n_number = n_bonds == 3 ? model.N_NPY : model.N_NAZA
        else
            Z_eff, n_number = atom.Z_eff, atom.n_number
        end
        
        return Atom(atom.symbol, atom.position, n_bonds, nz, atom.site_energy, atom.Hubbard_U, Z_eff, n_number)
    end
    
    n_electrons = calculate_n_electrons(atoms)
    
    return MolecularSystem(atoms, connectivity, n_electrons)
end

"""
    run_ppp_calculation(xyz_file::String)::Tuple{MolecularSystem,HuckelResult,SCFResult}

Run a complete PPP calculation for a molecule specified in an XYZ file.
"""
function run_ppp_calculation(xyz_file::String, model::AbstractModel)::Tuple{MolecularSystem,HuckelResult,SCFResult}
    system = read_geometry(xyz_file, model)
    display(system)
    
    Huckel_result = calculate_Huckel_Hamiltonian(system, model)
    display(Huckel_result)
  
    # TODO: Check for degeneracy
# @warn "HOMO/LUMO states show degeneracy! Results may be unreliable."
    
    SCF_result = run_SCF(system, Huckel_result, model)
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
    println(io, " #  Atom    X (Å)    Y (Å)   Z (Å)   Bonds  Site E (eV)  Hubbard U (eV)")
    println(io, "─"^65)  # Unicode box drawing character for line
    for (i, atom) in enumerate(system.atoms)
        @printf(io, "%2d  %2s   %7.3f  %7.3f  %7.3f    %d     %8.3f    %8.3f\n",
                i, atom.symbol, atom.position[1], atom.position[2], atom.position[3],
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
    n_sites = size(result.Hamiltonian, 1)
    n_occupied = result.n_occupied
 
    println(io, "Hückel Calculation Results")
    println(io, "========================")
    
    println(io, "\nHamiltonian Matrix (eV):")
    for i in 1:n_sites
        for j in 1:n_sites
            @printf(io, " %5.1f", result.Hamiltonian[i,j])
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
    n_sites = size(result.Hamiltonian, 1)
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

