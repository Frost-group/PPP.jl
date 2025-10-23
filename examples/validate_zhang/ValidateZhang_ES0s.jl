using Test
using PPP
using Printf
using StaticArrays
using CSV, DataFrames

#Calculate ground state energy for Ethene with Zhang 2011 parameters 

struct ES0_result
    rs::Vector{Float64}
    t_iis::Vector{Float64}
    t_ijs::Vector{Float64}
    gamma_iis::Vector{Float64}
    gamma_ijs::Vector{Float64}
    PPPs::Vector{Float64}
    Huckels::Vector{Float64}
    E_S0s::Vector{Float64}
end

@. E_S0(t_ii, t_ij, γ_ii, γ_ij) = 2t_ii + (γ_ii + γ_ij)/2 - sqrt(16t_ij^2 + (γ_ii - γ_ij)^2)/2

# E_S1(t_ii, γ_ii) = 2t_ii + γ_ii

# E_S2(t_ii, t_ij, γ_ii, γ_ij) = 2t_ii + (γ_ii + γ_ij)/2 + sqrt(16t_ij^2 + (γ_ii - γ_ij)^2)/2

# E_T1(t_ii, γ_ij) = 2t_ii + γ_ij

ethene_system = read_geometry("molecules/ethene.xyz", Zhang2011Model())

# Local helper function to deal with kruft of immutability in PPP
function update_system_atom_position(system::PPP.MolecularSystem, atom_idx::Int, new_position::SVector{3,Float64})
    new_atoms = copy(system.atoms)
    old_atom = new_atoms[atom_idx]
    new_atoms[atom_idx] = PPP.Atom(
        old_atom.symbol, 
        new_position, 
        old_atom.n_bonds, 
        old_atom.nz, 
        old_atom.site_energy, 
        old_atom.Hubbard_U, 
        old_atom.Z_eff, 
        old_atom.n_number
    )
    
    # Recalculate connectivity based on new positions
    new_positions = [atom.position for atom in new_atoms]
    new_connectivity = PPP.calculate_connectivity(new_positions, Zhang2011Model())
    
    return PPP.MolecularSystem(new_atoms, new_connectivity, system.n_electrons)
end

df = CSV.read("examples/validate_zhang/Zhang2011Ethene.dat", DataFrame;
              delim=' ',
              ignorerepeated=true,
              comment="#",
              header=false)

println(first(df, 5))

# First column of data file contains bond lengths (in Angstroms)
rs = df[:, 1] #Vector{Float64}

function vary_bond_lengths(system::PPP.MolecularSystem, rs::AbstractVector{Float64})::ES0_result
    t_iis = []
    t_ijs = []
    gamma_iis = []
    gamma_ijs = []
    PPPs = []
    Huckels = []
    E_S0s = []
    # rs = 1.30:0.01:1.60
    for r in rs
        new_position = system.atoms[1].position + SVector{3,Float64}(r, 0.0, 0.0)
        system = update_system_atom_position(system, 2, new_position)

        append!(t_iis, t_ii(Zhang2011Model(), system, 1))
        append!(t_ijs, t_ij(Zhang2011Model(), system, 1, 2))
        append!(gamma_iis, γ_ii(Zhang2011Model(), system, 1))
        append!(gamma_ijs, γ_ij(Zhang2011Model(), system, 1, 2))

        Huckel_result = PPP.calculate_Huckel_Hamiltonian(system, Zhang2011Model())
        SCF_result = PPP.run_SCF(system, Huckel_result, Zhang2011Model())

        append!(PPPs,SCF_result.total_energy)
        append!(Huckels,Huckel_result.total_energy)
    end

    E_S0s = E_S0(t_iis, t_ijs, gamma_iis, gamma_ijs)

    @info "Ground state energies from calculation:" E_S0s
    return ES0_result(rs, t_iis, t_ijs, gamma_iis, gamma_ijs, PPPs, Huckels, E_S0s)
end

ES0_ethene = vary_bond_lengths(ethene_system, rs)

df_ES0_ethene = DataFrame(
    rs = ES0_ethene.rs,
    t_iis = ES0_ethene.t_iis,
    t_ijs = ES0_ethene.t_ijs,
    gamma_iis = ES0_ethene.gamma_iis,
    gamma_ijs = ES0_ethene.gamma_ijs,
    Huckels = ES0_ethene.Huckels,
    PPPs = ES0_ethene.PPPs,
    E_S0s = ES0_ethene.E_S0s
)

CSV.write("examples/validate_zhang/ValidateZhang_ES0s_output.csv", df_ES0_ethene)

