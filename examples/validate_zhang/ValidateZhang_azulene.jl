using Test
using PPP
using Printf
using StaticArrays
using CSV, DataFrames

#Calculate ground state energy for Ethene with Zhang 2011 parameters, also run CIS calculations for each bond length

#NOT FINISHED!!!!!

struct Energy_results
    rs::Vector{Float64} #distance between C atoms
    t_iis::Vector{Float64}
    t_ijs::Vector{Float64}
    gamma_iis::Vector{Float64}
    gamma_ijs::Vector{Float64}
    PPPs::Vector{Float64} #total energy from SCF
    Huckels::Vector{Float64} #total energy from Huckel
    ES0s_PPP::Vector{Float64}
    ES1s_PPP::Vector{Float64}
    # ES2s_PPP::Vector{Float64} taken out as not used in Zhang calculations
    ET1s_PPP::Vector{Float64}
    ESs_CIS::Vector{Float64} #singlet CIS excitation energies
    ETs_CIS::Vector{Float64} #triplet CIS excitation energies
end

@. E_S0(t_ii, t_ij, γ_ii, γ_ij) = 2t_ii + (γ_ii + γ_ij)/2 - sqrt(16t_ij^2 + (γ_ii - γ_ij)^2)/2

@. E_S1(t_ii, γ_ii) = 2t_ii + γ_ii

@. E_S2(t_ii, t_ij, γ_ii, γ_ij) = 2t_ii + (γ_ii + γ_ij)/2 + sqrt(16t_ij^2 + (γ_ii - γ_ij)^2)/2 #not used

@. E_T1(t_ii, γ_ij) = 2t_ii + γ_ij

# Run one PPP calculation for CIS calculation
ethene_system, huckel_result, scf_result_ethene = PPP.run_ppp_calculation("molecules/ethene.xyz", Zhang2011Model())

singlet, triplet = run_cis_calculation(ethene_system, scf_result_ethene)

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

df_zhang_CASPT2 = CSV.read("examples/validate_zhang/Zhang2011Ethene.dat", DataFrame;
              delim=' ',
              ignorerepeated=true,
              comment="#",
              header=false)

# println(first(df, 5))

# First column of data file contains bond lengths (in Angstroms)
rs = df_zhang_CASPT2[:, 1] #Vector{Float64}
ES0s_Zhang = df_zhang_CASPT2[:, 2]
ES0s_Zhang .*= Ha2eV #Convert from Hartree to eV

function vary_bond_lengths(system::PPP.MolecularSystem, rs::AbstractVector{Float64})::Energy_results
    t_iis = []
    t_ijs = []
    gamma_iis = []
    gamma_ijs = []
    PPPs = []
    Huckels = []
    ES0s_PPP = []
    ES1s_PPP = []
    # E_S2s = []
    ET1s_PPP = []  
    ESs_CIS = []
    ETs_CIS = []
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

        singlet, triplet = run_cis_calculation(system, SCF_result)

        # FIXME! only works for Ethene because there's only one excitation energy from the CIS calculations

        append!(ESs_CIS, singlet.energies[1]) # vertical excitation energy
        append!(ETs_CIS, triplet.energies[1])

        append!(PPPs,SCF_result.total_energy)
        append!(Huckels,Huckel_result.total_energy)
    end

    ES0s_PPP = E_S0(t_iis, t_ijs, gamma_iis, gamma_ijs)
    ES1s_PPP = E_S1(t_iis, gamma_iis)
    # E_S2s = E_S2(t_iis, t_ijs, gamma_iis, gamma_ijs)
    ET1s_PPP = E_T1(t_iis, gamma_ijs)

    # @info "Ground state energies from calculation:" ES0s_PPP
    return Energy_results(rs, t_iis, t_ijs, gamma_iis, gamma_ijs, PPPs, Huckels, ES0s_PPP, ES1s_PPP, ET1s_PPP, ESs_CIS, ETs_CIS)
end

Energies_ethene = vary_bond_lengths(ethene_system, rs)

# Energy results: PPP parameters, ground state energies from PPP SCF and Huckel, energy levels calculated from PPP parameters, CIS excitation energies
df_Energies_ethene = DataFrame(
    rs = Energies_ethene.rs,
    t_iis = Energies_ethene.t_iis,
    t_ijs = Energies_ethene.t_ijs,
    gamma_iis = Energies_ethene.gamma_iis,
    gamma_ijs = Energies_ethene.gamma_ijs,
    Huckels = Energies_ethene.Huckels,
    PPPs = Energies_ethene.PPPs,
    ES0s_PPP = Energies_ethene.ES0s_PPP,
    ES1s_PPP = Energies_ethene.ES1s_PPP,
    ET1s_PPP = Energies_ethene.ET1s_PPP,
    ESs_CIS = Energies_ethene.ESs_CIS,
    ETs_CIS = Energies_ethene.ETs_CIS
)

CSV.write("examples/validate_zhang/ValidateZhang_CIS_energy_results.csv", df_Energies_ethene)