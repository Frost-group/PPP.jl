using Test
using PPP
using Printf
using StaticArrays
using CSV, DataFrames

#Calculate ground state energy for Ethene with Zhang 2011 parameters, also run CIS calculations for each bond length

#NOT FINISHED!!!!!!

# Run one PPP calculation for CIS calculation
azulene_system, huckel_result, scf_result_azulene = PPP.run_ppp_calculation("examples/validate_zhang/azulene/06-azulene_pbe1pbe-cc-pvtz.xyz", Zhang2011Model())

singlet, triplet = run_cis_calculation(azulene_system, scf_result_azulene)

println("First singlet excitation energy (eV): ", singlet.energies[1])
println("First triplet excitation energy (eV): ", triplet.energies[1])

ST_gap = singlet.energies[1] - triplet.energies[1]
println("Singlet-Triplet gap (eV): ", ST_gap)